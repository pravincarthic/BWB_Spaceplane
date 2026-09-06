"""Fusion 360 script: cyclic C4 microcavity groups placed on hand-picked upper-surface sketch points."""

import adsk.core
import adsk.fusion
# import adsk.cam
import csv
import math
import os
import tempfile
import time
import traceback

app = adsk.core.Application.get()
ui = app.userInterface

# Flip to -1 if a first test cut comes out inverted (cavity sitting proud of the skin).
DEPTH_SIGN = 1

DIALOG_TITLE = "Cavity Pattern C4"
SYMMETRY_LABEL = "C4"
DEFAULT_SKETCH_NAME = "c4_distr"
GROUP_CAP_PARAMETER = "cavityMaxGroups_C4"
GROUP_CAP_DEFAULT = 10.0
LOG_FOLDER = "CavityPattern_C4_Logs"
TOOL_BODY_PREFIX = "CavityTool_C4_"
PREVIEW_SKETCH_PREFIX = "CavityPreview_C4_"

MM_PER_CM = 10.0
DEG_PER_RAD = 180.0 / math.pi

DEFAULT_RING_RADIUS_MM = "15"
DEFAULT_PHASE_DEG = "0"
DEFAULT_MIN_SPACING_MM = "2"
DEFAULT_PENETRATION_MM = "0.05"
DEFAULT_DUP_DISTANCE_MM = "0.2"
DEFAULT_DUP_ANGLE_DEG = "2"
DEFAULT_SEARCH_TOLERANCE_MM = "50"
DEFAULT_SEAM_TOLERANCE_MM = "0.005"
DEFAULT_MAX_CAVITIES = "2000"
DEFAULT_GROUP_SPACING_MM = "100"

MIN_ALLOWED_AREA_MM2 = 100.0
MIN_ALLOWED_FACE_COUNT = 3
EXCLUSION_TOLERANCE_MM = 0.01
CURVE_JOIN_TOLERANCE_MM = 0.01
CENTRE_MERGE_FRACTION = 0.02

MANIFEST_COLUMNS = [
    "cavity_id", "group_id", "local_member_id", "symmetry_operation",
    "x", "y", "z", "normal_x", "normal_y", "normal_z",
    "streamwise_x", "streamwise_y", "streamwise_z", "rotation_angle",
    "accepted", "rejection_reason", "nearest_face_distance_mm",
]


def vec_add(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def vec_sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def vec_scale(a, s):
    return [a[0] * s, a[1] * s, a[2] * s]


def vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vec_cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def vec_length(a):
    return math.sqrt(vec_dot(a, a))


def vec_distance(a, b):
    return vec_length(vec_sub(a, b))


def vec_unit(a):
    n = vec_length(a)
    if n < 1.0e-12:
        return [0.0, 0.0, 0.0]
    return [a[0] / n, a[1] / n, a[2] / n]


# Placement transforms are composed as M = T * A * S: the local symmetry operation S acts on the
# seed first, then the surface-normal alignment rotation A, then the translation T onto the snapped
# point.  Everything stays in plain nested lists and only becomes a Matrix3D at the point of use.
def mat3_identity():
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def mat3_mul(a, b):
    out = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    for i in range(3):
        for j in range(3):
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]
    return out


def mat3_apply(m, v):
    return [m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]]


def mat3_determinant(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def mat3_rotation_z(angle):
    c = math.cos(angle)
    s = math.sin(angle)
    return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]


def mat3_axis_angle(axis, angle):
    u = vec_unit(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    t = 1.0 - c
    x, y, z = u[0], u[1], u[2]
    return [[t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c]]


def mat3_reflection(normal):
    """Reflection in the plane through the origin with the given unit normal (determinant -1)."""
    n = vec_unit(normal)
    x, y, z = n[0], n[1], n[2]
    return [[1.0 - 2.0 * x * x, -2.0 * x * y, -2.0 * x * z],
            [-2.0 * x * y, 1.0 - 2.0 * y * y, -2.0 * y * z],
            [-2.0 * x * z, -2.0 * y * z, 1.0 - 2.0 * z * z]]


def mat3_shortest_arc(source, target):
    """Minimal rotation carrying unit vector source onto unit vector target."""
    a = vec_unit(source)
    b = vec_unit(target)
    if vec_length(a) < 0.5 or vec_length(b) < 0.5:
        return mat3_identity()
    d = vec_dot(a, b)
    if d >= 1.0 - 1.0e-12:
        return mat3_identity()
    if d <= -1.0 + 1.0e-12:
        reference = [1.0, 0.0, 0.0]
        if abs(vec_dot(reference, a)) > 0.9:
            reference = [0.0, 1.0, 0.0]
        axis = vec_unit(vec_cross(a, reference))
        if vec_length(axis) < 0.5:
            axis = vec_unit(vec_cross(a, [0.0, 0.0, 1.0]))
        return mat3_axis_angle(axis, math.pi)
    v = vec_cross(a, b)
    k = [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
    kk = mat3_mul(k, k)
    f = 1.0 / (1.0 + d)
    out = mat3_identity()
    for i in range(3):
        for j in range(3):
            out[i][j] = out[i][j] + k[i][j] + kk[i][j] * f
    return out


def mat4_from_rotation_translation(rotation, translation):
    return [[rotation[0][0], rotation[0][1], rotation[0][2], translation[0]],
            [rotation[1][0], rotation[1][1], rotation[1][2], translation[1]],
            [rotation[2][0], rotation[2][1], rotation[2][2], translation[2]],
            [0.0, 0.0, 0.0, 1.0]]


def mat4_to_matrix3d(m):
    out = adsk.core.Matrix3D.create()
    for r in range(4):
        for c in range(4):
            out.setCell(r, c, m[r][c])
    return out


def to_point3d(v):
    return adsk.core.Point3D.create(v[0], v[1], v[2])


def from_point3d(p):
    return [p.x, p.y, p.z]


def from_vector3d(v):
    return [v.x, v.y, v.z]


class PlacementParameters:
    """User-configurable run settings; lengths in internal centimetres, angles in radians."""

    def __init__(self):
        self.sketch_name = DEFAULT_SKETCH_NAME
        self.group_spacing = 10.0
        self.ring_radius = 1.5
        self.phase_offset = 0.0
        self.include_centre = False
        self.min_edge_spacing = 0.2
        self.penetration_eps = 0.005
        self.duplicate_distance = 0.02
        self.duplicate_angle = 2.0 / DEG_PER_RAD
        self.search_tolerance = 5.0
        self.seam_tolerance = 0.0005
        self.exclusion_tolerance = EXCLUSION_TOLERANCE_MM / MM_PER_CM
        self.max_cavities = 2000
        self.max_groups = int(GROUP_CAP_DEFAULT)

    def describe(self):
        return [
            "sketch_name: " + str(self.sketch_name),
            "group_separation_mm: " + format_number(self.group_spacing * MM_PER_CM),
            "ring_radius_mm: " + format_number(self.ring_radius * MM_PER_CM),
            "phase_offset_deg: " + format_number(self.phase_offset * DEG_PER_RAD),
            "include_centre_cavity: " + str(self.include_centre),
            "min_edge_to_edge_spacing_mm: " + format_number(self.min_edge_spacing * MM_PER_CM),
            "penetration_epsilon_mm: " + format_number(self.penetration_eps * MM_PER_CM),
            "duplicate_distance_tolerance_mm: " + format_number(self.duplicate_distance * MM_PER_CM),
            "duplicate_angle_tolerance_deg: " + format_number(self.duplicate_angle * DEG_PER_RAD),
            "surface_search_tolerance_mm: " + format_number(self.search_tolerance * MM_PER_CM),
            "seam_tolerance_mm: " + format_number(self.seam_tolerance * MM_PER_CM),
            "max_cavity_count: " + str(self.max_cavities),
            "max_group_count_parameter: " + GROUP_CAP_PARAMETER + " = " + str(self.max_groups),
            "depth_sign: " + str(DEPTH_SIGN),
        ]


class SurfaceFrame:
    """Orthonormal tangent frame at a point snapped onto an allowed upper-surface face."""

    def __init__(self, origin, normal, u_axis, v_axis, face, nearest_distance, snap_mode):
        self.origin = origin
        self.normal = normal
        self.u_axis = u_axis
        self.v_axis = v_axis
        self.face = face
        self.nearest_distance = nearest_distance
        self.snap_mode = snap_mode

    def to_world(self, offset):
        return vec_add(self.origin,
                       vec_add(vec_add(vec_scale(self.u_axis, offset[0]),
                                       vec_scale(self.v_axis, offset[1])),
                               vec_scale(self.normal, offset[2])))


class SymmetryTransform:
    """One member of the local symmetry orbit: a 3x3 operation plus its in-frame offset."""

    def __init__(self, member_id, station_id, label, matrix3, offset):
        self.member_id = member_id
        self.station_id = station_id
        self.label = label
        self.matrix3 = matrix3
        self.offset = offset
        self.determinant = mat3_determinant(matrix3)

    def is_mirror(self):
        return self.determinant < 0.0


class CavityPlacement:
    """One candidate cavity: where it landed, how it is oriented, and whether it survived."""

    def __init__(self, cavity_id, group_id, symmetry):
        self.cavity_id = cavity_id
        self.group_id = group_id
        self.member_id = symmetry.member_id
        self.station_id = symmetry.station_id
        self.symmetry_operation = symmetry.label
        self.determinant = symmetry.determinant
        self.point = [0.0, 0.0, 0.0]
        self.normal = [0.0, 0.0, 0.0]
        self.u_axis = [0.0, 0.0, 0.0]
        self.depth_axis = [0.0, 0.0, 0.0]
        self.rotation_angle = 0.0
        self.transform4 = None
        self.accepted = False
        self.rejection_reason = ""
        self.nearest_face_distance = -1.0

    def reject(self, reason, nearest_distance=None):
        self.accepted = False
        self.rejection_reason = reason
        if nearest_distance is not None:
            self.nearest_face_distance = nearest_distance

    def manifest_row(self):
        distance_mm = ""
        if self.nearest_face_distance >= 0.0:
            distance_mm = format_number(self.nearest_face_distance * MM_PER_CM)
        return [
            self.cavity_id, self.group_id, self.member_id, self.symmetry_operation,
            format_number(self.point[0] * MM_PER_CM),
            format_number(self.point[1] * MM_PER_CM),
            format_number(self.point[2] * MM_PER_CM),
            format_number(self.normal[0]), format_number(self.normal[1]), format_number(self.normal[2]),
            format_number(self.u_axis[0]), format_number(self.u_axis[1]), format_number(self.u_axis[2]),
            format_number(self.rotation_angle * DEG_PER_RAD),
            "yes" if self.accepted else "no",
            self.rejection_reason,
            distance_mm,
        ]


class PlacementResult:
    """Aggregate outcome of one generation pass."""

    def __init__(self):
        self.placements = []
        self.groups_evaluated = 0
        self.duplicates_removed = 0
        self.cancelled = False
        self.centre_source = ""
        self.timings = {}
        self.notes = []

    def accepted(self):
        return [p for p in self.placements if p.accepted]

    def rejection_counts(self):
        counts = {}
        for p in self.placements:
            if p.accepted:
                continue
            reason = p.rejection_reason or "unspecified"
            counts[reason] = counts.get(reason, 0) + 1
        return counts

    def nearest_distance_summary(self):
        values = [p.nearest_face_distance for p in self.placements if p.nearest_face_distance >= 0.0]
        if not values:
            return "no surface distances recorded"
        return ("min " + format_number(min(values) * MM_PER_CM) + " mm, max "
                + format_number(max(values) * MM_PER_CM) + " mm")


class SpatialHash:
    """Uniform bucket grid so near-neighbour tests never degrade into all-pairs comparisons."""

    def __init__(self, cell_size):
        self.cell = max(float(cell_size), 1.0e-6)
        self.buckets = {}

    def key_of(self, point):
        c = self.cell
        return (int(math.floor(point[0] / c)), int(math.floor(point[1] / c)), int(math.floor(point[2] / c)))

    def insert(self, point, payload):
        self.buckets.setdefault(self.key_of(point), []).append(payload)

    def neighbours(self, point):
        kx, ky, kz = self.key_of(point)
        found = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.buckets.get((kx + dx, ky + dy, kz + dz))
                    if bucket:
                        found.extend(bucket)
        return found


def format_number(value):
    return "{0:.6f}".format(float(value))


RUN_LOG_LINES = []


def log_line(text):
    """Record one run-log line and mirror it to the Fusion Text Commands window."""
    RUN_LOG_LINES.append(str(text))
    try:
        app.log(str(text))
    except:  #pylint:disable=bare-except
        pass


def log_folder_path():
    """Resolve a writable log folder, falling back past Documents to the temp folder."""
    candidates = [os.path.join(os.path.expanduser("~"), "Documents", LOG_FOLDER),
                  os.path.join(os.path.expanduser("~"), LOG_FOLDER),
                  os.path.join(tempfile.gettempdir(), LOG_FOLDER)]
    for folder in candidates:
        try:
            os.makedirs(folder, exist_ok=True)
            probe = os.path.join(folder, "write_test.tmp")
            with open(probe, "w") as handle:
                handle.write("ok")
            os.remove(probe)
            return folder
        except:  #pylint:disable=bare-except
            continue
    return None


def write_run_log(log_path):
    """Write the accumulated run log; returns the path actually written, or None."""
    if log_path is None:
        return None
    try:
        with open(log_path, "w") as handle:
            handle.write("\n".join(RUN_LOG_LINES) + "\n")
        return log_path
    except:  #pylint:disable=bare-except
        return None


def ask_text(prompt, default):
    text, cancelled = ui.inputBox(prompt, DIALOG_TITLE, default)
    if cancelled:
        return None
    return text


def ask_yes_no(prompt):
    answer = ui.messageBox(prompt, DIALOG_TITLE,
                           adsk.core.MessageBoxButtonTypes.YesNoButtonType,
                           adsk.core.MessageBoxIconTypes.QuestionIconType)
    return answer == adsk.core.DialogResults.DialogYes


def ask_length(units, prompt, default_text, allow_zero=True):
    """Evaluate a user length through the units manager; the returned value is in internal cm."""
    while True:
        text = ask_text(prompt, default_text)
        if text is None:
            return None
        try:
            value = units.evaluateExpression(text, "mm")
        except:  #pylint:disable=bare-except
            ui.messageBox("Could not read '" + str(text) + "' as a length in mm.", DIALOG_TITLE)
            continue
        if value < 0.0 or (value == 0.0 and not allow_zero):
            ui.messageBox("Value must be greater than zero.", DIALOG_TITLE)
            continue
        return value


def ask_angle(units, prompt, default_text):
    """Evaluate a user angle through the units manager; the returned value is in radians."""
    while True:
        text = ask_text(prompt, default_text)
        if text is None:
            return None
        try:
            return units.evaluateExpression(text, "deg")
        except:  #pylint:disable=bare-except
            ui.messageBox("Could not read '" + str(text) + "' as an angle in degrees.", DIALOG_TITLE)


def ask_count(prompt, default_text, minimum=1):
    while True:
        text = ask_text(prompt, default_text)
        if text is None:
            return None
        try:
            value = int(round(float(text)))
        except:  #pylint:disable=bare-except
            ui.messageBox("Could not read '" + str(text) + "' as a whole number.", DIALOG_TITLE)
            continue
        if value < minimum:
            ui.messageBox("Value must be at least " + str(minimum) + ".", DIALOG_TITLE)
            continue
        return value


def ensure_group_cap_parameter(design):
    """Read cavityMaxGroups_C4, creating it at its default on first run."""
    try:
        parameter = design.userParameters.itemByName(GROUP_CAP_PARAMETER)
    except:  #pylint:disable=bare-except
        parameter = None
    if parameter is None:
        try:
            value_input = adsk.core.ValueInput.createByReal(GROUP_CAP_DEFAULT)
            parameter = design.userParameters.add(
                GROUP_CAP_PARAMETER, value_input, "",
                "Maximum number of C4 cavity groups accepted per run")
        except:  #pylint:disable=bare-except
            parameter = None
    if parameter is None:
        ui.messageBox("Could not create the user parameter " + GROUP_CAP_PARAMETER
                      + ". Enter the group cap for this run instead.", DIALOG_TITLE)
        return ask_count("Maximum number of groups for this run:", str(int(GROUP_CAP_DEFAULT)))
    try:
        return max(1, int(round(parameter.value)))
    except:  #pylint:disable=bare-except
        return int(GROUP_CAP_DEFAULT)


def collect_user_inputs(design):
    """Prompt for every run setting; returns None if the user cancels."""
    units = design.unitsManager
    params = PlacementParameters()

    sketch_name = ask_text("Name of the distribution sketch holding the group centres:",
                           DEFAULT_SKETCH_NAME)
    if sketch_name is None:
        return None
    params.sketch_name = sketch_name.strip()

    value = ask_length(units,
                       "Separation between cavity GROUPS along the distribution curve."
                       " Enter mm, or append m for metres (for example 250 or 0.25 m):",
                       DEFAULT_GROUP_SPACING_MM, allow_zero=False)
    if value is None:
        return None
    params.group_spacing = value

    value = ask_length(units, "Ring radius in mm:", DEFAULT_RING_RADIUS_MM, allow_zero=False)
    if value is None:
        return None
    params.ring_radius = value

    if params.group_spacing < 2.0 * params.ring_radius:
        message = ("Group separation is "
                   + format_number(params.group_spacing * MM_PER_CM)
                   + " mm but each group spans about "
                   + format_number(2.0 * params.ring_radius * MM_PER_CM)
                   + " mm across its ring, so neighbouring groups will interleave along the"
                   + " curve.\n\nContinue anyway?")
        if not ask_yes_no(message):
            return None

    value = ask_angle(units, "Angular phase offset in degrees:", DEFAULT_PHASE_DEG)
    if value is None:
        return None
    params.phase_offset = value

    params.include_centre = ask_yes_no("Include an extra cavity at the centre of each group?")

    value = ask_length(units, "Minimum edge-to-edge spacing in mm:", DEFAULT_MIN_SPACING_MM)
    if value is None:
        return None
    params.min_edge_spacing = value

    value = ask_length(units, "Penetration epsilon in mm (tool lift above the skin):",
                       DEFAULT_PENETRATION_MM)
    if value is None:
        return None
    params.penetration_eps = value

    value = ask_length(units, "Duplicate distance tolerance in mm:", DEFAULT_DUP_DISTANCE_MM)
    if value is None:
        return None
    params.duplicate_distance = value

    value = ask_angle(units, "Duplicate angle tolerance in degrees:", DEFAULT_DUP_ANGLE_DEG)
    if value is None:
        return None
    params.duplicate_angle = abs(value)

    value = ask_length(units, "Surface search tolerance in mm (how rough the sketch points may be):",
                       DEFAULT_SEARCH_TOLERANCE_MM, allow_zero=False)
    if value is None:
        return None
    params.search_tolerance = value

    value = ask_length(units, "Seam tolerance in mm (trim-boundary slack between adjacent faces):",
                       DEFAULT_SEAM_TOLERANCE_MM)
    if value is None:
        return None
    params.seam_tolerance = value

    count = ask_count("Maximum cavity count (safety ceiling):", DEFAULT_MAX_CAVITIES)
    if count is None:
        return None
    params.max_cavities = count

    cap = ensure_group_cap_parameter(design)
    if cap is None:
        return None
    params.max_groups = cap
    return params


def select_body(prompt):
    try:
        selection = ui.selectEntity(prompt, "SolidBodies")
    except:  #pylint:disable=bare-except
        return None
    if selection is None:
        return None
    return adsk.fusion.BRepBody.cast(selection.entity)


def collect_faces(prompt):
    """Repeated face picking, de-duplicated by entity token; Esc finishes the list."""
    faces = []
    tokens = set()
    duplicates = 0
    while True:
        try:
            selection = ui.selectEntity(prompt, "SolidFaces")
        except:  #pylint:disable=bare-except
            break
        if selection is None:
            break
        face = adsk.fusion.BRepFace.cast(selection.entity)
        if face is None:
            continue
        try:
            token = face.entityToken
        except:  #pylint:disable=bare-except
            token = "index_" + str(len(faces))
        if token in tokens:
            duplicates = duplicates + 1
            continue
        tokens.add(token)
        faces.append(face)
    return faces, duplicates


def total_face_area_mm2(faces):
    total = 0.0
    for face in faces:
        try:
            total = total + face.area * MM_PER_CM * MM_PER_CM
        except:  #pylint:disable=bare-except
            continue
    return total


def validate_model_entities(design):
    """Collect and sanity-check the seed body, target body and allowed / excluded faces."""
    root = design.rootComponent

    seed_body = select_body("Select the seed cavity tool body, then click OK")
    if seed_body is None:
        ui.messageBox("No seed cavity body was selected.", DIALOG_TITLE)
        return None
    if not seed_body.isSolid:
        ui.messageBox("The seed cavity must be a solid body.", DIALOG_TITLE)
        return None

    target_body = select_body("Select the spaceplane body to cut, then click OK")
    if target_body is None:
        ui.messageBox("No target body was selected.", DIALOG_TITLE)
        return None
    if target_body == seed_body:
        ui.messageBox("The target body and the seed cavity body must be different bodies.",
                      DIALOG_TITLE)
        return None

    allowed_faces, allowed_duplicates = collect_faces(
        "Select each allowed upper-surface face, press Esc when finished")
    if not allowed_faces:
        ui.messageBox("No allowed upper-surface faces were selected.", DIALOG_TITLE)
        return None

    allowed_area = total_face_area_mm2(allowed_faces)
    if allowed_area < MIN_ALLOWED_AREA_MM2:
        message = ("The selected allowed faces total only "
                   + format_number(allowed_area) + " mm2. That usually means a fillet, seam or"
                   + " sliver face was picked instead of a skin panel, and every group centre"
                   + " will collapse onto the same point.\n\nContinue anyway?")
        if not ask_yes_no(message):
            return None
    if len(allowed_faces) < MIN_ALLOWED_FACE_COUNT:
        message = ("Only " + str(len(allowed_faces)) + " allowed face(s) selected. A curved upper"
                   + " surface is normally split into several trimmed faces, so parts of the patch"
                   + " may be unreachable.\n\nContinue anyway?")
        if not ask_yes_no(message):
            return None

    excluded_faces = []
    excluded_duplicates = 0
    if ask_yes_no("Select excluded faces as well (keep-out panels)?"):
        excluded_faces, excluded_duplicates = collect_faces(
            "Select each excluded face, press Esc when finished")

    excluded_tokens = set()
    for face in excluded_faces:
        try:
            excluded_tokens.add(face.entityToken)
        except:  #pylint:disable=bare-except
            continue
    overlap = 0
    filtered = []
    for face in allowed_faces:
        try:
            token = face.entityToken
        except:  #pylint:disable=bare-except
            token = None
        if token is not None and token in excluded_tokens:
            overlap = overlap + 1
            continue
        filtered.append(face)
    if not filtered:
        ui.messageBox("Every allowed face was also marked as excluded.", DIALOG_TITLE)
        return None

    return {
        "root": root,
        "seed_body": seed_body,
        "target_body": target_body,
        "seed_name": seed_body.name,
        "target_name": target_body.name,
        "allowed_faces": filtered,
        "excluded_faces": excluded_faces,
        "allowed_area_mm2": allowed_area,
        "allowed_duplicates": allowed_duplicates,
        "excluded_duplicates": excluded_duplicates,
        "allowed_excluded_overlap": overlap,
    }


def seed_anchor_point(seed_body):
    """Topmost seed vertex along the local depth axis; ties broken deterministically."""
    best = None
    for vertex in seed_body.vertices:
        p = vertex.geometry
        candidate = (p.z * DEPTH_SIGN, -p.x, -p.y)
        if best is None or candidate > best[0]:
            best = (candidate, [p.x, p.y, p.z])
    if best is None:
        raise ValueError("The seed cavity body has no vertices to anchor against.")
    return best[1]


def seed_footprint_radius(seed_body, anchor):
    """Largest in-plane distance from the anchor axis, used for edge-to-edge spacing."""
    radius = 0.0
    for vertex in seed_body.vertices:
        p = vertex.geometry
        distance = math.hypot(p.x - anchor[0], p.y - anchor[1])
        if distance > radius:
            radius = distance
    return radius


def build_local_symmetry_orbit(params):
    """C4 orbit: four pure rotations of the seed about the group normal, plus an optional centre copy."""
    seed_offset = [params.ring_radius * math.cos(params.phase_offset),
                   params.ring_radius * math.sin(params.phase_offset),
                   0.0]
    members = []
    if params.include_centre:
        members.append(SymmetryTransform(0, 0, "centre", mat3_identity(), [0.0, 0.0, 0.0]))
    for k in range(4):
        operation = mat3_rotation_z(k * math.pi / 2.0)
        offset = mat3_apply(operation, seed_offset)
        members.append(SymmetryTransform(len(members), k + 1,
                                         "C4_rot_" + str(k * 90) + "deg", operation, offset))
    return members


def find_sketch(design, name):
    root = design.rootComponent
    sketch = root.sketches.itemByName(name)
    if sketch is not None:
        return sketch
    for component in design.allComponents:
        sketch = component.sketches.itemByName(name)
        if sketch is not None:
            return sketch
    return None


def curve_segments(sketch):
    """World-space evaluator, arc length and endpoints for every curve in the distribution sketch."""
    segments = []
    curves = sketch.sketchCurves
    for index in range(curves.count):
        curve = curves.item(index)
        if curve is None or not curve.isValid:
            continue
        try:
            evaluator = curve.worldGeometry.evaluator
            ok, start_param, end_param = evaluator.getParameterExtents()
            if not ok:
                continue
            ok, length = evaluator.getLengthAtParameter(start_param, end_param)
            if not ok or length <= 0.0:
                continue
            ok, start_point = evaluator.getPointAtParameter(start_param)
            if not ok:
                continue
            ok, end_point = evaluator.getPointAtParameter(end_param)
            if not ok:
                continue
        except:  #pylint:disable=bare-except
            continue
        segments.append({"evaluator": evaluator,
                         "start_param": start_param,
                         "end_param": end_param,
                         "length": length,
                         "start": from_point3d(start_point),
                         "end": from_point3d(end_point)})
    return segments


def build_curve_chains(segments, join_tolerance):
    """Order curves head to tail so a multi-segment line is walked as one continuous path."""
    used = set()
    chains = []

    def touches(a, b):
        return vec_distance(a, b) <= join_tolerance

    def tail_of(entry):
        segment = segments[entry[0]]
        return segment["start"] if entry[1] else segment["end"]

    def head_of(entry):
        segment = segments[entry[0]]
        return segment["end"] if entry[1] else segment["start"]

    def find_unused(point):
        for index in range(len(segments)):
            if index in used:
                continue
            segment = segments[index]
            if touches(segment["start"], point) or touches(segment["end"], point):
                return index
        return None

    for seed in range(len(segments)):
        if seed in used:
            continue
        used.add(seed)
        chain = [(seed, False)]
        while True:
            tail = tail_of(chain[-1])
            index = find_unused(tail)
            if index is None:
                break
            used.add(index)
            chain.append((index, not touches(segments[index]["start"], tail)))
        while True:
            head = head_of(chain[0])
            index = find_unused(head)
            if index is None:
                break
            used.add(index)
            chain.insert(0, (index, touches(segments[index]["start"], head)))
        chains.append(chain)
    return chains


def point_at_arc_length(segment, local_length, flipped):
    """Evaluate one curve at an arc length from the chain head, honouring traversal direction."""
    distance = segment["length"] - local_length if flipped else local_length
    evaluator = segment["evaluator"]
    if evaluator is not None:
        try:
            ok, parameter = evaluator.getParameterAtLength(segment["start_param"], distance)
            if ok:
                ok, point = evaluator.getPointAtParameter(parameter)
                if ok:
                    return from_point3d(point)
        except:  #pylint:disable=bare-except
            pass
    if segment["length"] <= 0.0:
        return list(segment["start"])
    fraction = max(0.0, min(1.0, distance / segment["length"]))
    return vec_add(segment["start"],
                   vec_scale(vec_sub(segment["end"], segment["start"]), fraction))


def sample_chain(segments, chain, spacing):
    """Group centres spaced by arc length along one chain, starting at its head."""
    total = 0.0
    for index, flipped in chain:
        total = total + segments[index]["length"]
    if total <= 0.0 or spacing <= 0.0:
        return [], 0.0
    epsilon = max(total, spacing) * 1.0e-9
    count = int(math.floor((total + epsilon) / spacing)) + 1
    points = []
    cursor = 0
    travelled = 0.0
    for step in range(count):
        target = step * spacing
        while (cursor < len(chain) - 1
               and target > travelled + segments[chain[cursor][0]]["length"]):
            travelled = travelled + segments[chain[cursor][0]]["length"]
            cursor = cursor + 1
        index, flipped = chain[cursor]
        segment = segments[index]
        local = max(0.0, min(segment["length"], target - travelled))
        point = point_at_arc_length(segment, local, flipped)
        if point is not None:
            points.append(point)
    return points, total


def merge_coincident_centres(points, tolerance):
    """Drop samples landing on the same spot, such as a shared vertex or a closed loop wrap."""
    grid = SpatialHash(max(tolerance, 1.0e-6))
    kept = []
    for point in points:
        duplicate = False
        for other in grid.neighbours(point):
            if vec_distance(point, other) <= tolerance:
                duplicate = True
                break
        if duplicate:
            continue
        grid.insert(point, point)
        kept.append(point)
    return kept


def standalone_sketch_points(sketch):
    """Sketch points owned by no curve, so hand-placed markers still count as group centres."""
    points = []
    for index in range(sketch.sketchPoints.count):
        point = sketch.sketchPoints.item(index)
        if point is None or not point.isValid:
            continue
        local = point.geometry
        if index == 0 and abs(local.x) < 1.0e-9 and abs(local.y) < 1.0e-9 and abs(local.z) < 1.0e-9:
            continue
        try:
            connected = point.connectedEntities
            if connected is not None and connected.count > 0:
                continue
        except:  #pylint:disable=bare-except
            pass
        try:
            points.append(from_point3d(point.worldGeometry))
        except:  #pylint:disable=bare-except
            continue
    return points


def build_global_group_centres(design, params):
    """Group centres walked along the distribution curves at the user-set group separation."""
    sketch = find_sketch(design, params.sketch_name)
    if sketch is None:
        ui.messageBox("No sketch named '" + params.sketch_name + "' was found in this document.",
                      DIALOG_TITLE)
        return None

    join_tolerance = CURVE_JOIN_TOLERANCE_MM / MM_PER_CM
    segments = curve_segments(sketch)
    chains = build_curve_chains(segments, join_tolerance)

    curve_centres = []
    total_length = 0.0
    for chain in chains:
        sampled, length = sample_chain(segments, chain, params.group_spacing)
        curve_centres.extend(sampled)
        total_length = total_length + length

    merge_tolerance = max(join_tolerance, params.group_spacing * CENTRE_MERGE_FRACTION)
    sampled_count = len(curve_centres)
    curve_centres = merge_coincident_centres(curve_centres, merge_tolerance)
    merged = sampled_count - len(curve_centres)

    loose_points = standalone_sketch_points(sketch)
    centres = curve_centres + loose_points

    if not centres:
        ui.messageBox("Sketch '" + params.sketch_name + "' holds no curve to distribute groups"
                      " along and no standalone sketch point.", DIALOG_TITLE)
        return None

    source = ("Group centres: " + str(len(curve_centres)) + " walked along "
              + str(len(chains)) + " curve chain(s) totalling "
              + format_number(total_length * MM_PER_CM) + " mm at "
              + format_number(params.group_spacing * MM_PER_CM) + " mm separation, plus "
              + str(len(loose_points)) + " standalone sketch point(s).")

    if len(centres) > params.max_groups:
        message = (source + "\n\nThat is " + str(len(centres)) + " groups in total, more than the"
                   " cap of " + str(params.max_groups) + " set by the user parameter "
                   + GROUP_CAP_PARAMETER + ".\n\nEdit that parameter in Modify, Change Parameters"
                   " to change the cap permanently, or raise the group separation to thin the"
                   " groups out.\n\nProcess all " + str(len(centres)) + " group centres now?")
        if not ask_yes_no(message):
            return None

    return {"centres": centres,
            "curve_count": len(segments),
            "chain_count": len(chains),
            "total_length": total_length,
            "curve_centres": len(curve_centres),
            "loose_points": len(loose_points),
            "merged": merged,
            "source": source,
            "sketch": sketch}


def snap_point_to_faces(point, faces, params):
    """Snap a world point onto the nearest allowed face.

    Returns (surface_point, normal, face, nearest_distance, snap_mode); the nearest raw distance is
    always reported, even when nothing is close enough to accept."""
    nearest_distance = -1.0
    best_on_face = None
    best_seam = None
    world = to_point3d(point)

    for face in faces:
        if not face.isValid:
            continue
        try:
            evaluator = face.evaluator
        except:  #pylint:disable=bare-except
            continue
        if evaluator is None:
            continue
        try:
            ok, parameter = evaluator.getParameterAtPoint(world)
        except:  #pylint:disable=bare-except
            continue
        if not ok:
            continue
        try:
            ok, surface_point = evaluator.getPointAtParameter(parameter)
        except:  #pylint:disable=bare-except
            continue
        if not ok:
            continue

        candidate = from_point3d(surface_point)
        distance = vec_distance(point, candidate)
        if nearest_distance < 0.0 or distance < nearest_distance:
            nearest_distance = distance

        try:
            on_face = evaluator.isParameterOnFace(parameter)
        except:  #pylint:disable=bare-except
            on_face = False

        record = (distance, face, parameter, candidate)
        if on_face:
            if distance <= params.search_tolerance and (best_on_face is None or distance < best_on_face[0]):
                best_on_face = record
        elif distance <= params.seam_tolerance and (best_seam is None or distance < best_seam[0]):
            best_seam = record

    chosen = best_on_face
    mode = "on_face"
    if chosen is None:
        chosen = best_seam
        mode = "seam_tolerance"
    if chosen is None:
        return None, None, None, nearest_distance, "none"

    distance, face, parameter, candidate = chosen
    try:
        ok, raw_normal = face.evaluator.getNormalAtParameter(parameter)
    except:  #pylint:disable=bare-except
        ok = False
        raw_normal = None
    if not ok or raw_normal is None:
        return None, None, face, nearest_distance, "normal_failed"

    normal = vec_unit(from_vector3d(raw_normal))
    if vec_length(normal) < 0.5:
        return None, None, face, nearest_distance, "normal_failed"
    try:
        if face.isParamReversed:
            normal = vec_scale(normal, -1.0)
    except:  #pylint:disable=bare-except
        pass
    return candidate, normal, face, nearest_distance, mode


def calculate_surface_frame(origin, normal, face, nearest_distance, snap_mode):
    """Right-handed tangent frame, falling back when the primary reference is parallel to the normal."""
    n = vec_unit(normal)
    reference = [1.0, 0.0, 0.0]
    if abs(vec_dot(reference, n)) > 0.9:
        reference = [0.0, 1.0, 0.0]
    if abs(vec_dot(reference, n)) > 0.9:
        reference = [0.0, 0.0, 1.0]
    u = vec_sub(reference, vec_scale(n, vec_dot(reference, n)))
    if vec_length(u) < 1.0e-9:
        u = vec_sub([0.0, 0.0, 1.0], vec_scale(n, vec_dot([0.0, 0.0, 1.0], n)))
    u = vec_unit(u)
    if vec_length(u) < 0.5:
        return None
    v = vec_unit(vec_cross(n, u))
    if vec_length(v) < 0.5:
        return None
    return SurfaceFrame(origin, n, u, v, face, nearest_distance, snap_mode)


def project_candidate_to_surface(point, faces, params):
    """Snap one candidate point and build its own surface frame."""
    origin, normal, face, nearest_distance, mode = snap_point_to_faces(point, faces, params)
    if origin is None or normal is None:
        return None, nearest_distance, mode
    frame = calculate_surface_frame(origin, normal, face, nearest_distance, mode)
    if frame is None:
        return None, nearest_distance, "degenerate_frame"
    return frame, nearest_distance, mode


def orient_cavity_tool(frame, symmetry, anchor, params):
    """Compose M = T * A * S for one member and return it with its rotation and arc angle."""
    depth_local = [0.0, 0.0, float(DEPTH_SIGN)]
    align = mat3_shortest_arc(depth_local, frame.normal)
    rotation = mat3_mul(align, symmetry.matrix3)
    target = vec_add(frame.origin, vec_scale(frame.normal, params.penetration_eps))
    translation = vec_sub(target, mat3_apply(rotation, anchor))
    transform = mat4_from_rotation_translation(rotation, translation)
    arc = math.acos(max(-1.0, min(1.0, vec_dot(depth_local, frame.normal))))
    depth_axis = vec_unit(mat3_apply(rotation, depth_local))
    return transform, depth_axis, arc, target


def hits_excluded_face(point, excluded_faces, params):
    """True when the snapped point sits on one of the excluded keep-out faces."""
    if not excluded_faces:
        return False
    world = to_point3d(point)
    for face in excluded_faces:
        if not face.isValid:
            continue
        try:
            evaluator = face.evaluator
            ok, parameter = evaluator.getParameterAtPoint(world)
            if not ok:
                continue
            ok, surface_point = evaluator.getPointAtParameter(parameter)
            if not ok:
                continue
            if vec_distance(point, from_point3d(surface_point)) > params.exclusion_tolerance:
                continue
            if evaluator.isParameterOnFace(parameter):
                return True
        except:  #pylint:disable=bare-except
            continue
    return False


def apply_exclusion_rules(placement, excluded_faces, spacing_hash, spacing_radius, params):
    """Keep-out faces and minimum edge-to-edge spacing; returns a reason string or None."""
    if hits_excluded_face(placement.point, excluded_faces, params):
        return "on_excluded_face"
    for other in spacing_hash.neighbours(placement.point):
        if other.group_id == placement.group_id and other.station_id == placement.station_id:
            continue
        if vec_distance(placement.point, other.point) < spacing_radius:
            return "spacing_violation"
    return None


def remove_duplicate_placements(placements, params):
    """Flag co-located, co-aligned placements of matching handedness using a spatial hash grid."""
    grid = SpatialHash(max(params.duplicate_distance, 1.0e-6))
    cosine_limit = math.cos(min(math.pi, max(0.0, params.duplicate_angle)))
    removed = 0
    for placement in placements:
        if not placement.accepted:
            continue
        is_duplicate = False
        for other in grid.neighbours(placement.point):
            if vec_distance(placement.point, other.point) > params.duplicate_distance:
                continue
            if (placement.determinant < 0.0) != (other.determinant < 0.0):
                continue
            if vec_dot(placement.depth_axis, other.depth_axis) < cosine_limit:
                continue
            is_duplicate = True
            break
        if is_duplicate:
            placement.reject("duplicate")
            removed = removed + 1
        else:
            grid.insert(placement.point, placement)
    return removed


def generate_all_placements(entities, params, orbit, anchor, seed_radius, centres):
    """Dry-run pass: snap, orient, filter and de-duplicate every candidate cavity."""
    result = PlacementResult()
    started = time.time()
    allowed_faces = entities["allowed_faces"]
    excluded_faces = entities["excluded_faces"]
    spacing_radius = 2.0 * seed_radius + params.min_edge_spacing
    spacing_hash = SpatialHash(max(spacing_radius, 1.0e-6))

    progress = ui.createProgressDialog()
    progress.isBackgroundTranslucent = False
    progress.show(DIALOG_TITLE, "Evaluating group %v of %m", 0, len(centres), 0)

    cavity_id = 0
    accepted_count = 0
    try:
        for group_id, centre in enumerate(centres):
            if progress.wasCancelled:
                result.cancelled = True
                break
            result.groups_evaluated = result.groups_evaluated + 1

            group_frame, group_distance, group_mode = project_candidate_to_surface(
                centre, allowed_faces, params)

            for symmetry in orbit:
                placement = CavityPlacement(cavity_id, group_id, symmetry)
                cavity_id = cavity_id + 1
                result.placements.append(placement)

                if group_frame is None:
                    placement.point = list(centre)
                    placement.reject("group_centre_off_surface_" + group_mode, group_distance)
                    continue

                candidate = group_frame.to_world(symmetry.offset)
                placement.point = list(candidate)

                member_frame, member_distance, member_mode = project_candidate_to_surface(
                    candidate, allowed_faces, params)
                if member_frame is None:
                    placement.reject("member_off_surface_" + member_mode, member_distance)
                    continue

                if accepted_count >= params.max_cavities:
                    placement.reject("max_cavity_count_reached", member_distance)
                    continue

                try:
                    transform, depth_axis, arc, target = orient_cavity_tool(
                        member_frame, symmetry, anchor, params)
                except:  #pylint:disable=bare-except
                    placement.reject("orientation_failed", member_distance)
                    continue

                placement.point = list(member_frame.origin)
                placement.normal = list(member_frame.normal)
                placement.u_axis = list(member_frame.u_axis)
                placement.depth_axis = depth_axis
                placement.rotation_angle = arc
                placement.transform4 = transform
                placement.nearest_face_distance = member_distance

                reason = apply_exclusion_rules(placement, excluded_faces, spacing_hash,
                                               spacing_radius, params)
                if reason is not None:
                    placement.reject(reason, member_distance)
                    continue

                placement.accepted = True
                accepted_count = accepted_count + 1
                spacing_hash.insert(placement.point, placement)

            progress.progressValue = group_id + 1
    finally:
        progress.hideDialog()

    result.duplicates_removed = remove_duplicate_placements(result.placements, params)
    result.timings["generate_seconds"] = time.time() - started
    return result


def summary_lines(result, entities, centre_count):
    accepted = result.accepted()
    lines = ["Candidates evaluated: " + str(len(result.placements)),
             "Group centres read: " + str(centre_count),
             "Groups evaluated: " + str(result.groups_evaluated),
             "Accepted placements: " + str(len(accepted)),
             "Duplicates removed: " + str(result.duplicates_removed)]
    counts = result.rejection_counts()
    other = sum(v for k, v in counts.items() if k != "duplicate")
    lines.append("Other rejections: " + str(other))
    for reason in sorted(counts.keys()):
        if reason == "duplicate":
            continue
        lines.append("  " + reason + ": " + str(counts[reason]))
    lines.append("Nearest surface distance across all candidates: " + result.nearest_distance_summary())
    if result.centre_source:
        lines.append(result.centre_source)
    lines.append("Allowed faces: " + str(len(entities["allowed_faces"]))
                 + " covering " + format_number(entities["allowed_area_mm2"]) + " mm2")
    if result.cancelled:
        lines.append("Run was cancelled from the progress dialog before all groups were evaluated.")
    return lines


def preview_placements(entities, result, params, centre_count, stamp):
    """Build the preview sketch and ask for confirmation before any boolean work."""
    root = entities["root"]
    accepted = result.accepted()
    sketch = None
    if not accepted:
        lines = summary_lines(result, entities, centre_count)
        lines.append("")
        lines.append("Nothing to cut, so no preview sketch was created. If the nearest distances"
                     " above are a few mm, raise the surface search tolerance. If they run to"
                     " hundreds of mm or metres, the wrong faces were selected.")
        ui.messageBox("\n".join(lines), DIALOG_TITLE)
        return False, None
    try:
        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.name = PREVIEW_SKETCH_PREFIX + stamp
        sketch.isComputeDeferred = True
        for placement in accepted:
            sketch.sketchPoints.add(sketch.modelToSketchSpace(to_point3d(placement.point)))
        sketch.isComputeDeferred = False
    except:  #pylint:disable=bare-except
        if sketch is not None:
            try:
                sketch.isComputeDeferred = False
            except:  #pylint:disable=bare-except
                pass
        log_line("Preview sketch failed: " + traceback.format_exc())

    lines = summary_lines(result, entities, centre_count)
    lines.append("")
    lines.append("Preview sketch created: " + (sketch.name if sketch is not None else "not created"))
    lines.append("")
    lines.append("The cut copies the selected seed cavity once per placement and moves each copy"
                 " onto its snapped point, which adds about two features per cavity to the"
                 " timeline plus one combine.")
    lines.append("")
    lines.append("Cut these " + str(len(accepted)) + " cavities into "
                 + entities["target_name"] + " now?")
    proceed = ask_yes_no("\n".join(lines))
    return proceed, sketch


def copy_seed_body(root, source_body, known_names):
    """Copy a selected workspace body in place; returns the new body or None."""
    feature = root.features.copyPasteBodies.add(source_body)
    body = None
    try:
        if feature.bodies.count > 0:
            body = feature.bodies.item(0)
    except:  #pylint:disable=bare-except
        body = None
    if body is None:
        for candidate in root.bRepBodies:
            if candidate.name not in known_names:
                body = candidate
                break
    return body


def move_body(root, body, matrix4):
    """Free-move a workspace body by an explicit 4x4."""
    collection = adsk.core.ObjectCollection.create()
    collection.add(body)
    matrix = mat4_to_matrix3d(matrix4)
    move_features = root.features.moveFeatures
    try:
        move_input = move_features.createInput2(collection)
        move_input.defineAsFreeMove(matrix)
        move_features.add(move_input)
        return True
    except:  #pylint:disable=bare-except
        pass
    try:
        move_features.add(move_features.createInput(collection, matrix))
        return True
    except:  #pylint:disable=bare-except
        log_line("Move feature failed: " + traceback.format_exc())
        return False


def build_cavity_tool_bodies(root, entities, result, stamp):
    """One workspace copy of the selected seed cavity per placement, moved onto its snapped point."""
    seed_body = entities["seed_body"]
    accepted = result.accepted()
    known_names = set()
    for body in root.bRepBodies:
        known_names.add(body.name)
    tool_names = []

    progress = ui.createProgressDialog()
    progress.isBackgroundTranslucent = False
    progress.show(DIALOG_TITLE, "Placing cavity tool body %v of %m", 0, len(accepted), 0)
    try:
        for index, placement in enumerate(accepted):
            if progress.wasCancelled:
                for remaining in accepted[index:]:
                    remaining.reject("cancelled_during_cut")
                break
            try:
                body = copy_seed_body(root, seed_body, known_names)
            except:  #pylint:disable=bare-except
                log_line("Body copy failed: " + traceback.format_exc())
                body = None
            if body is None:
                placement.reject("body_copy_failed")
                continue
            tool_name = TOOL_BODY_PREFIX + stamp + "_" + str(placement.cavity_id)
            try:
                body.name = tool_name
            except:  #pylint:disable=bare-except
                tool_name = body.name
            known_names.add(tool_name)
            if not move_body(root, body, placement.transform4):
                placement.reject("transform_failed")
                try:
                    body.deleteMe()
                except:  #pylint:disable=bare-except
                    pass
                continue
            tool_names.append(tool_name)
            progress.progressValue = index + 1
    finally:
        progress.hideDialog()
    adsk.doEvents()
    return tool_names


def refetch_body(root, name):
    """Look a workspace body up by name rather than trusting a stale handle."""
    body = root.bRepBodies.itemByName(name)
    if body is None or not body.isValid:
        return None
    return body


def execute_boolean_cuts(design, entities, result, stamp):
    """Copy the selected seed cavity once per placement, then exactly one combine cut."""
    started = time.time()
    root = entities["root"]

    if not entities["seed_body"].isValid:
        seed_body = refetch_body(root, entities["seed_name"])
        if seed_body is None:
            return False, "The selected seed cavity body is no longer valid."
        entities["seed_body"] = seed_body
    if not entities["target_body"].isValid:
        target_body = refetch_body(root, entities["target_name"])
        if target_body is None:
            return False, "The selected target body is no longer valid."
        entities["target_body"] = target_body

    tool_names = build_cavity_tool_bodies(root, entities, result, stamp)
    if not tool_names:
        return False, "No cavity tool body could be placed in the workspace."

    tools = adsk.core.ObjectCollection.create()
    lost = 0
    for name in tool_names:
        body = refetch_body(root, name)
        if body is None:
            lost = lost + 1
            continue
        tools.add(body)
    if tools.count == 0:
        return False, "Every placed cavity tool body was lost before the cut."

    target_body = refetch_body(root, entities["target_name"])
    if target_body is None:
        target_body = entities["target_body"]
    if target_body is None or not target_body.isValid:
        return False, "The target body was lost before the cut."

    try:
        combines = root.features.combineFeatures
        combine_input = combines.createInput(target_body, tools)
        combine_input.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
        combine_input.isKeepToolBodies = False
        combine_input.isNewComponent = False
        combines.add(combine_input)
        adsk.doEvents()
    except:  #pylint:disable=bare-except
        log_line("Combine cut failed: " + traceback.format_exc())
        return False, ("The combine cut failed. " + str(tools.count) + " cavity tool bodies named "
                       + TOOL_BODY_PREFIX + stamp + " have been left in the workspace.")

    result.timings["cut_seconds"] = time.time() - started
    message = "Cut " + str(tools.count) + " cavities into " + entities["target_name"] + "."
    if lost:
        message = message + " " + str(lost) + " tool bodies were lost before the cut."
    return True, message


def export_placement_manifest(result, params, entities, centre_count, folder, stamp):
    """Write the CSV manifest and record the full run summary in the run log."""
    manifest_path = None
    if folder is not None:
        manifest_path = os.path.join(folder, "cavity_manifest_" + stamp + ".csv")
        try:
            with open(manifest_path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(MANIFEST_COLUMNS)
                for placement in result.placements:
                    writer.writerow(placement.manifest_row())
        except:  #pylint:disable=bare-except
            log_line("Manifest could not be written: " + traceback.format_exc())
            manifest_path = None

    log_line("")
    log_line("Selected entities")
    log_line("  seed cavity body: " + entities["seed_name"])
    log_line("  target body: " + entities["target_name"])
    log_line("  allowed faces: " + str(len(entities["allowed_faces"])))
    log_line("  allowed face area: " + format_number(entities["allowed_area_mm2"]) + " mm2")
    log_line("  excluded faces: " + str(len(entities["excluded_faces"])))
    log_line("  duplicate allowed picks ignored: " + str(entities["allowed_duplicates"]))
    log_line("  duplicate excluded picks ignored: " + str(entities["excluded_duplicates"]))
    log_line("  allowed faces dropped for also being excluded: "
             + str(entities["allowed_excluded_overlap"]))
    log_line("")
    log_line("Counts")
    for entry in summary_lines(result, entities, centre_count):
        log_line("  " + entry)
    log_line("")
    log_line("Timing")
    for key in sorted(result.timings.keys()):
        log_line("  " + key + ": " + format_number(result.timings[key]))
    log_line("")
    log_line("Notes")
    log_line("  Manifest columns streamwise_x, streamwise_y and streamwise_z carry the local"
             " tangent reference axis of each surface frame. This script derives no flow"
             " direction and forces no cavity to face any global direction.")
    log_line("  Duplicate detection compares position, depth-axis angle and handedness, so a"
             " mirrored member co-located with its rotation partner is kept rather than"
             " silently discarded.")
    log_line("  Positions in the manifest are in mm; Fusion stores lengths internally in cm.")
    log_line("  Manifest file: " + str(manifest_path))
    return manifest_path


def report_incomplete_run(outcome, log_path):
    """Say why the run stopped and exactly where to read the log for it."""
    lines = ["Cavity pattern " + SYMMETRY_LABEL + " did not complete.",
             "",
             "Reason: " + str(outcome),
             ""]
    if log_path:
        lines.append("A log of this run was written to:")
        lines.append(log_path)
    else:
        lines.append("No log file could be written anywhere. Open the Text Commands window from"
                     " the Fusion View menu to read the full traceback.")
    try:
        ui.messageBox("\n".join(lines), DIALOG_TITLE)
    except:  #pylint:disable=bare-except
        pass


def execute_run(stamp, folder, log_path):
    """Do the work; returns 'completed' or a short reason the run stopped early."""
    design = adsk.fusion.Design.cast(app.activeProduct)
    if design is None:
        ui.messageBox("Open a Fusion design before running this script.", DIALOG_TITLE)
        return "no active Fusion design"

    try:
        log_line("document: " + str(app.activeDocument.name))
    except:  #pylint:disable=bare-except
        log_line("document: unavailable")
    if design.designType == adsk.fusion.DesignTypes.ParametricDesignType:
        log_line("design mode: parametric")
    else:
        log_line("design mode: direct modelling")
    log_line("Cut architecture: the seed cavity and the target body are both selected from the"
             " workspace; each placement is a copy of the selected seed moved into position, and"
             " all of them are consumed by a single combine cut.")

    log_line("")
    log_line("Stage: collecting settings")
    params = collect_user_inputs(design)
    if params is None:
        return "cancelled while collecting settings"
    log_line("Parameters")
    for entry in params.describe():
        log_line("  " + entry)

    log_line("")
    log_line("Stage: selecting bodies and faces")
    entities = validate_model_entities(design)
    if entities is None:
        return "body or face selection was cancelled or rejected"
    log_line("  seed cavity body: " + entities["seed_name"])
    log_line("  target body: " + entities["target_name"])
    log_line("  allowed faces: " + str(len(entities["allowed_faces"]))
             + " covering " + format_number(entities["allowed_area_mm2"]) + " mm2")
    log_line("  excluded faces: " + str(len(entities["excluded_faces"])))

    log_line("")
    log_line("Stage: measuring the seed cavity")
    try:
        anchor = seed_anchor_point(entities["seed_body"])
    except:  #pylint:disable=bare-except
        log_line("Anchor search failed: " + traceback.format_exc())
        ui.messageBox("Could not find the seed cavity anchor vertex.", DIALOG_TITLE)
        return "the seed cavity anchor vertex could not be found"
    seed_radius = seed_footprint_radius(entities["seed_body"], anchor)
    log_line("  seed anchor vertex in mm: "
             + format_number(anchor[0] * MM_PER_CM) + ", "
             + format_number(anchor[1] * MM_PER_CM) + ", "
             + format_number(anchor[2] * MM_PER_CM))
    log_line("  seed footprint radius about the anchor: "
             + format_number(seed_radius * MM_PER_CM) + " mm")
    entities["anchor"] = anchor
    entities["seed_radius"] = seed_radius

    orbit = build_local_symmetry_orbit(params)
    log_line("  local symmetry members per group: " + str(len(orbit)))
    log_line("  mirror members per group: " + str(sum(1 for s in orbit if s.is_mirror())))

    log_line("")
    log_line("Stage: reading the distribution sketch")
    centre_data = build_global_group_centres(design, params)
    if centre_data is None:
        return "the distribution sketch was missing, empty, or the group cap was declined"
    centres = centre_data["centres"]
    log_line("  " + centre_data["source"])
    log_line("  distribution sketch curves: " + str(centre_data["curve_count"])
             + " curve(s) joined into " + str(centre_data["chain_count"]) + " chain(s)")
    if centre_data["merged"]:
        log_line("  merged " + str(centre_data["merged"]) + " coincident centre(s) at shared"
                 " vertices or closed-loop wraps")
    log_line("  group separation is arc length along the curve before snapping; on-surface"
             " spacing differs slightly wherever the curve sits off the skin")

    log_line("")
    log_line("Stage: generating placements")
    result = generate_all_placements(entities, params, orbit, anchor, seed_radius, centres)
    result.centre_source = centre_data["source"]

    proceed, preview_sketch = preview_placements(entities, result, params, len(centres), stamp)

    cut_message = "No boolean cut was performed."
    if proceed:
        log_line("")
        log_line("Stage: cutting")
        if preview_sketch is not None and preview_sketch.isValid:
            try:
                preview_sketch.deleteMe()
            except:  #pylint:disable=bare-except
                log_line("  the preview sketch could not be deleted automatically")
        ok, cut_message = execute_boolean_cuts(design, entities, result, stamp)
        log_line("  boolean cut: " + ("succeeded" if ok else "failed") + ". " + cut_message)
    elif result.accepted():
        log_line("User declined the cut; the preview sketch and manifest were kept.")
    else:
        log_line("No placement was accepted, so no cut was attempted.")

    manifest_path = export_placement_manifest(result, params, entities, len(centres),
                                              folder, stamp)

    final = summary_lines(result, entities, len(centres))
    final.append("")
    final.append(cut_message)
    final.append("")
    final.append("Manifest and log folder:")
    final.append(str(folder))
    if manifest_path:
        final.append(os.path.basename(manifest_path))
    if log_path:
        final.append(os.path.basename(log_path))
    ui.messageBox("\n".join(final), DIALOG_TITLE)
    return "completed"


def run(context: str):
    """Script entry point."""
    stamp = time.strftime("%Y%m%d_%H%M%S")
    del RUN_LOG_LINES[:]
    folder = log_folder_path()
    log_path = None
    if folder is not None:
        log_path = os.path.join(folder, "cavity_run_" + stamp + ".log")
    outcome = "did not complete"
    try:
        log_line("Cavity pattern run log, local symmetry " + SYMMETRY_LABEL)
        log_line("timestamp: " + stamp)
        log_line("log file: " + str(log_path))
        outcome = execute_run(stamp, folder, log_path)
    except:  #pylint:disable=bare-except
        outcome = "failed with an exception"
        log_line("")
        log_line("TRACEBACK")
        log_line(traceback.format_exc())
        app.log(f"Failed:\n{traceback.format_exc()}")
    finally:
        log_line("")
        log_line("Run outcome: " + str(outcome))
        written = write_run_log(log_path)
        if outcome != "completed":
            report_incomplete_run(outcome, written)


def stop(context: str):
    """Script teardown."""
    pass
