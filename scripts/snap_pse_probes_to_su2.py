"""
snap_pse_probes_to_su2.py

Reads an SU2 native mesh, extracts the wall marker surface, classifies nodes as
windward or leeward using surface normals, snaps PSE stations to real surface
nodes, and emits a CUSTOM_OUTPUTS block for SU2 v8 (Harrier).

Usage:
    python snap_pse_probes_to_su2.py mesh.su2 [marker] [--offset 0.002]
                                              [--tol 0.5] [--out probes.cfg]
"""

import argparse
import sys

import numpy as np
from scipy.spatial import cKDTree

# ----------------------------------------------------------------------
# Target stations (from OpenFOAM system/probesPSE)
# ----------------------------------------------------------------------
X_STATIONS = [
    0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00,
    8.00, 10.00, 12.00, 14.00, 17.00, 20.00, 23.00, 26.00, 29.00, 32.00,
    35.00, 38.00, 41.00, 44.00, 47.00,
]

SPANWISE_STATIONS = [
    (20.0, -12.0), (20.0, -8.0), (20.0, -4.0), (20.0, 0.0),
    (20.0, 4.0), (20.0, 8.0), (20.0, 12.0),
    (35.0, -12.0), (35.0, -8.0), (35.0, -4.0), (35.0, 0.0),
    (35.0, 4.0), (35.0, 8.0), (35.0, 12.0),
]

FIELDS = ["PRESSURE", "TEMPERATURE", "DENSITY"]

# VTK element type -> node count, for surface elements
SURF_TYPE_NODES = {3: 2, 5: 3, 9: 4}


# ----------------------------------------------------------------------
# Mesh parsing
# ----------------------------------------------------------------------
def _value_after_equals(line):
    return line.split("=", 1)[1].strip()


def read_su2_surface(mesh_filename, wall_marker_name):
    """Return (all_vertices, faces) where faces is a list of node-index lists."""
    vertices = None
    faces = []
    ndime = 3
    marker_found = False

    with open(mesh_filename, "r") as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("NDIME"):
                ndime = int(_value_after_equals(stripped))
                if ndime != 3:
                    raise ValueError("Mesh is %dD; this script expects 3D" % ndime)
                continue

            if stripped.startswith("NPOIN"):
                # NPOIN can carry extra tokens, e.g. "NPOIN= 1000 800"
                n_points = int(_value_after_equals(stripped).split()[0])
                print("Global mesh points: %d. Loading coordinates." % n_points)
                pts = np.empty((n_points, 3), dtype=np.float64)
                for i in range(n_points):
                    parts = f.readline().split()
                    if len(parts) < 3:
                        raise ValueError(
                            "Truncated point line at index %d" % i)
                    # Trailing point index, if present, is ignored
                    pts[i] = (float(parts[0]), float(parts[1]), float(parts[2]))
                vertices = pts
                continue

            if stripped.startswith("MARKER_TAG"):
                tag = _value_after_equals(stripped)
                is_target = tag == wall_marker_name
                # MARKER_ELEMS must be the next non-blank line
                elems_line = f.readline()
                while elems_line.strip() == "":
                    elems_line = f.readline()
                n_elems = int(_value_after_equals(elems_line).split()[0])

                if not is_target:
                    for _ in range(n_elems):
                        f.readline()
                    continue

                marker_found = True
                print("Marker %s: %d surface elements." % (tag, n_elems))
                for i in range(n_elems):
                    parts = f.readline().split()
                    etype = int(parts[0])
                    if etype not in SURF_TYPE_NODES:
                        raise ValueError(
                            "Unsupported surface element type %d on marker %s"
                            % (etype, tag))
                    nn = SURF_TYPE_NODES[etype]
                    if len(parts) < nn + 1:
                        raise ValueError(
                            "Truncated element line %d on marker %s" % (i, tag))
                    # Slice by node count so a trailing element index is dropped
                    faces.append([int(p) for p in parts[1:nn + 1]])

    if vertices is None:
        raise ValueError("No NPOIN block found in %s" % mesh_filename)
    if not marker_found:
        raise ValueError(
            "Marker %s not found in %s" % (wall_marker_name, mesh_filename))
    if not faces:
        raise ValueError("Marker %s has no elements" % wall_marker_name)

    return vertices, faces


# ----------------------------------------------------------------------
# Surface classification
# ----------------------------------------------------------------------
def nodal_normals(vertices, faces):
    """Area-weighted nodal normals over the marker surface.

    Returns (node_ids, coords, normals) for nodes used by the marker only.
    """
    used = sorted({n for face in faces for n in face})
    remap = {g: i for i, g in enumerate(used)}
    coords = vertices[used]
    normals = np.zeros_like(coords)

    for face in faces:
        if len(face) < 3:
            continue
        local = [remap[n] for n in face]
        # Fan-triangulate; cross product magnitude gives 2x triangle area
        p0 = coords[local[0]]
        for k in range(1, len(local) - 1):
            p1 = coords[local[k]]
            p2 = coords[local[k + 1]]
            n = np.cross(p1 - p0, p2 - p0)
            for li in (local[0], local[k], local[k + 1]):
                normals[li] += n

    mag = np.linalg.norm(normals, axis=1)
    good = mag > 0
    normals[good] /= mag[good][:, None]
    return np.array(used), coords, normals


def split_surfaces(coords, normals):
    """Split into windward (normal points down) and leeward (points up)."""
    nz = normals[:, 2]
    wind = nz < -1e-6
    leew = nz > 1e-6
    ambiguous = int(np.count_nonzero(~(wind | leew)))
    if ambiguous:
        print("Note: %d nodes have near-horizontal normals "
              "(side or trailing edge); excluded from both surfaces."
              % ambiguous)
    return wind, leew


def report_geometry(coords):
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    print("Wall bounding box:")
    print("  x: %.4f to %.4f" % (lo[0], hi[0]))
    print("  y: %.4f to %.4f" % (lo[1], hi[1]))
    print("  z: %.4f to %.4f" % (lo[2], hi[2]))

    x_needed = max(X_STATIONS)
    if x_needed > hi[0] or min(X_STATIONS) < lo[0]:
        print("WARNING: requested x stations fall outside the wall x range. "
              "Check mesh units and origin.")

    half_model = lo[1] > -1e-3 or hi[1] < 1e-3
    if half_model:
        print("WARNING: wall spans only one side of y=0. "
              "Mirrored spanwise stations will be dropped.")
    return lo, hi, half_model


# ----------------------------------------------------------------------
# Probe construction
# ----------------------------------------------------------------------
def snap(tree, subset_coords, subset_normals, target_xy, offset):
    dist, idx = tree.query(target_xy)
    pt = subset_coords[idx]
    if offset != 0.0:
        # Move inward, opposite the outward normal, off the boundary face
        pt = pt - offset * subset_normals[idx]
    return pt, dist


def build_probe_definitions(coords, normals, offset, tol, half_model):
    wind_mask, leew_mask = split_surfaces(coords, normals)
    if not wind_mask.any() or not leew_mask.any():
        raise ValueError("Failed to identify both windward and leeward nodes")

    wind_c, wind_n = coords[wind_mask], normals[wind_mask]
    leew_c, leew_n = coords[leew_mask], normals[leew_mask]
    print("Windward nodes: %d. Leeward nodes: %d."
          % (len(wind_c), len(leew_c)))

    tree_wind = cKDTree(wind_c[:, :2])
    tree_leew = cKDTree(leew_c[:, :2])

    probes = []
    rejected = []

    def add(tag, tree, c, n, target):
        pt, dist = snap(tree, c, n, target, offset)
        if dist > tol:
            rejected.append((tag, target, dist))
            return
        probes.append((tag, pt, dist))

    for i, x in enumerate(X_STATIONS):
        add("WIND_%02d" % i, tree_wind, wind_c, wind_n, [x, 0.0])
    for i, x in enumerate(X_STATIONS):
        add("LEEW_%02d" % i, tree_leew, leew_c, leew_n, [x, 0.0])
    for i, (x, y) in enumerate(SPANWISE_STATIONS):
        if half_model and y < 0.0:
            continue
        add("SPAN_%02d" % i, tree_wind, wind_c, wind_n, [x, y])

    print("\nSnap distances in the x-y plane:")
    for tag, pt, dist in probes:
        print("  %-10s -> %10.4f %10.4f %10.4f   d = %.4f"
              % (tag, pt[0], pt[1], pt[2], dist))

    if rejected:
        print("\nDROPPED %d station(s) exceeding tolerance %.4f:"
              % (len(rejected), tol))
        for tag, target, dist in rejected:
            print("  %-10s target (%.4f, %.4f) nearest node d = %.4f"
                  % (tag, target[0], target[1], dist))

    return probes


# ----------------------------------------------------------------------
# Config generation
# ----------------------------------------------------------------------
def generate_config_snippet(probes):
    entries = []
    for tag, pt, _ in probes:
        for field in FIELDS:
            name = "%s_%s" % (field[0], tag)
            entries.append("  %s : Probe{%s}[%.6f, %.6f, %.6f]"
                           % (name, field, pt[0], pt[1], pt[2]))

    # Separator is ";" plus a line continuation; the final line takes neither
    body = ";\\\n".join(entries)

    lines = []
    lines.append("% ---------------------- CUSTOM PROBES ----------------------- %")
    lines.append("CUSTOM_OUTPUTS= '\\")
    lines.append(body + "'")
    lines.append("")
    lines.append("HISTORY_OUTPUT= ( ITER, RMS_RES, AERO_COEFF, HEAT, CUSTOM )")
    return "\n".join(lines), len(entries)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("marker", nargs="?", default="Solid_Walls")
    ap.add_argument("--offset", type=float, default=0.0,
                    help="inward offset along the nodal normal, in mesh units")
    ap.add_argument("--tol", type=float, default=0.5,
                    help="max allowed x-y snap distance before a station is dropped")
    ap.add_argument("--out", default=None, help="write the config block to a file")
    args = ap.parse_args()

    print("Reading SU2 mesh: %s" % args.mesh)
    vertices, faces = read_su2_surface(args.mesh, args.marker)
    _, coords, normals = nodal_normals(vertices, faces)
    print("Unique wall nodes: %d." % len(coords))

    _, _, half_model = report_geometry(coords)
    probes = build_probe_definitions(
        coords, normals, args.offset, args.tol, half_model)

    snippet, n_entries = generate_config_snippet(probes)
    print("\nGenerated %d custom outputs from %d probes."
          % (n_entries, len(probes)))

    if args.out:
        with open(args.out, "w") as f:
            f.write(snippet + "\n")
        print("Written to %s" % args.out)
    else:
        print()
        print(snippet)


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as exc:
        print("Error: file not found: %s" % exc.filename)
        sys.exit(1)
    except ValueError as exc:
        print("Error: %s" % exc)
        sys.exit(1)
