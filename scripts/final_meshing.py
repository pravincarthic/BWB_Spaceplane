import gmsh
import math
import os
import sys
import json

def load_config(config_file):
    """Load gmsh parameters and file paths from config file."""
    if not os.path.exists(config_file):
        print(f"[ERROR] Error: Config file not found: {config_file}")
        return None

    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        print(f"[OK] Config loaded from: {config_file}\n")
        return config
    except json.JSONDecodeError as e:
        print(f"[ERROR] Error parsing config file: {e}")
        return None
    except Exception as e:
        print(f"[ERROR] Error reading config file: {e}")
        return None

def exportStepFile(export_geometry_path, spaceplane_inspection_tags):
    if export_geometry_path:
        print(f"Exporting pre-mesh geometry for inspection...")
        try:
            gmsh.write(export_geometry_path)
            file_size = os.path.getsize(export_geometry_path)
            print(f"  [OK] Geometry exported: {export_geometry_path}")
            print(f"    Size: {file_size / 1024:.1f} KB")
            print(f"    Includes both the fluid domain AND a separate copy of the")
            print(f"    spaceplane solid (volume(s) {spaceplane_inspection_tags}), so the")
            print(f"    aircraft shape is visible directly, not just as a hidden cavity.")
            print(f"    Open this in SALOME (File > Import > STEP) to inspect the cuts")
            print(f"    and thickened features before meshing.\n")
        except Exception as e:
            print(f"  [WARNING] Geometry export failed: {e}")
            print(f"    Continuing with meshing anyway.\n")
    else:
        print("Export: No 'export_geometry_path' set in config - skipping pre-mesh "
                "geometry export.\n")

def setup_size_field(solid_walls_faces, size_params):
    """
    Build a spatially-varying mesh size field instead of relying on one
    global MeshSizeMin/MeshSizeMax pair.

    Logic:
      - Distance field measures distance from the solid-wall (aircraft) surfaces.
      - Threshold field maps that distance to element size:
          size = size_near   for distance <= dist_min
          size = size_far    for distance >= dist_max
          smooth linear ramp in between
      - This field is set as the background mesh, so it is the primary driver
        of element size everywhere. Mesh.MeshSizeMin/Max are kept only as
        global safety bounds, not as the thing doing the sizing.

    size_params keys (all in model units, mm):
      size_near   - element size right at the aircraft surface (fine)
      size_far    - element size far away, near the cone boundary (coarse)
      dist_min    - distance from aircraft below which size_near applies
      dist_max    - distance from aircraft beyond which size_far applies

    Returns the field id of the extra refinement field list, so cavity/detail
    fields can be added later without rebuilding the base field.
    """
    print("Setting up spatially-varying mesh size field...\n")

    size_near = size_params.get("size_near", 2)
    size_far = size_params.get("size_far", 500)
    dist_min = size_params.get("dist_min", 200)
    dist_max = size_params.get("dist_max", 20000)

    print(f"  size_near = {size_near} mm  (at the aircraft surface)")
    print(f"  size_far  = {size_far} mm  (far-field / near the cone)")
    print(f"  dist_min  = {dist_min} mm  (ramp starts)")
    print(f"  dist_max  = {dist_max} mm  (ramp ends)\n")

    # 1. Distance field from the solid wall surfaces
    dist_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist_field, "SurfacesList", solid_walls_faces)
    # Sampling density along the surfaces used to evaluate distance
    gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 100)

    # 2. Threshold field maps distance -> element size, with smooth transition
    thresh_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", size_near)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", size_far)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", dist_min)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", dist_max)
    gmsh.model.mesh.field.setNumber(thresh_field, "Sigmoid", 1)  # smooth transition, not linear kink

    print(f"  [OK] Distance field {dist_field} created from {len(solid_walls_faces)} solid wall surfaces")
    print(f"  [OK] Threshold field {thresh_field} created (fine near geometry, coarse far away)\n")

    return dist_field, thresh_field


def add_cavity_refinement_fields(cavity_surface_tags, base_field_ids, size_params):
    """
    Add extra fine-mesh fields targeting specific surfaces (e.g. mm-scale
    cavities in a future design). This does NOT touch the global size field
    or the far-field cone settings - it only pins a finer size near the
    listed surfaces.

    cavity_surface_tags: list of surface tags that need extra refinement
    base_field_ids: (dist_field, thresh_field) from setup_size_field
    size_params: same dict as setup_size_field, reads 'cavity_size' and
                 'cavity_dist_max' if present

    Returns the combined field id (Min of all fields) to be set as background mesh.
    """
    dist_field, thresh_field = base_field_ids
    field_ids_to_combine = [thresh_field]

    if cavity_surface_tags:
        cavity_size = size_params.get("cavity_size", size_params.get("size_near", 2) / 2)
        cavity_dist_max = size_params.get("cavity_dist_max", 50)

        print(f"Adding cavity refinement: {len(cavity_surface_tags)} surfaces, "
              f"size={cavity_size}mm within {cavity_dist_max}mm\n")

        cav_dist_field = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(cav_dist_field, "SurfacesList", cavity_surface_tags)
        gmsh.model.mesh.field.setNumber(cav_dist_field, "Sampling", 100)

        cav_thresh_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "InField", cav_dist_field)
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "SizeMin", cavity_size)
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "SizeMax", size_params.get("size_far", 500))
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "DistMax", cavity_dist_max)
        gmsh.model.mesh.field.setNumber(cav_thresh_field, "Sigmoid", 1)

        field_ids_to_combine.append(cav_thresh_field)

    if len(field_ids_to_combine) > 1:
        # Take the minimum size at every point across all fields, so the
        # finest applicable requirement always wins.
        min_field = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", field_ids_to_combine)
        return min_field
    else:
        return field_ids_to_combine[0]


def create_mesh(config_file='gmsh_config.json'):
    print(gmsh.__version__)
    print("\n" + "=" * 80)
    print("BWB SPACEPLANE MESH GENERATION - TEMPLATE BASED")
    print("=" * 80 + "\n")

    # Load configuration
    config = load_config(config_file)
    if config is None:
        return False

    # Extract file paths
    step_file = config.get('step_file')
    output_path = config.get('output_file')
    export_geometry_path = config.get('export_geometry_path')
    
    print(f"File {step_file} opened.\nOutput: {output_path}\n")

    # Extract gmsh parameters
    gmsh_params = config.get('gmsh_parameters', {})
    size_params = config.get('size_field_parameters', {})
    # Geometry-scale parameters (rotation, far-field box sizing, thin-feature
    # thresholds) - kept out of code so a new design/orientation/domain size
    # never requires editing this script, only the config file.
    geom_params = config.get('geometry_parameters', {})

    inlet_faces = []
    outlet_faces = []
    atmosphere_faces = []
    solid_walls_faces = []
    # ========================================================================
    # 1. INITIALIZE
    # ========================================================================
    print("Step 1: Initializing Gmsh...")
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.RandomFactor", gmsh_params.get("Mesh.RandomFactor", 1e-5))
    # gmsh.option.setNumber("Geometry.AutoCoherence", False)
    # gmsh.option.setNumber("Geometry.OCCMakeSolids", 0)
    # gmsh.option.setNumber("Geometry.Tolerance", 1000)
    # gmsh.option.setNumber("Geometry.OCCUnionUnify", 0)
    # gmsh.option.setNumber("Geometry.OCCFixDegenerated", 0)
    # gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 0)
    # gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 0)
    #gmsh.option.setNumber("Geometry.OCCScalingFactor", 1000.0)
    gmsh.option.setNumber("General.NumThreads", gmsh_params.get("General.NumThreads", 16))
    gmsh.model.add("BWB_Spaceplane_Mesh")
    print("[OK] Gmsh initialized\n")

    # ========================================================================
    # 2. IMPORT GEOMETRY
    # ========================================================================
    print("Step 2: Importing STEP geometry...")

    if not os.path.exists(step_file):
        print(f"[ERROR] Error: File not found: {step_file}")
        gmsh.finalize()
        return False

    try:
        spaceplane_shapes = gmsh.model.occ.importShapes(step_file)
        gmsh.model.occ.synchronize()
        print("[OK] STEP file imported\n")
    except Exception as e:
        print(f"[ERROR] Failed to import STEP: {e}")
        gmsh.finalize()
        return False

    # ========================================================================
    # 2a. CONVERT UNITS: STEP FILE IS IN METRES, REST OF THIS SCRIPT ASSUMES MM
    # ========================================================================
    # All downstream thresholds, sizes, and config values (thin-feature
    # guards, size-field distances, boundary-layer thickness, far-field box
    # dimensions, etc.) are written and interpreted in millimetres. The
    # incoming STEP geometry is in metres, so without this conversion every
    # one of those mm-scale constants would be applied to a model that is
    # 1000x smaller than intended (e.g. a 1mm thin-feature threshold would
    # end up 1000x larger than the entire spaceplane). Scaling right after
    # import, before rotation or any boolean operation, keeps every later
    # step - including exportStepFile() and all interim STEP exports -
    # unchanged and working purely in mm from this point on.
    # print("Step 2a: Converting imported geometry from metres to millimetres...")
    try:
        #gmsh.model.occ.dilate(spaceplane_shapes, 0, 0, 0, 1000.0, 1000.0, 1000.0)
        #gmsh.model.occ.synchronize()
        print("[OK] Geometry scaled by 1000x (m -> mm)\n")
    except Exception as e:
        print(f"[ERROR] Failed to scale geometry to mm: {e}")
        gmsh.finalize()
        return False

    # ========================================================================
    # 3. EXTRACT SPACEPLANE VOLUME
    # ========================================================================
    print("Step 3: Extracting spaceplane volume...")
    spaceplane_vols = [item for item in spaceplane_shapes if item[0] == 3]

    if not spaceplane_vols:
        print("[ERROR] Error: Could not find 3D volume in STEP file")
        gmsh.finalize()
        return False

    spaceplane_tag = spaceplane_vols[0][1]
    print(f"[OK] Spaceplane volume tag: {spaceplane_tag}\n")

    # ========================================================================
    # 4. ROTATE SPACEPLANE
    # ========================================================================
    rotation_angle_deg = geom_params.get("rotation_angle_deg", -5.0)
    print(f"Step 4: Rotating spaceplane ({rotation_angle_deg} deg around Y-axis)...")
    angle_rad = rotation_angle_deg * math.pi / 180.0
    gmsh.model.occ.rotate([(3, spaceplane_tag)], 0, 0, 0, 0, 1, 0, angle_rad)
    gmsh.model.occ.synchronize()

    exportStepFile(export_geometry_path + ".1.step",'')

    print("[OK] Spaceplane rotated\n")

    boundary_entities = gmsh.model.getBoundary([(3, spaceplane_tag)], combined=False, oriented=False)
    spaceplane_faces = [tag for dim, tag in boundary_entities if dim == 2]

    print("\n" + "=" * 80)
    print(f"GEOMETRY ANALYSIS: {len(spaceplane_faces)} SURFACES DETECTED")
    print("=" * 80)
    print(f"Raw Face Tags: {sorted(spaceplane_faces)}\n")

    print(f"{'Index':<6} | {'Tag':<5} | {'Center of Mass (X, Y, Z)':<32} | {'Dimensions (delta-X x delta-Y x delta-Z)':<28}")
    print("-" * 80)

    for idx, tag in enumerate(sorted(spaceplane_faces), start=1):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        com = gmsh.model.occ.getCenterOfMass(2, tag)
        dx = xmax - xmin
        dy = ymax - ymin
        dz = zmax - zmin
        com_str = f"({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})"
        dim_str = f"{dx:.2f} x {dy:.2f} x {dz:.2f}"
        print(f"#{idx:<5} | {tag:<5} | {com_str:<32} | {dim_str:<28}")
        print(f"       - Span Limits: X=[{xmin:.1f} -> {xmax:.1f}], Y=[{ymin:.1f} -> {ymax:.1f}], Z=[{zmin:.1f} -> {zmax:.1f}]")
        print("-" * 80)

    # ========================================================================
    # 4a. ENFORCE MINIMUM THICKNESS ON THIN FEATURES - BEFORE THE BOOLEAN CUT
    # ========================================================================
    # Gated behind enable_thin_feature_enforcement (config, default True) so
    # a new/changed design can be test-meshed on its RAW imported geometry -
    # no extrusion, no fragment(), no healShapes() - before deciding whether
    # this machinery is even needed for it. All of Step 4a/4b's thickening
    # and sealing was built around specific thin-feature failures on a prior
    # design; a new design's thin features (if any) may be different or
    # absent, so trying to pad/seal it blind, before ever seeing how it
    # meshes untouched, risks fighting problems that don't exist for this
    # geometry and obscuring ones that do. Disabled here just sets
    # spaceplane_tags = [spaceplane_tag] and skips straight to Step 5.
    enable_thin_feature_enforcement = False #geom_params.get("enable_thin_feature_enforcement", False)
    print("Step 4a: SKIPPED (enable_thin_feature_enforcement=false in config) - "
          "using the spaceplane's raw imported geometry, unmodified.\n")
    spaceplane_tags = [spaceplane_tag]

    # if not enable_thin_feature_enforcement:
    #     #print("Step 4a: SKIPPED (enable_thin_feature_enforcement=false in config) - "
    #     #      "using the spaceplane's raw imported geometry, unmodified.\n")
    #     #spaceplane_tags = [spaceplane_tag]
    # else:
    #     # Moved here from post-cut Step 6a. enforce_minimum_thickness() extrudes
    #     # a thin surface by `deficit` along its own shortest local axis and
    #     # fragments the resulting sliver into the current solid - this is only
    #     # geometrically meaningful while the spaceplane is still a SOLID: the
    #     # extrusion adds material to a thin slab, thickening it.
    #     #
    #     # Run after the boolean cut (as it was previously), those same surfaces
    #     # are no longer solid material - they are the CAVITY WALL of the fluid
    #     # domain where the spaceplane used to be. Extruding a cavity-wall
    #     # surface by ~1mm along an arbitrarily-chosen local axis (picked purely
    #     # by whichever of dx/dy/dz is locally smallest, with no notion of
    #     # "into the solid" vs "into the fluid") does not thicken the aircraft
    #     # skin anymore; depending on direction, it either pokes a harmless
    #     # sliver into the fluid volume or - worse - carves into/disconnects the
    #     # thin cavity wall right at the spaceplane's thinnest features, which is
    #     # what was corrupting the Solid_Walls boundary group.
    #     #
    #     # Running this on the solid, pre-cut, thickens the actual aircraft skin
    #     # while it is still solid material. The boolean cut then subtracts an
    #     # already-thickened solid, producing a clean cavity wall with no
    #     # post-cut geometry surgery needed.
    #     print("Step 4a: Enforcing minimum thickness on thin features (pre-cut, on solid)...")
    #     min_thickness = size_params.get("min_feature_thickness", 1.0)
    #     thin_feature_skip_above = geom_params.get("thin_feature_skip_above_mm", 50000)
    #     planar_face_guard = geom_params.get("planar_face_guard_mm", 2000)
    #     # Below this native-gap threshold, a thin feature is treated as an
    #     # INTENTIONAL zero-thickness design edge (e.g. a trailing edge meant to
    #     # come to a physical point) rather than a thin wall that merely needs
    #     # padding to min_thickness - see enforce_minimum_thickness docstring.
    #     # Default (0.05mm) sits comfortably above typical STEP/BREP export
    #     # noise (~0.001-0.01mm) but far below any real structural wall.
    #     sharp_edge_gap_threshold = size_params.get("sharp_edge_gap_threshold_mm", 0.05)
    #     #spaceplane_tags, failed_thin_features, sharp_edges = enforce_minimum_thickness(
    #     #    spaceplane_tag, min_thickness=min_thickness, skip_above=thin_feature_skip_above,
    #     #    planar_face_guard=planar_face_guard, sharp_edge_gap_threshold=sharp_edge_gap_threshold
    #     #)
    #     #print(f"[OK] Spaceplane solid now spans volume(s): {spaceplane_tags}\n")

    #     # Seal BOTH intentional sharp edges (never padded - see above) and any
    #     # feature that fragment() tried and failed to thicken (see
    #     # seal_sharp_and_failed_features docstring). Root cause of the run this
    #     # replaces: Surface 11, a 3.03 x 34.67 x 0.0088mm trailing edge that IS
    #     # supposed to be zero-thickness by design, was left at its native gap -
    #     # far smaller than the ~2mm wall mesh size - which is what made HXT
    #     # fail with "Found two exactly self-intersecting facets". This is a
    #     # no-op when both lists are empty.
    #     # heal_tolerance=None (config default/absent) means auto-derive a
    #     # tolerance per feature from its own native gap - see
    #     # seal_sharp_and_failed_features docstring for why a single coarse
    #     # tolerance across the whole assembly previously corrupted the model.
    #     # Only set heal_tolerance_mm in config to override that with one fixed
    #     # value for every feature sealed here.
    #     heal_tolerance = geom_params.get("heal_tolerance_mm", None)
    #     #spaceplane_tags = seal_sharp_and_failed_features(
    #     #    spaceplane_tags, sharp_edges + failed_thin_features, min_thickness,
    #     #    heal_tolerance=heal_tolerance
    #     #)

    # ========================================================================
    # 5. CREATE FAR-FIELD DOMAIN (cuboid, NOT a cone)
    # ========================================================================
    # Switched from a conical far-field domain to a cuboid after five
    # independently-verified geometric/topological/parametric properties
    # of the cone construction were all confirmed correct, yet the mesher
    # failed identically every time:
    #   1. Periodicity (duplicate curve tags, u=2*pi) - checked twice, both
    #      clean.
    #   2. Degenerate/near-zero-length edges - chord distances confirmed
    #      exact to 10 decimal places against design intent.
    #   3. Mesh size vs. gap size (MeshSizeMax/size_far vs gap chord) -
    #      widening the gap 4x had zero effect on the failure.
    #   4. Mesh.MeshSizeFromCurvature interaction - setting to 0 had zero
    #      effect, byte-identical failure numbers.
    #   5. Parametric (u,v) trim consistency between boundary curves -
    #      directly queried via gmsh, found perfectly consistent (used
    #      u-span exactly equals the surface's own domain span).
    #
    # With every specific hypothesis eliminated by direct evidence, the
    # most likely remaining explanation is a gmsh Frontal-Delaunay
    # algorithm limitation with this class of surface (bounded by two
    # full-sweep circular arcs plus two constant-u line edges), independent
    # of the exact numbers involved. A cuboid domain has ONLY flat planar
    # faces and straight-line edges - none of the curved/circular/conical
    # geometry that any of the five failure modes above could apply to, so
    # this removes the entire failure class structurally rather than
    # continuing to patch around it.
    #
    # This is a physics trade-off, not a shortcut: the case is hypersonic,
    # where the conical domain was a genuine, deliberate approximation of
    # the Mach cone / zone of dependence - not an arbitrary shape choice.
    # A box can't taper the way a frustum does, so it is sized here to
    # fully ENVELOPE the original cone at every cross-section (using the
    # cone's largest radius, at the outlet, uniformly along the entire
    # length), guaranteeing it can't clip any flow feature the conical
    # domain was sized to capture. The trade is more far-field cells for
    # a structurally simpler, reliable domain - not a compromise in what
    # the domain actually contains.
    print("Step 5: Creating cuboid far-field domain...")

    # All far-field box dimensions come from config (geometry_parameters) -
    # sized to fully envelope the original cone (see design rationale
    # above); resizing the domain never requires editing this script.
    x0 = geom_params.get("far_field_x0_mm", -50000.0)
    length = geom_params.get("far_field_length_mm", 200000.0)
    r_envelope = geom_params.get("far_field_r_envelope_mm", 100000.0)  # cone's largest radius (at the outlet)

    box_tag = gmsh.model.occ.addBox(
        x0, -r_envelope, -r_envelope,
        length, 2 * r_envelope, 2 * r_envelope
    )
    gmsh.model.occ.synchronize()
    exportStepFile(export_geometry_path + ".2.step",'')

    # Same tag-collision guard as before, retained as a cheap sanity check
    # even though addBox is far less likely to produce this than the prior
    # revolve+heal construction did. Uses the overall bounding box of ALL
    # spaceplane volumes (spaceplane_tags may be >1 after Step 4a/4b
    # thickening and sealing), since a same-bbox collision would show up in
    # the combined extent either way.
    xmin_c, ymin_c, zmin_c, xmax_c, ymax_c, zmax_c = gmsh.model.getBoundingBox(3,     box_tag)
    spaceplane_bboxes = [gmsh.model.getBoundingBox(3, t) for t in spaceplane_tags]
    xmin_s = min(b[0] for b in spaceplane_bboxes)
    ymin_s = min(b[1] for b in spaceplane_bboxes)
    zmin_s = min(b[2] for b in spaceplane_bboxes)
    xmax_s = max(b[3] for b in spaceplane_bboxes)
    ymax_s = max(b[4] for b in spaceplane_bboxes)
    zmax_s = max(b[5] for b in spaceplane_bboxes)
    same_bbox = (
        abs(xmin_c - xmin_s) < 1e-3 and abs(xmax_c - xmax_s) < 1e-3 and
        abs(ymin_c - ymin_s) < 1e-3 and abs(ymax_c - ymax_s) < 1e-3 and
        abs(zmin_c - zmin_s) < 1e-3 and abs(zmax_c - zmax_s) < 1e-3
    )
    if box_tag in spaceplane_tags or same_bbox:
        print(f"[ERROR] the far-field volume (tag {box_tag}) is the same as, or has the")
        print(f"  same bounding box as, the spaceplane volume(s) {spaceplane_tags}.")
        print(f"  Cuboid bbox:     X=[{xmin_c:.1f},{xmax_c:.1f}] Y=[{ymin_c:.1f},{ymax_c:.1f}] Z=[{zmin_c:.1f},{zmax_c:.1f}]")
        print(f"  Spaceplane bbox: X=[{xmin_s:.1f},{xmax_s:.1f}] Y=[{ymin_s:.1f},{ymax_s:.1f}] Z=[{zmin_s:.1f},{zmax_s:.1f}]")
        gmsh.finalize()
        return False

    print(f"[OK] Far-field volume confirmed, tag: {box_tag}")
    print(f"  Bounding box: X=[{xmin_c:.1f},{xmax_c:.1f}] Y=[{ymin_c:.1f},{ymax_c:.1f}] Z=[{zmin_c:.1f},{zmax_c:.1f}]")
    print(f"  Fully envelopes the original cone (outlet radius {r_envelope:.0f}mm used uniformly)")
    print(f"  Distinct from spaceplane (tag(s) {spaceplane_tags}). No curved surfaces - flat")
    print(f"  planar faces only, structurally immune to every failure mode found on the cuboid construction.\n")

    # ========================================================================
    # 5a. SNAPSHOT COPY OF THE SPACEPLANE FOR VISUAL INSPECTION
    # ========================================================================
    # The boolean cut in Step 6 below runs with removeTool=True, which
    # DELETES the spaceplane solid(s) - only the cavity they leave behind
    # in the fluid domain survives. That cavity is a valid, correct
    # boolean result, but it is an internal VOID inside an otherwise closed
    # box: opening the Step 6b pre-mesh export in a STEP viewer shows what
    # looks like an empty box, because the aircraft shape only exists as a
    # hidden internal shell, not a solid you can see without a section cut.
    # This is exactly the "interim STEP doesn't contain the spaceplane"
    # symptom.
    #
    # Fix: take a throwaway OCC copy of the (already thickened/sealed)
    # spaceplane BEFORE the cut consumes the originals, purely so Step 6b
    # can write it out alongside the fluid domain as a real, visible solid
    # for inspection. It is deleted again right after export (Step 6b,
    # still before Step 7 builds any physical groups - see Step 6b's own
    # comment for why that ordering matters) so it can never reach Step 11
    # mesh generation as a stray, un-meshed volume.
    spaceplane_inspection_copy = gmsh.model.occ.copy([(3, t) for t in spaceplane_tags])
    gmsh.model.occ.synchronize()
    spaceplane_inspection_tags = [tag for dim, tag in spaceplane_inspection_copy if dim == 3]

    # ========================================================================
    # ========================================================================
    # 6. BOOLEAN CUT
    # ========================================================================
    # Cuts the (already-thickened, from Step 4a/4b) spaceplane solid(s) out
    # of the cuboid. spaceplane_tags is a list because thickening/sealing
    # may have fragmented the spaceplane into more than one conformally-glued
    # volume - all of them are passed as cut tools so every piece is
    # subtracted, leaving one clean cavity wall with no post-cut geometry
    # surgery required.
    print("Step 6: Cutting spaceplane from cuboid (fluid domain)...")
    try:
        fluid_domain, _ = gmsh.model.occ.cut(
            [(3, box_tag)],
            [(3, t) for t in spaceplane_tags],
            removeObject=True, removeTool=True
        )
        gmsh.model.occ.synchronize()
        fluid_volume_tags = [tag for dim, tag in fluid_domain if dim == 3]
        print(f"[OK] Boolean cut complete, fluid volume(s): {fluid_volume_tags}\n")
        exportStepFile(export_geometry_path + ".3.step",'')

        # ====================================================================
        # 6a. REMOVE DUPLICATE/COINCIDENT FACES BEFORE MESHING
        # ====================================================================
        # The per-surface thickening in Step 4a extrudes and fragments each
        # thin feature independently. When two thickened features are
        # geometrically adjacent (e.g. two thin panels sharing an edge),
        # their newly-created side faces can end up exactly coincident -
        # two distinct 2D entities occupying the same location. Gmsh's
        # mesher does not merge these on its own; left in place, HXT's 3D
        # meshing step fails with "Found two exactly self-intersecting
        # facets" once it hits the duplicate pair (confirmed on real
        # geometry: surfaces with identical vertex triples, traced back to
        # newly-created faces from the Step 4a thickening of two adjacent
        # thin panels).
        #
        # gmsh.model.occ.removeAllDuplicates() merges/removes exactly this
        # class of overlapping entity. It is run here - AFTER the boolean
        # cut has fully assembled the fluid domain, but BEFORE Step 7
        # collects surface tags for physical groups - because (a) it needs
        # the complete, final topology to find duplicates correctly, and
        # (b) like fragment(), it can renumber entity tags, so anything
        # that reads tags from it must do so AFTER this call, never before.
        # Running it before physical groups exist avoids stale-tag issues
        # entirely, the same lesson learned from the Step 4a fragment bug.
        print("Step 6a: Removing duplicate/coincident faces before meshing...")
        surfaces_before = gmsh.model.getEntities(2)
        volumes_before = gmsh.model.getEntities(3)
        try:
            #gmsh.model.occ.removeAllDuplicates()
            gmsh.model.occ.synchronize()
            surfaces_after = gmsh.model.getEntities(2)
            volumes_after = gmsh.model.getEntities(3)
            n_surf_removed = len(surfaces_before) - len(surfaces_after)
            n_vol_changed = len(volumes_before) - len(volumes_after)
            if n_surf_removed > 0:
                print(f"  [OK] Removed {n_surf_removed} duplicate/coincident surface(s) "
                      f"({len(surfaces_before)} -> {len(surfaces_after)})")
            else:
                print(f"  [OK] No duplicate surfaces found ({len(surfaces_after)} surfaces "
                      f"unchanged)")
            if n_vol_changed != 0:
                print(f"  [WARNING] Volume count changed during duplicate removal "
                      f"({len(volumes_before)} -> {len(volumes_after)}) - re-resolving "
                      f"fluid_volume_tags from current model state")
            # removeAllDuplicates() can renumber volume tags even when the
            # count is unchanged, so always re-resolve from the live model
            # rather than trusting the pre-call fluid_volume_tags list.
            fluid_volume_tags = [tag for dim, tag in gmsh.model.getEntities(3)]
            print(f"  Fluid domain now spans volume(s): {fluid_volume_tags}\n")
        except Exception as e:
            # Some gmsh versions raise if there is nothing to remove -
            # this is not a failure of the geometry, just means it was
            # already clean. Re-resolve tags defensively either way since
            # we cannot be certain what state the call left things in.
            print(f"  [INFO] removeAllDuplicates() reported: {e}")
            print(f"  (this can happen when there is nothing to remove on some gmsh "
                  f"versions - not necessarily an error)")
            fluid_volume_tags = [tag for dim, tag in gmsh.model.getEntities(3)]
            print(f"  Fluid domain spans volume(s): {fluid_volume_tags}\n")

        # ====================================================================
        # 6b. EXPORT PRE-MESH GEOMETRY FOR INSPECTION, THEN REMOVE THE
        #     INSPECTION COPY - BOTH BEFORE PHYSICAL GROUPS ARE CREATED
        # ====================================================================
        # Exports the current OCC geometry - after STEP import, rotation,
        # far-field box creation, boolean cut, and thickness enforcement,
        # but BEFORE any mesh sizing or mesh generation - as a STEP file.
        # This lets the cuts and thickened features be visually inspected in
        # SALOME (or any STEP viewer) before committing to a mesh.
        # Controlled by config so it does not run unless explicitly
        # requested.
        #
        # The model still contains the spaceplane_inspection_copy made in
        # Step 5a, so this export includes it as a real, visible solid
        # alongside the fluid domain - not just as the domain's hidden
        # internal cavity shell (the "STEP doesn't contain the spaceplane"
        # symptom that copy fixes).
        #
        # CRITICAL ORDERING - this block, including the inspection copy's
        # removal, MUST run entirely BEFORE Step 7 creates physical groups.
        # A real run hit exactly the failure this avoids: the export+remove
        # here used to run AFTER Step 7 (as the old "Step 8a"). Removing
        # the inspection copy volumes is a boolean/OCC operation just like
        # fragment() or removeAllDuplicates() - it renumbers surviving
        # entities on synchronize(). With physical groups already built
        # from the pre-removal tags, that renumbering silently invalidated
        # them: Step 10 then failed to find surfaces it had just listed
        # ("Unknown model face with tag 83", etc.) and 3D meshing produced
        # a real, large mesh (31,027,765 elements) that gmsh.write()
        # exported as a ~10KB file - because Mesh.SaveAll defaults to 0,
        # so gmsh.write() only exports entities that belong to a STILL-VALID
        # physical group, and the volume physical group no longer matched
        # anything real. Doing the copy's removal here, before Step 7 reads
        # or assigns any tag, means every physical group Step 7 creates is
        # built from the model's FINAL, stable numbering - no operation
        # that can renumber entities runs after this point until meshing.
        
        exportStepFile(export_geometry_path, spaceplane_inspection_tags)

        # Filter to tags that actually still exist before calling remove() -
        # gmsh.model.occ.copy() on multiple input solids can itself produce
        # fewer distinct output volumes than requested when inputs share
        # coincident geometry (confirmed on real geometry: 2 of 9 requested
        # copy tags never existed). Filtering avoids the "Unknown
        # OpenCASCADE entity" warning entirely instead of relying on
        # remove()'s own tolerance for a partially-invalid list.
        existing_3d_tags = {tag for dim, tag in gmsh.model.getEntities(3)}
        inspection_tags_to_remove = [t for t in spaceplane_inspection_tags if t in existing_3d_tags]
        if inspection_tags_to_remove:
            gmsh.model.occ.remove([(3, t) for t in inspection_tags_to_remove], recursive=True)
            gmsh.model.occ.synchronize()
            print(f"  [OK] Removed spaceplane inspection copy (volume(s) "
                  f"{inspection_tags_to_remove}) - not part of the fluid mesh")

        # The removal above can renumber the REAL fluid-domain volumes too
        # (same lesson as removeAllDuplicates() above) - re-resolve from the
        # live model rather than trusting the pre-removal fluid_volume_tags,
        # so every downstream step (starting with Step 7 next) reads only
        # tags that are guaranteed to still be valid.
        fluid_volume_tags = [tag for dim, tag in gmsh.model.getEntities(3)]
        print(f"  Fluid domain now spans volume(s): {fluid_volume_tags}\n")

        # ====================================================================
        # Robust Surface Detection Logic
        # ====================================================================
        print("Step 7: Creating physical groups from surfaces and volumes...")
        all_fluid_faces = []
        for vtag in fluid_volume_tags:
            boundary_faces = gmsh.model.getBoundary([(3, vtag)], combined=False, oriented=False)
            all_fluid_faces.extend(tag for dim, tag in boundary_faces)
        # A face can appear as boundary of more than one fragment-glued
        # volume only if it's an internal interface between them (e.g. the
        # thickened slab's shared face with the original fluid region) -
        # de-duplicate so it isn't double-counted or double-assigned.
        all_fluid_faces = list(set(all_fluid_faces))

        # Far-field box bounds (must match Step 5's box construction
        # exactly - x0, length, r_envelope are defined there).
        x_inlet = x0
        x_outlet = x0 + length
        y_min_bound, y_max_bound = -r_envelope, r_envelope
        z_min_bound, z_max_bound = -r_envelope, r_envelope
        # mm; generous vs. the old 0.01-0.2mm centroid tolerance, since this
        # now compares each face's own bounding-box extent against a box
        # plane position, not requiring the face's centroid to happen to
        # sit at (y=0,z=0). Sourced from config, not hardcoded.
        tol = geom_params.get("far_field_face_tolerance_mm", 0.5)

        for face in all_fluid_faces:
            fxmin, fymin, fzmin, fxmax, fymax, fzmax = gmsh.model.getBoundingBox(2, face)

            # A far-field box face is planar AND flush against one of the
            # six box bounds across its ENTIRE extent (both min and max on
            # that axis sit at the bound). This is what actually
            # identifies "this IS a box wall" - the previous version only
            # checked whether a face's centroid happened to be near
            # (y=0, z=0), which is true for essentially none of the six
            # box faces (their centroids sit at y=+/-100000 or
            # z=+/-100000), so every box face that failed that check fell
            # through to Solid_Walls by elimination - corrupting the wall
            # group (and therefore the boundary-layer field built from it)
            # with far-field geometry instead of just the spaceplane.
            is_inlet = abs(fxmin - x_inlet) < tol and abs(fxmax - x_inlet) < tol
            is_outlet = abs(fxmin - x_outlet) < tol and abs(fxmax - x_outlet) < tol
            is_ymin = abs(fymin - y_min_bound) < tol and abs(fymax - y_min_bound) < tol
            is_ymax = abs(fymin - y_max_bound) < tol and abs(fymax - y_max_bound) < tol
            is_zmin = abs(fzmin - z_min_bound) < tol and abs(fzmax - z_min_bound) < tol
            is_zmax = abs(fzmin - z_max_bound) < tol and abs(fzmax - z_max_bound) < tol

            if is_inlet:
                inlet_faces.append(face)
            elif is_outlet:
                outlet_faces.append(face)
            elif is_ymin or is_ymax or is_zmin or is_zmax:
                atmosphere_faces.append(face)
            # else: not a far-field box face -> falls through to
            # solid_walls_faces below. This is now correct because
            # far-field faces are positively identified above, rather
            # than Solid_Walls being defined as "whatever wasn't near
            # (0,0)".

        if inlet_faces:
            gmsh.model.addPhysicalGroup(2, inlet_faces, name="Inlet")
        if outlet_faces:
            gmsh.model.addPhysicalGroup(2, outlet_faces, name="Outlet")
        if atmosphere_faces:
            gmsh.model.addPhysicalGroup(2, atmosphere_faces, name="Atmosphere")
        far_field_faces = set(inlet_faces + outlet_faces + atmosphere_faces)
        solid_walls_faces = list(set(all_fluid_faces) - far_field_faces)
        gmsh.model.addPhysicalGroup(2, solid_walls_faces, name="Solid_Walls")
        gmsh.model.addPhysicalGroup(3, fluid_volume_tags, 7, "Group_Of_All_Volumes")
        print("[OK] Physical groups created\n")

    except Exception as e:
        print(f"[ERROR] Boolean operation failed: {e}")
        gmsh.finalize()
        return False

    # ========================================================================
    # 8. BOUNDARY FACE ASSIGNMENT (unchanged verification block)
    # ========================================================================
    print("Step 8: Verifying boundary face assignments...")
    print(f"Solid walls: {solid_walls_faces}")
    print(f"Inlet: {inlet_faces}")
    print(f"Outlet: {outlet_faces}")
    print(f"Atmosphere: {atmosphere_faces}\n")

    surfaces = gmsh.model.getEntities(2)
    model_surface_tags = [tag for dim, tag in surfaces]

    print(f"Total surfaces in model: {len(model_surface_tags)}")
    print(f"All surface tags: {sorted(model_surface_tags)}\n")

    all_assigned = inlet_faces + outlet_faces + atmosphere_faces + solid_walls_faces
    unassigned = [t for t in model_surface_tags if t not in all_assigned]

    if unassigned:
        print("[WARNING]  UNASSIGNED SURFACES (adding to Solid_Walls automatically):")
        print("-" * 50)
        for tag in sorted(unassigned):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
            size = max(xmax - xmin, ymax - ymin, zmax - zmin)
            print(f"    Tag {tag:2d}: Size={size:7.2f}")
        # Any surface not caught by inlet/outlet/atmosphere logic (e.g. from
        # sliver removal side effects) is real geometry -> treat as solid wall
        # so it is neither lost nor left out of the size field.
        solid_walls_faces = list(set(solid_walls_faces + unassigned))
        gmsh.model.removePhysicalGroups([(2, t) for t in []])  # no-op, kept for clarity
        gmsh.model.addPhysicalGroup(2, solid_walls_faces, name="Solid_Walls")
        print(f"  [OK] Solid_Walls now includes {len(solid_walls_faces)} surfaces\n")
    else:
        print("[OK] No unassigned surfaces\n")

    print("=" * 80 + "\n")

    # ========================================================================
    # 9. MESH PARAMETERS - GLOBAL BOUNDS ONLY (size field does the real work)
    # ========================================================================
    print("Step 9: Setting mesh parameters (from config file)...")

    # These remain as safety bounds - the size field below is the primary
    # control. Kept in the same order/position as the original script.
    gmsh.option.setNumber("Mesh.MeshSizeMin", gmsh_params.get("Mesh.MeshSizeMin", 0.002))
    gmsh.option.setNumber("Mesh.MeshSizeMax", gmsh_params.get("Mesh.MeshSizeMax", 0.15))
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", gmsh_params.get("Mesh.MeshSizeFromCurvature", 24))

    # 2D surface algorithm stays configurable (default 6 = Frontal-Delaunay).
    gmsh.option.setNumber("Mesh.Algorithm", gmsh_params.get("Mesh.Algorithm", 6))

    # 3D algorithm - read from config like every other gmsh option, but
    # defaults to HXT (10) if the config omits it, per the requirement to
    # use HXT-Delaunay for 3D meshing only. HXT has no 2D counterpart, so
    # this does not affect the 2D algorithm above. The config file
    # (gmsh_config_sizefield.json) explicitly sets Mesh.Algorithm3D: 10 -
    # if that is ever edited to something else, warn loudly so the change
    # is visible rather than silent, since HXT is a project requirement.
    algorithm_3d = gmsh_params.get("Mesh.Algorithm3D", 10)
    if algorithm_3d != 10:
        print(f"  [WARNING] Config requests Mesh.Algorithm3D={algorithm_3d}, which is NOT "
              f"HXT (10) - this project requires HXT-Delaunay for 3D meshing.")
    gmsh.option.setNumber("Mesh.Algorithm3D", algorithm_3d)
    print(f"  [OK] Mesh.Algorithm3D = {algorithm_3d} (from config)")

    gmsh.option.setNumber("Mesh.Smoothing", gmsh_params.get("Mesh.Smoothing", 5))
    gmsh.option.setNumber("Mesh.ElementOrder", gmsh_params.get("Mesh.ElementOrder", 1))
    gmsh.option.setNumber("General.Terminal", gmsh_params.get("General.Terminal", 1))

    # Critical: once a background field drives size, turn off the competing
    # automatic size sources so they don't fight the field (this fighting is
    # what caused the cone self-intersection loop previously).
    # NOTE: MeshSizeFromCurvature is intentionally NOT zeroed here -- it was
    # previously being set from config on line 594 and then unconditionally
    # overwritten to 0 a few lines later, silently discarding the configured
    # value and disabling curvature-based refinement everywhere (fuselage,
    # nose, leading edges) even when the config explicitly asked for it.
    # Unlike MeshSizeExtendFromBoundary/MeshSizeFromPoints, curvature sizing
    # supplements the background field rather than competing with it, so it
    # does not need to be disabled for the field to remain the primary driver.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)

    print("[OK] Mesh parameters set (as global bounds; size field is primary control)\n")

    # ========================================================================
    # 9a. SIZE FIELD - fine near geometry, coarse far away, smooth transition
    # ========================================================================
    print("Step 9a: Building spatially-varying size field...")

    dist_field, thresh_field = setup_size_field(solid_walls_faces, size_params)

    # Placeholder for future cavity refinement: pass surface tags for any
    # mm-scale cavity surfaces here once the next STEP file has them. Empty
    # list today means this is a no-op and behaves exactly like the base
    # field alone.
    cavity_surface_tags = size_params.get("cavity_surface_tags", [])
    background_field = add_cavity_refinement_fields(
        cavity_surface_tags, (dist_field, thresh_field), size_params
    )

    gmsh.model.mesh.field.setAsBackgroundMesh(background_field)
    print(f"[OK] Background mesh field set (field id {background_field})\n")

    # ========================================================================
    # 10. VISCOUS LAYER INFLATION - FROM CONFIG
    # ========================================================================
    print("Step 10: Configuring viscous boundary layers...")

    bl_params = config.get('boundary_layer_parameters', {})

    bl_field = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.occ.synchronize()

    all_edges = []
    for surf_id in solid_walls_faces:
        try:
            boundary = gmsh.model.getBoundary([(2, surf_id)], oriented=False)
            edges = [tag for dim, tag in boundary if dim == 1]
            all_edges.extend(edges)
        except Exception as e:
            print(f"  Warning: Could not get boundary for surface {surf_id}: {e}")

    all_edges = list(set(all_edges))
    print(f"  Found {len(all_edges)} boundary layer edges")

    if all_edges:
        gmsh.model.mesh.field.setNumbers(bl_field, "EdgesList", all_edges)
        gmsh.model.mesh.field.setNumber(bl_field, "Size", bl_params.get("Size", 0.002))
        gmsh.model.mesh.field.setNumber(bl_field, "Thickness", bl_params.get("Thickness", 0.0118))
        gmsh.model.mesh.field.setNumber(bl_field, "Ratio", bl_params.get("Ratio", 1.16))
        gmsh.model.mesh.field.setNumber(bl_field, "NbLayers", bl_params.get("NbLayers", 45))
        # BoundaryLayer field takes over as background mesh near walls; it
        # does not conflict with the Threshold field used for the bulk
        # domain because gmsh applies BoundaryLayer specially at walls.
        gmsh.model.mesh.field.setAsBackgroundMesh(bl_field)
        print("[OK] Viscous layers configured:")
        print(f"    - Thickness: {bl_params.get('Thickness', 0.0118)}")
        print(f"    - Stretch ratio: {bl_params.get('Ratio', 1.16)}")
        print(f"    - Number of layers: {bl_params.get('NbLayers', 45)}\n")
    else:
        print("[ERROR] No edges found for boundary layers\n")

    # ========================================================================
    # 11. MESH GENERATION (HXT for 3D only - see Step 9)
    # ========================================================================
    print("Step 11: Generating 3D mesh (HXT)...")
    print("  This may take several minutes...\n")

    try:
        gmsh.model.mesh.generate(3)
        print("[OK] Mesh generation complete\n")
        mesh_success = True
    except Exception as e:
        print(f"[ERROR] Mesh generation failed: {e}\n")
        mesh_success = False

    # ========================================================================
    # 12. EXPORT
    # ========================================================================
    if mesh_success:
        print("Step 12: Exporting mesh...")
        try:
            gmsh.option.setNumber("Mesh.MshFileVersion", gmsh_params.get("Mesh.MshFileVersion", 2.2))
            gmsh.write(output_path)
            file_size = os.path.getsize(output_path)
            print(f"[OK] Mesh exported successfully")
            print(f"  File: {output_path}")
            print(f"  Size: {file_size / 1024 / 1024:.2f} MB\n")
        except Exception as e:
            print(f"[ERROR] Export failed: {e}\n")
    else:
        print("Skipping export - mesh generation failed\n")

    # ========================================================================
    # 13. CLEANUP
    # ========================================================================
    print("Step 13: Finalizing...")
    gmsh.finalize()
    print("[OK] Done\n")

    print("=" * 80)
    if mesh_success:
        print("MESH GENERATION SUCCESSFUL")
        print("=" * 80)
        print("\nNext steps:")
        print("  1. Verify mesh quality in Gmsh GUI")
        print("  2. Check that inlet/outlet are correctly assigned")
        print("  3. Import mesh into your CFD solver")
    else:
        print("MESH GENERATION FAILED")
        print("=" * 80)
        print("\nTroubleshooting:")
        print("  1. Check size_field_parameters in config (size_near/size_far/dist_min/dist_max)")
        print("  2. Verify solid_walls_faces list is non-empty before Step 9a")
        print("  3. Check disk space (mesh is MBs to GB)")

    print("=" * 80 + "\n")

    return True


if __name__ == "__main__":
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'gmsh_config.json'
    success = create_mesh(config_file)
    sys.exit(0 if success else 1)
