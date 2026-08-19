import gmsh
import math
import os
import sys
import json

def load_config(config_file):
    """Load gmsh parameters and file paths from config file."""
    if not os.path.exists(config_file):
        print(f"✗ Error: Config file not found: {config_file}")
        return None
    
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        print(f"✓ Config loaded from: {config_file}\n")
        return config
    except json.JSONDecodeError as e:
        print(f"✗ Error parsing config file: {e}")
        return None
    except Exception as e:
        print(f"✗ Error reading config file: {e}")
        return None

def create_mesh(config_file='gmsh_config.json'):
    print("\n" + "="*80)
    print("BWB SPACEPLANE MESH GENERATION - TEMPLATE BASED")
    print("="*80 + "\n")

    # Load configuration
    config = load_config(config_file)
    if config is None:
        return False
    
    # Extract file paths
    step_file = config.get('step_file')
    output_path = config.get('output_file')
    
    # Extract gmsh parameters
    gmsh_params = config.get('gmsh_parameters', {})
    
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
    gmsh.model.add("BWB_Spaceplane_Mesh")
    gmsh.option.setNumber("General.NumThreads", gmsh_params.get("General.NumThreads", 28))
    print("✓ Gmsh initialized\n")
    
    # ========================================================================
    # 2. IMPORT GEOMETRY
    # ========================================================================
    print("Step 2: Importing STEP geometry...")
    
    if not os.path.exists(step_file):
        print(f"✗ Error: File not found: {step_file}")
        gmsh.finalize()
        return False
    
    try:
        spaceplane_shapes = gmsh.model.occ.importShapes(step_file)
        gmsh.model.occ.synchronize()
        print("✓ STEP file imported\n")
    except Exception as e:
        print(f"✗ Failed to import STEP: {e}")
        gmsh.finalize()
        return False
    
    # ========================================================================
    # 3. EXTRACT SPACEPLANE VOLUME
    # ========================================================================
    print("Step 3: Extracting spaceplane volume...")
    spaceplane_vols = [item for item in spaceplane_shapes if item[0] == 3]
    
    if not spaceplane_vols:
        print("✗ Error: Could not find 3D volume in STEP file")
        gmsh.finalize()
        return False
    
    spaceplane_tag = spaceplane_vols[0][1]
    print(f"✓ Spaceplane volume tag: {spaceplane_tag}\n")
    
    # ========================================================================
    # 4. ROTATE SPACEPLANE
    # ========================================================================
    print("Step 4: Rotating spaceplane (-5° around Z-axis)...")
    angle_rad = -5.0 * math.pi / 180.0
    gmsh.model.occ.rotate([(3, spaceplane_tag)], 0, 0, 0, 0, 0, 1, angle_rad)
    gmsh.model.occ.synchronize()

    boundary_entities = gmsh.model.getBoundary([(3, spaceplane_tag)], combined=False, oriented=False)
    spaceplane_faces = [tag for dim, tag in boundary_entities if dim == 2]
    
    print("\n" + "="*80)
    print(f"GEOMETRY ANALYSIS: {len(spaceplane_faces)} SURFACES DETECTED")
    print("="*80)
    print(f"Raw Face Tags: {sorted(spaceplane_faces)}\n")
    
    # 6. Print Coordinates and Geometry Details
    print(f"{'Index':<6} | {'Tag':<5} | {'Center of Mass (X, Y, Z)':<32} | {'Dimensions (ΔX x ΔY x ΔZ)':<28}")
    print("-" * 80)
    
    for idx, tag in enumerate(sorted(spaceplane_faces), start=1):
        # Get absolute bounding box
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        
        # Get center of mass
        com = gmsh.model.occ.getCenterOfMass(2, tag)
        
        # Calculate span lengths
        dx = xmax - xmin
        dy = ymax - ymin
        dz = zmax - zmin
        
        com_str = f"({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})"
        dim_str = f"{dx:.2f} x {dy:.2f} x {dz:.2f}"
        
        print(f"#{idx:<5} | {tag:<5} | {com_str:<32} | {dim_str:<28}")
        print(f"       └─ Span Limits: X=[{xmin:.1f} → {xmax:.1f}], Y=[{ymin:.1f} → {ymax:.1f}], Z=[{zmin:.1f} → {zmax:.1f}]")
        print("-" * 80)

    print("✓ Spaceplane rotated\n")
    
    # ========================================================================
    # 5. CREATE FAR-FIELD CONE
    # ========================================================================
    print("Step 5: Creating conical far-field domain...")
    cone_tag = gmsh.model.occ.addCone(-50.0*1000, 0.0, 0.0, 200.0*1000, 0.0, 0.0, 40.0*1000, 100.0*1000)
    gmsh.model.occ.synchronize()
    print(f"✓ Cone created with tag: {cone_tag}\n")
    
    # ========================================================================
    # 6. BOOLEAN CUT
    # ========================================================================
    print("Step 6: Cutting spaceplane from cone (fluid domain)...")
    try:
        fluid_domain, _ = gmsh.model.occ.cut([(3, cone_tag)], [(3, spaceplane_tag)],
                                             removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        fluid_volume_tag = fluid_domain[0][1]
        print(f"✓ Boolean cut complete, fluid volume tag: {fluid_volume_tag}\n")

        # --------------------------------------------------------------------
        # Robust Surface Detection Logic
        # --------------------------------------------------------------------
        # Get every face on the boundary of our newly generated fluid domain
        print("Step 7: Creating physical groups from surfaces and volumes...")
        boundary_faces = gmsh.model.getBoundary([(3, fluid_volume_tag)], combined=False, oriented=False)
        all_fluid_faces = [tag for dim, tag in boundary_faces]

        # Loop through all fluid faces and query their precise spatial properties
        for face in all_fluid_faces:
            # Check center of mass (the average physical location of the surface)
            com = gmsh.model.occ.getCenterOfMass(2, face)
            
            # 1. Catch Inlet: Face center is right at X = -50
            if abs(com[0] - (-50.0*1000)) < 0.2:
                inlet_faces.append(face)
                
            # 2. Catch Outlet: Face center is right at X = 150
            elif abs(com[0] - 150.0*1000) < 0.2:
                outlet_faces.append(face)
                
            # 3. Catch Atmosphere: Center sits on the central X-axis line, but in the middle
            elif -50.0*1000 < com[0] < 150.0*1000:
                # The curved side of the cone balances perfectly along Y=0 and Z=0
                if abs(com[1]) < 0.01 and abs(com[2]) < 0.01:
                    atmosphere_faces.append(face)

        # Assign to Gmsh Physical Groups immediately
        if inlet_faces:
            gmsh.model.addPhysicalGroup(2, inlet_faces, name="Inlet")
        if outlet_faces:
            gmsh.model.addPhysicalGroup(2, outlet_faces, name="Outlet")
        if atmosphere_faces:
            gmsh.model.addPhysicalGroup(2, atmosphere_faces, name="Atmosphere")
        far_field_faces = set(inlet_faces + outlet_faces + atmosphere_faces)
        solid_walls_faces = list(set(all_fluid_faces) - far_field_faces)
        gmsh.model.addPhysicalGroup(2, solid_walls_faces, name="Solid_Walls")
        gmsh.model.addPhysicalGroup(3, [fluid_volume_tag], 7, "Group_Of_All_Volumes")
        print("✓ Physical groups created\n")

    except Exception as e:
        print(f"✗ Boolean operation failed: {e}")
        gmsh.finalize()
        return False
    
    # ========================================================================
    # 8. BOUNDARY FACE ASSIGNMENT
    # ========================================================================
    print("Step 8: Verifying boundary face assignments...")
    print(f"Solid walls: {solid_walls_faces}")
    print(f"Inlet: {inlet_faces}")
    print(f"Outlet: {outlet_faces}")
    print(f"Atmosphere: {atmosphere_faces}\n")

    print("Verifying face assignments:\n")
    
    # Get all surface tags in model
    surfaces = gmsh.model.getEntities(2)
    model_surface_tags = [tag for dim, tag in surfaces]
    
    print(f"Total surfaces in model: {len(model_surface_tags)}")
    print(f"All surface tags: {sorted(model_surface_tags)}\n")
    
    # Check each assignment
    print("="*80)
    all_assigned = inlet_faces + outlet_faces + atmosphere_faces + solid_walls_faces
    all_invalid = []
    
    if solid_walls_faces:
        print("SOLID WALLS:")
        print("-" * 50)
        for tag in solid_walls_faces:
            if tag in model_surface_tags:
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
                size = max(xmax-xmin, ymax-ymin, zmax-zmin)
                print(f"  ✓ Tag {tag:2d}: Size={size:7.2f}, Y range: {ymin:7.2f}→{ymax:7.2f}, Z range: {zmin:7.2f}→{zmax:7.2f}")
            else:
                print(f"  ✗ Tag {tag:2d}: DOES NOT EXIST IN MODEL ⚠️")
                all_invalid.append(tag)
    
    if inlet_faces:
        print("\nINLET:")
        print("-" * 50)
        for tag in inlet_faces:
            if tag in model_surface_tags:
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
                size = max(xmax-xmin, ymax-ymin, zmax-zmin)
                print(f"  ✓ Tag {tag:2d}: Size={size:7.2f}, Y range: {ymin:7.2f}→{ymax:7.2f}, Z range: {zmin:7.2f}→{zmax:7.2f}")
            else:
                print(f"  ✗ Tag {tag:2d}: DOES NOT EXIST IN MODEL ⚠️")
                all_invalid.append(tag)
    
    if outlet_faces:
        print("\nOUTLET:")
        print("-" * 50)
        for tag in outlet_faces:
            if tag in model_surface_tags:
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
                size = max(xmax-xmin, ymax-ymin, zmax-zmin)
                print(f"  ✓ Tag {tag:2d}: Size={size:7.2f}, Y range: {ymin:7.2f}→{ymax:7.2f}, Z range: {zmin:7.2f}→{zmax:7.2f}")
            else:
                print(f"  ✗ Tag {tag:2d}: DOES NOT EXIST IN MODEL ⚠️")
                all_invalid.append(tag)
    
    if atmosphere_faces:
        print("\nATMOSPHERE:")
        print("-" * 50)
        for tag in atmosphere_faces:
            if tag in model_surface_tags:
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
                size = max(xmax-xmin, ymax-ymin, zmax-zmin)
                print(f"  ✓ Tag {tag:2d}: Size={size:7.2f}, Y range: {ymin:7.2f}→{ymax:7.2f}, Z range: {zmin:7.2f}→{zmax:7.2f}")
            else:
                print(f"  ✗ Tag {tag:2d}: DOES NOT EXIST IN MODEL ⚠️")
                all_invalid.append(tag)
    else:
        print("\nATMOSPHERE: None assigned (OK for STEP-only geometry)")
    
    # Check for unassigned surfaces
    unassigned = [t for t in model_surface_tags if t not in all_assigned]
    if unassigned:
        print("\n⚠️  UNASSIGNED SURFACES:")
        print("-" * 50)
        print(f"  These surfaces are in the model but not assigned to any boundary:")
        for tag in sorted(unassigned):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
            size = max(xmax-xmin, ymax-ymin, zmax-zmin)
            print(f"    Tag {tag:2d}: Size={size:7.2f}, Y range: {ymin:7.2f}→{ymax:7.2f}, Z range: {zmin:7.2f}→{zmax:7.2f}")
    
    # ========================================================================
    # FINAL CHECK
    # ========================================================================
    print("\n" + "="*80)
    if all_invalid:
        print("❌ CONFIGURATION HAS ERRORS")
        print("="*80)
        print(f"\nInvalid face tags: {all_invalid}")
        print("\nYou CANNOT proceed with meshing until these are fixed.")
        print("\nSolutions:")
        print("  1. Run find_face_tags_interactive.py to find correct tags")
        print("  2. Update your meshing script with correct tags")
        print("  3. Re-run this verification script")
    else:
        print("✅ CONFIGURATION VALID - READY FOR MESHING")
        print("="*80)
        print(f"\nTotal faces assigned: {len(all_assigned)}")
        print(f"  - Solid walls: {len(solid_walls_faces)}")
        print(f"  - Inlet: {len(inlet_faces)}")
        print(f"  - Outlet: {len(outlet_faces)}")
        print(f"  - Atmosphere: {len(atmosphere_faces)}")
        
        if unassigned:
            print(f"\n⚠️  Note: {len(unassigned)} unassigned surface(s)")
            print("   If this is intentional (e.g., internal surfaces), that's OK")
            print("   Otherwise, add them to one of the boundary groups")
        
        print("\nYou can now run your meshing script safely.")
    print("\n" + "="*80 + "\n")

    # ========================================================================
    # 9. MESH PARAMETERS - FROM CONFIG FILE
    # ========================================================================
    print("Step 9: Setting mesh parameters (from config file)...")
    
    # Element size
    gmsh.option.setNumber("Mesh.MeshSizeMin", gmsh_params.get("Mesh.MeshSizeMin", 0.002))
    gmsh.option.setNumber("Mesh.MeshSizeMax", gmsh_params.get("Mesh.MeshSizeMax", 0.15))
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", gmsh_params.get("Mesh.MeshSizeFromCurvature", 24))
    
    # Algorithms
    gmsh.option.setNumber("Mesh.Algorithm", gmsh_params.get("Mesh.Algorithm", 6))
    gmsh.option.setNumber("Mesh.Algorithm3D", gmsh_params.get("Mesh.Algorithm3D", 10))
    
    # Quality
    gmsh.option.setNumber("Mesh.Smoothing", gmsh_params.get("Mesh.Smoothing", 5))
    gmsh.option.setNumber("Mesh.ElementOrder", gmsh_params.get("Mesh.ElementOrder", 1))
    gmsh.option.setNumber("General.Terminal", gmsh_params.get("General.Terminal", 1))
    
    print("✓ Mesh parameters set\n")
    
    # ========================================================================
    # 10. VISCOUS LAYER INFLATION - FROM CONFIG
    # ========================================================================
    print("Step 10: Configuring viscous boundary layers...")
    
    bl_params = config.get('boundary_layer_parameters', {})
    
    bl_field = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.occ.synchronize()

    # Get edges from solid wall surfaces
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
        gmsh.model.mesh.field.setAsBackgroundMesh(bl_field)
        print("✓ Viscous layers configured:")
        print(f"    - Thickness: {bl_params.get('Thickness', 0.0118)}")
        print(f"    - Stretch ratio: {bl_params.get('Ratio', 1.16)}")
        print(f"    - Number of layers: {bl_params.get('NbLayers', 45)}\n")
    else:
        print("✗ No edges found for boundary layers\n")
    
    # ========================================================================
    # 11. MESH GENERATION
    # ========================================================================
    print("Step 11: Generating 3D mesh...")
    print("  This may take several minutes...\n")
    
    try:
        gmsh.model.mesh.generate(3)
        print("✓ Mesh generation complete\n")
        mesh_success = True
    except Exception as e:
        print(f"✗ Mesh generation failed: {e}\n")
        mesh_success = False
    
    # ========================================================================
    # 12. EXPORT
    # ========================================================================
    if mesh_success:
        print("Step 12: Exporting mesh...")
        
        try:
            # Use Gmsh 2.2 format (legacy, same as existing mesh)
            gmsh.option.setNumber("Mesh.MshFileVersion", 22)
            gmsh.write(output_path)
            
            file_size = os.path.getsize(output_path)
            print(f"✓ Mesh exported successfully")
            print(f"  File: {output_path}")
            print(f"  Size: {file_size/1024/1024:.2f} MB\n")
            
        except Exception as e:
            print(f"✗ Export failed: {e}\n")
    else:
        print("Skipping export - mesh generation failed\n")

    # ========================================================================
    # 13. CLEANUP
    # ========================================================================
    print("Step 13: Finalizing...")
    gmsh.finalize()
    print("✓ Done\n")

    print("="*80)
    if mesh_success:
        print("MESH GENERATION SUCCESSFUL")
        print("="*80)
        print("\nNext steps:")
        print("  1. Verify mesh quality in Gmsh GUI")
        print("  2. Check that inlet/outlet are correctly assigned")
        print("  3. Import mesh into your CFD solver")
    else:
        print("MESH GENERATION FAILED")
        print("="*80)
        print("\nTroubleshooting:")
        print("  1. Check that STEP file exists and is valid")
        print("  2. Verify face tags are correct (run inspect_geometry_15surfaces.py)")
        print("  3. Check disk space (mesh is ~300+ MB)")

    print("="*80 + "\n")

    return True  # Return True if configuration is valid


if __name__ == "__main__":
    # You can optionally pass a custom config file path
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'gmsh_config.json'
    success = create_mesh(config_file)
    sys.exit(0 if success else 1)
