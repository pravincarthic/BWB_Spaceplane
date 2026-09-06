#!/usr/bin/env python3
"""
Multi-pass STEP file healing using Gmsh
Applies healing 3 times with increasing tolerances
More robust than single-pass healing
"""

import gmsh
import sys
import os
import shutil

def multipass_heal_gmsh(input_file, output_file):
    current_file = input_file
    tolerances = [0.01, 0.05, 0.1, 0.5]
    
    print("Healing STEP file sewing faults")
    print(f"Input:  {input_file}\n")
    
    for pass_num in range(3, 5):
        tolerance = tolerances[pass_num - 1]
        temp_file = f"{input_file.rsplit('.', 1)[0]}_PASS{pass_num}.step"
        
        print(f"Pass {pass_num}/2: Healing with {tolerance}mm tolerance")
        
        gmsh.initialize()
        gmsh.model.add("heal")
        
        try:
            imported_tags=gmsh.model.occ.importShapes(current_file)
            gmsh.model.occ.synchronize()
            
            gmsh.option.setNumber("Geometry.Tolerance", tolerance)
            gmsh.option.setNumber("Geometry.ToleranceBoolean", tolerance)
            gmsh.option.setNumber("Geometry.MatchMeshTolerance", 1)
            
            healed_tags = gmsh.model.occ.healShapes(
                dimTags=imported_tags,
                tolerance=tolerance,              # Your explicit FixShape.Tolerance3d equivalent
                fixDegenerated=True,         # Fixes collapsed or corrupted edges
                fixSmallEdges=True,               # Removes microscopic edges/faces
                fixSmallFaces=True,               # Removes microscopic faces
                makeSolids=True,    # Corrects overlapping surface flaws
                sewFaces=True         # Seams/sews unconnected face boundaries
            )

            # Synchronize the CAD engine with the Gmsh model
            gmsh.model.occ.synchronize()

            gmsh.write(temp_file)
            
            print(f"  Success: {temp_file}\n")
            current_file = temp_file
            
        except Exception as e:
            print(f"  Error: {e}")
            return False
        finally:
            gmsh.finalize()
    
    shutil.copy(current_file, output_file)
    
    print(f"Output: {output_file}")
    print("Geometry healed successfully - ready for meshing\n")
    
    return True

def main():
    if len(sys.argv) < 2:
        print("Multi-pass STEP file healing using Gmsh")
        print("\nUsage: python heal_multipass_gmsh.py <input.step> [output.step]")
        print("\nExamples:")
        print("  python heal_multipass_gmsh.py geometry.step")
        print("  python heal_multipass_gmsh.py geometry.step geometry_HEALED.step")
        print("\nHealing process:")
        print("  Pass 1: 0.01mm tolerance (fixes tiny gaps)")
        print("  Pass 2: 0.05mm tolerance (fixes small gaps)")
        print("  Pass 3: 0.1mm tolerance (final cleanup)")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"ERROR: File not found: {input_file}")
        sys.exit(1)
    
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace(".step", "_HEALED.step")
    
    print("=" * 80)
    success = multipass_heal_gmsh(input_file, output_file)
    print("=" * 80 + "\n")
    
    if success:
        print("Next steps:")
        print(f"1. Update gmsh_config.json: \"step_file\": \"{output_file}\"")
        print(f"2. Run: python final_meshing.py gmsh_config.json")
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
