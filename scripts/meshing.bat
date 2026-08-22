cd "D:\BWB Spaceplane\Meshing Scripts\V2_cleaned"
d:
set PYTHONUNBUFFERED=1

python3 final_meshing_sizefield_updated.py gmsh_config_sizefield.json > log.log 2>&1
