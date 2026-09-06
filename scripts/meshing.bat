cd "D:\BWB_git\BWB_Spaceplane\scripts"
d:
set PYTHONUNBUFFERED=1

python3 final_meshing.py gmsh_config.json > log.log 2>&1
