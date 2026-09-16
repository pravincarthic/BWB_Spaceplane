cd "D:\BWB_git\BWB_Spaceplane\scripts"
d:
set PYTHONUNBUFFERED=1

python3 final_meshing.py gmsh_config.json > ..\Output\PlainBWB_mesh\PlainBWB_log.log 2>&1
