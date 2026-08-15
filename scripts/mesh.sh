#!/bin/bash
export OMP_NUM_THREADS=28
nohup env PYTHONUNBUFFERED=1 /home/pravin/SALOME-9.16.0/salome -b ./Gmsh_cloud_run.py > ./output_meshing.log 2>&1 &
