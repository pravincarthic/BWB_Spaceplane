#!/bin/bash
export OMP_NUM_THREADS=64
nohup /home/pravin/SALOME-9.16.0/salome -b ./Gmsh_cloud_run.py > ./meshing_output.log 2>&1 &
