mkdir -p $SCRATCH/$USER/OpenFOAM
cd $SCRATCH/$USER/OpenFOAM

tar -xzf $SCRATCH/$USER/OpenFOAM-v1706.tgz
tar -xzf $SCRATCH/$USER/ThirdParty-v1706.tgz
tar -xzf $SCRATCH/$USER/hypersonicfoam.tar.gz -C $SCRATCH/$USER/OpenFOAM/hypersonicfoam

#nano $SCRATCH/$USER/OpenFOAM/OpenFOAM-v1706/etc/bashrc
#change: foamInstall=$SCRATCH/$USER/OpenFOAM
#change:WM_COMPILER=Gcc
#       WM_MPLIB=INTELMPI

