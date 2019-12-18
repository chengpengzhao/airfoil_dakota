#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory
#------------------------------------------------------------------------------
#generate blocMeshDict
python3 PARSEC.py;

caseName="airfoil_run";
runpath="./cases/${caseName}";
if [ -d ${runpath} ];then
echo "Delete files......\n\n"
rm -rf ${runpath}
fi

foamCloneCase airfoil_template/ ${runpath};
cp blockMeshDict ${runpath}/system;
cd cases/airfoil_run;
blockMesh;
checkMesh;
renumberMesh -overwrite;
#paraFoam 
decomposePar;
mpirun -np 4 simpleFoam -parallel | tee solve.log;
reconstructPar -constant;
pyFoamPlotWatcher.py --progress solve.log;
