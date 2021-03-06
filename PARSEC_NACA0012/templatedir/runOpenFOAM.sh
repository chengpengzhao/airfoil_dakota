#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory
#------------------------------------------------------------------------------
#generate blocMeshDict
python3 PARSEC_0012.py;

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
#checkMesh;
renumberMesh -overwrite;
#paraFoam 
decomposePar;
mpirun --allow-run-as-root -np 12 simpleFoam -parallel | tee solve.log;
reconstructPar -constant;
#pyFoamPlotWatcher.py --progress solve.log;
rm -r processor*;
mv 2000 result;
rm -rf *[0-9]*;
cd ..;
cd ..;
python3 readClCd.py;
