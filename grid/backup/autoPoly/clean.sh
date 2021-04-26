#!/bin/bash
#preprocessing stuff
rm log.*
rm *.obj
rm -r VTK
rm -rf constant/polyMesh/*
#run time stuff
rm -rf 0.* 1* 2* 3* 4* 5* sets parameters.log pyresult.dat processor*
rm autoPoly.foam
