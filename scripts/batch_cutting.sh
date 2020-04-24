#!/bin/bash

EXE=../build-Desktop-Release/volumetric_cutter

for D in `find . -type d`
do
   OBJ=`find $D/*_splitted.obj`
   LOOP=`find $D/*_loop.txt`
   echo 
   echo ------------------------------------
   echo Processing directory $D
   echo Object $OBJ
   echo Loop set $LOOP
   echo ------------------------------------
   echo 
   $EXE $OBJ $LOOP -batch-mode $D > $D/out_log.txt   
done