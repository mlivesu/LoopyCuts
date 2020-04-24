#!/bin/bash

EXE=

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
   $EXE $OBJ $LOOP -batch-mode $D > $D/log.txt   
done