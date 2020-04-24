# Surface Field Driven Strongly Hex-Dominant Meshing

This repository contains the meshing software developed as part of the publication

**LoopyCuts: Practical Feature-Preserving Block Decomposition for Strongly Hex-Dominant Meshing**<br>
[Marco Livesu](http://pers.ge.imati.cnr.it/livesu/), 
[Nico Pietroni](http://vcg.isti.cnr.it/~pietroni/), 
[Enrico Puppo](https://www.disi.unige.it/person/PuppoE/), 
[Alla Sheffer](http://www.cs.ubc.ca/%7Esheffa/), 
[Paolo Cignoni](http://vcg.isti.cnr.it/~cignoni/)<br>
ACM Transactions on Graphics (SIGGRAPH 2020)


## Step #1: Generation of cutting loops
NICO TODO

## Step #2: Cutting and Hex-dominant meshing
This part reads the loop data generated at the previous step, and outputs a hex-dominant mesh. The software depends on [Qt](https://www.qt.io/download) for the GUI, [Cinolib](https://github.com/mlivesu/cinolib) for geometry processing, and [Tetgen](http://wias-berlin.de/software/tetgen/) for tetrahedralization. All dependencies must be installed beforehand, and properly referred to in the project file `volumetric_cutter/volumetric_cutter.pro`. One configured, compiling the project should as easy as opening a terminal in the same folder and typing
```
qmake .
make -j4
```

### Usage
The program can be used either with a GUI, or by command line (useful to batch run entire datasets of models).
```
./volumetric_cutter <mesh> <loops> [ -batch-mode <output_folder> ]
```

