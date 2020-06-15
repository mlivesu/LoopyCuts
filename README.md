# Surface Field Driven Strongly Hex-Dominant Meshing

<p align="center"><img src="LoopyCuts.jpg" width="700"></p>

This repository contains the meshing software developed as part of the publication

**LoopyCuts: Practical Feature-Preserving Block Decomposition for Strongly Hex-Dominant Meshing**<br>
[Marco Livesu](http://pers.ge.imati.cnr.it/livesu/), 
[Nico Pietroni](http://www.nicopietroni.com),
[Enrico Puppo](https://www.disi.unige.it/person/PuppoE/), 
[Alla Sheffer](http://www.cs.ubc.ca/%7Esheffa/), 
[Paolo Cignoni](http://vcg.isti.cnr.it/~cignoni/)<br>
ACM Transactions on Graphics (SIGGRAPH 2020)<br>
[PDF](http://pers.ge.imati.cnr.it/livesu/papers/LPPSC20/LPPSC20.pdf)

## Code and external dependencies
To install LoopyCuts on your machine clone this repo in recursive mode with
```
git clone --recursive https://github.com/mlivesu/LoopyCuts.git
```
The software has a number of external dependencies, some of which must be downloaded and installed in your system separately
* [Qt](https://www.qt.io/download): download and install it
* [Tetgen](http://wias-berlin.de/software/tetgen/): install it on your machine and verify that references in `volumetric_cutter/volumetric_cutter.pro` are correct
* [CinoLib](https://github.com/mlivesu/cinolib): already included in the project as a submodule, no action required
* [VcgLib](http://vcg.isti.cnr.it/vcglib/): already included in the project as a submodule, no action required
* [AntTweakBar1.16](http://anttweakbar.sourceforge.net): already included in the project


## Step 1: Generation of cutting loops
This part reads an input mesh, the sharp features and the field. The software depends on [Qt](https://www.qt.io/download) for the GUI, [vcgLibrary](http://vcg.isti.cnr.it/vcglib/) for geometry processing, and [AntTweakBar1.16](http://anttweakbar.sourceforge.net). To compile it, open a terminal in `loop_distributor` and type

```
qmake .
make -j4
```
The program can be used either with a GUI, or by command line (useful to batch run entire datasets of models).
```
./loop_distributor <mesh> [batch]
```
The mesh can be either an obj or a ply.
It requires to have in the same folder a .rosy and a .sharp file (with the same name of the mesh file)

##.rosy file format
```
fn                   // number of faces of the mesh
4                    // directions of the field (always 4 for a cross-field)
x0 y0 z0             // xyz directions of one vector of the cross field of the first face
...
xn yn zn             // xyz directions of one vector of the cross field of the n-th face
```

##.sharp file format
```
sn                   // number of sharp features
t0 f0 e0             // for each sharp edge: first integer is 0 if the edge is concave 1 if convex then the face and the index of the sharp edge
...
tn fn en             // nth sharp edge
```

## Output files

_splitted mesh file  :  the mesh traced and splitted
_loop file           :  the mesh traced and splitted


## _loop file format
```
ln                  // number of loops
for each loop:
REGULAR|CONCAVE|CONVEX      //the kind of loop
Closed|Open                 //if the loop has been succesfully closed or not
Cross OK|FAIL               //if the loop has the right number of cross for topology of the cut
en                          //number of edges of the loop
indexF indexE 0|1           //the index of the face, the edge and if is on a sharp feature or not
```

## Step 2: Cutting and Hex-dominant meshing
This part reads the refined mesh and loops generated at the previous step, and outputs a hex-dominant mesh. The software depends on [Qt](https://www.qt.io/download) for the GUI, [CinoLib](https://github.com/mlivesu/cinolib) for geometry processing, and [Tetgen](http://wias-berlin.de/software/tetgen/) for tetrahedralization. To compile it, open a terminal in `volumetric_cutter` and type
```
qmake .
make -j4
```
The program can be used either with a GUI, or by command line (useful to batch run entire datasets of models). 
```
./volumetric_cutter <mesh> <loops> [ -batch-mode <output_folder> ]
```
We recommend using the command line version, because it is much faster. The `scripts` folder contains a useful bash script for processing large collections of shapes with a single call.

## Output Format
Although almost entirely composed of hexahedra, our output meshes may contain arbitrary polyhedra which cannot be ecnoded in popular volumetric mesh formats such as `.mesh` and `.vtk`. All our outputs are therefore encoded using the `.hedra` format, which is structured as follows
```
nv nf np             // number of vertices, faces and polyhedra, respectively
x0 y0 z0             // xyz coordinates of the 1st point 
x1 y1 z1             // xyz coordinates of the 2nd point
...                  // 
f0 v1 v2 ... vf1     // f1: number of vertices of the 1st face, followed by the (CCW ordered) list of vertices
f1 v1 v2 ... vf2     // f2: number of vertices of the 2nd face, followed by the (CCW ordered) list of vertices
...                  //
p1  f1 -f2 ...  fp1  // p1: number of faces of the 1st polyhedron, followed by the list of faces
p2 -f1  f2 ... -fp2  // p2: number of faces of the 2nd polyhedron, followed by the list of faces
...                  // (references with negative numbers (e.g. -f) denote that face |f| is seen CW by the polyhedron)
```
These meshes can be visualized using [CinoLib](https://github.com/mlivesu/cinolib) (see e.g. example [#06](https://github.com/mlivesu/cinolib/tree/master/examples/06_base_app_polyhedralmesh)). Note that in case the output is a full hexahedral mesh a `.mesh` file will be also produced. Such a file can be visually inspected directly on browser connecting to [HexaLab](https://www.hexalab.net).

## Acknowldegment
If you use LoopyCuts, please consider citing the associated scientific paper using the following 
BibTeX entry:

```bibtex
@article{LoopyCuts2020,
  title   = {LoopyCuts: Practical Feature-Preserving Block Decomposition for Strongly Hex-Dominant Meshing},
  author  = {Livesu, Marco and Pietroni, Nico and Puppo, Enrico and Sheffer, Alla and Cignoni, Paolo},
  journal = {ACM Transactions on Graphics},
  year    = {2020},
  volume  = {39},
  number  = {4},
  doi     = {(10.1145/3386569.3392472)}}
```

