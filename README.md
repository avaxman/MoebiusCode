# MoebiusCode

This code serves as a Demo for the paper [Conformal Mesh Deformations with Möbius Transformations](https://www.staff.science.uu.nl/~vaxma001/Conformal_Mesh_Deformations_with_Mobius_Transformations.pdf) by Vaxman *et al.* from SIGGRAPH 2015. The demo uses the optimization package, and specialized traits from the [libhedra](https://github.com/avaxman/libhedra) package, which is in turn based on [libigl](http://libigl.github.io/libigl/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It implements the following features from the paper:

1. 2D deformation and interpolation using complex numbers.
2. 3D deformation using quaternions.
3. Exact metric conformal (2D and 3D) and intersection-angle preserving (2D).
4. Working with any polyhedral (non-triangular) meshes, where every face undergoes a single M\"obius transformation

See paper for details on the mathematical formulation of the above.

##Demo

To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/MoebiusCode.git
```

to compile the demo, enter the following commands in a terminal:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

You can use `cmake-gui` in windows to create a visual studio project. The dependencies are already recusrisvely downloaded and built.

Upon running, the demo will ask to load a mesh. It automatically infers if the mesh is 2D (if the $z$ coordinate is negligible) of 3D, and then would run the respective algorithms. The demo offers the following possibilities:

![Demo Screen](MoebiusCodeDemoScreen.png)

The editing options are:

| Editing Options                     | Description                                                                         |
| :----------------------- | :---------------------------------------------------------------------------------- |
| `Original`            | Shows the original mesh for comparison. It cannot be edited in this mode.                                     |
| `Deformation`               | Shows the currently deforming mesh. The handles will be spheres, where the green sphere is the currently edited one.|
| `Interpolation`              | Shows the interpolated mesh at the current time frame.|



