---
icon: lucide/rocket
---

# Getting Started with CabanaPD

## Requirements

CabanaPD has the following dependencies:

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.11+    | Yes      | Build system
|Cabana     | edb0c7cd | Yes      | Performance portable particle algorithms
|GTest      | 1.10+    | No       | Unit test framework

Cabana must be built with the following in order to work with CabanaPD:

|Cabana Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|CMake      | 3.16+   | Yes      | Build system
|MPI        | GPU-Aware if CUDA/HIP enabled | Yes | Message Passing Interface
|Kokkos     | 4.1.0+  | Yes      | Performance portable on-node parallelism
|HDF5       | master  | No       | Particle output
|SILO       | master  | No       | Particle output

The underlying parallel programming models are available on most systems, as is
CMake. Those must be installed first, if not available. Kokkos and Cabana are
available on some systems or can be installed with [`spack`](https://spack.readthedocs.io/en/latest/getting_started.html):

```
spack install cabana@master+grid+hdf5
```

Alternatively, Kokkos can be built locally, followed by Cabana:
https://github.com/ECP-CoPA/Cabana/wiki/1-Build-Instructions

Build instructions are available for both CPU and GPU. Note that Cabana must be
compiled with MPI and the Grid sub-package.


## Installation

### Obtaining CabanaPD

Clone the master branch:

```
git clone https://github.com/ORNL/CabanaPD.git
```

### Build and install
#### CPU Build

After building Kokkos and Cabana for CPU, the following script will build and install CabanaPD:

```
#Change directory as needed
export CABANA_DIR=$HOME/Cabana/build/install

cd ./CabanaPD
mkdir build
cd build
pwd
cmake \
    -D CMAKE_PREFIX_PATH="$CABANA_DIR" \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CabanaPD_ENABLE_TESTING=ON \
    -D CabanaPD_ENABLE_EXAMPLES=ON \
    .. ;
make install
```

#### CUDA Build

[After building Kokkos and Cabana for Cuda](https://github.com/ECP-CoPA/Cabana/wiki/Build-CUDA):

The CUDA build script is identical to that above, but again note that Kokkos
must be compiled with the CUDA backend.

Note that the same compiler should be used for Kokkos, Cabana, and CabanaPD.

#### HIP Build

[After building Kokkos and Cabana for HIP](https://github.com/ECP-CoPA/Cabana/wiki/Build-HIP-and-SYCL#HIP):

The HIP build script is identical to that above, except that `hipcc` compiler
must be used:

```
-D CMAKE_CXX_COMPILER=hipcc
```

Note that `hipcc` should be used for Kokkos, Cabana, and CabanaPD.


## Usage

### Examples

CabanaPD includes many examples;
see details on running individual examples in [examples](user/examples.md).

New examples can be created by using any of the current cases as a template.
Examples can be built by updating the CabanaPD CMake configuration in the
script above with:

```
-D CabanaPD_ENABLE_EXAMPLES=ON
```

Once built and installed, the CabanaPD examples in the `examples/` directory can be run. Timing and energy
information is written to file and particle output is written if enabled within Cabana in formats that can be [visualized](#visualizing-with-paraview).

Most inputs are specified in the example JSON files; some inputs are set directly in the example `.cpp` files. See more details on inputs in [inputs](API/input_files.md).


## Visualizing with Paraview

As mentioned above, the simulation results can be visualized with Paraview or similar applications.

### How to Install

The installation instructions can be found [here](https://www.paraview.org/download/). Ensure you select the appropriate version based on your operating system.

### Importing Files

Once Paraview is installed, the following simulation output file group should be imported to view the results: `particles_..xmf` for HDF5 or `particles_..silo` for SILO.

If shown the option to select a reader type, select `XDMF Reader` in the "Open Data With ..." window for HDF5.

### Viewing Results

Below are some basic guidelines for how to perform the initial steps in order to view and analyze the results. A more in-depth tutorial for Paraview can be found [here](https://docs.paraview.org/en/latest/Tutorials/SelfDirectedTutorial/index.html).

1. Select `Apply` in the lower left-hand Properties window. This will load your simulation data.

2. In the Properties window, under Representation, `Surface` will be selected by default as the geometry representation. Change this to `Point Gaussian`.

3. Different output fields can be selected within the Coloring menu below Representation.

4. To control the size of the visualized points, scroll down within the Properties window until the Point Gaussian menu and choose a value for Gaussian Radius.
