# Hadrons [![Teamcity status](http://ci.cliath.ph.ed.ac.uk/app/rest/builds/aggregated/strob:(buildType:(affectedProject(id:GridBasedSoftware_Hadrons)),branch:name:develop)/statusIcon.svg)](http://ci.cliath.ph.ed.ac.uk/project.html?projectId=GridBasedSoftware_Hadrons&tab=projectOverview)
_[Grid](https://github.com/paboyle/Grid)-based workflow management system for
lattice field theory simulations_

## Install
Download and compile the Grid library and install it. Please refer to the
instructions from the [Grid repository](https://github.com/paboyle/Grid). Using
the `develop` branch of Grid is recommended.

Then Hadrons can be downloaded and built using

``` bash
git clone https://github.com/aportelli/Hadrons.git
cd Hadrons
./bootstrap.sh
mkdir build; cd build
../configure --with-grid=<dir>
make -j<N>
```
`<dir>` is the installation prefix ogf Grid and `<N>` is the number of parallel
build tasks. All the compilation flags used for compiling Grid will be reused to
compile Hadrons. You can extend these flags or change the compiler by modifying
the `CXXFLAGS` and `CXX` environment variables.

## Run
The main Hadrons executables are in the `utilities` directory, examples can be
found in the `tests` directory, and can be built using `make tests`.
