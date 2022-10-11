<img src="graphics/hadrons-icon-title.png" height=100px>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4063666.svg)](https://doi.org/10.5281/zenodo.4063666)  
_[Grid](https://github.com/paboyle/Grid)-based workflow management system for
lattice field theory simulations_  

---

__If Grid & Hadrons played an important role in the production of scientific results leading to a peer-reviewed publication, we would be grateful if you consider citing the GitHub repository in your paper, and/or invite some contributors for authorship if relevant (especially PhD students & postdoctoral researchers).__

To generate a BibTeX citation, please see the [Zenodo page](https://doi.org/10.5281/zenodo.4063666).

Documentation (work in progress): https://aportelli.github.io/Hadrons-doc/.

## Install
Download and compile the Grid library and install it. Please refer to the
instructions from the [Grid repository](https://github.com/paboyle/Grid). Using
the `develop` branch of Grid is recommended.

Hadrons can be downloaded and built using

``` bash
git clone https://github.com/aportelli/Hadrons.git
cd Hadrons
./bootstrap.sh
mkdir build; cd build
../configure --with-grid=<dir>
make -j<N>
```
`<dir>` is the installation prefix of Grid and `<N>` is the number of parallel
build tasks. All the compilation flags used for compiling Grid will be reused to
compile Hadrons. You can extend these flags or change the compiler by modifying
the `CXXFLAGS` and `CXX` environment variables.

## Run
The main Hadrons executables are in the `utilities` directory, examples can be
found in the `tests` directory, and can be built using `make tests`.
