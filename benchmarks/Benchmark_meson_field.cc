/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: benchmarks/Benchmark_meson_field.cc

Copyright (C) 2015-2018

Author: Felix Erben <felix.erben@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/qcd/utils/A2Autils.h>

using namespace Grid;

const int VDIM = 20; // length of each vector
const int NMOM = 27; // number of momenta
const int NGAM = 16; // number of gamma matrices

typedef typename DomainWallFermionR::ComplexField ComplexField;
typedef typename DomainWallFermionR::FermionField FermionField;

int main(int argc, char *argv[])
{
  // initialization
  Grid_init(&argc, &argv);
  std::cout << GridLogMessage << "Grid initialized" << std::endl;

  // Lattice and rng setup 
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(4, vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian    grid(latt_size,simd_layout,mpi_layout);
  int Nt = GridDefaultLatt()[Tp];
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&grid);
  pRNG.SeedFixedIntegers(seeds);

  // MesonField lhs and rhs vectors
  std::vector<FermionField> phi(VDIM,&grid);
  std::cout << GridLogMessage << "Initialising random meson field" << std::endl;
  for (unsigned int i = 0; i < VDIM; ++i){
    random(pRNG,phi[i]);
  }
  std::cout << GridLogMessage << "Meson field initialised" << std::endl;

  // Gamma matrices used in the contraction
  std::vector<Gamma::Algebra> Gmu(NGAM); 
  for (unsigned int j = 0; j < NGAM; ++j)
  {
    Gmu[j]=Gamma::Algebra::GammaT;
  }

  // momentum phases e^{ipx}
  std::vector<std::vector<double>> momenta (NMOM);
  for (unsigned int j = 0; j < NMOM; ++j)
  {
    momenta[j]= {1.,1.,1.};
  }

  std::cout << GridLogMessage << "Meson fields will be created for " << Gmu.size() << " Gamma matrices and " << momenta.size() << " momenta." << std::endl;

  std::cout << GridLogMessage << "Computing complex phases" << std::endl;
  std::vector<ComplexField> phases(momenta.size(),&grid);
  ComplexField coor(&grid);
  Complex Ci(0.0,1.0);
  for (unsigned int j = 0; j < momenta.size(); ++j)
  {
    phases[j] = Zero();
    for(unsigned int mu = 0; mu < momenta[j].size(); mu++)
    {
      LatticeCoordinate(coor, mu);
      phases[j] = phases[j] + momenta[j][mu]/GridDefaultLatt()[mu]*coor;
    }
    phases[j] = exp((Real)(2*M_PI)*Ci*phases[j]);
  }
  std::cout << GridLogMessage << "Computing complex phases done." << std::endl;

  Eigen::Tensor<ComplexD,5, Eigen::RowMajor> Mpp(momenta.size(),Gmu.size(),Nt,VDIM,VDIM);

  // timer
  double start,stop;

  //execute meson field routine
  start = usecond();
  A2Autils<WilsonImplR>::MesonField(Mpp,&phi[0],&phi[0],Gmu,phases,Tp);
  stop = usecond();
  std::cout << GridLogMessage << "M(phi,phi) created, execution time " << stop-start << " us" << std::endl;
  start = usecond();

  std::string FileName = "Meson_Field";
#ifdef HAVE_HDF5
  using Default_Reader = Grid::Hdf5Reader;
  using Default_Writer = Grid::Hdf5Writer;
  FileName.append(".h5");
#else
  using Default_Reader = Grid::BinaryReader;
  using Default_Writer = Grid::BinaryWriter;
  FileName.append(".bin");
#endif

  Default_Writer w(FileName);
  write(w,"phi_phi",Mpp);

  // epilogue
  std::cout << GridLogMessage << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
