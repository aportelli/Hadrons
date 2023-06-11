/*
 * LoadA2AMatrixDiskVector.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_
#define Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadA2AMatrixDiskVector                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadA2AMatrixDiskVectorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadA2AMatrixDiskVectorPar,
                                    std::string,  file,
                                    std::string,  dataset,
                                    std::string,  diskVectorDir,
                                    int,  cacheSize);
};

template <typename FImpl>
class TLoadA2AMatrixDiskVector: public Module<LoadA2AMatrixDiskVectorPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    // constructor
    TLoadA2AMatrixDiskVector(const std::string name);
    // destructor
    virtual ~TLoadA2AMatrixDiskVector(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadA2AMatrixDiskVector, TLoadA2AMatrixDiskVector<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadA2AMatrixDiskVector implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadA2AMatrixDiskVector<FImpl>::TLoadA2AMatrixDiskVector(const std::string name)
: Module<LoadA2AMatrixDiskVectorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadA2AMatrixDiskVector<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadA2AMatrixDiskVector<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AMatrixDiskVector<FImpl>::setup(void)
{
    int Ls = 1;

    std::string dvDir = par().diskVectorDir;
    std::string dataset = par().dataset;
    std::string dvFile = dvDir + "/" + getName() + "." + std::to_string(vm().getTrajectory());

    int nt = env().getDim(Tp);
    int cacheSize = par().cacheSize;
    bool clean = true;
    GridBase *grid = envGetGrid(FermionField);

    envCreate(EigenDiskVector<ComplexD>, getName(), Ls, dvFile, nt, cacheSize, clean, grid);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AMatrixDiskVector<FImpl>::execute(void)
{
    int nt = env().getDim(Tp);
    std::string file  = par().file;
    std::string dataset  = par().dataset;
    GridBase *grid = envGetGrid(FermionField);

    auto &mesonFieldDV = envGet(EigenDiskVector<ComplexD>, getName());

    int traj = vm().getTrajectory();
    tokenReplace(file, "traj", traj);
    LOG(Message) << "-- Loading '" << file << "'-- " << std::endl;
    double t;
    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> mfIO(file, dataset, nt);
    mfIO.load(mesonFieldDV, &t, grid);
    LOG(Message) << "Read " << mfIO.getSize() << " bytes in " << t << " usec, " << mfIO.getSize() / t * 1.0e6 / 1024 / 1024 << " MB/s" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadA2AMatrixDiskVector_hpp_
