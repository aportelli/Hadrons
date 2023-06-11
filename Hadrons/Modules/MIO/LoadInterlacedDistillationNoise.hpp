/*
 * LoadInterlacedDistillationNoise.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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
#ifndef Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_
#define Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadInterlacedDistillationNoise                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadInterlacedDistillationNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadInterlacedDistillationNoisePar,
                                    unsigned int, ti,
                                    unsigned int, li,
                                    unsigned int, si,
                                    unsigned int, nNoise,
                                    std::string, lapEigenPack,
                                    std::string, fileName,);
};

template <typename FImpl>
class TLoadInterlacedDistillationNoise: public Module<LoadInterlacedDistillationNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename DistillationNoise<FImpl>::Index Index;
public:
    // constructor
    TLoadInterlacedDistillationNoise(const std::string name);
    // destructor
    virtual ~TLoadInterlacedDistillationNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadInterlacedDistillationNoise, TLoadInterlacedDistillationNoise<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadInterlacedDistillationNoise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadInterlacedDistillationNoise<FImpl>::TLoadInterlacedDistillationNoise(const std::string name)
: Module<LoadInterlacedDistillationNoisePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadInterlacedDistillationNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack};
    
    return in;
}

template <typename FImpl>
DependencyMap TLoadInterlacedDistillationNoise<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().lapEigenPack, getName()});

    return dep;
}

template <typename FImpl>
std::vector<std::string> TLoadInterlacedDistillationNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedDistillationNoise<FImpl>::setup(void)
{
    auto          &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                                    par().lapEigenPack);
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);

    envCreateDerived(DistillationNoise<FImpl>, InterlacedDistillationNoise<FImpl>,
                     getName(), 1, g, g3d, epack, par().ti, par().li, par().si, 
                     par().nNoise);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedDistillationNoise<FImpl>::execute(void)
{
    LOG(Message) << "Loading interlaced distillation noise with (ti, li, si) = (" 
                 << par().ti << ", " << par().li << ", " << par().si << ")"
                 << std::endl;

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());
    noise.load(par().fileName, "InterlacedDistillation", vm().getTrajectory());
    noise.dumpDilutionMap();
    auto hash = noise.generateHash();
    LOG(Message) << "Noise hit hashes : " << std::endl;
    for(auto& h: hash)
    {
        LOG(Message) << h << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_
