/*
 * ExactDistillation.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MNoise_ExactDistillation_hpp_
#define Hadrons_MNoise_ExactDistillation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Interlaced distillation noise module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class ExactDistillationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDistillationPar,
                                    std::string, lapEigenPack);
};

template <typename FImpl>
class TExactDistillation: public Module<ExactDistillationPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TExactDistillation(const std::string name);
    // destructor
    virtual ~TExactDistillation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDistillation, TExactDistillation<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZExactDistillation, TExactDistillation<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TExactDistillation implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExactDistillation<FImpl>::TExactDistillation(const std::string name)
: Module<ExactDistillationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExactDistillation<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack};
    
    return in;
}

template <typename FImpl>
DependencyMap TExactDistillation<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().lapEigenPack, getName()});

    return dep;
}

template <typename FImpl>
std::vector<std::string> TExactDistillation<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactDistillation<FImpl>::setup(void)
{
    auto          &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                                  par().lapEigenPack);
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);

    envCreateDerived(DistillationNoise<FImpl>, ExactDistillationPolicy<FImpl>,
                     getName(), 1, g, g3d, epack);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactDistillation<FImpl>::execute(void)
{
    LOG(Message) << "Generating exact distillation policy" << std::endl;

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                ExactDistillationPolicy<FImpl>, getName());

    noise.dumpDilutionMap();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_ExactDistillation_hpp_
