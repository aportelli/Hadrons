/*
 * Loop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MContraction_Loop_hpp_
#define Hadrons_MContraction_Loop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Noise loop propagator
 -----------------------------
 * loop_x = q_x * adj(eta_x)
 
 * options:
 - q = Result of inversion on noise source.
 - eta = noise source.

 */

/******************************************************************************
 *                         Loop                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class LoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoopPar,
                                    std::string, q,
                                    std::string, eta);
};

template <typename FImpl>
class TLoop: public Module<LoopPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoop(const std::string name);
    // destructor
    virtual ~TLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Loop, TLoop<FIMPL>, MContraction);

/******************************************************************************
 *                 TLoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoop<FImpl>::TLoop(const std::string name)
: Module<LoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().eta};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoop<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoop<FImpl>::execute(void)
{
    auto &loop = envGet(PropagatorField, getName());
    auto &q    = envGet(PropagatorField, par().q);
    auto &eta  = envGet(PropagatorField, par().eta);
    loop = q*adj(eta);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Loop_hpp_
