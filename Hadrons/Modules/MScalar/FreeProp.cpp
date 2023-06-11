/*
 * FreeProp.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
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
#include <Hadrons/Modules/MScalar/FreeProp.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>
#include <Hadrons/Serialization.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                        TFreeProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TFreeProp::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TFreeProp::getInput(void)
{
    return {par().source};
}

std::vector<std::string> TFreeProp::getOutput(void)
{
    return {getName(), getName()+"_sliceSum"};
}

// setup ///////////////////////////////////////////////////////////////////////
void TFreeProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    
    freePropDone_ = env().hasCreatedObject(freeMomPropName_);
    envCacheLat(ScalarField, freeMomPropName_);
    envCreateLat(ScalarField, getName());
    envCreate(HadronsSerializable, getName()+"_sliceSum", 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
void TFreeProp::execute(void)
{
    auto &freeMomProp = envGet(ScalarField, freeMomPropName_);
    auto &prop        = envGet(ScalarField, getName());
    auto &source      = envGet(ScalarField, par().source);

    if (!freePropDone_)
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        SIMPL::MomentumSpacePropagator(freeMomProp, par().mass);
    }
    LOG(Message) << "Computing free scalar propagator..." << std::endl;
    SIMPL::FreePropagator(source, prop, freeMomProp);
    
    std::vector<TComplex> buf;
    std::vector<Complex>  result;
    
    sliceSum(prop, buf, Tp);
    result.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result[t] = TensorRemove(buf[t]);
    }
    envGet(HadronsSerializable, getName()+"_sliceSum") = result;
    saveResult(par().output, "freeprop", result);
}
