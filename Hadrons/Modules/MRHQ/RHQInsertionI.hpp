/*
 * RHQInsertionI.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Alessandro Barone <barone1618@gmail.com>
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

#ifndef Hadrons_MRHQ_RHQInsertionI_hpp_
#define Hadrons_MRHQ_RHQInsertionI_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              RHQInsertionI                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQInsertionIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIPar,
                                    std::string,    q,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge);
};

// See https://arxiv.org/abs/1501.05373 equations 11 and 12 for the operators
// implemented in this module.
template <typename FImpl, typename GImpl>
class TRHQInsertionI: public Module<RHQInsertionIPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQInsertionI(const std::string name);
    // destructor
    virtual ~TRHQInsertionI(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionI, ARG(TRHQInsertionI<FIMPL, GIMPL>), MRHQ);
MODULE_REGISTER_TMP(RHQInsertionII, ARG(TRHQInsertionI<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                            RHQInsertionI                                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionI<FImpl, GImpl>::TRHQInsertionI(const std::string name)
: Module<RHQInsertionIPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionI<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionI<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionI<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());//, 1, env().getDim(Tp));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionI<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Applying Improvement term I with index " << par().index
                 << " and gamma5=" << par().gamma5 
                 << " to '" << par().q 
                 << std::endl;

    if (par().gamma5 != Gamma::Algebra::Gamma5 && par().gamma5 != Gamma::Algebra::Identity)
    {
        HADRONS_ERROR(Argument, "gamma5 must be either 'Gamma5' or 'Identity'."); 
    }
    Gamma g5(par().gamma5);

    if (par().index < 0 || par().index>3)
    {
        HADRONS_ERROR(Argument, "Index must be in {0, 1, 2, 3}."); 
    }
    const auto &index = par().index;

    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    const auto internal_gauge = peekLorentz(gaugefield, index);
    PropagatorField insertion = g5*(GImpl::CovShiftForward(internal_gauge,index,field) - GImpl::CovShiftBackward(internal_gauge,index,field));
    
    auto &out = envGet(PropagatorField, getName());
    out = insertion;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionI_hpp_
