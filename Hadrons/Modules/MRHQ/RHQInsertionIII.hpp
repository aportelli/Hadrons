/*
 * RHQInsertionIII.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MRHQ_RHQInsertionIII_hpp_
#define Hadrons_MRHQ_RHQInsertionIII_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                             RHQInsertionIII                                *
 ******************************************************************************/
GRID_SERIALIZABLE_ENUM(OpIIIFlag, undef, Chroma, 0, LeftRight, 1);

BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQInsertionIIIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIIIPar,
                                    std::string,    q,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge,
                                    OpIIIFlag,      flag);
};

// See https://arxiv.org/abs/1501.05373 equation 13 for the (Chroma convention) 
// operator implemented in this module.
template <typename FImpl, typename GImpl>
class TRHQInsertionIII: public Module<RHQInsertionIIIPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQInsertionIII(const std::string name);
    // destructor
    virtual ~TRHQInsertionIII(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionIII, ARG(TRHQInsertionIII<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                           RHQInsertionIII                                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionIII<FImpl, GImpl>::TRHQInsertionIII(const std::string name)
: Module<RHQInsertionIIIPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIII<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIII<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIII<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIII<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Applying Improvement term III with index " << par().index
                 << " and gamma5=" << par().gamma5 
                 << " to '" << par().q 
                 << "' with flag '" << par().flag << "'"
                 << std::endl;
    
    if (par().gamma5 != Gamma::Algebra::Gamma5 && par().gamma5 != Gamma::Algebra::Identity)
    {
        HADRONS_ERROR(Argument, "gamma5 must be either 'Gamma5' or 'Identity'."); 
    }
    Gamma g5(par().gamma5);
    
    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    const auto gauge_x = peekLorentz(gaugefield, 0);
    const auto gauge_y = peekLorentz(gaugefield, 1);
    const auto gauge_z = peekLorentz(gaugefield, 2);

    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);

    const PropagatorField Dx = GImpl::CovShiftForward(gauge_x,0,field) - GImpl::CovShiftBackward(gauge_x,0,field);
    const PropagatorField Dy = GImpl::CovShiftForward(gauge_y,1,field) - GImpl::CovShiftBackward(gauge_y,1,field);
    const PropagatorField Dz = GImpl::CovShiftForward(gauge_z,2,field) - GImpl::CovShiftBackward(gauge_z,2,field);

    Gamma::Algebra gi; 
    switch(par().index){
        case 0:
            gi = Gamma::Algebra::GammaX;
            break;
        case 1:
            gi = Gamma::Algebra::GammaY;
            break;
        case 2:
            gi = Gamma::Algebra::GammaZ;
            break;
        case 3:
            gi = Gamma::Algebra::GammaT;
            break;
        default:
            HADRONS_ERROR(Argument, "Index must be in {0, 1, 2, 3}."); 
    }
    
    // "Flag" flips between the conventions used in Chroma and a reformulation
    // using Left and Right derivatives
    auto &out = envGet(PropagatorField, getName());
    if (par().flag == OpIIIFlag::Chroma)
    {     
        PropagatorField insertion =
            gi*g5*gx * Dx
          + gi*g5*gy * Dy
          + gi*g5*gz * Dz;
        
        out = insertion;
    }
    else if (par().flag == OpIIIFlag::LeftRight)
    {        
        PropagatorField insertion = 
            gi*gx*g5 * Dx - gx*gi*g5 * Dx
          + gi*gy*g5 * Dy - gy*gi*g5 * Dy
          + gi*gz*g5 * Dz - gz*gi*g5 * Dz;

        out = 0.5*insertion;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionIII_hpp_
