/*
 * RHQInsertionIII.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Alessandro Barone <ab1n19@soton.ac.uk>
 * Author: Alessandro Barone <barone1618@gmail.com>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Matthew Black <matthewkblack@protonmail.com>
 * Author: Matthew Black <mbr-phys@protonmail.com>
 * Author: RChrHill <75032435+RChrHill@users.noreply.github.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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
BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQInsertionIIIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIIIPar,
                                    std::string,    q,
                                    std::string,    index1,
                                    std::string,    index2,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge);
};

// See https://arxiv.org/abs/1501.05373 equation 13 for the 
// operator implemented in this module.
// To convert to the charge conjugation "invariant" basis:
// sJ = +1 if gamma5 = Identity; sJ = -1 if gamma5 = Gamma5
// if either index1 or index2 empty:
//     RHQIII_mu -> sJ * RHQIII_mu - delta_{mu,i} * RHQI(D_i)
// if both index1 and index2 used:
//     RHQIII_{mu,nu} -> sJ * RHQIII_{mu,nu} - 2 * ( delta_{nu,i} * RHQI(gamma_mu,D_i) - delta_{mu,i} * RHQI(gamma_nu,D_i) )
// also see ale-barone/feature/RHQImprTensor
template <typename FImpl, typename GImpl>
class TRHQInsertionIII: public Module<RHQInsertionIIIPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
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
    envTmpLat(ColourMatrixField, "gauge_x");
    envTmpLat(ColourMatrixField, "gauge_y");
    envTmpLat(ColourMatrixField, "gauge_z");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIII<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Applying Improvement term III with (index1, index2)=" 
                 << "(" << par().index1 << ", " << par().index2 << ")"
                 << " and gamma5=" << par().gamma5 
                 << " to '" << par().q 
                 << std::endl;
    
    if (par().gamma5 != Gamma::Algebra::Gamma5 && par().gamma5 != Gamma::Algebra::Identity)
    {
        HADRONS_ERROR(Argument, "gamma5 must be either 'Gamma5' or 'Identity'."); 
    }
    Gamma g5(par().gamma5);

    int index1;
    int index2;
    bool hasIndex1 = false;
    bool hasIndex2 = false;
    if (!par().index1.empty())
    {
        index1 = std::stoi(par().index1);
        hasIndex1 = true;
    }
    if (!par().index2.empty())
    {
        index2 = std::stoi(par().index2);
        hasIndex2 = true;
    }
    if (!hasIndex1 && !hasIndex2)
    {
        HADRONS_ERROR(Argument, "index1 and index2 cannot be both empty."); 
    }
    
    auto &field = envGet(PropagatorField, par().q);
    const auto &gaugefield = envGet(GaugeField, par().gauge);
    envGetTmp(ColourMatrixField, gauge_x);
    envGetTmp(ColourMatrixField, gauge_y);
    envGetTmp(ColourMatrixField, gauge_z);
    gauge_x = peekLorentz(gaugefield, 0);
    gauge_y = peekLorentz(gaugefield, 1);
    gauge_z = peekLorentz(gaugefield, 2);

    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);

    const PropagatorField Dx = GImpl::CovShiftForward(gauge_x,0,field) - GImpl::CovShiftBackward(gauge_x,0,field);
    const PropagatorField Dy = GImpl::CovShiftForward(gauge_y,1,field) - GImpl::CovShiftBackward(gauge_y,1,field);
    const PropagatorField Dz = GImpl::CovShiftForward(gauge_z,2,field) - GImpl::CovShiftBackward(gauge_z,2,field);

    Gamma::Algebra gi;
    if (hasIndex1 && hasIndex2)
    {
        switch(index1){
            case 0:
                switch(index2)
                    {
                        case 0:
                            gi = 0.*Gamma::Algebra::Identity;
                            break;
                        case 1:
                            gi = Gamma::Algebra::SigmaXY;
                            break;
                        case 2:
                            gi = Gamma::Algebra::SigmaXZ;
                            break;
                        case 3:
                            gi = Gamma::Algebra::SigmaXT;
                            break;
                        default:
                            HADRONS_ERROR(Argument, "index2 must be in {0, 1, 2, 3}."); 
                    }
                break;
            case 1:
                switch(index2)
                    {
                        case 0:
                            gi = Gamma::Algebra::MinusSigmaXY;
                            break;
                        case 1:
                            gi = 0.*Gamma::Algebra::Identity;
                            break;
                        case 2:
                            gi = Gamma::Algebra::SigmaYZ;
                            break;
                        case 3:
                            gi = Gamma::Algebra::SigmaYT;
                            break;
                        default:
                            HADRONS_ERROR(Argument, "index2 must be in {0, 1, 2, 3}."); 
                    }
                break;
            case 2:
                switch(index2)
                    {
                        case 0:
                            gi = Gamma::Algebra::MinusSigmaXZ;
                            break;
                        case 1:
                            gi = Gamma::Algebra::MinusSigmaYZ;
                            break;
                        case 2:
                            gi = 0.*Gamma::Algebra::Identity;
                            break;
                        case 3:
                            gi = Gamma::Algebra::SigmaZT;
                            break;
                        default:
                            HADRONS_ERROR(Argument, "index2 must be in {0, 1, 2, 3}."); 
                    }
                break;
            case 3:
                switch(index2)
                    {
                        case 0:
                            gi = Gamma::Algebra::MinusSigmaXT;
                            break;
                        case 1:
                            gi = Gamma::Algebra::MinusSigmaYT;
                            break;
                        case 2:
                            gi = Gamma::Algebra::MinusSigmaZT;
                            break;
                        case 3:
                            gi = 0.*Gamma::Algebra::Identity;
                            break;
                        default:
                            HADRONS_ERROR(Argument, "index2 must be in {0, 1, 2, 3}."); 
                    } 
                break;
            default:
                HADRONS_ERROR(Argument, "index1 must be in {0, 1, 2, 3}."); 
        }
    }
    else if (hasIndex1 && !hasIndex2)
    {
        switch(index1)
            {
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
                    HADRONS_ERROR(Argument, "index1 must be in {0, 1, 2, 3}."); 
            }
    }
    else if (!hasIndex1 && hasIndex2)
    {
        switch(index2)
            {
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
                    HADRONS_ERROR(Argument, "index2 must be in {0, 1, 2, 3}."); 
            }
    }
        
    auto &out = envGet(PropagatorField, getName());
    PropagatorField insertion =
        gi*g5*gx * Dx
      + gi*g5*gy * Dy
      + gi*g5*gz * Dz;
    
    out = insertion;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionIII_hpp_
