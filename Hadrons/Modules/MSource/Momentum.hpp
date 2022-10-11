/*
 * Momentum.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef Hadrons_Momentum_hpp_
#define Hadrons_Momentum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/* 
Plane Wave source
-----------------
src_x = e^i2pi/L * p *position
*/

/******************************************************************************
 *                          Plane Wave source                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class MomentumPar: Serializable
{
    public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumPar,
                                    std::string, mom);
};

template <typename FImpl>
class TMomentum: public Module<MomentumPar>
{
    public:
    FERM_TYPE_ALIASES(FImpl,);
    public:
    // constructor
    TMomentum(const std::string name);
    // destructor
    virtual ~TMomentum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Momentum, TMomentum<FIMPL>, MSource);

/******************************************************************************
*                       TMomentum template implementation                     *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMomentum<FImpl>::TMomentum(const std::string name)
: Module<MomentumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMomentum<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    return in;
}

template <typename FImpl>
std::vector<std::string> TMomentum<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomentum<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmpLat(ComplexField, "C");
    envTmpLat(ComplexField, "xMu");
}

//execution//////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomentum<FImpl>::execute(void)
{
    LOG(Message) << "Generating planewave momentum source with momentum " << par().mom << std::endl;
    PropagatorField        &src = envGet(PropagatorField, getName());
    std::vector<Real>      p = strToVec<Real>(par().mom);;
    Coordinate                  latt_size = GridDefaultLatt();
    Complex                Ci(0.0,1.0);
    envGetTmp(ComplexField, C);
    envGetTmp(ComplexField, xMu);

    src = Zero();
    C = Zero();

    for(int mu = 0; mu < Nd; mu++){
        Real TwoPiL =  M_PI * 2.0 / latt_size[mu];
        LatticeCoordinate(xMu,mu);
        C = C + (TwoPiL * p[mu]) * xMu;
    }
    C = exp(C*Ci);
    src = src + C;
    LOG(Message) << "source created" << std::endl;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Momentum_hpp_
