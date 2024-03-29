/*
 * Wall.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

#ifndef Hadrons_MSource_WallSource_hpp_
#define Hadrons_MSource_WallSource_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Wall source
 -----------------------------
 * src_x = delta(x_3 - tW) * exp(i x.mom)
 
 * options:
 - tW: source timeslice (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class WallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WallPar,
                                    unsigned int, tW,
                                    std::string, mom);
};

template <typename FImpl>
class TWall: public Module<WallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWall(const std::string name);
    // destructor
    virtual ~TWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(Wall, TWall<FIMPL>, MSource);

/******************************************************************************
 *                 TWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWall<FImpl>::TWall(const std::string name)
: Module<WallPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::execute(void)
{    
    LOG(Message) << "Generating wall source at t = " << par().tW 
                 << " with momentum " << par().mom << std::endl;
    
    auto  &src = envGet(PropagatorField, getName());
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for(unsigned int mu = 0; mu < p.size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }

    src = 1.;
    src = where((t == par().tW), src*ph, 0.*src);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_WallSource_hpp_
