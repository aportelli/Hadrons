/*
 * SeqGammaWall.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
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

#ifndef Hadrons_MSource_SeqGammaWall_hpp_
#define Hadrons_MSource_SeqGammaWall_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential wall source
 -----------------------------
 * Compute SitePropagator on each timeslice as
     wallsum[t] = \sum_x theta(t - x_3) q(x)
 * Then construct source as sum of the wall on each timeslice
     src(x) = sum_{t \in \set{t_A ... tB}} ( gamma * wallsum[t] * exp(i x.mom) )

 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqGammaWallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqGammaWallPar,
                                    std::string,    q,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    Gamma::Algebra, gamma,
                                    std::string,    mom);
};

template <typename FImpl>
class TSeqGammaWall: public Module<SeqGammaWallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqGammaWall(const std::string name);
    // destructor
    virtual ~TSeqGammaWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeSource(PropagatorField &src, const PropagatorField &q);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(SeqGammaWall, TSeqGammaWall<FIMPL>, MSource);
MODULE_REGISTER_TMP(ZSeqGammaWall, TSeqGammaWall<ZFIMPL>, MSource);

/******************************************************************************
 *                         TSeqGammaWall implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqGammaWall<FImpl>::TSeqGammaWall(const std::string name)
: Module<SeqGammaWallPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqGammaWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqGammaWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGammaWall<FImpl>::setup(void)
{
    if (envHasType(PropagatorField, par().q))
    {
        envCreateLat(PropagatorField, getName());
    }
    else if (envHasType(std::vector<PropagatorField>, par().q))
    {
        auto &q = envGet(std::vector<PropagatorField>, par().q);

        envCreate(std::vector<PropagatorField>, getName(), 1, q.size(),
                envGetGrid(PropagatorField));
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().q 
                          + "' has an incompatible type ("
                          + env().getObjectType(par().q)
                          + ")", env().getObjectAddress(par().q))
    }
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(PropagatorField, "wallTmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGammaWall<FImpl>::makeSource(PropagatorField &src,
                                  const PropagatorField &q)
{
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Gamma g(par().gamma);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }
    envGetTmp(PropagatorField, wallTmp);
    SlicedPropagator qSliced;
    sliceSum(q, qSliced, Tp);
    src = Zero();
    for(unsigned int loop_t = par().tA; loop_t <= par().tB; ++loop_t)
    {
        wallTmp = g * qSliced[loop_t];
        src = src + where((t == loop_t), ph*wallTmp, 0.*wallTmp);
    }
}

template <typename FImpl>
void TSeqGammaWall<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating " << par().gamma
                     << " sequential source(s) at t= " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating " << par().gamma
                     << " sequential source(s) for "
                     << par().tA << " <= t <= " << par().tB << std::endl;
    }

    if (envHasType(PropagatorField, par().q))
    {
        auto  &src = envGet(PropagatorField, getName()); 
        auto  &q   = envGet(PropagatorField, par().q);

        LOG(Message) << "Using propagator '" << par().q << "'" << std::endl;
        makeSource(src, q);
    }
    else
    {
        auto  &src = envGet(std::vector<PropagatorField>, getName()); 
        auto  &q   = envGet(std::vector<PropagatorField>, par().q);

        for (unsigned int i = 0; i < q.size(); ++i)
        {
            LOG(Message) << "Using element " << i << " of propagator vector '" 
                         << par().q << "'" << std::endl;
            makeSource(src[i], q[i]);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqGammaWall_hpp_
