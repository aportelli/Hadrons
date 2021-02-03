/*
 * SeqGamma.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
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

#ifndef Hadrons_MSource_SeqGamma_Region_hpp_
#define Hadrons_MSource_SeqGamma_Region_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source region
 -----------------------------
 * src_x = q_x * theta(x \in region) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - LowerLeft: bottom left corner of region (Nd x integer, eg "0 0 0 0")
 - RegionSize: size of region (Nd x integer, eg "1 1 1 1")
 - gamma: gamma product to insert (Gamma::Algebra)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

// This doesn't need to be specialised by template type
namespace SeqGammaRegionHelper
{
    inline std::vector<int> ErrorCheck(const std::string &Value, unsigned int Nd, const std::string &FieldName)
    {
        std::vector<int> vi = strToVec<int>(Value);
        if (vi.size() != Nd)
        {
            HADRONS_ERROR(Size,FieldName+" \""+Value+"\" should have "+std::to_string(Nd)+" dimensions");
        }
        return vi;
    }
};

class SeqGammaRegionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqGammaRegionPar,
                                    std::string,    q,
                                    std::string,    LowerLeft,
                                    std::string,    RegionSize,
                                    Gamma::Algebra, gamma,
                                    std::string,    mom);
};

template <typename FImpl>
class TSeqGammaRegion: public Module<SeqGammaRegionPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqGammaRegion(const std::string name);
    // destructor
    virtual ~TSeqGammaRegion(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

    void makeSource(PropagatorField &src, const PropagatorField &q);

    Coordinate        LowerLeft;
    Coordinate        RegionSize;
    Coordinate        Momentum;
    bool              bPropVec;
    bool              hasPhase_{false};
    const std::string momphName_;
    const std::string qMaskName_;
};

MODULE_REGISTER_TMP(SeqGammaRegion,  TSeqGammaRegion<FIMPL>,  MSource);
MODULE_REGISTER_TMP(ZSeqGammaRegion, TSeqGammaRegion<ZFIMPL>, MSource);

/******************************************************************************
 *                         TSeqGamma implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqGammaRegion<FImpl>::TSeqGammaRegion(const std::string name)
: Module<SeqGammaRegionPar>(name)
, momphName_ (name + "_momph")
, qMaskName_ (name + "_qMask")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqGammaRegion<FImpl>::getInput(void)
{
    return {par().q};
}

template <typename FImpl>
std::vector<std::string> TSeqGammaRegion<FImpl>::getOutput(void)
{
  return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGammaRegion<FImpl>::setup(void)
{
    if (envHasType(PropagatorField, par().q))
    {
        bPropVec = false;
        envCreateLat(PropagatorField, getName());
    }
    else if (envHasType(std::vector<PropagatorField>, par().q))
    {
        bPropVec = true;
        auto &q = envGet(std::vector<PropagatorField>, par().q);
        envCreate(std::vector<PropagatorField>, getName(), 1, q.size(), envGetGrid(PropagatorField));
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().q
                          + "' has an incompatible type ("
                          + env().getObjectType(par().q)
                          + ")", env().getObjectAddress(par().q))
    }
    // Validate parameters - fail early
    LowerLeft  = SeqGammaRegionHelper::ErrorCheck(par().LowerLeft,  env().getNd(), "LowerLeft");
    RegionSize = SeqGammaRegionHelper::ErrorCheck(par().RegionSize, env().getNd(), "RegionSize");
    Momentum   = SeqGammaRegionHelper::ErrorCheck(par().mom       , env().getNd(), "mom");
    Gamma g(par().gamma);
    // Create temporaries
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(PropagatorField, qMaskName_);
    envTmpLat(LatticeComplex, "coor");
    // The region never changes, so only need to zero outside the region once
    auto &qMask = envGet(PropagatorField, qMaskName_);
    qMask = Zero();
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGammaRegion<FImpl>::makeSource(PropagatorField &src, const PropagatorField &q)
{
    auto &ph    = envGet(LatticeComplex, momphName_);
    auto &qMask = envGet(PropagatorField, qMaskName_);
    Gamma g(par().gamma);

    // Get the specified region of the propagator. The rest is already zero
    localCopyRegion(q, qMask, LowerLeft, LowerLeft, RegionSize);

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
        hasPhase_ = true;
    }
    src = ph*(g*qMask);
}

template <typename FImpl>
void TSeqGammaRegion<FImpl>::execute(void)
{
    LOG(Message) << "Generating " << par().gamma
                 << " sequential source region. LowerLeft=" << par().LowerLeft
                 << " . RegionSize=" << par().RegionSize
                 << " . mom=" << par().mom
                 << std::endl;

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

#endif // Hadrons_MSource_SeqGamma_Region_hpp_
