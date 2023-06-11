/*
 * SeqGammaRegion.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
               Regions do not wrap around the edge of the lattice
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
    template <typename T>
    inline std::vector<T> ErrorCheck(const std::string &Value, unsigned int Nd, const std::string &FieldName)
    {
        std::vector<T> vi = strToVec<T>(Value);
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
    using LatSInt = Lattice<iScalar<vInteger>>;
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
    const std::string coorName_;
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
, coorName_ (name + "_coor")
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
    if (env().getNd() !=4 )
    {
        // TODO: make sure the where clause works on the correct number of dimensions
        HADRONS_ERROR(Size,"Expected 4 dimensions, but have "+std::to_string(env().getNd()));
    }
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
    LowerLeft  = SeqGammaRegionHelper::ErrorCheck<int>(par().LowerLeft,  env().getNd(), "LowerLeft");
    RegionSize = SeqGammaRegionHelper::ErrorCheck<int>(par().RegionSize, env().getNd(), "RegionSize");
    Momentum   = SeqGammaRegionHelper::ErrorCheck<int>(par().mom       , env().getNd(), "mom");
    Gamma g(par().gamma);
    // Create temporaries
    envCacheLat(LatticeComplex, momphName_);
    envCache(std::vector<LatSInt>, coorName_, 1, env().getNd(), envGetGrid(LatticeComplex)); // coords for where clause
    envTmpLat(LatticeComplex, "pcoor"); // This temporary is used for momentum
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqGammaRegion<FImpl>::makeSource(PropagatorField &src, const PropagatorField &q)
{
    auto &ph   = envGet(LatticeComplex, momphName_);
    auto &coor = envGet(std::vector<LatSInt>, coorName_);
    Gamma g(par().gamma);

    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, pcoor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(pcoor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*pcoor;
            LatticeCoordinate(coor[mu], mu);
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        hasPhase_ = true;
    }
    src = where(    (coor[0] >= LowerLeft[0]) and (coor[0] < LowerLeft[0] + RegionSize[0])
                and (coor[1] >= LowerLeft[1]) and (coor[1] < LowerLeft[1] + RegionSize[1])
                and (coor[2] >= LowerLeft[2]) and (coor[2] < LowerLeft[2] + RegionSize[2])
                and (coor[3] >= LowerLeft[3]) and (coor[3] < LowerLeft[3] + RegionSize[3]),
                ph*(g*q), 0.*q);
}

template <typename FImpl>
void TSeqGammaRegion<FImpl>::execute(void)
{
    LOG(Message) << "Generating " << par().gamma
                 << " sequential source region: LowerLeft=" << par().LowerLeft
                 << ", RegionSize=" << par().RegionSize
                 << ", mom=" << par().mom
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
