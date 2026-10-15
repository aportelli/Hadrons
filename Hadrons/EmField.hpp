/*
 * EmField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_EmField_hpp_
#define Hadrons_EmField_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Class to generate stochatic EM fields                         *
 ******************************************************************************/
// gauge enum
GRID_SERIALIZABLE_ENUM(QedGauge, undef, feynman, 0, coulomb, 1);

// main class
template <typename VType>
class TEmFieldGenerator
{
public:
    // type aliases
    typedef iVector<iScalar<iScalar<VType>>, Nd>  vSiteField;
    typedef iScalar<iScalar<iScalar<VType>>>      vSiteScalar;
    typedef typename vSiteScalar::scalar_object   SiteScalar;
    typedef typename SiteScalar::scalar_type      ComplexType;
    typedef Lattice<vSiteField>                   GaugeField;
    typedef Lattice<vSiteScalar>                  ScalarField;
    typedef std::function<void(GaugeField &)>     TransformFn;
public:
    // constructor
    TEmFieldGenerator(GridBase *g);
    // static functions to construct useful lattices
    static void makeSpatialNorm(LatticeInteger &out);
    static void makeKHat(std::vector<ScalarField> &out);
    static void makeKHatSquared(ScalarField &out);
    // spatial transverse projection in momentum space (Coulomb gauge projection)
    static void transverseProjectSpatial(GaugeField &out);
    // get gauge transformation out of enum
    static TransformFn getGaugeTranform(const QedGauge gauge);
    // field generation using given weights
    void operator()(GaugeField &out, GridParallelRNG &rng, const ScalarField &feynmanProp, 
                    TransformFn momSpaceTransform = nullptr);
    // QED_L weights
    void makeFeynmanPropQedL(ScalarField &feynmanProp, std::vector<double> improvement = {});
    // QED_TL weights
    void makeFeynmanPropQedTL(ScalarField &feynmanProp);
    // QED_Zeta weights
    void makeFeynmanPropQedZeta(ScalarField &feynmanProp, const double zeta);
    // QED_M weights
    void makeFeynmanPropQedM(ScalarField &feynmanProp, const double mass);
private:
    GridBase       *g_;
    unsigned int   nd_;
    LatticeInteger spNrm_;
    ScalarField    latOne_, kHatSquared_;
    SiteScalar     one_, z_;
    Coordinate     zm_;
};

// default alias
typedef TEmFieldGenerator<vComplex> EmFieldGenerator;

/******************************************************************************
 *                 TEmFieldGenerator implementation                           *
 ******************************************************************************/
// constructor ////////////////////////////////////////////////////////////////
template <typename VType>
TEmFieldGenerator<VType>::TEmFieldGenerator(GridBase *g)
: g_(g), nd_(g->Nd()), spNrm_(g), latOne_(g), kHatSquared_(g), zm_(g->Nd(), 0)
{
    one_    = ComplexType(1., 0.);
    z_      = ComplexType(0., 0.);
    latOne_ = ComplexType(1., 0.);
}

// static functions to construct useful lattices //////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::makeSpatialNorm(LatticeInteger &out)
{
    GridBase *g = out.Grid();
    LatticeInteger coor(g);
    auto l = g->FullDimensions();

    out = Zero();
    for(int mu = 0; mu < g->Nd() - 1; mu++)
    {
      LatticeCoordinate(coor, mu);
      coor  = where(coor < Integer(l[mu]/2), coor, coor - Integer(l[mu]));
      out = out + coor*coor;
    }
}

template <typename VType>
void TEmFieldGenerator<VType>::makeKHat(std::vector<ScalarField> &out)
{
    GridBase *g = out.front().Grid();
    const unsigned int nd = g->Nd();
    auto l = g->FullDimensions();
    ComplexType ci(0., 1.);

    out.resize(nd, g);
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      Real piL = M_PI/l[mu];

      LatticeCoordinate(out[mu], mu);
      out[mu] = exp(piL*ci*out[mu])*2.*sin(piL*out[mu]);
    }
}

template <typename VType>
void TEmFieldGenerator<VType>::makeKHatSquared(ScalarField &out)
{
    GridBase *g = out.Grid();
    const unsigned int nd = g->Nd();
    std::vector<ScalarField> khat(nd, g);
    
    out = Zero();
    makeKHat(khat);
    for(int mu = 0; mu < nd; mu++)
    {
      out = out + khat[mu]*conjugate(khat[mu]);
    }
}

// spatial transverse projection in momentum space (Coulomb gauge projection)
template <typename VType>
void TEmFieldGenerator<VType>::transverseProjectSpatial(GaugeField &out)
{
    GridBase                 *g = out.Grid();
    const unsigned int       nd = g->Nd();
    ScalarField              invKHat(g), cst(g), spdiv(g);
    LatticeInteger           spNrm(g);
    std::vector<ScalarField> khat(nd, g), a(nd, g), aProj(nd, g);

    invKHat = Zero();
    makeSpatialNorm(spNrm);
    makeKHat(khat);
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      a[mu] = peekLorentz(out, mu);
      if (mu < nd - 1)
      {
        invKHat += khat[mu]*conjugate(khat[mu]);
      }
    }
    cst     = ComplexType(1., 0.);
    invKHat = where(spNrm == Integer(0), cst, invKHat);
    invKHat = cst/invKHat;
    cst     = Zero();
    invKHat = where(spNrm == Integer(0), cst, invKHat);
    spdiv   = Zero();
    for (unsigned int nu = 0; nu < nd - 1; ++nu)
    {
      spdiv += conjugate(khat[nu])*a[nu];
    }
    spdiv *= invKHat;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      aProj[mu] = a[mu] - khat[mu]*spdiv;
      pokeLorentz(out, aProj[mu], mu);
    }
}

// get gauge transformation out of enum ///////////////////////////////////////
template <typename VType>
typename TEmFieldGenerator<VType>::TransformFn 
TEmFieldGenerator<VType>::getGaugeTranform(const QedGauge gauge)
{
    switch (gauge)
    {
    case QedGauge::feynman:
      return nullptr;
      break;
    case QedGauge::coulomb:
      return TransformFn(&TEmFieldGenerator<VType>::transverseProjectSpatial);
      break;
    default:
      HADRONS_ERROR(Definition, "invalid gauge")
      break;
    }
}

// field generation using given weights ///////////////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::operator()(GaugeField &out, GridParallelRNG &rng, 
                                          const ScalarField &feynmanProp, TransformFn momSpaceTransform)
{
    ScalarField        r(g_), sqrtW(g_);
    GaugeField         aTilde(g_);
    FFT                fft(dynamic_cast<GridCartesian *>(g_));
    double             vol = g_->gSites();

    sqrtW = sqrt(vol)*sqrt(feynmanProp);
    for(unsigned int mu = 0; mu < nd_; mu++)
    {
      gaussian(rng, r);
      r = sqrtW*r;
      pokeLorentz(aTilde, r, mu);
    }
    if (momSpaceTransform != nullptr)
    {
        momSpaceTransform(aTilde);
    }
    fft.FFT_all_dim(out, aTilde, FFT::backward);
    out = real(out);
}

// QED_L weights //////////////////////////////////////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::makeFeynmanPropQedL(ScalarField &feynmanProp, 
                                                   std::vector<double> improvement)
{
    makeKHatSquared(feynmanProp);
    pokeSite(one_, feynmanProp, zm_);
    feynmanProp = latOne_/feynmanProp;
    pokeSite(z_, feynmanProp, zm_);
    makeSpatialNorm(spNrm_);
    feynmanProp = where(spNrm_ == Integer(0), 0.*feynmanProp, feynmanProp);
    for(int i = 0; i < improvement.size(); i++)
    {
      Real f = improvement[i] + 1;
      feynmanProp = where(spNrm_ == Integer(i + 1), f*feynmanProp, feynmanProp);
    }
}

// QED_TL weights /////////////////////////////////////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::makeFeynmanPropQedTL(ScalarField &feynmanProp)
{
    makeKHatSquared(feynmanProp);
    pokeSite(one_, feynmanProp, zm_);
    feynmanProp = latOne_/feynmanProp;
    pokeSite(z_, feynmanProp, zm_);
}

// QED_Zeta weights ///////////////////////////////////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::makeFeynmanPropQedZeta(ScalarField &feynmanProp, const double zeta)
{
    auto l = g_->FullDimensions()[0];
    ComplexType zl(zeta*l, 0.), zm = 1./(zl*zl);

    makeKHatSquared(feynmanProp);
    makeSpatialNorm(spNrm_);
    feynmanProp = where(spNrm_ == Integer(0), feynmanProp + zm*latOne_, feynmanProp);
    feynmanProp = latOne_/feynmanProp;
}

// QED_M weights //////////////////////////////////////////////////////////////
template <typename VType>
void TEmFieldGenerator<VType>::makeFeynmanPropQedM(ScalarField &feynmanProp, const double m)
{
    makeKHatSquared(feynmanProp);
    feynmanProp += m*m*latOne_;
    feynmanProp = latOne_/feynmanProp;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_EmField_hpp_
