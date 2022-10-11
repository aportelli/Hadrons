/*
 * NPRUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Felix Erben <ferben@ed.ac.uk>
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

#ifndef Hadrons_MNPR_NPRUtils_hpp_
#define Hadrons_MNPR_NPRUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MNPR)

template <typename FImpl>
class NPRUtils
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    static void tensorProd(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b);
    static void tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b);
    // covariant derivative
    static void dslash(PropagatorField &in, const PropagatorField &out,
        const GaugeField &Umu);
    static void phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut);
    static void dot(ComplexField &pDotX, std::vector<Real> p);

};

// Tensor product of two PropagatorFields (Lattice Spin Colour Matrices in many FImpls)
template <typename FImpl>
void NPRUtils<FImpl>::tensorProd(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b)
{
    autoView(lret_v, lret, CpuWrite);
    autoView(a_v, a, CpuRead);
    autoView(b_v, b, CpuRead);

    thread_for( site, lret_v.size(), {
        vTComplex left;
        for(int si=0; si < Ns; ++si)
	{
        for(int sj=0; sj < Ns; ++sj)
	{
            for (int ci=0; ci < Nc; ++ci)
	    {
            for (int cj=0; cj < Nc; ++cj)
	    {
                left()()() = a_v[site]()(si,sj)(ci,cj);
                lret_v[site]()(si,sj)(ci,cj)=left()*b_v[site]();
            }}
        }}
    });
}

// Tensor product on a single site only
template <typename FImpl>
void NPRUtils<FImpl>::tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b)
{
    for(int si=0; si < Ns; ++si)
    {
    for(int sj=0; sj < Ns; ++sj)
    {
        for (int ci=0; ci < Nc; ++ci)
	{
        for (int cj=0; cj < Nc; ++cj)
	{
            const ComplexD val = TensorRemove(a()(si,sj)(ci,cj));
            lret()(si,sj)(ci,cj) = val * b();
        }}
    }}
}

// Computes gamma^mu D_mu for the given input field. Currently uses the
// symmetric derivative, though this could change in the future.
template <typename FImpl>
void NPRUtils<FImpl>::dslash(PropagatorField &out, const PropagatorField &in,
        const GaugeField &Umu)
{
    assert(&out != &in);
    out = Zero();
    PropagatorField tmp(Umu.Grid());
    typename FImpl::GaugeLinkField U(Umu.Grid());
    for (int mu = 0; mu < Nd; mu++)
    {
        // Overall formula:
        // tmp(x) = U_\mu(x) in(x + \hat{\mu}) - U_\mu^\dag(x - \hat{\mu}) in(x - \hat{\mu})
        U = peekLorentz(Umu, mu);
        tmp = FImpl::CovShiftForward(U, mu, in);
        tmp = tmp - FImpl::CovShiftBackward(U, mu, in);

        Gamma gamma_mu = Gamma::gmu[mu];
        out += gamma_mu * tmp;
    }
    out = 0.5 * out;
}


//// Compute phases for phasing propagators
// bilinearPhase = exp(-i (pIn - pOut) \cdot x)
template <typename FImpl>
void NPRUtils<FImpl>::phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut)
{
    bilinearPhase = Zero();
    ComplexField coordinate(bilinearPhase.Grid());
    Coordinate                  latt_size = GridDefaultLatt();
    for (int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;

        bilinearPhase += coordinate * (pIn[mu] - pOut[mu]);
    }
    Complex Ci = Complex(0.0, 1.0);
    bilinearPhase = exp(-Ci * bilinearPhase);
}


// pDotX = p \cdot x
template <typename FImpl>
void NPRUtils<FImpl>::dot(ComplexField &pDotX, std::vector<Real> p)
{
    ComplexField coordinate(pDotX.Grid());
    Coordinate                  latt_size = GridDefaultLatt();
    pDotX = Zero();
    for (int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;
        pDotX += coordinate * p[mu];
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
