/*
 * ChargedProp.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
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
#include <Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>
#include <Hadrons/Serialization.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                     TChargedProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TChargedProp::TChargedProp(const std::string name)
: Module<ChargedPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TChargedProp::getInput(void)
{
    return {par().source, par().emField};
}

std::vector<std::string> TChargedProp::getOutput(void)
{
    return {getName(), getName()+"_0", getName()+"_Q", getName()+"_Sun", getName()+"_Tad", getName()+"_projections"};
}

// setup ///////////////////////////////////////////////////////////////////////
void TChargedProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    phaseName_.clear();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        phaseName_.push_back("_shiftphase_" + std::to_string(mu));
    }
    GFSrcName_ = getName() + "_DinvSrc";
	prop0Name_ = getName() + "_0";
    propQName_ = getName() + "_Q";
    propSunName_ = getName() + "_Sun";
    propTadName_ = getName() + "_Tad";
    fftName_   = getName() + "_fft";

    freeMomPropDone_ = env().hasCreatedObject(freeMomPropName_);
    GFSrcDone_       = env().hasCreatedObject(GFSrcName_);
    phasesDone_      = env().hasCreatedObject(phaseName_[0]);
	prop0Done_		 = env().hasCreatedObject(prop0Name_);
    envCacheLat(ScalarField, freeMomPropName_);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        envCacheLat(ScalarField, phaseName_[mu]);
    }
    envCacheLat(ScalarField, GFSrcName_);
	envCacheLat(ScalarField, prop0Name_);
    envCreateLat(ScalarField, getName());
    envCreateLat(ScalarField, propQName_);
    envCreateLat(ScalarField, propSunName_);
    envCreateLat(ScalarField, propTadName_);
    envCreate(HadronsSerializable, getName()+"_projections", 1, 0);
    envTmpLat(ScalarField, "buf");
    envTmpLat(ScalarField, "result");
    envTmpLat(ScalarField, "Amu");
    envCache(FFT, fftName_, 1, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
void TChargedProp::execute(void)
{
    // CACHING ANALYTIC EXPRESSIONS
    makeCaches();

    // PROPAGATOR CALCULATION
    LOG(Message) << "Computing charged scalar propagator"
                 << " (mass= " << par().mass
                 << ", charge= " << par().charge << ")..." << std::endl;
    
    auto   &prop    = envGet(ScalarField, getName());
	auto   &prop0   = envGet(ScalarField, prop0Name_);
	auto   &propQ   = envGet(ScalarField, propQName_);
	auto   &propSun = envGet(ScalarField, propSunName_);
	auto   &propTad = envGet(ScalarField, propTadName_);
    auto   &GFSrc   = envGet(ScalarField, GFSrcName_);
    auto   &G       = envGet(ScalarField, freeMomPropName_);
    auto   &fft     = envGet(FFT, fftName_);
    double q        = par().charge;
    envGetTmp(ScalarField, buf);

    // -G*momD1*G*F*Src (momD1 = F*D1*Finv)
    propQ = GFSrc;
    momD1(propQ, fft);
    propQ = -G*propQ;
    propSun = -propQ;
    fft.FFT_dim(propQ, propQ, env().getNd()-1, FFT::backward);

    // G*momD1*G*momD1*G*F*Src (here buf = G*momD1*G*F*Src)
    momD1(propSun, fft);
    propSun = G*propSun;
    fft.FFT_dim(propSun, propSun, env().getNd()-1, FFT::backward);

    // -G*momD2*G*F*Src (momD2 = F*D2*Finv)
    propTad = GFSrc;
    momD2(propTad, fft);
    propTad = -G*propTad;
    fft.FFT_dim(propTad, propTad, env().getNd()-1, FFT::backward);
    
    // full charged scalar propagator
    fft.FFT_dim(buf, GFSrc, env().getNd()-1, FFT::backward);
    prop = buf + q*propQ + q*q*propSun + q*q*propTad;

    // output selected momenta
    Result result;
    TComplex            site;
    std::vector<int>    siteCoor;

    LOG(Message) << "Saving momentum-projected propagator to file '"
                 << resultFilename(par().output) << "' and object '"
                 << getName()+"_projections" << "'..." << std::endl;
    result.projection.resize(par().outputMom.size());
    result.lattice_size = env().getGrid()->FullDimensions().toVector();
    result.mass = par().mass;
    result.charge = q;
    auto nd = env().getNd();
    auto nt = env().getGrid()->FullDimensions()[nd-1];
    siteCoor.resize(nd);
    for (unsigned int i_p = 0; i_p < par().outputMom.size(); ++i_p)
    {
        result.projection[i_p].momentum = strToVec<int>(par().outputMom[i_p]);

        LOG(Message) << "Calculating (" << par().outputMom[i_p]
                     << ") momentum projection" << std::endl;

        result.projection[i_p].corr_0.resize(nt);
        result.projection[i_p].corr.resize(nt);
        result.projection[i_p].corr_Q.resize(nt);
        result.projection[i_p].corr_Sun.resize(nt);
        result.projection[i_p].corr_Tad.resize(nt);

        for (unsigned int j = 0; j < nd-1; ++j)
        {
            siteCoor[j] = result.projection[i_p].momentum[j];
        }

        for (unsigned int t = 0; t < result.projection[i_p].corr.size(); ++t)
        {
            siteCoor[nd-1] = t;
            peekSite(site, prop, siteCoor);
            result.projection[i_p].corr[t]=TensorRemove(site);
            peekSite(site, buf, siteCoor);
            result.projection[i_p].corr_0[t]=TensorRemove(site);
            peekSite(site, propQ, siteCoor);
            result.projection[i_p].corr_Q[t]=TensorRemove(site);
            peekSite(site, propSun, siteCoor);
            result.projection[i_p].corr_Sun[t]=TensorRemove(site);
            peekSite(site, propTad, siteCoor);
            result.projection[i_p].corr_Tad[t]=TensorRemove(site);
        }
    }
    saveResult(par().output, "prop", result);
    envGet(HadronsSerializable, getName()+"_projections") = result;

    std::vector<int> mask(nd,1);
    mask[nd-1] = 0;
    fft.FFT_dim_mask(prop, prop, mask, FFT::backward);
    fft.FFT_dim_mask(propQ, propQ, mask, FFT::backward);
    fft.FFT_dim_mask(propSun, propSun, mask, FFT::backward);
    fft.FFT_dim_mask(propTad, propTad, mask, FFT::backward);
}

void TChargedProp::makeCaches(void)
{
    auto &freeMomProp = envGet(ScalarField, freeMomPropName_);
    auto &GFSrc       = envGet(ScalarField, GFSrcName_);
	auto &prop0		  = envGet(ScalarField, prop0Name_);
    auto &fft         = envGet(FFT, fftName_);

    if (!freeMomPropDone_)
    {
        LOG(Message) << "Caching momentum-space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        SIMPL::MomentumSpacePropagator(freeMomProp, par().mass);
    }
    if (!GFSrcDone_)
    {   
        auto &source = envGet(ScalarField, par().source);
        
        LOG(Message) << "Caching G*F*src..." << std::endl;
        fft.FFT_all_dim(GFSrc, source, FFT::forward);
        GFSrc = freeMomProp*GFSrc;
    }
	if (!prop0Done_)
	{
		LOG(Message) << "Caching position-space free scalar propagator..."
                     << std::endl;
		fft.FFT_all_dim(prop0, GFSrc, FFT::backward);
	}
    if (!phasesDone_)
    {
        auto l = env().getGrid()->FullDimensions();
        Complex          ci(0.0,1.0);
        
        LOG(Message) << "Caching shift phases..." << std::endl;
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            Real twoPiL = M_PI*2./l[mu];
            auto &phmu  = envGet(ScalarField, phaseName_[mu]);
            
            LatticeCoordinate(phmu, mu);
            phmu = exp(ci*twoPiL*phmu);
            phase_.push_back(&phmu);
        }
    }
    else
    {
        phase_.clear();
        for (unsigned int mu = 0; mu < env().getNd(); ++mu)
        {
            phase_.push_back(env().getObject<ScalarField>(phaseName_[mu]));
        }
    }
}

void TChargedProp::momD1(ScalarField &s, FFT &fft)
{
    auto        &A = envGet(EmField, par().emField);
    Complex     ci(0.0,1.0);

    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, result);
    envGetTmp(ScalarField, Amu);

    result = Zero();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = Amu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result - ci*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = Amu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + ci*adj(*phase_[mu])*buf;
    }

    s = result;
}

void TChargedProp::momD2(ScalarField &s, FFT &fft)
{
    auto &A = envGet(EmField, par().emField);

    envGetTmp(ScalarField, buf);
    envGetTmp(ScalarField, result);
    envGetTmp(ScalarField, Amu);

    result = Zero();
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);
        buf = (*phase_[mu])*s;
        fft.FFT_all_dim(buf, buf, FFT::backward);
        buf = Amu*Amu*buf;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*buf;
    }
    fft.FFT_all_dim(s, s, FFT::backward);
    for (unsigned int mu = 0; mu < env().getNd(); ++mu)
    {
        Amu = peekLorentz(A, mu);        
        buf = Amu*Amu*s;
        fft.FFT_all_dim(buf, buf, FFT::forward);
        result = result + .5*adj(*phase_[mu])*buf;
    }

    s = result;
}
