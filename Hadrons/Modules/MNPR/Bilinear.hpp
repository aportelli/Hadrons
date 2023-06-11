/*
 * Bilinear.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
 * Author: felixerben <46817371+felixerben@users.noreply.github.com>
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

#ifndef Hadrons_Bilinear_hpp_
#define Hadrons_Bilinear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TBilinear                                       *
        Performs bilinear contractions of the type tr[g5*adj(qOut')*g5*G*qIn']
        Suitable for non exceptional momenta in Rome-Southampton NPR

Compute the bilinear vertex needed for the NPR.
V(G) = sum_x  [ g5 * adj(S'(x,p2)) * g5 * G * S'(x,p1) ]_{si,sj,ci,cj}
G is one of the 16 gamma vertices [I,gmu,g5,g5gmu,sig(mu,nu)]

        * G
       / \
    p1/   \p2
     /     \
    /       \

Returns a spin-colour matrix, with indices si,sj, ci,cj

Conventions:
p1 - incoming momenta
p2 - outgoing momenta
**************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class BilinearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BilinearPar,
                                    std::string,    qIn,
                                    std::string,    qOut,
                                    std::string,    pIn,
                                    std::string,    pOut,
                                    std::string,    output);
};

template <typename FImpl>
class TBilinear: public Module<BilinearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gamma,
                                        std::string,  pIn,
                                        std::string,  pOut);
    };
    typedef Correlator<Metadata, SpinColourMatrix> Result;
public:
    // constructor
    TBilinear(const std::string name);
    // destructor
    virtual ~TBilinear(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Bilinear, ARG(TBilinear<FIMPL>), MNPR);

/******************************************************************************
 *                           TBilinear implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBilinear<FImpl>::TBilinear(const std::string name)
: Module<BilinearPar>(name)
{}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBilinear<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for Bilinear" << std::endl;

    envTmpLat(PropagatorField, "qIn_phased");
    envTmpLat(PropagatorField, "qOut_phased");
    envTmpLat(PropagatorField, "lret");
    envTmpLat(ComplexField, "pDotXIn");
    envTmpLat(ComplexField, "pDotXOut");
    envTmpLat(ComplexField, "xMu");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBilinear<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().qIn, par().qOut};

    return input;
}

template <typename FImpl>
std::vector<std::string> TBilinear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TBilinear<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

template <typename FImpl>
void TBilinear<FImpl>::execute(void)
{

    LOG(Message) << "Computing bilinear contractions '" << getName() << "' using"
                 << " propagators '" << par().qIn << "' and '" << par().qOut << "'"
                 << std::endl;

    // Propagators
    auto  &qIn    = envGet(PropagatorField, par().qIn);
    auto  &qOut   = envGet(PropagatorField, par().qOut);
    envGetTmp(PropagatorField, qIn_phased);
    envGetTmp(PropagatorField, qOut_phased);
    envGetTmp(PropagatorField, lret);
    envGetTmp(ComplexField, pDotXIn);
    envGetTmp(ComplexField, pDotXOut);
    envGetTmp(ComplexField, xMu);

    // momentum on legs
    std::vector<Real>           pIn  = strToVec<Real>(par().pIn),
	                        pOut = strToVec<Real>(par().pOut);
    Coordinate                  latt_size = GridDefaultLatt();
    Gamma                       g5(Gamma::Algebra::Gamma5);
    Complex                     Ci(0.0,1.0);
    std::vector<Result>         result;
    Result                      r;

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    NPRUtils<FImpl>::dot(pDotXIn,pIn);
    qIn_phased  = qIn  * exp(-Ci * pDotXIn);
    NPRUtils<FImpl>::dot(pDotXOut,pOut);
    qOut_phased = qOut * exp(-Ci * pDotXOut);

    r.info.pIn  = par().pIn;
    r.info.pOut = par().pOut;
    for (auto &G: Gamma::gall)
    {
        r.info.gamma = G.g;
        lret = g5 * adj(qOut_phased) * g5 * G * qIn_phased;
        r.corr.push_back( (1.0 / volume) * sum_large(lret) );
        result.push_back(r);
        r.corr.erase(r.corr.begin());
    }

    //////////////////////////////////////////////////
    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;
    saveResult(par().output, "Bilinear", result);
    auto& out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Bilinear_hpp_
