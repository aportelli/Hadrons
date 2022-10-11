/*
 * ExternalLeg.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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

#ifndef Hadrons_ExternalLeg_hpp_
#define Hadrons_ExternalLeg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TExternalLeg                                       *
        Computes the propagator from a momentum source with
        appropriate phase corrections. Output is a SpinColourMatrix.
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class ExternalLegPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExternalLegPar,
                                    std::string,    qIn,
                                    std::string,    pIn,
                                    std::string,    output);
};

template <typename FImpl>
class TExternalLeg: public Module<ExternalLegPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string,  pIn);
    };
    typedef Correlator<Metadata, SpinColourMatrix> Result;
public:
    // constructor
    TExternalLeg(const std::string name);
    // destructor
    virtual ~TExternalLeg(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExternalLeg, ARG(TExternalLeg<FIMPL>), MNPR);

/******************************************************************************
 *                           TExternalLeg implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExternalLeg<FImpl>::TExternalLeg(const std::string name)
: Module<ExternalLegPar>(name)
{}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExternalLeg<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for ExternalLeg" << std::endl;

    envTmpLat(PropagatorField, "qIn_phased");
    envTmpLat(ComplexField, "pDotXIn");
    envTmpLat(ComplexField, "xMu");
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExternalLeg<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().qIn};

    return input;
}

template <typename FImpl>
std::vector<std::string> TExternalLeg<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
std::vector<std::string> TExternalLeg<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExternalLeg<FImpl>::execute(void)
{
    LOG(Message) << "Computing propagator '" << getName() << "' using"
                 << " momentum '" << par().pIn << "'"
                 << std::endl;
    auto                &qIn    = envGet(PropagatorField, par().qIn);
    envGetTmp(PropagatorField, qIn_phased);
    envGetTmp(ComplexField, pDotXIn);
    envGetTmp(ComplexField, xMu);
    std::vector<Real>   pIn  = strToVec<Real>(par().pIn);
    Coordinate          latt_size = GridDefaultLatt();
    Gamma               g5(Gamma::Algebra::Gamma5);
    Complex             Ci(0.0,1.0);
    Result              r;


    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    NPRUtils<FImpl>::dot(pDotXIn,pIn);
    qIn_phased = qIn * exp(-Ci * pDotXIn); // phase corrections

    r.info.pIn  = par().pIn;
    r.corr.push_back( (1.0 / volume) * sum(qIn_phased) );

    saveResult(par().output, "ExternalLeg", r);
    LOG(Message) << "Complete. Writing results to " << par().output << std:: endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ExternalLeg_hpp_
