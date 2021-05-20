/*
 * Bilinear.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
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

#ifndef Hadrons_Bilinear_hpp_
#define Hadrons_Bilinear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

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
                                        std::string,  description,
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
    envTmpLat(ComplexField, "pDotXIn");
    envTmpLat(ComplexField, "pDotXOut");
    envTmpLat(ComplexField, "xMu");
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
    std::vector<std::string> out = {};
    
    return out;
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
    LatticeSpinColourMatrix qIn_phased(env().getGrid()), qOut_phased(env().getGrid());
    envGetTmp(ComplexField, pDotXIn);
    envGetTmp(ComplexField, pDotXOut);
    envGetTmp(ComplexField, xMu);

    // momentum on legs
    //TODO: Do we want to check the momentum input format? Not done in MSink::Point, so probably ok like this.
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

    pDotXIn=Zero();
    pDotXOut=Zero();
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0 / latt_size[mu];
        LatticeCoordinate(xMu,mu);
        pDotXIn  = pDotXIn  + (TwoPiL * pIn[mu])  * xMu;
        pDotXOut = pDotXOut + (TwoPiL * pOut[mu]) * xMu;
    }
    qIn_phased  = qIn  * exp(-Ci * pDotXIn); //phase corrections
    qOut_phased = qOut * exp(-Ci * pDotXOut);
    
    r.info.pIn  = par().pIn; // Redundant to write these into every group
    r.info.pOut = par().pOut; // Redundant to write these into every group
    for (auto &G: Gamma::gall)
    {
    	r.info.description = Gamma::name[G.g]; // The change from Gamma::Algebra to string causes all strings to have the same length
                                               // which leads to trailing spaces in the string. Is there an easy way to avoid this?
    	r.corr.push_back( (1.0 / volume) * sum(g5 * adj(qOut_phased) * g5 * G * qIn_phased) );
        result.push_back(r);
    	//This is all still quite hacky - we probably want to think about the output format a little more!
    	r.corr.erase(r.corr.begin());
    }

    // Also write the propagators to the outfile
    r.info.pIn  = par().pIn;
    r.info.pOut = " ";
    r.info.description = "qIn";
    r.corr.push_back( (1.0 / volume) * sum(qIn_phased) );
    result.push_back(r);
    r.corr.erase(r.corr.begin());

    r.info.pIn  = " ";
    r.info.pOut = par().pOut;
    r.info.description = "qOut";
    r.corr.push_back( (1.0 / volume) * sum(g5 * adj(qOut_phased) * g5) );
    result.push_back(r);
    r.corr.erase(r.corr.begin());

    //////////////////////////////////////////////////
    saveResult(par().output, "bilinear", result);
    LOG(Message) << "Complete. Writing results to " << par().output << std:: endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Bilinear_hpp_
