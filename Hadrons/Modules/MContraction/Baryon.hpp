/*
 * Baryon.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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

#ifndef Hadrons_MContraction_Baryon_hpp_
#define Hadrons_MContraction_Baryon_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Baryon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaAB;
typedef std::pair<GammaAB, GammaAB> GammaABPair;

class BaryonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, quarks,
                                    std::string, shuffle,
                                    std::string, sinkq1,
                                    std::string, sinkq2,
                                    std::string, sinkq3,
                                    bool, sim_sink,
                                    std::string, gammas,
                                    int, parity,
                                    std::string, output);
};

template <typename FImpl>
class TBaryon: public Module<BaryonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaA_left,
                                        Gamma::Algebra, gammaB_left,
                                        Gamma::Algebra, gammaA_right,
                                        Gamma::Algebra, gammaB_right,
                                        std::string, quarksR,
                                        std::string, quarksL,
                                        std::string, shuffle,
                                        int, parity);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TBaryon(const std::string name);
    // destructor
    virtual ~TBaryon(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void parseGammaString(std::vector<GammaABPair> &gammaList);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Which gamma algebra was specified
    Gamma::Algebra  al;
};

MODULE_REGISTER_TMP(Baryon, ARG(TBaryon<FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBaryon<FImpl>::TBaryon(const std::string name)
: Module<BaryonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBaryon<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3, par().sinkq1, par().sinkq2, par().sinkq3};
    
    return input;
}

template <typename FImpl>
std::vector<std::string> TBaryon<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}
template <typename FImpl>
std::vector<std::string> TBaryon<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

template <typename FImpl>
void TBaryon<FImpl>::parseGammaString(std::vector<GammaABPair> &gammaList)
{
    gammaList.clear();
    
    std::string gammaString = par().gammas;
    //Shorthands for standard baryon operators
    gammaString = regex_replace(gammaString, std::regex("j12"),"(Identity SigmaXZ)");
    gammaString = regex_replace(gammaString, std::regex("j32X"),"(Identity MinusGammaZGamma5)");
    gammaString = regex_replace(gammaString, std::regex("j32Y"),"(Identity GammaT)");
    gammaString = regex_replace(gammaString, std::regex("j32Z"),"(Identity GammaXGamma5)");
    //Shorthands for less common baryon operators
    gammaString = regex_replace(gammaString, std::regex("j12_alt1"),"(Gamma5 MinusSigmaYT)");
    gammaString = regex_replace(gammaString, std::regex("j12_alt2"),"(Identity GammaYGamma5)");
    
    //A single gamma matrix 
    std::regex rex_g("([0-9a-zA-Z]+)");
    //The full string we expect
    std::regex rex("( *\\(( *\\(([0-9a-zA-Z]+) +([0-9a-zA-Z]+) *\\)){2} *\\) *)+");
    std::smatch sm;
    std::regex_match(gammaString, sm, rex);
    assert(sm[0].matched && "invalid gamma structure.");

    auto gamma_begin = std::sregex_iterator(gammaString.begin(), gammaString.end(), rex_g);
    auto gamma_end = std::sregex_iterator();

    int nGamma = std::distance(gamma_begin, gamma_end); 
    //couldn't find out how to count the size in the iterator, other than looping through it...
  /*  int nGamma=0;
    for (std::sregex_iterator i = gamma_begin; i != gamma_end; ++i) {
	nGamma++;
    }
*/   
    gammaList.resize(nGamma/4);
    std::vector<std::string> gS;
    gS.resize(nGamma);
    //even more ugly workarounds here...
    int iG=0;
    for (std::sregex_iterator i = gamma_begin; i != gamma_end; ++i) {
        std::smatch match = *i;                                                 
        gS[iG] = match.str(); 
	iG++;
    }
    for (int i = 0; i < gammaList.size(); i++){
	std::vector<Gamma::Algebra> gS1 = strToVec<Gamma::Algebra>(gS[4*i]);
	std::vector<Gamma::Algebra> gS2 = strToVec<Gamma::Algebra>(gS[4*i+1]);
	std::vector<Gamma::Algebra> gS3 = strToVec<Gamma::Algebra>(gS[4*i+2]);
	std::vector<Gamma::Algebra> gS4 = strToVec<Gamma::Algebra>(gS[4*i+3]);
        gammaList[i].first.first=gS1[0];
        gammaList[i].first.second=gS2[0];
        gammaList[i].second.first=gS3[0];
        gammaList[i].second.second=gS4[0];
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon<FImpl>::execute(void)
{
    // Check shuffle is a permutation of "123"
    assert(par().shuffle.size()==3 && "shuffle parameter must be 3 characters long");
    std::string shuffle_tmp = par().shuffle;
    std::sort(shuffle_tmp.begin(), shuffle_tmp.end());
    assert(shuffle_tmp == "123" && "shuffle parameter must be a permulation of 123");

    std::vector<int> shuffle = { std::stoi( par().shuffle.substr(0,1) ) -1,
                                 std::stoi( par().shuffle.substr(1,1) ) -1,
                                 std::stoi( par().shuffle.substr(2,1) ) -1 };


    assert(par().quarks.size()==3 && "quark-structure must consist of 3 quarks");
    std::string quarksR = par().quarks;
    std::string quarksL = quarksR;

    // Shuffle quark flavours
    quarksL[0] = quarksR[shuffle[0]];
    quarksL[1] = quarksR[shuffle[1]];
    quarksL[2] = quarksR[shuffle[2]];

    std::vector<std::string> props(3, "");
    props[0] = par().q1;
    props[1] = par().q2;
    props[2] = par().q3;

    // Shuffle propagators
    std::vector<std::string> propsL(props);
    propsL[0] = props[shuffle[0]];
    propsL[1] = props[shuffle[1]];
    propsL[2] = props[shuffle[2]];

    if (par().sim_sink)
        assert(par().sinkq1==par().sinkq2 && par().sinkq2==par().sinkq3 && "when sim_sink is true all three sinks must be the same");

    assert(par().parity == 1 || par().parity == -1 && "parity must be 1 or -1");

    std::vector<GammaABPair> gammaList;
    parseGammaString(gammaList);


    LOG(Message) << "Computing baryon contractions '" << getName() << "'" << std::endl;
    LOG(Message) << "  using shuffle " << shuffle << " and parity " << par().parity << std::endl;
    LOG(Message) << "  using quarksL (" << quarksL << ") with left propagators (" << propsL[0] << ", " << propsL[1] << ", and " << propsL[2] << ")" << std::endl;
    LOG(Message) << "  using quarksR (" << quarksR << ") ";
    if (par().sim_sink)
        LOG(Message) << "with simultaneous sink " << par().sinkq1 << std::endl;
    else 
        LOG(Message) << "with sinks (" << par().sinkq1 << ", " << par().sinkq2 << ", and " << par().sinkq3 << ")" << std::endl;

    for (int iG = 0; iG < gammaList.size(); iG++)
        LOG(Message) << "    with (Gamma^A,Gamma^B)_left = ( " << gammaList[iG].first.first << " , " << gammaList[iG].first.second << "') and (Gamma^A,Gamma^B)_right = ( " << gammaList[iG].second.first << " , " << gammaList[iG].second.second << ")" << std::endl;
    
    envGetTmp(LatticeComplex, c);
    int nt = env().getDim(Tp);

    std::vector<Result> result;
    Result              r;
    r.info.parity  = par().parity;

    r.info.quarksR  = quarksR;
    r.info.quarksL  = quarksL;
    r.info.shuffle = par().shuffle;
        
    bool wick_contractions[6];
    BaryonUtils<FIMPL>::Wick_Contractions(quarksL,quarksR,wick_contractions);
    
    PropagatorField &q1  = envGet(PropagatorField, propsL[0]);
    PropagatorField &q2  = envGet(PropagatorField, propsL[1]);
    PropagatorField &q3  = envGet(PropagatorField, propsL[2]);

    if (par().sim_sink) { 
        for (unsigned int i = 0; i < gammaList.size(); ++i)
        {
            std::vector<TComplex> buf;

            r.info.gammaA_left = gammaList[i].first.first;
            r.info.gammaB_left = gammaList[i].first.second;
            r.info.gammaA_right = gammaList[i].second.first;
            r.info.gammaB_right = gammaList[i].second.second;

            Gamma gAl(gammaList[i].first.first);
            Gamma gBl(gammaList[i].first.second);
            Gamma gAr(gammaList[i].second.first);
            Gamma gBr(gammaList[i].second.second);
        
            std::string ns = vm().getModuleNamespace(env().getObjectModule(par().sinkq1));
            if (ns == "MSource")
            {
                c=Zero();
                BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,
                                                    gAl,gBl,gAr,gBr,
                                                    wick_contractions,
                                                    par().parity,
                                                    c);

                PropagatorField &sink = envGet(PropagatorField, par().sinkq1);
                auto test = closure(trace(sink*c));     
                sliceSum(test, buf, Tp); 
            }
            else if (ns == "MSink")
            {
                c=Zero();
                BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,
                                                    gAl,gBl,gAr,gBr,
                                                    wick_contractions,
                                                    par().parity,
                                                    c);

                SinkFnScalar &sink = envGet(SinkFnScalar, par().sinkq1);
                buf = sink(c);
            } 
            r.corr.clear();
            for (unsigned int t = 0; t < buf.size(); ++t)
            {
                r.corr.push_back(TensorRemove(buf[t]));
            }
            result.push_back(r);
        }
    } else {
        const int epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};

        SinkFn* sinkFn[3];
        sinkFn[0] = &envGet(SinkFn, par().sinkq1);
        sinkFn[1] = &envGet(SinkFn, par().sinkq2);
        sinkFn[2] = &envGet(SinkFn, par().sinkq3);
        
        SlicedPropagator q1_slice[6];
        SlicedPropagator q2_slice[6];
        SlicedPropagator q3_slice[6];

        // Compute all sinked propagators
        for (int ie=0; ie < 6 ; ie++) {
            q1_slice[ie] = (*sinkFn[epsilon[ie][0]])(q1);
            q2_slice[ie] = (*sinkFn[epsilon[ie][1]])(q2);
            q3_slice[ie] = (*sinkFn[epsilon[ie][2]])(q3);
        }

        std::vector<TComplex> buf;
        for (int iG = 0; iG < gammaList.size(); iG++) {
            buf = std::vector<TComplex>(nt, Zero());

            Gamma gAl(gammaList[iG].first.first);
            Gamma gBl(gammaList[iG].first.second);
            Gamma gAr(gammaList[iG].second.first);
            Gamma gBr(gammaList[iG].second.second);

            for (int ie=0; ie < 6 ; ie++) {
                if (wick_contractions[ie]) {
                    // Perform the current contraction only
                    bool wc[6];
                    for (int i=0; i<6; i++)
                        wc[i] = (i == ie);

                    BaryonUtils<FIMPL>::ContractBaryons_Sliced( q1_slice[ie],q2_slice[ie],q3_slice[ie],
                                                                gAl,gBl,gAr,gBr,
                                                                wc,
                                                                par().parity,
                                                                nt,
                                                                buf);
                }
            }
            r.info.gammaA_left = gammaList[iG].first.first;
            r.info.gammaB_left = gammaList[iG].first.second;
            r.info.gammaA_right = gammaList[iG].second.first;
            r.info.gammaB_right = gammaList[iG].second.second;

            r.corr.clear();
            for (int t = 0; t < nt; t++) {
                r.corr.push_back(TensorRemove(buf[t]));
            }
            result.push_back(r);
        }
        
    }

    saveResult(par().output, "baryon", result);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_hpp_
