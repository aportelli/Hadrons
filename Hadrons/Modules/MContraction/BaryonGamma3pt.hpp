/*
 * BaryonGamma3pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk.com>
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

#ifndef Hadrons_MContraction_BaryonGamma3pt_hpp_
#define Hadrons_MContraction_BaryonGamma3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/* 
 * 3pt contraction with gamma matrix insertion.
 *
 * Schematic:
 *
 * Wickgroup 1:
 *             qL1          qR1
 *        /---->------*------<---\
 *       /          gamma         \
 *      /                          \
 *   i * ------------->------------ * f
 *      \            qL2           /
 *       \                        /
 *        \----------->----------/
 *                   qL3
 *
 * Wickgroup 2:
 *                   qL1     
 *        /----------->----------\
 *       /                        \
 *      /      qL2          qR2    \
 *   i * ------>------*------<----- * f
 *      \           gamma          /
 *       \                        /
 *        \----------->----------/
 *                   qL3
 *
 * Wickgroup 3:
 *                   qL1 
 *        /----------->----------\
 *       /                        \
 *      /             qL2          \
 *   i * ------------->------------ * f
 *      \                          /
 *       \          gamma         /
 *        \---->------*------<---/
 *             qL3          qR3
 *
 *  options:
 *   - quarksL: ordered triplet of single character quark types
 *   - quarksR: ordered triplet of single character quark types
 *   - quarksJ: ordered double of single character quark types
 *   - qL1: propagator with source at i
 *   - qL2: propagator with source at i
 *   - qL3: propagator with source at i
 *   - qR1: propagator with source at f
 *   - qR2: propagator with source at f
 *   - qR3: propagator with source at f
 *   - sinkq1: sink corresponding to qR1 source
 *   - sinkq2: sink corresponding to qR2 source
 *   - sinkq3: sink corresponding to qR3 source
 *   - gammaLR: list or pairs of pairs of gamma matricies
 *             e.g. "((Identity MinusGammaZGamma5) (Identity GammaT))"
 *   - gammaJ: gamma matrices to insert
 *             (space-separated strings e.g. "GammaT GammaX GammaY")
 *   - mom: 3-momentum componenets
 *             (space-separated strings e.g. "0 0 0")
 *   - tSnk: sink position for spectator propagators
 */

/******************************************************************************
 *                               BaryonGamma3pt                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


#if (!defined(GRID_HIP))
typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaAB;
typedef std::pair<GammaAB, GammaAB> GammaABPair;

class BaryonGamma3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonGamma3ptPar,
                                    std::string,  quarksL,
                                    std::string,  quarksR,
                                    std::string,  quarksJ,
                                    std::string,  qL1,
                                    std::string,  qL2,
                                    std::string,  qL3,
                                    std::string,  qR1,
                                    std::string,  qR2,
                                    std::string,  qR3,
                                    std::string,  sink1,
                                    std::string,  sink2,
                                    std::string,  sink3,
                                    std::string,  gammaLR,
                                    std::string,  gammaJ,
                                    std::string,  mom,
                                    unsigned int, tf,
                                    std::string,  output);
};

template <typename FImpl>
class TBaryonGamma3pt: public Module<BaryonGamma3ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaJ,
                                        Gamma::Algebra, gammaAi,
                                        Gamma::Algebra, gammaBi,
                                        Gamma::Algebra, gammaAf,
                                        Gamma::Algebra, gammaBf,
                                        std::string,  quarksL,
                                        std::string,  quarksR,
                                        std::string,  quarksJ,
                                        std::string,  mom);
    };
    typedef Correlator<Metadata, SpinMatrix> Result;
public:
    // constructor
    TBaryonGamma3pt(const std::string name);
    // destructor
    virtual ~TBaryonGamma3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void parseGammaJString(std::string gamma, std::vector<Gamma::Algebra> &gammaList);
    virtual void parseGammaLRString(std::string gammas, std::vector<GammaABPair> &gammaList);

    std::vector<Real> mom_;
    std::string momphName_;
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BaryonGamma3pt, ARG(TBaryonGamma3pt<FIMPL>), MContraction);

/******************************************************************************
 *                       TBaryonGamma3pt implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBaryonGamma3pt<FImpl>::TBaryonGamma3pt(const std::string name)
: Module<BaryonGamma3ptPar>(name),
  momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBaryonGamma3pt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qL1, par().qL2, par().qL3, par().sink1, par().sink2, par().sink3};

    std::string qR = {par().quarksR[0], par().quarksR[1], par().quarksR[2]};
    char qJR = par().quarksJ[1];

    if (qR[0] == qJR)
        in.push_back(par().qR1);
    if (qR[1] == qJR)
        in.push_back(par().qR2);
    if (qR[2] == qJR)
        in.push_back(par().qR3);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TBaryonGamma3pt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}
template <typename FImpl>
std::vector<std::string> TBaryonGamma3pt<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryonGamma3pt<FImpl>::setup(void)
{
    auto parse_vector = [](const std::string &vec, int dim,
            const std::string &desc)
    {
        std::vector<Real> res = strToVec<Real>(vec);
        if(res.size() != dim) {
            HADRONS_ERROR(Size, desc + " has "
                    + std::to_string(res.size()) + " instead of "
                    + std::to_string(dim) + " components");
        }
        return res;
    };
    mom_      = parse_vector(par().mom, env().getNd()-1, "momentum");

    envTmpLat(SpinMatrixField, "c");
    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
}

template <typename FImpl>
void TBaryonGamma3pt<FImpl>::parseGammaJString(std::string gamma, std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at the current insertion.
    if (gamma.compare("all") == 0) {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2) {
            gammaList.push_back((Gamma::Algebra)i);
        }
    } else {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(gamma);
    } 
}

template <typename FImpl>
void TBaryonGamma3pt<FImpl>::parseGammaLRString(std::string gammas, std::vector<GammaABPair> &gammaList)
{
    gammaList.clear();
    
    std::string gammaString = gammas;
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
    for (int i = 0; i < gammaList.size(); i++) {
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


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryonGamma3pt<FImpl>::execute(void)
{
    std::string qL = {par().quarksL[0], par().quarksL[1], par().quarksL[2]};
    std::string qR = {par().quarksR[0], par().quarksR[1], par().quarksR[2]};
    char qLJ = par().quarksJ[0];
    char qRJ = par().quarksJ[1];
    
    std::vector<GammaABPair>    gammaLRList;
    parseGammaLRString(par().gammaLR, gammaLRList);

    std::vector<Gamma::Algebra> gammaJList;
    parseGammaJString(par().gammaJ, gammaJList);


    LOG(Message) << "Computing baryon 3pt contractions '" << getName() << "'" << std::endl;
    LOG(Message) << "using quarksL (" << qL << "), quarksR (" << qR << ") and quarksJ (" << qLJ << " -> " << qRJ << ")" << std::endl;
    LOG(Message) << "using propagatorsL '" << par().qL1 << "', '" << par().qL2 << "' and '" << par().qL3 << "'" << std::endl;
    LOG(Message) << "using propagatorsR '" << par().qR1 << "', '" << par().qR2 << "' and '" << par().qR3 << "'" << std::endl;
    if (qR[0] != qRJ && par().qR1 != "")
        LOG(Message) << "  note: qR1 was specified as '" << par().qR1 << "' but is not used" << std::endl;
    if (qR[1] != qRJ && par().qR2 != "")
        LOG(Message) << "  note: qR2 was specified as '" << par().qR2 << "' but is not used" << std::endl;
    if (qR[2] != qRJ && par().qR3 != "")
        LOG(Message) << "  note: qR3 was specified as '" << par().qR3 << "' but is not used" << std::endl;
    LOG(Message) << "using sinks '" << par().sink1 << "', '" << par().sink2 << "' and '" << par().sink3 << "'" << std::endl;
    for (int iG = 0; iG < gammaLRList.size(); iG++)
        LOG(Message) << "  with (Gamma^A,Gamma^B)_left = ( " << gammaLRList[iG].first.first << " , " << gammaLRList[iG].first.second << "') and (Gamma^A,Gamma^B)_right = ( " << gammaLRList[iG].second.first << " , " << gammaLRList[iG].second.second << ")" << std::endl; 
    for (int iG = 0; iG < gammaJList.size(); iG++)
        LOG(Message) << "  inserting gammaJ = (" << gammaJList[iG] << ")" << std::endl;
    LOG(Message) << "with momentum " << mom_ << std::endl;


    int epsilon[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    bool wick_groups[3];
    bool wick_contractions[3][6];

    wick_groups[0] = (qL[0] == qLJ)? true : false ;
    wick_groups[1] = (qL[1] == qLJ)? true : false ;
    wick_groups[2] = (qL[2] == qLJ)? true : false ;

    for (int i=0; i<6; i++) 
    {
        wick_contractions[0][i] = ( wick_groups[0]
                                    && qRJ   == qR[epsilon[i][0]] 
                                    && qL[1] == qR[epsilon[i][1]] 
                                    && qL[2] == qR[epsilon[i][2]] )? true : false ;
        wick_contractions[1][i] = ( wick_groups[1]
                                    && qL[0] == qR[epsilon[i][0]] 
                                    && qRJ   == qR[epsilon[i][1]] 
                                    && qL[2] == qR[epsilon[i][2]] )? true : false ;
        wick_contractions[2][i] = ( wick_groups[2]
                                    && qL[0] == qR[epsilon[i][0]] 
                                    && qL[1] == qR[epsilon[i][1]] 
                                    && qRJ   == qR[epsilon[i][2]] )? true : false ;
    }

    int n_contractions = 0;
    for (int i=0; i<3; i++) 
    {
        if (wick_groups[i]) 
        {
            for (int j=0; j<6; j++) 
            {
                LOG(Message) << "Wick group " << i+1 << " : contraction " << j+1 << " : " << ( (wick_contractions[i][j]) ? "true" : "false" )  << std::endl;
                n_contractions += (wick_contractions[i][j])? 1 : 0 ;
            }
        }
    }
    if (n_contractions == 0)
        LOG(Error) << "No contractions are performed. Check the quark flavors" << std::endl; // TODO: Should this error be fatal?


    SinkFn* sink[3];
    sink[0] = &envGet(SinkFn, par().sink1);
    sink[1] = &envGet(SinkFn, par().sink2);
    sink[2] = &envGet(SinkFn, par().sink3);

    PropagatorField &propQL1  = envGet(PropagatorField, par().qL1);
    PropagatorField &propQL2  = envGet(PropagatorField, par().qL2);
    PropagatorField &propQL3  = envGet(PropagatorField, par().qL3);

    PropagatorField* propQR[3];
    propQR[0]  = ( qR[0] == qRJ ) ? &envGet(PropagatorField, par().qR1) : nullptr;
    propQR[1]  = ( qR[1] == qRJ ) ? &envGet(PropagatorField, par().qR2) : nullptr;
    propQR[2]  = ( qR[2] == qRJ ) ? &envGet(PropagatorField, par().qR3) : nullptr;

    SlicedPropagator propQL1_slice[6];
    SlicedPropagator propQL2_slice[6];
    SlicedPropagator propQL3_slice[6];

    for (int ie=0; ie<6; ie++) 
    {
        propQL1_slice[ie] = (*sink[epsilon[ie][0]])(propQL1);
        propQL2_slice[ie] = (*sink[epsilon[ie][1]])(propQL2);
        propQL3_slice[ie] = (*sink[epsilon[ie][2]])(propQL3);
    }


    int nt = env().getDim(Tp);
    envGetTmp(SpinMatrixField, c);
    std::vector<SpinMatrix> buf;

    std::vector<Result>         result;
    Result                      r;

    r.info.quarksL = par().quarksL;
    r.info.quarksR = par().quarksR;
    r.info.quarksJ = par().quarksJ;

    r.info.mom = par().mom;

    for (auto &GJ:  gammaJList)  
    {
        for (auto &GLR: gammaLRList) 
        {
            LOG(Message) << "Contracting for gammaJ  = " << GJ << std::endl;
            LOG(Message) << "                gammaLR = " << GLR << std::endl;

            r.info.gammaJ  = GJ;
            r.info.gammaAi = GLR.first.first;
            r.info.gammaBi = GLR.first.second;
            r.info.gammaAf = GLR.second.first;
            r.info.gammaBf = GLR.second.second;

            const Gamma GAL(GLR.first.first);
            const Gamma GBL(GLR.first.second);
            const Gamma GAR(GLR.second.first);
            const Gamma GBR(GLR.second.second);

            c=Zero();

            for (int ie=0; ie<6; ie++) 
            {
                if (wick_groups[0]) 
                {
                    if (wick_contractions[0][ie]) 
                    {
                        auto propQ2_spec = propQL2_slice[ie][par().tf];
                        auto propQ3_spec = propQL3_slice[ie][par().tf];

                        BaryonUtils<FIMPL>::BaryonGamma3pt(propQL1, propQ2_spec, propQ3_spec, *propQR[epsilon[ie][0]],
                                                         1, ie+1, GJ, GLR.first.second, GLR.second.second, c);
                    }
                }
                if (wick_groups[1]) 
                {
                    if (wick_contractions[1][ie]) 
                    {
                        auto propQ1_spec = propQL1_slice[ie][par().tf];
                        auto propQ3_spec = propQL3_slice[ie][par().tf];

                        BaryonUtils<FIMPL>::BaryonGamma3pt(propQL2, propQ1_spec, propQ3_spec, *propQR[epsilon[ie][1]],
                                                         2, ie+1, GJ, GLR.first.second, GLR.second.second, c);
                    }
                }
                if (wick_groups[2]) 
                {
                    if (wick_contractions[2][ie]) 
                    {
                        auto propQ1_spec = propQL1_slice[ie][par().tf];
                        auto propQ2_spec = propQL2_slice[ie][par().tf];

                        BaryonUtils<FIMPL>::BaryonGamma3pt(propQL3, propQ1_spec, propQ2_spec, *propQR[epsilon[ie][2]],
                                                         3, ie+1, GJ, GLR.first.second, GLR.second.second, c);
                    }
                }
            }
        
            auto &ph = envGet(LatticeComplex, momphName_);
        
            if (mom_[0] != 0 || mom_[1] != 0 || mom_[2] != 0) 
            {
                LOG(Message) << "Adding momentum phase " << mom_ << std::endl;

                Complex           i(0.0,1.0);

                envGetTmp(LatticeComplex, coor);
                ph = Zero();
                for(unsigned int mu = 0; mu < 3; mu++)
                {
                    LatticeCoordinate(coor, mu);
                    ph = ph + (mom_[mu]/env().getDim(mu))*coor;
                }
                ph = exp((Real)(2*M_PI)*i*ph);
                c = ph*c;
            }


            sliceSum(c,buf,Tp);
            r.corr.clear();
            for (unsigned int t = 0; t < buf.size(); ++t) 
            {
                buf[t]() = GAR * buf[t]() * GAL;
                r.corr.push_back( buf[t] );
            }
            result.push_back(r);
        }
    }

    saveResult(par().output, "baryongamma3pt", result);
}

#endif
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_BaryonGamma3pt_hpp_
