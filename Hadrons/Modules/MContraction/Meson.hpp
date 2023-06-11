/*
 * Meson.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
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

#ifndef Hadrons_MContraction_Meson_hpp_
#define Hadrons_MContraction_Meson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Meson contractions
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string)
 - q2: input propagator 2 (string)
 - gammas: gamma products to insert at sink & source, pairs of gamma matrices 
           (space-separated strings) in round brackets (i.e. (g_sink g_src)),
           in a sequence (e.g. "(Gamma5 Gamma5)(Gamma5 GammaT)").

           Special values: "all" - perform all possible contractions.
 - sink: module to compute the sink to use in contraction (string).
*/

/******************************************************************************
 *                                TMeson                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class MesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
                                    std::string, sink,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TMeson: public Module<MesonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TMeson(const std::string name);
    // destructor
    virtual ~TMeson(void) {};
    // parse arguments
    virtual int parseGammaString(std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> &gammaMap);
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Meson, ARG(TMeson<FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                           TMeson implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TMeson<FImpl1, FImpl2>::TMeson(const std::string name)
: Module<MesonPar>(name)
{}

// parse arguments /////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
int TMeson<FImpl1, FImpl2>::parseGammaString(std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> &gammaMap)
{
    // Determine gamma matrices to insert at source/sink.
    int nGammaPairs=0;
    if (par().gammas.compare("all") == 0)
    {
        for (Gamma gammaA: Gamma::gall) {
            std::vector<Gamma::Algebra> gMapTmp;
            for (Gamma gammaB: Gamma::gall) {      
                gMapTmp.push_back(gammaB.g);
                nGammaPairs++;
            }
            gammaMap[gammaA.g]=gMapTmp;
        }
    }
    else
    {
        // Parse individual contractions from input string.
        std::vector<GammaPair> tmp;
        tmp = strToVec<GammaPair>(par().gammas);
        // assert that only legal Gamma matrices are used
        for (unsigned int j = 0; j < tmp.size(); j++)
        {
            if( (tmp[j].first==Gamma::Algebra::undef) ||  (tmp[j].second==Gamma::Algebra::undef) )
            {
                HADRONS_ERROR(Argument, "Wrong Argument for Gamma matrices. " + par().gammas); 
            }
        }
        for (Gamma gammaA: Gamma::gall) {
            std::vector<Gamma::Algebra> gMapTmp;
            for (unsigned int j = 0; j < tmp.size(); j++)
            {
                if(tmp[j].first==gammaA.g)
                {
                    gMapTmp.push_back(tmp[j].second);
                    nGammaPairs++;
                }
            }
            gammaMap[gammaA.g]=gMapTmp;
        }
    }

    return nGammaPairs; 
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2};

    if (!par().sink.empty())
    {
        input.push_back(par().sink);
    }
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TMeson<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envTmpLat(LatticePropagator, "q1Gq2");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
#define mesonConnected(q1, q2, gSnk, gSrc) \
(g5*(gSnk))*(q1)*(adj(gSrc)*g5)*adj(q2)

#define mesonConnected1(q1, q2, gSnk) \
adj(q2)*(g5*(gSnk))*(q1)

#define mesonConnected2(q, gSrc) \
(q)*(adj(gSrc)*g5)

template <typename FImpl1, typename FImpl2>
void TMeson<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;
    
    std::vector<TComplex>  buf;
    std::vector<Result>    result;
    Gamma                  g5(Gamma::Algebra::Gamma5);
    int                    nt = env().getDim(Tp);

    std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> gammaMap;
    int nGammas = parseGammaString(gammaMap);
    result.resize(nGammas);
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].corr.resize(nt);
    }
    if (envHasType(SlicedPropagator1, par().q1) and
        envHasType(SlicedPropagator2, par().q2))
    {
        auto &q1 = envGet(SlicedPropagator1, par().q1);
        auto &q2 = envGet(SlicedPropagator2, par().q2);
        
        LOG(Message) << "(propagator already sinked)" << std::endl;
        unsigned int i = 0;
        for(auto &ss: gammaMap)
        {
            Gamma::Algebra gammaSink = ss.first;
            Gamma gSnk(gammaSink);
            for (Gamma::Algebra &gammaSource: ss.second)
            {
                Gamma gSrc(gammaSource);
            
                startTimer("mesonConnected");
                for (unsigned int t = 0; t < nt; ++t)
                {
                    result[i].corr[t] = TensorRemove(trace(mesonConnected(q1[t], q2[t], gSnk, gSrc)));
                }
                stopTimer("mesonConnected");
                result[i].gamma_snk = gSnk.g;
                result[i].gamma_src = gSrc.g;
                i++;
            }
        }
    }
    else
    {
        auto &q1 = envGet(PropagatorField1, par().q1);
        auto &q2 = envGet(PropagatorField2, par().q2);
        
        envGetTmp(LatticeComplex, c);
        envGetTmp(LatticePropagator, q1Gq2);
        if (par().sink.empty())
        {
            HADRONS_ERROR(Definition, "no sink provided");
        }
        LOG(Message) << "(using sink '" << par().sink << "')" << std::endl;
        unsigned int i = 0;
        for(auto &ss: gammaMap)
        {
            Gamma::Algebra gammaSink = ss.first;
            Gamma gSnk(gammaSink);
            startTimer("mesonConnectedSnk");
            q1Gq2 = mesonConnected1(q1,q2,gSnk);
            stopTimer("mesonConnectedSnk");
            for (Gamma::Algebra &gammaSource: ss.second)
            {
                Gamma gSrc(gammaSource);
                std::string ns;
                    
                ns = vm().getModuleNamespace(env().getObjectModule(par().sink));
                if (ns == "MSource")
                {
                    PropagatorField1 &sink = envGet(PropagatorField1, par().sink);
                    
                    startTimer("mesonConnected");
                    c = trace(mesonConnected2(q1Gq2, gSrc)*sink);
                    stopTimer("mesonConnected");
                    startTimer("sliceSum");
                    sliceSum(c, buf, Tp);
                    stopTimer("sliceSum");
                }
                else if (ns == "MSink")
                {
                    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);
                    
                    startTimer("mesonConnected");
                    c   = trace(mesonConnected2(q1Gq2, gSrc));
                    stopTimer("mesonConnected");
                    startTimer("sliceSum");
                    buf = sink(c);
                    stopTimer("sliceSum");
                }
                for (unsigned int t = 0; t < buf.size(); ++t)
                {
                    result[i].corr[t] = TensorRemove(buf[t]);
                }
                result[i].gamma_snk = gSnk.g;
                result[i].gamma_src = gSrc.g;
                i++;
            }
        }
    }
    startTimer("I/O");
    saveResult(par().output, "meson", result);
    stopTimer("I/O");
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Meson_hpp_
