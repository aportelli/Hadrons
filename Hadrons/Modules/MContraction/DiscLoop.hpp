/*
 * DiscLoop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

#ifndef Hadrons_MContraction_DiscLoop_hpp_
#define Hadrons_MContraction_DiscLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Disconnected loop contractions
 ------------------------------
 
 * options:
 - q_loop: input propagator (string)
 - gammas: gammas: gamma matrices to insert
           (space-separated strings e.g. "GammaT GammaX GammaY") 

           Special values: "all" - perform all possible contractions.
 - mom:    momentum insertion, vector of space-separated int sequence
           (e.g {"0 0 0", "1 0 0", "0 2 0"})
*/

/******************************************************************************
 *                                DiscLoop                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DiscLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DiscLoopPar,
                                    std::string,              q_loop,
                                    std::string,              gammas,
                                    std::vector<std::string>, mom,
                                    std::string,              output);
};

template <typename FImpl>
class TDiscLoop: public Module<DiscLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef std::vector<SitePropagator> SlicedOp;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<int>, mom,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TDiscLoop(const std::string name);
    // destructor
    virtual ~TDiscLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<std::vector<int>> mom_;
};

MODULE_REGISTER_TMP(DiscLoop, TDiscLoop<FIMPL>, MContraction);

/******************************************************************************
 *                       TDiscLoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDiscLoop<FImpl>::TDiscLoop(const std::string name)
: Module<DiscLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDiscLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q_loop};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDiscLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TDiscLoop<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiscLoop<FImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();
    mom_.resize(par().mom.size());
    for (unsigned int i = 0; i < mom_.size(); ++i)
    {
        mom_[i] = strToVec<int>(par().mom[i]);
        if (mom_[i].size() != nd - 1)
        {
            HADRONS_ERROR(Size, "momentum number of components different from " 
                               + std::to_string(nd-1));
        }
        for (unsigned int j = 0; j < nd - 1; ++j)
        {
            mom_[i][j] = (mom_[i][j] + env().getDim(j)) % env().getDim(j);
        }
    }
    envTmpLat(PropagatorField, "ftBuf");
    envTmpLat(PropagatorField, "op");
}

template <typename FImpl>
void TDiscLoop<FImpl>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gammas);
    } 
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiscLoop<FImpl>::execute(void)
{
    LOG(Message) << "Computing disconnected loop contraction '" << getName() 
                 << "' using '" << par().q_loop << "' with " << par().gammas 
                 << " insertion and momentum " << par().mom
                 << "." << std::endl;

    const unsigned int                 nt      = env().getDim(Tp);
    const unsigned int                 nd      = env().getDim().size();
    const unsigned int                 nmom    = mom_.size();
    auto                               &q_loop = envGet(PropagatorField, par().q_loop);
    std::vector<Gamma::Algebra>        gammaList;
    SitePropagator                     buf;
    std::vector<std::vector<SlicedOp>> slicedOp;
    std::vector<std::vector<Result>>   result;
    FFT                                fft(envGetGrid(PropagatorField));
    std::vector<int>                   dMask(nd, 1);

    dMask[nd - 1] = 0;
    envGetTmp(PropagatorField, ftBuf);
    envGetTmp(PropagatorField, op);
    parseGammaString(gammaList);
    const unsigned int ngam = gammaList.size();
    result.resize(ngam);
    for (unsigned int g = 0; g < ngam; ++g)
    {
        result[g].resize(nmom);
        for (unsigned int m = 0; m < nmom; ++m)
        {
            result[g][m].gamma = gammaList[g];
            result[g][m].mom   = mom_[m];
            result[g][m].corr.resize(nt);
        }
    }

    slicedOp.resize(ngam);
    for (unsigned int g = 0; g < ngam; ++g)
    {
        Gamma gamma(gammaList[g]);
        op = gamma*q_loop;
        fft.FFT_dim_mask(ftBuf, op, dMask, FFT::forward);
        slicedOp[g].resize(nmom);
        for (unsigned int m = 0; m < nmom; ++m)
        {
            auto qt = mom_[m];
            qt.resize(nd);
            slicedOp[g][m].resize(nt);
            for (unsigned int t = 0; t < nt; ++t)
            {
                qt[nd - 1] = t;
                peekSite(buf, ftBuf, qt);
                slicedOp[g][m][t] = buf;
            }
        }
    }
    for (unsigned int g = 0; g < ngam; ++g)
    for (unsigned int m = 0; m < nmom; ++m)
    {
        for (unsigned int t = 0; t < nt; ++t)
        {
            result[g][m].corr[t] = TensorRemove(trace(slicedOp[g][m][t]));
        }
    }
    saveResult(par().output, "disc", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DiscLoop_hpp_
