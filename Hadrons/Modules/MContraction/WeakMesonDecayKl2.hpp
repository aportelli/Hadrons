/*
 * WeakMesonDecayKl2.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Andrew Zhen Ning Yong <andrew.yong@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
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

#ifndef Hadrons_MContraction_WeakMesonDecayKl2_hpp_
#define Hadrons_MContraction_WeakMesonDecayKl2_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
* Kl2 contraction
* -----------------------------
*
* contraction for Kl2 decay, including the lepton
*
* 	trace(q1*adj(q2)*g5*gL[mu]) * (gL[mu] * lepton)_{a,b}
*
* with open spinor indices (a,b) for the lepton part
*
*             q1                  lepton
*        /------------\       /------------
*       /              \     /
*      /                \H_W/
* g_5 *                  * * 
*      \                /
*       \              / 
*        \____________/
*             q2
*
* * options:
* - q1: input propagator 1 (string)
* - q2: input propagator 2 (string)
* - lepton: input lepton (string)
*/

/******************************************************************************
 *                               TWeakMesonDecayKl2                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WeakMesonDecayKl2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakMesonDecayKl2Par,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, lepton,
				                    std::string, output);
};

template <typename FImpl, typename LImpl>
class TWeakMesonDecayKl2: public Module<WeakMesonDecayKl2Par>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(LImpl, Lepton);
    typedef typename SpinMatrixField::vector_object::scalar_object SpinMatrix;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<SpinMatrix>, corr);
    };
public:
    // constructor
    TWeakMesonDecayKl2(const std::string name);
    // destructor
    virtual ~TWeakMesonDecayKl2(void) {};
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

MODULE_REGISTER_TMP(WeakMesonDecayKl2, ARG(TWeakMesonDecayKl2<FIMPL, LIMPL>), MContraction);

/******************************************************************************
 *                           TWeakMesonDecayKl2 implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename LImpl>
TWeakMesonDecayKl2<FImpl, LImpl>::TWeakMesonDecayKl2(const std::string name)
: Module<WeakMesonDecayKl2Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename LImpl>
std::vector<std::string> TWeakMesonDecayKl2<FImpl, LImpl>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().lepton};
    
    return input;
}

template <typename FImpl, typename LImpl>
std::vector<std::string> TWeakMesonDecayKl2<FImpl, LImpl>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

template <typename FImpl, typename LImpl>
std::vector<std::string> TWeakMesonDecayKl2<FImpl, LImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ////////////////////////////////////////////////////////////////////////
template <typename FImpl, typename LImpl>
void TWeakMesonDecayKl2<FImpl, LImpl>::setup(void)
{
    envTmpLat(ComplexField, "c");
    envTmpLat(PropagatorFieldLepton, "res_buf");
    envTmpLat(SpinMatrixField, "buf");
    envCreate(HadronsSerializable, getName(), 1, 0);
    envTmpLat(LatticeComplex, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename LImpl>
void TWeakMesonDecayKl2<FImpl, LImpl>::execute(void)
{
    LOG(Message) << "Computing QED Kl2 contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "' and"
		         << "lepton '"  << par().lepton << "'" << std::endl;

    Gamma                   g5(Gamma::Algebra::Gamma5);
    int                     nt = env().getDim(Tp);
    std::vector<SpinMatrix> res_summed;
    Result                  r;

    auto &q1     = envGet(PropagatorField, par().q1);
    auto &q2     = envGet(PropagatorField, par().q2);
    auto &lepton = envGet(PropagatorFieldLepton, par().lepton);
    envGetTmp(SpinMatrixField, buf);
    envGetTmp(ComplexField, c);
    envGetTmp(PropagatorFieldLepton, res_buf); res_buf = Zero();
    
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        c = Zero();
        //hadronic part: trace(q1*adj(q2)*g5*gL[mu]) 
        c = trace(q1*adj(q2)*g5*GammaL(Gamma::gmu[mu]));
        //multiply lepton part
        res_buf += c * (GammaL(Gamma::gmu[mu]) * lepton);
    }
    buf = peekColour(res_buf, 0, 0);
    sliceSum(buf, r.corr, Tp);
    saveResult(par().output, "weakdecay", r);
    auto &out = envGet(HadronsSerializable, getName());
    out = r;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakMesonDecayKl2_hpp_
