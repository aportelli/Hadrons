/*
 * MixedPrecisionRBPrecCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_
#define Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Mixed precision schur red-black preconditioned CG             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MixedPrecisionRBPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedPrecisionRBPrecCGPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    double      , residual,
                                    std::string , innerGuesser,
                                    std::string , outerGuesser);
};

template <typename FImplInner, typename FImplOuter>
class TMixedPrecisionRBPrecCG: public Module<MixedPrecisionRBPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatInner, FermionFieldInner> SchurFMatInner;
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatOuter, FermionFieldOuter> SchurFMatOuter;
private:
    template <typename Field>
    class OperatorFunctionWrapper: public OperatorFunction<Field>
    {
    public:
        using OperatorFunction<Field>::operator();
        OperatorFunctionWrapper(LinearFunction<Field> &fn): fn_(fn) {};
        virtual ~OperatorFunctionWrapper(void) = default;
        virtual void operator()(LinearOperatorBase<Field> &op, 
                                const Field &in, Field &out)
        {
            fn_(in, out);
        }
    private:
        LinearFunction<Field> &fn_;
    };
public:
    // constructor
    TMixedPrecisionRBPrecCG(const std::string name);
    // destructor
    virtual ~TMixedPrecisionRBPrecCG(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MixedPrecisionRBPrecCG, ARG(TMixedPrecisionRBPrecCG<FIMPLF, FIMPLD>), MSolver);
MODULE_REGISTER_TMP(ZMixedPrecisionRBPrecCG, ARG(TMixedPrecisionRBPrecCG<ZFIMPLF, ZFIMPLD>), MSolver);

/******************************************************************************
 *                 TMixedPrecisionRBPrecCG implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::TMixedPrecisionRBPrecCG(const std::string name)
: Module<MixedPrecisionRBPrecCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
std::vector<std::string> TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::getInput(void)
{
    std::vector<std::string> in = {par().innerAction, par().outerAction};
    
    if (!par().innerGuesser.empty())
    {
        in.push_back(par().innerGuesser);
    }
    if (!par().outerGuesser.empty())
    {
        in.push_back(par().outerGuesser);
    }

    return in;
}

template <typename FImplInner, typename FImplOuter>
std::vector<std::string> TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}

template <typename FImplInner, typename FImplOuter>
DependencyMap TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().innerAction, getName()});
    dep.insert({par().innerAction, getName() + "_subtract"});
    dep.insert({par().outerAction, getName()});
    dep.insert({par().outerAction, getName() + "_subtract"});
    if (!par().innerGuesser.empty())
    {
        dep.insert({par().innerGuesser, getName(),             });
        dep.insert({par().innerGuesser, getName() + "_subtract"});
    }
    if (!par().outerGuesser.empty())
    {
        dep.insert({par().outerGuesser, getName(),             });
        dep.insert({par().outerGuesser, getName() + "_subtract"});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                                                                   \
typedef typename FermionFieldInner::vector_type VTypeInner;                                           \
ZeroGuesser<FermionFieldInner> iguesserDefault;                                                       \
ZeroGuesser<FermionFieldOuter> oguesserDefault;                                                       \
LinearFunction<FermionFieldInner> &iguesser = (iguesserPt == nullptr) ? iguesserDefault : *iguesserPt;\
LinearFunction<FermionFieldOuter> &oguesser = (oguesserPt == nullptr) ? oguesserDefault : *oguesserPt;\
SchurFMatInner simat(imat);                                                                           \
SchurFMatOuter somat(omat);                                                                           \
MixedPrecisionConjugateGradient<FermionFieldOuter, FermionFieldInner>                                 \
    mpcg(par().residual, par().maxInnerIteration,                                                     \
         par().maxOuterIteration,                                                                     \
         getGrid<FermionFieldInner>(true, Ls),                                                        \
         simat, somat);                                                                               \
mpcg.useGuesser(iguesser);                                                                            \
OperatorFunctionWrapper<FermionFieldOuter> wmpcg(mpcg);                                               \
HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter> schurSolver(wmpcg);                                    \
schurSolver.subtractGuess(subGuess);                                                                  \
schurSolver(omat, source, sol, oguesser);

template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned mixed-precision "
                 << "CG for inner/outer action '" << par().innerAction 
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration " 
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;

    auto                              Ls          = env().getObjectLs(par().innerAction);
    auto                              &imat       = envGet(FMatInner, par().innerAction);
    auto                              &omat       = envGet(FMatOuter, par().outerAction);
    LinearFunction<FermionFieldInner> *iguesserPt = nullptr; 
    LinearFunction<FermionFieldOuter> *oguesserPt = nullptr;

    if (!par().innerGuesser.empty())
    {
        iguesserPt = &envGet(LinearFunction<FermionFieldInner>, par().innerGuesser);
    }
    if (!par().outerGuesser.empty())
    {
        oguesserPt = &envGet(LinearFunction<FermionFieldOuter>, par().outerGuesser);
    }
    auto makeSolver = [&imat, &omat, iguesserPt, oguesserPt, Ls, this](bool subGuess)
    {
        return [&imat, &omat, iguesserPt, oguesserPt, subGuess, Ls, this]
            (FermionFieldOuter &sol, const FermionFieldOuter &source) 
        {
            SOLVER_BODY;
        };
    };
    auto makeVecSolver = [&imat, &omat, iguesserPt, oguesserPt, Ls, this](bool subGuess)
    {
        return [&imat, &omat, iguesserPt, oguesserPt, subGuess, Ls, this]
            (std::vector<FermionFieldOuter> &sol, const std::vector<FermionFieldOuter> &source) 
        {
            SOLVER_BODY;
        };
    };
    auto solver    = makeSolver(false);
    auto vecSolver = makeVecSolver(false);
    envCreate(Solver, getName(), Ls, solver, vecSolver, omat);
    auto solver_subtract    = makeSolver(true);
    auto vecSolver_subtract = makeVecSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, vecSolver_subtract, omat);
}

#undef SOLVER_BODY

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecCG<FImplInner, FImplOuter>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MixedPrecisionRBPrecCG_hpp_
