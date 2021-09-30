/*
 * RBPrecCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSolver_RBPrecCG_hpp_
#define Hadrons_MSolver_RBPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned CG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecCGPar ,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , guesser);
};

template <typename FImpl, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TRBPrecCG: public Module<RBPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRBPrecCG(const std::string name);
    // destructor
    virtual ~TRBPrecCG(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::multimap<std::string, std::string> getObjectDependencies(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RBPrecCG, ARG(TRBPrecCG<FIMPL>), MSolver);
MODULE_REGISTER_TMP(ZRBPrecCG, ARG(TRBPrecCG<ZFIMPL>), MSolver);

/******************************************************************************
 *                      TRBPrecCG template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TRBPrecCG<FImpl, nBasis>::TRBPrecCG(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCG<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    if (!par().guesser.empty())
    {
        in.push_back(par().guesser);
    }

    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCG<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};
    
    return out;
}

template <typename FImpl, int nBasis>
std::multimap<std::string, std::string> TRBPrecCG<FImpl, nBasis>::getObjectDependencies(void)
{
    std::multimap<std::string, std::string> dep;

    dep.insert({par().action, getName()});
    dep.insert({par().action, getName() + "_subtract"});
    if (!par().guesser.empty())
    {
        dep.insert({par().guesser, getName(),             });
        dep.insert({par().guesser, getName() + "_subtract"});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
// C++11 does not support template lambdas so it is easier
// to make a macro with the solver body
#define SOLVER_BODY                                       \
ConjugateGradient<FermionField> cg(par().residual,        \
                                   par().maxIteration);   \
HADRONS_DEFAULT_SCHUR_SOLVE<FermionField> schurSolver(cg);\
schurSolver.subtractGuess(subGuess);                      \
schurSolver(mat, source, sol, guesser);

template <typename FImpl, int nBasis>
void TRBPrecCG<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto Ls       = env().getObjectLs(par().action);
    auto &mat     = envGet(FMat, par().action);
    auto &guesser = envGet(LinearFunction<FermionField>, par().guesser);

    auto makeSolver = [&mat, &guesser, this](bool subGuess) {
        return [&mat, &guesser, subGuess, this]
        (FermionField &sol, const FermionField &source) 
        {
            SOLVER_BODY;
        };
    };
    auto makeVecSolver = [&mat, &guesser, this](bool subGuess) {
        return [&mat, &guesser, subGuess, this]
        (std::vector<FermionField> &sol, const std::vector<FermionField> &source) 
        {
            SOLVER_BODY;
        };
    };
    auto solver    = makeSolver(false);
    auto vecSolver = makeVecSolver(false);
    envCreate(Solver, getName(), Ls, solver, vecSolver, mat);
    auto solver_subtract    = makeSolver(true);
    auto vecSolver_subtract = makeVecSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, vecSolver_subtract, mat);
}

#undef SOLVER_BODY

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecCG<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecCG_hpp_
