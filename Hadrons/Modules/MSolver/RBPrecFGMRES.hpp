/*
 * RBPrecFGMRES.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MSolver_RBPrecFGMRES_hpp_
#define Hadrons_MSolver_RBPrecFGMRES_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  Schur red-black preconditioned FGMRES                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecFGMRESPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecFGMRESPar,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    unsigned int, restartLength,
                                    double      , residual,
                                    std::string , guesser);
};

template <typename FImpl>
class TRBPrecFGMRES: public Module<RBPrecFGMRESPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    TRBPrecFGMRES(const std::string name);
    virtual ~TRBPrecFGMRES(void) = default;
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RBPrecFGMRES, ARG(TRBPrecFGMRES<FIMPL>), MSolver);
MODULE_REGISTER_TMP(ZRBPrecFGMRES, ARG(TRBPrecFGMRES<ZFIMPL>), MSolver);

/******************************************************************************
 *                    TRBPrecFGMRES implementation                            *
 ******************************************************************************/
template <typename FImpl>
TRBPrecFGMRES<FImpl>::TRBPrecFGMRES(const std::string name)
: Module<RBPrecFGMRESPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TRBPrecFGMRES<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action};

    if (!par().guesser.empty())
    {
        in.push_back(par().guesser);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TRBPrecFGMRES<FImpl>::getOutput(void)
{
    return {getName(), getName() + "_subtract"};
}

template <typename FImpl>
DependencyMap TRBPrecFGMRES<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName()});
    dep.insert({par().action, getName() + "_subtract"});
    if (!par().guesser.empty())
    {
        dep.insert({par().guesser, getName()});
        dep.insert({par().guesser, getName() + "_subtract"});
    }

    return dep;
}

#define SOLVER_BODY                                                                          \
ZeroGuesser<FermionField>    defaultGuesser;                                                 \
LinearFunction<FermionField> &guesser = (guesserPt == nullptr) ? defaultGuesser : *guesserPt;\
TrivialPrecon<FermionField> identityPreconditioner;                                          \
FlexibleGeneralisedMinimalResidual<FermionField> fgmres(                                    \
    par().residual, par().maxIteration, identityPreconditioner,                              \
    par().restartLength);                                                                    \
HADRONS_DEFAULT_NON_HERMITIAN_SCHUR_SOLVE<FermionField> schurSolver(fgmres);                 \
schurSolver.subtractGuess(subGuess);                                                         \
schurSolver(mat, source, sol, guesser);

template <typename FImpl>
void TRBPrecFGMRES<FImpl>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }
    if (par().restartLength == 0)
    {
        HADRONS_ERROR(Argument, "zero FGMRES restart length");
    }
    if (par().restartLength > par().maxIteration)
    {
        HADRONS_ERROR(Argument, "FGMRES restart length exceeds maximum iteration");
    }
    if (!std::isfinite(par().residual) || (par().residual <= 0.0))
    {
        HADRONS_ERROR(Argument, "residual must be finite and positive");
    }

    LOG(Message) << "Setting up non-Hermitian Schur red-black FGMRES for action '"
                 << par().action << "' with residual " << par().residual
                 << ", maximum iteration " << par().maxIteration
                 << ", restart length " << par().restartLength << std::endl;

    auto                         Ls          = env().getObjectLs(par().action);
    auto                         &mat        = envGet(FMat, par().action);
    LinearFunction<FermionField> *guesserPt = nullptr;

    if (!par().guesser.empty())
    {
        guesserPt = &envGet(LinearFunction<FermionField>, par().guesser);
    }
    auto makeSolver = [&mat, guesserPt, this](bool subGuess)
    {
        return [&mat, guesserPt, subGuess, this]
        (FermionField &sol, const FermionField &source)
        {
            SOLVER_BODY;
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solverSubtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solverSubtract, mat);
}

#undef SOLVER_BODY

template <typename FImpl>
void TRBPrecFGMRES<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecFGMRES_hpp_
