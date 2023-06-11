/*
 * CGNE.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MSolver_CGNE_hpp_
#define Hadrons_MSolver_CGNE_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        CG on normal equations (CGNE)                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class CGNEPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CGNEPar,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , guesser);
};

template <typename FImpl>
class TCGNE: public Module<CGNEPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCGNE(const std::string name);
    // destructor
    virtual ~TCGNE(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CGNE, TCGNE<FIMPL>, MSolver);

/******************************************************************************
 *                           TCGNE implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCGNE<FImpl>::TCGNE(const std::string name)
: Module<CGNEPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCGNE<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    if (!par().guesser.empty())
    {
        in.push_back(par().guesser);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TCGNE<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl>
DependencyMap TCGNE<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

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
template <typename FImpl>
void TCGNE<FImpl>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up normal equation CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;
    
    auto Ls                                 = env().getObjectLs(par().action);
    auto &mat                               = envGet(FMat, par().action);
    LinearFunction<FermionField> *guesserPt = nullptr;
    
    if (!par().guesser.empty())
    {
        guesserPt = &envGet(LinearFunction<FermionField>, par().guesser);
    }
    auto makeSolver = [&mat, guesserPt, this](bool subGuess) 
    {
        return [&mat, guesserPt, subGuess, this](FermionField &sol,
                                                 const FermionField &source) 
        {
            GridBase                                *g = sol.Grid();
            FermionField                            guess(g), tmp(g);
            MdagMLinearOperator<FMat, FermionField> hermOp(mat);
            ConjugateGradient<FermionField>         cg(par().residual,
                                                       par().maxIteration);
            ZeroGuesser<FermionField>               defaultGuesser;

            guess = sol;
            mat.Mdag(source, tmp);
            if (guesserPt != nullptr)
            {
                (*guesserPt)(tmp, sol);
            }
            else
            {
                defaultGuesser(tmp, sol);
            }
            cg(hermOp, tmp, sol);
            if (subGuess)
            {
                sol -= guess;
            }
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCGNE<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_CGNE_hpp_
