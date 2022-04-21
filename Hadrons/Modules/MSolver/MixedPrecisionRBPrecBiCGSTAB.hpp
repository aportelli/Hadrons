/*
 * MixedPrecisionRBPrecBiCGSTAB.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
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
#ifndef Hadrons_MSolver_MixedPrecisionRBPrecBiCGSTAB_hpp_
#define Hadrons_MSolver_MixedPrecisionRBPrecBiCGSTAB_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Mixed precision schur red-black preconditioned BiCGSTAB             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MixedPrecisionRBPrecBiCGSTABPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedPrecisionRBPrecBiCGSTABPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    double      , residual,
                                    std::string , eigenPack);
};

template <typename FImplInner, typename FImplOuter, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TMixedPrecisionRBPrecBiCGSTAB: public Module<MixedPrecisionRBPrecBiCGSTABPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    typedef HADRONS_DEFAULT_NON_HERMITIAN_SCHUR_OP<FMatInner, FermionFieldInner> SchurFMatInner;
    typedef HADRONS_DEFAULT_NON_HERMITIAN_SCHUR_OP<FMatOuter, FermionFieldOuter> SchurFMatOuter;
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
    TMixedPrecisionRBPrecBiCGSTAB(const std::string name);
    // destructor
    virtual ~TMixedPrecisionRBPrecBiCGSTAB(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MixedPrecisionRBPrecBiCGSTAB, ARG(TMixedPrecisionRBPrecBiCGSTAB<FIMPLF, FIMPLD>), MSolver);
MODULE_REGISTER_TMP(ZMixedPrecisionRBPrecBiCGSTAB, ARG(TMixedPrecisionRBPrecBiCGSTAB<ZFIMPLF, ZFIMPLD>), MSolver);

/******************************************************************************
 *                 TMixedPrecisionRBPrecBiCGSTAB implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::TMixedPrecisionRBPrecBiCGSTAB(const std::string name)
: Module<MixedPrecisionRBPrecBiCGSTABPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::getReference(void)
{
    std::vector<std::string> ref = {par().innerAction, par().outerAction};

    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned mixed-precision "
                 << "BiCGSTAB for inner/outer action '" << par().innerAction
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration "
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;

    auto Ls        = env().getObjectLs(par().innerAction);
    auto &imat     = envGet(FMatInner, par().innerAction);
    auto &omat     = envGet(FMatOuter, par().outerAction);

    auto guesserPt64 = makeGuesser<FImplOuter, nBasis>("");
    auto guesserPt32 = makeGuesser<FImplInner, nBasis>("");

    try
    {
        guesserPt64 = makeGuesser<FImplOuter, nBasis>(par().eigenPack);
    }
    catch (Exceptions::ObjectType &e)
    {
        guesserPt32 = makeGuesser<FImplInner, nBasis>(par().eigenPack);
    }

    auto makeSolver = [&imat, &omat, guesserPt32, guesserPt64, Ls, this](bool subGuess)
    {
        return [&imat, &omat, guesserPt32, guesserPt64, subGuess, Ls, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source)
        {
            SchurFMatInner simat(imat);
            SchurFMatOuter somat(omat);
            MixedPrecisionBiCGSTAB<FermionFieldOuter, FermionFieldInner>
                mpcg(par().residual, par().maxInnerIteration,
                     par().maxOuterIteration,
                     getGrid<FermionFieldInner>(true, Ls),
                     simat, somat);
                mpcg.useGuesser(*guesserPt32);
            OperatorFunctionWrapper<FermionFieldOuter> wmpcg(mpcg);
            HADRONS_DEFAULT_NON_HERMITIAN_SCHUR_SOLVE<FermionFieldOuter> schurSolver(wmpcg);
            schurSolver.subtractGuess(subGuess);
            schurSolver(omat, source, sol, *guesserPt64);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, omat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMixedPrecisionRBPrecBiCGSTAB<FImplInner, FImplOuter, nBasis>
::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MixedPrecisionRBPrecBiCGSTAB_hpp_
