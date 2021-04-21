/*
 * MADWFCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
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
#ifndef Hadrons_MSolver_MADWFCG_hpp_
#define Hadrons_MSolver_MADWFCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              MADWF CG Solver                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MADWFCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MADWFCGPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    unsigned int, maxPVIteration,
                                    double      , innerResidual,
                                    double      , outerResidual,
                                    std::string , eigenPack);
};

template <typename FImplInner, typename FImplOuter, int nBasis>
class TMADWFCG: public Module<MADWFCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
public:
    // constructor
    TMADWFCG(const std::string name);
    // destructor
    virtual ~TMADWFCG(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    struct CGincreaseTol;
};

MODULE_REGISTER_TMP(ZMADWFCG,          ARG(TMADWFCG<ZFIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZMADWFMixedPrecCG, ARG(TMADWFCG<ZFIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP( MADWFCG,          ARG(TMADWFCG< FIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                        TMADWFCG implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
TMADWFCG<FImplInner, FImplOuter, nBasis>
::TMADWFCG(const std::string name)
: Module<MADWFCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, nBasis>
::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, nBasis>
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
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, nBasis>
::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}


template <typename FImplInner, typename FImplOuter, int nBasis>
struct TMADWFCG<FImplInner, FImplOuter, nBasis>
::CGincreaseTol : public MADWFinnerIterCallbackBase {
    ConjugateGradient<FermionFieldInner> &cg_inner;  
    RealD outer_resid;

    CGincreaseTol(ConjugateGradient<FermionFieldInner> &cg_inner,
    RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}

    void operator()(const RealD current_resid){
        LOG(Message) << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
        while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;

        //cg_inner.Tolerance = outer_resid/current_resid;
        LOG(Message) << cg_inner.Tolerance << std::endl;
    }
};


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMADWFCG<FImplInner, FImplOuter, nBasis>
::setup(void)
{
    LOG(Message) << "Setting up MADWF solver " << std::endl;
    LOG(Message) << "with inner/outer action  '"  << par().innerAction   << "'/'" << par().outerAction   << "'" << std::endl;
    LOG(Message) << "     inner/outer residual "  << par().innerResidual <<  "/"  << par().outerResidual        << std::endl;
    LOG(Message) << "     maximum inner/outer/PV iterations " << par().maxInnerIteration <<  "/"  << par().maxOuterIteration <<  "/"  << par().maxPVIteration << std::endl;
    if (!par().eigenPack.empty())
        LOG(Message) << "     eigenpack '" << par().eigenPack << "'" << std::endl;

    auto Ls_outer  = env().getObjectLs(par().outerAction);
    auto &omat     = envGet(FMatOuter, par().outerAction);
    auto guesserPt = makeGuesser<FImplInner, nBasis>(par().eigenPack);

    CayleyFermion5D<FImplOuter> &D_outer = envGetDerived(FMatOuter, CayleyFermion5D<FImplOuter>, par().outerAction);
    CayleyFermion5D<FImplInner> &D_inner = envGetDerived(FMatInner, CayleyFermion5D<FImplInner>, par().innerAction);

    auto makeSolver = [&D_outer, &D_inner, guesserPt, this] (bool subGuess)
    {
        return [&D_outer, &D_inner, guesserPt, subGuess, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source)
        {
            if (subGuess) {
                HADRONS_ERROR(Implementation, "MADWF solver with subtracted guess is not implemented!");
            }

            ConjugateGradient<FermionFieldOuter> CG_PV(par().outerResidual, par().maxPVIteration);
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter> Schur_PV(CG_PV);
            typedef PauliVillarsSolverRBprec<FermionFieldOuter, HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter>> PVtype;
            PVtype PV_outer(Schur_PV);

            ConjugateGradient<FermionFieldInner> CG_inner(par().innerResidual, par().maxInnerIteration, 0);
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldInner> SchurSolver_inner(CG_inner);

            CGincreaseTol update(CG_inner, par().outerResidual);

            MADWF<CayleyFermion5D<FImplOuter>, CayleyFermion5D<FImplInner>,
                  PVtype, HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldInner>, 
                  LinearFunction<FermionFieldInner> >  
                madwf(D_outer, D_inner,
                      PV_outer, SchurSolver_inner,
                      *guesserPt,
                      par().outerResidual, par().maxOuterIteration,
                      &update);

            madwf(source, sol);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls_outer, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls_outer, solver_subtract, omat);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMADWFCG<FImplInner, FImplOuter, nBasis>
::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MADWFCG_hpp_
