/*
 * ImplicitlyRestartedLanczos.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MSolver_ImplicitlyRestartedLanczos_hpp_
#define Hadrons_MSolver_ImplicitlyRestartedLanczos_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Implicitly Restarted Lanczos module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class ImplicitlyRestartedLanczosPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ImplicitlyRestartedLanczosPar,
                                    LanczosParams, lanczosParams,
                                    std::string,   op,
                                    std::string,   output,
                                    bool,          redBlack,
                                    bool,          multiFile);
};

template <typename Field, typename FieldIo = Field>
class TImplicitlyRestartedLanczos: public Module<ImplicitlyRestartedLanczosPar>
{
public:
    typedef BaseEigenPack<Field>      BasePack;
    typedef EigenPack<Field, FieldIo> Pack;
    typedef LinearOperatorBase<Field> Op;
public:
    // constructor
    TImplicitlyRestartedLanczos(const std::string name);
    // destructor
    virtual ~TImplicitlyRestartedLanczos(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionImplicitlyRestartedLanczos, TImplicitlyRestartedLanczos<FIMPL::FermionField>, MSolver);
MODULE_REGISTER_TMP(FermionImplicitlyRestartedLanczosIo32, ARG(TImplicitlyRestartedLanczos<FIMPL::FermionField, FIMPLF::FermionField>), MSolver);

/******************************************************************************
 *                 TImplicitlyRestartedLanczos implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
TImplicitlyRestartedLanczos<Field, FieldIo>::TImplicitlyRestartedLanczos(const std::string name)
: Module<ImplicitlyRestartedLanczosPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
std::vector<std::string> TImplicitlyRestartedLanczos<Field, FieldIo>::getInput(void)
{
    std::vector<std::string> in = {par().op};
    
    return in;
}

template <typename Field, typename FieldIo>
std::vector<std::string> TImplicitlyRestartedLanczos<Field, FieldIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TImplicitlyRestartedLanczos<Field, FieldIo>::setup(void)
{
    LOG(Message) << "Setting up implicitly restarted Lanczos eigensolver for"
                 << " operator '" << par().op << "' (" << par().lanczosParams.Nstop
                 << " eigenvectors)..." << std::endl;

    GridBase     *grid = nullptr, *gridIo = nullptr;
    unsigned int Ls = env().getObjectLs(par().op);
    auto &op = envGet(Op, par().op);

    grid = getGrid<Field>(par().redBlack, Ls);
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = getGrid<FieldIo>(par().redBlack, Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), Ls, 
                     par().lanczosParams.Nm, grid, gridIo);
    envTmp(Chebyshev<Field>, "cheby", Ls, par().lanczosParams.Cheby);
    envGetTmp(Chebyshev<Field>, cheby);
    envTmp(FunctionHermOp<Field>, "chebyOp", Ls, cheby, op);
    envGetTmp(FunctionHermOp<Field>, chebyOp);
    envTmp(PlainHermOp<Field>, "hermOp", Ls, op);
    envGetTmp(PlainHermOp<Field>, hermOp);
    envTmp(ImplicitlyRestartedLanczos<Field>, "irl", Ls, chebyOp, hermOp, 
        par().lanczosParams.Nstop, par().lanczosParams.Nk, par().lanczosParams.Nm,
        par().lanczosParams.resid, par().lanczosParams.MaxIt, par().lanczosParams.betastp, 
        par().lanczosParams.MinRes);
    envTmp(Field, "src", Ls, grid);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TImplicitlyRestartedLanczos<Field, FieldIo>::execute(void)
{
    int          nconv;
    auto         &epack = envGetDerived(BasePack, Pack, getName());
    GridBase     *grid = nullptr;
    [[maybe_unused]] GridBase *gridIo = nullptr;
    unsigned int Ls = env().getObjectLs(par().op);
    
    envGetTmp(ImplicitlyRestartedLanczos<Field>, irl);
    envGetTmp(Field, src);

    grid = getGrid<Field>(par().redBlack, Ls);
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = getGrid<FieldIo>(par().redBlack, Ls);
    }
    gaussian(rng4d(), src);
    if (par().redBlack)
    {
        src.Checkerboard() = Odd;
    }
    irl.calc(epack.eval, epack.evec, src, nconv, false);
    epack.eval.resize(par().lanczosParams.Nstop);
    epack.evec.resize(par().lanczosParams.Nstop, grid);
    epack.record.operatorXml = vm().getModule(env().getObjectModule(par().op))->parString();
    epack.record.solverXml   = parString();
    if (!par().output.empty())
    {
        epack.write(par().output, par().multiFile, vm().getTrajectory());
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_ImplicitlyRestartedLanczos_hpp_
