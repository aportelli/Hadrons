/*
 * DeflationSrcCast.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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
#ifndef Hadrons_MGuesser_DeflationSrcCast_hpp_
#define Hadrons_MGuesser_DeflationSrcCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Source cast deflation module                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class DeflationSrcCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DeflationSrcCastPar,
                                    std::string, guesser);
};

template <typename FImplOuter, typename FImplInner>
class TDeflationSrcCast: public Module<DeflationSrcCastPar>
{
public:
    typedef typename FImplOuter::FermionField FieldOuter;
    typedef typename FImplInner::FermionField FieldInner;
public:
    // constructor
    TDeflationSrcCast(const std::string name);
    // destructor
    virtual ~TDeflationSrcCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DeflationSrcCastEPackF, ARG(TDeflationSrcCast<FIMPL, FIMPLF>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImplOuter, typename FImplInner>
class DeflationSrcCastGuesser: public LinearFunction<typename FImplOuter::FermionField>
{
public:
    typedef typename FImplOuter::FermionField FieldOuter;
    typedef typename FImplInner::FermionField FieldInner;

    DeflationSrcCastGuesser(LinearFunction<FieldInner>& guesser, GridBase* gridInner)
    : guesser_(guesser) , gridInner_(gridInner)
    {};

    virtual void operator() (const FieldOuter &in, FieldOuter &out)
    {
        std::vector<FieldOuter> inVec  = {in};
        std::vector<FieldOuter> outVec = {out};

        (*this)(inVec,outVec);

        out = outVec[0];
    }

    virtual void operator() (const std::vector<FieldOuter> &in, std::vector<FieldOuter> &out)
    {
        assert(in.size() == out.size());

        unsigned int sourceSize = out.size();

        for (auto &v: out)
        {
            v = Zero();
        }

        double cast_t = 0.;
        double guess_t = 0.;

        std::vector<FieldInner> inCast( sourceSize, FieldInner(gridInner_) );
        std::vector<FieldInner> outCast(sourceSize, FieldInner(gridInner_) );

        for (unsigned int i = 0; i < sourceSize; ++i) {
            cast_t -= usecond();
            precisionChange(inCast[i],in[i]);
            cast_t += usecond();

            outCast[i] = Zero();
        }

        guess_t -= usecond();
        guesser_(inCast, outCast);
        guess_t += usecond();

        for (unsigned int i = 0; i < sourceSize; ++i) {
            cast_t -= usecond();
            precisionChange(out[i],outCast[i]);
            cast_t += usecond();
        }

        LOG(Message) << "DeflationSrcCastGuesser: Total source casting time " << cast_t/1.e6  << " s" << std::endl;
        LOG(Message) << "DeflationSrcCastGuesser: Total inner guesser time  " << guess_t/1.e6 << " s" <<  std::endl;
    }

private:
    LinearFunction<FieldInner>& guesser_;
    GridBase* gridInner_;
};


/******************************************************************************
 *                     TDeflationSrcCast implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplOuter, typename FImplInner>
TDeflationSrcCast<FImplOuter,FImplInner>::TDeflationSrcCast(const std::string name)
: Module<DeflationSrcCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplOuter, typename FImplInner>
std::vector<std::string> TDeflationSrcCast<FImplOuter,FImplInner>::getInput(void)
{
    std::vector<std::string> in = {par().guesser};
    
    return in;
}

template <typename FImplOuter, typename FImplInner>
std::vector<std::string> TDeflationSrcCast<FImplOuter,FImplInner>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImplOuter, typename FImplInner>
DependencyMap TDeflationSrcCast<FImplOuter,FImplInner>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().guesser, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplOuter, typename FImplInner>
void TDeflationSrcCast<FImplOuter,FImplInner>::setup(void)
{
    LinearFunction<FieldInner>& guesser = envGet(LinearFunction<FieldInner>, par().guesser);

    envCreateDerived(LinearFunction<FieldOuter>, ARG(DeflationSrcCastGuesser<FImplOuter,FImplInner>),
                     getName(), env().getObjectLs(par().guesser), guesser, getGrid<FieldInner>(true, env().getObjectLs(par().guesser)));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplOuter, typename FImplInner>
void TDeflationSrcCast<FImplOuter,FImplInner>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_DeflationSrcCast_hpp_
