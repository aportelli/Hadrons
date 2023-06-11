/*
 * ExactDeflation.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGuesser_ExactDeflation_hpp_
#define Hadrons_MGuesser_ExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Exact deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class ExactDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDeflationPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename EPack>
class TExactDeflation: public Module<ExactDeflationPar>
{
public:
    typedef typename EPack::Field Field;
public:
    // constructor
    TExactDeflation(const std::string name);
    // destructor
    virtual ~TExactDeflation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDeflation, TExactDeflation<BaseFermionEigenPack<FIMPL>>, MGuesser);
MODULE_REGISTER_TMP(ExactDeflationF, TExactDeflation<BaseFermionEigenPack<FIMPLF>>, MGuesser);


/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template<class Field>
class ExactDeflatedGuesser: public LinearFunction<Field> {
private:
  const std::vector<Field> &evec;
  const std::vector<RealD> &eval;
  const unsigned int       epackSize;

public:
    using LinearFunction<Field>::operator();

    ExactDeflatedGuesser(const std::vector<Field>& evec_,const std::vector<RealD>& eval_)
    : ExactDeflatedGuesser(evec_, eval_, evec_.size())
    {}

    ExactDeflatedGuesser(const std::vector<Field>& evec_, const std::vector<RealD>& eval_, unsigned int epackSize_)
    : evec(evec_), eval(eval_), epackSize(epackSize_)
    {
        assert(evec.size()==eval.size());
        assert(epackSize <= evec.size());
    } 

    virtual void operator()(const Field &src,Field &guess) {
        std::vector<Field> srcVec   = {src};
        std::vector<Field> guessVec = {guess};

        (*this)(srcVec,guessVec);

        guess = guessVec[0];
    }

    virtual void operator() (const std::vector<Field> &src, std::vector<Field> &guess)
    {
        assert(src.size() == guess.size());

        unsigned int sourceSize = src.size();

        unsigned int evBatchSize     = epackSize;
        unsigned int sourceBatchSize = sourceSize;
        // These choices of evBatchSize and sourceBatchSize make the loops
        // below trivial. Left machinery in place in case we want to change
        // this in the future.

        for (auto &v: guess)
            v = Zero();

        double time_axpy = 0.;

        for (unsigned int bv = 0; bv < epackSize;  bv += evBatchSize) {
        for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize)
        {
            unsigned int evBlockSize     = std::min(epackSize - bv, evBatchSize);
            unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize);

            time_axpy -= usecond();
            BatchDeflationUtils::projAccumulate(src, guess, evec, eval, 
                                                bv, bv + evBlockSize, 
                                                bs, bs + sourceBlockSize);
            time_axpy += usecond();
        }}

        for (int k=0; k<sourceSize; k++)
            guess[k].Checkerboard() = src[k].Checkerboard();

        LOG(Message) << "ExactDeflatedGuesser: Total axpy time " << time_axpy/1.e6 << " s" <<  std::endl;
    }
};


/******************************************************************************
 *                 TExactDeflation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TExactDeflation<EPack>::TExactDeflation(const std::string name)
: Module<ExactDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TExactDeflation<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename EPack>
std::vector<std::string> TExactDeflation<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename EPack>
DependencyMap TExactDeflation<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflation<EPack>::setup(void)
{
    LOG(Message) << "Setting exact deflation guesser with eigenpack '"
                 << par().eigenPack << "' (" 
                 << par().size << " modes)" << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ExactDeflatedGuesser<Field>, getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflation<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_ExactDeflation_hpp_
