/*
 * SparseSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
#ifndef Hadrons_MNoise_SparseSpinColorDiagonal_hpp_
#define Hadrons_MNoise_SparseSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SparseSpinColorDiagonal                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class SparseSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SparseSpinColorDiagonalPar,
                                    unsigned int, nsrc,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TSparseSpinColorDiagonal: public Module<SparseSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSparseSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TSparseSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SparseSpinColorDiagonal, TSparseSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZSparseSpinColorDiagonal, TSparseSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TSparseSpinColorDiagonal implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSparseSpinColorDiagonal<FImpl>::TSparseSpinColorDiagonal(const std::string name)
: Module<SparseSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSparseSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSparseSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     SparseNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, par().nsparse);    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    LOG(Message) << "Generating sparse spin-color diagonal noise with" 
                 << " nsrc = " << par().nsrc
                 << " and nSparse = " << par().nsparse << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_SparseSpinColorDiagonal_hpp_
