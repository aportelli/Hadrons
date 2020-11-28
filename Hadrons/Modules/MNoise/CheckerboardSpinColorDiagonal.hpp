/*
 * CheckerboardSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_
#define Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         CheckerboardSpinColorDiagonal                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class CheckerboardSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CheckerboardSpinColorDiagonalPar,
                                    unsigned int, nsrc,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TCheckerboardSpinColorDiagonal: public Module<CheckerboardSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCheckerboardSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TCheckerboardSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CheckerboardSpinColorDiagonal, TCheckerboardSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZCheckerboardSpinColorDiagonal, TCheckerboardSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TCheckerboardSpinColorDiagonal implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCheckerboardSpinColorDiagonal<FImpl>::TCheckerboardSpinColorDiagonal(const std::string name)
: Module<CheckerboardSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCheckerboardSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TCheckerboardSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCheckerboardSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     CheckerboardNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, par().nsparse);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCheckerboardSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    LOG(Message) << "Generating checkerboard spin-color diagonal noise with" 
                 << " nsrc = " << par().nsrc
                 << " and nSparse = " << par().nsparse << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_
