/*
 * Random.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_MSource_Random_hpp_
#define Hadrons_MSource_Random_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                          Random source                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

template <typename FImpl>
class TRandom: public Module<NoPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandom(const std::string name);
    // destructor
    virtual ~TRandom(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Random,       TRandom<FIMPL>,        MSource);

/******************************************************************************
 *                       TRandom implementation                               *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandom<FImpl>::TRandom(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandom<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandom<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandom<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmpLat(LatticeComplex, "eta");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandom<FImpl>::execute(void)
{
    auto    &src = envGet(PropagatorField, getName());

    envGetTmp(LatticeComplex, eta);
    random(rng4d(), eta);
    src = 1.;
    src = src*eta;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Random_hpp_
