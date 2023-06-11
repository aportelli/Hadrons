/*
 * JacobiSmear.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
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
#ifndef Hadrons_MSource_JacobiSmear_hpp_
#define Hadrons_MSource_JacobiSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         JacobiSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class JacobiSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(JacobiSmearPar,
                                    std::string, gauge,
                                    double, width,
                                    int, iterations,
                                    int, orthog,
                                    std::string, source);
};

template <typename FImpl>
class TJacobiSmear: public Module<JacobiSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TJacobiSmear(const std::string name);
    // destructor
    virtual ~TJacobiSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(JacobiSmear, TJacobiSmear<FIMPL>, MSource);

/******************************************************************************
 *                 TJacobiSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TJacobiSmear<FImpl>::TJacobiSmear(const std::string name)
: Module<JacobiSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmp(std::vector<GaugeMat>, "Umu", 1, 4, envGetGrid(LatticeColourMatrix));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::execute(void)
{
    auto &out = envGet(PropagatorField, getName());
    auto &src = envGet(PropagatorField, par().source);
    auto &U = envGet(GaugeField, par().gauge);
    envGetTmp(std::vector<GaugeMat>, Umu);
    for(int mu=0; mu<4; mu++)
    {
       Umu.at(mu)=peekLorentz(U,mu);
    }
    CovariantSmearing<FImpl> covsmear;
    out=src;
    startTimer("Jacobi iteration");
    covsmear.GaussianSmear(Umu, out, par().width, par().iterations, par().orthog);
    stopTimer("Jacobi iteration");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_JacobiSmear_hpp_
