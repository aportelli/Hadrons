/*
 * LinearComb.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@desy.de>
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
#ifndef Hadrons_MUtilities_LinearComb_hpp_
#define Hadrons_MUtilities_LinearComb_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LinearComb                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class LinearCombPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LinearCombPar,
                                    std::vector<std::string>, fields,
                                    std::vector<double>,      coeffs);
};

template <typename FImpl>
class TLinearComb: public Module<LinearCombPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLinearComb(const std::string name);
    // destructor
    virtual ~TLinearComb(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LinearComb, TLinearComb<FIMPL>, MUtilities);

/******************************************************************************
 *                 TLinearComb implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLinearComb<FImpl>::TLinearComb(const std::string name)
: Module<LinearCombPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLinearComb<FImpl>::getInput(void)
{
    std::vector<std::string> in = par().fields;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLinearComb<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLinearComb<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLinearComb<FImpl>::execute(void)
{
    assert(par().fields.size() == par().coeffs.size());

    LOG(Message) << par().fields << std::endl;
    LOG(Message) << par().coeffs << std::endl;

    auto &res  = envGet(PropagatorField, getName());

    for (int i=0; i<par().fields.size(); i++) {
        auto &q = envGet(PropagatorField, par().fields[i]);
        if (i==0)
            res = par().coeffs[i]*q;
        else
            res += par().coeffs[i]*q;
        LOG(Message) << "Field[" << i << "] : Norm = " << norm2(q) << std::endl;  
    }
    LOG(Message) << std::endl;
    LOG(Message) << "LinearComb : Norm = " << norm2(res) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_LinearComb_hpp_