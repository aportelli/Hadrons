/*
 * TrMag.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
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
#ifndef Hadrons_MScalarSUN_TrMag_hpp_
#define Hadrons_MScalarSUN_TrMag_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Trace of powers of the magnetisation                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TrMagPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrMagPar,
                                    std::string,  field,
                                    unsigned int, maxPow,
                                    std::string,  output);
};

class TrMagResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrMagResult,
                                    std::string, op,
                                    Real,        value);
};

template <typename SImpl>
class TTrMag: public Module<TrMagPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TTrMag(const std::string name);
    // destructor
    virtual ~TTrMag(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TrMagSU2, TTrMag<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TrMagSU3, TTrMag<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TrMagSU4, TTrMag<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TrMagSU5, TTrMag<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TrMagSU6, TTrMag<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                         TTrMag implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTrMag<SImpl>::TTrMag(const std::string name)
: Module<TrMagPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTrMag<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TTrMag<SImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrMag<SImpl>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTrMag<SImpl>::execute(void)
{
    LOG(Message) << "Computing tr(mag^n) for n even up to " << par().maxPow
                 << std::endl;

    std::vector<TrMagResult> result;
    auto                     &phi = envGet(Field, par().field);

    auto m2 = sum(phi);
    auto mn = m2;

    m2 = -m2*m2;
    mn = 1.;
    for (unsigned int n = 2; n <= par().maxPow; n += 2)
    {
        TrMagResult r;

        mn = mn*m2;
        r.op    = "tr(mag^" + std::to_string(n) + ")";
        r.value = TensorRemove(trace(mn)).real();
        result.push_back(r);
    }
    saveResult(par().output, "trmag", result);
    envGet(HadronsSerializable, getName()) = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TrMag_hpp_
