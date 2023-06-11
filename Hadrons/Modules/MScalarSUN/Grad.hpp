/*
 * Grad.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MScalarSUN_Grad_hpp_
#define Hadrons_MScalarSUN_Grad_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Gradient of a complex field                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class GradPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GradPar,
                                    std::string, op,
                                    DiffType,    type,
                                    std::string, output);
};

class GradResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GradResult,
                                    DiffType,              type,
                                    std::vector<Complex>,  value);
};

template <typename SImpl>
class TGrad: public Module<GradPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    TGrad(const std::string name);
    // destructor
    virtual ~TGrad(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GradSU2, TGrad<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(GradSU3, TGrad<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(GradSU4, TGrad<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(GradSU5, TGrad<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(GradSU6, TGrad<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                         TGrad implementation                               *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TGrad<SImpl>::TGrad(const std::string name)
: Module<GradPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TGrad<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().op};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TGrad<SImpl>::getOutput(void)
{
    std::vector<std::string> out;
    const auto               nd = env().getNd();

    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        out.push_back(varName(getName(), mu));
    }
    out.push_back(getName());

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TGrad<SImpl>::setup(void)
{
    const auto nd = env().getNd();

    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        envCreateLat(ComplexField, varName(getName(), mu));
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TGrad<SImpl>::execute(void)
{
    LOG(Message) << "Computing the " << par().type << " gradient of '"
                 << par().op << "'" << std::endl;

    const unsigned int nd = env().getNd();
    GradResult         result;
    auto               &op = envGet(ComplexField, par().op);

    if (!par().output.empty())
    {
        result.type = par().type;
        result.value.resize(nd);
    }
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        auto &der = envGet(ComplexField, varName(getName(), mu));

        dmu(der, op, mu, par().type);
        if (!par().output.empty())
        {
            result.value[mu] = TensorRemove(sum(der));
        }
    }

    saveResult(par().output, "grad", result);
    envGet(HadronsSerializable, getName()) = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Grad_hpp_
