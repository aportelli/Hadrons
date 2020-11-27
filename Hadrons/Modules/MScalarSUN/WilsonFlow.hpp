/*
 * WilsonFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Joseph K. L. Lee <joseph.lee@ed.ac.uk>
 * Author: Joseph Lee <joseph.lee@ed.ac.uk>
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
#ifndef Hadrons_MScalarSUN_WilsonFlow_hpp_
#define Hadrons_MScalarSUN_WilsonFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WilsonFlow                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class WilsonFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonFlowPar,
                                    std::string, field,
                                    double, flowtime);
};

class ScalarActionParameters: Serializable 
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarActionParameters,
                                    double, mass_squared,
                                    double, lambda,
                                    double, g);
};

template <typename SImpl>
class TWilsonFlow: public Module<WilsonFlowPar>
{
public:
    typedef typename SImpl::Field                    Field;
    typedef typename SImpl::ComplexField             ComplexField;
    typedef typename SImpl::Group                    Group;
    typedef typename SImpl::SiteField::scalar_object Site;
public:
    // constructor
    TWilsonFlow(const std::string name);
    // destructor
    virtual ~TWilsonFlow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFlowSU2, TWilsonFlow<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(WilsonFlowSU3, TWilsonFlow<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(WilsonFlowSU4, TWilsonFlow<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(WilsonFlowSU5, TWilsonFlow<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(WilsonFlowSU6, TWilsonFlow<ScalarNxNAdjImplR<6>>, MScalarSUN);
/******************************************************************************
 *                 TWilsonFlow implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TWilsonFlow<SImpl>::TWilsonFlow(const std::string name)
: Module<WilsonFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TWilsonFlow<SImpl>::getInput(void)
{
    std::vector<std::string> in= {par().field};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TWilsonFlow<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TWilsonFlow<SImpl>::setup(void)
{
    envTmpLat(Field, "phift");
    envTmpLat(ComplexField, "flowfactor");
    envCreateLat(Field, getName());    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TWilsonFlow<SImpl>::execute(void)
{
    auto    &phi              = envGet(Field, par().field);
    auto    &flowed_phi       = envGet(Field, getName());
    envGetTmp(Field, phift);
    envGetTmp(ComplexField, flowfactor);
    FFT     fft(envGetGrid(Field));

    
    fft.FFT_all_dim(phift, phi, FFT::forward);
    SImpl::MomentaSquare(flowfactor);
    LOG(Message) << "     flowtime = " << par().flowtime << std::endl;
    flowfactor = exp(-par().flowtime*flowfactor);
    phift *= flowfactor;
    fft.FFT_all_dim(flowed_phi, phift, FFT::backward);    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_WilsonFlow_hpp_
