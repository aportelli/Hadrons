/*
 * RandomGaugeTransformation.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MGauge_RandomGaugeTransformation_hpp_
#define Hadrons_MGauge_RandomGaugeTransformation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Generate a random gauge transformation                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

template <typename GImpl>
class TRandomGaugeTransformation: public Module<NoPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRandomGaugeTransformation(const std::string name);
    // destructor
    virtual ~TRandomGaugeTransformation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RandomGaugeTransformation, TRandomGaugeTransformation<GIMPL>, MGauge);

/******************************************************************************
 *                 TRandomGaugeTransformation implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TRandomGaugeTransformation<GImpl>::TRandomGaugeTransformation(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TRandomGaugeTransformation<GImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TRandomGaugeTransformation<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TRandomGaugeTransformation<GImpl>::setup(void)
{
    envCreateLat(GaugeLinkField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TRandomGaugeTransformation<GImpl>::execute(void)
{
    LOG(Message) << "Generating random gauge transformation" << std::endl;
    
    auto &g = envGet(GaugeLinkField, getName());

    Group::LieRandomize(rng4d(), g, 1.);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_RandomGaugeTransformation_hpp_
