/*
 * LoadOpenQcd.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
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
#ifndef Hadrons_MIO_LoadOpenQcd_hpp_
#define Hadrons_MIO_LoadOpenQcd_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Load an OpenQcd configuration

 file    Namestem of the files to read in. In contrast to the LoadNersc
         module, the separator between the namestem and the configuration
         number is 'n' instead of '.' compliant with the openQCD convention.         
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadOpenQcdPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadOpenQcdPar,
                                    std::string, file);
};

template <typename GImpl>
class TLoadOpenQcd: public Module<LoadOpenQcdPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TLoadOpenQcd(const std::string name);
    // destructor
    virtual ~TLoadOpenQcd(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadOpenQcd,  TLoadOpenQcd<GIMPL>,  MIO);

/******************************************************************************
*                       TLoadOpenQcd implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadOpenQcd<GImpl>::TLoadOpenQcd(const std::string name)
: Module<LoadOpenQcdPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadOpenQcd<GImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TLoadOpenQcd<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadOpenQcd<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadOpenQcd<GImpl>::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "n"
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading OpenQcd configuration from file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(GaugeField, getName());
    OpenQcdIO::readConfiguration(U, header, fileName);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadOpenQcd_hpp_
