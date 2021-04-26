/*
 * SaveNersc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
#ifndef Hadrons_MIO_SaveNersc_hpp_
#define Hadrons_MIO_SaveNersc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Save a NERSC configuration

 gauge        Name of the gauge field object to write
 file         Name of the file to write the gauge field to
 ensebleLabel Label of the ensemble. Recommended this includes
              a suffix identifying this as gauge-fixed and which gauge
 ******************************************************************************/


BEGIN_MODULE_NAMESPACE(MIO)

class SaveNerscPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveNerscPar,
                                    std::string, gauge,
                                    std::string, fileStem,
                                    std::string, ensembleLabel);
};

template <typename GImpl>
class TSaveNersc: public Module<SaveNerscPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TSaveNersc(const std::string name);
    // destructor
    virtual ~TSaveNersc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveNersc,  TSaveNersc<GIMPL>,  MIO);

/******************************************************************************
*                       TSaveNersc implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TSaveNersc<GImpl>::TSaveNersc(const std::string name)
: Module<SaveNerscPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TSaveNersc<GImpl>::getInput(void)
{
  return { par().gauge };
}

template <typename GImpl>
std::vector<std::string> TSaveNersc<GImpl>::getOutput(void)
{
    return {};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveNersc<GImpl>::setup(void)
{
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveNersc<GImpl>::execute(void)
{
    std::string fileName = par().fileStem + "." + std::to_string(vm().getTrajectory());
    LOG(Message) << "Saving NERSC configuration to file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(GaugeField, par().gauge);
    makeFileDir(fileName, U.Grid());
    NerscIO::writeConfiguration(U, fileName, par().ensembleLabel);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveNersc_hpp_
