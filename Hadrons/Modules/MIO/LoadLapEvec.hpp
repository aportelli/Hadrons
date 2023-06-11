/*
 * LoadLapEvec.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 *  Author: Felix Erben <felix.erben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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
#ifndef Hadrons_MIO_LoadLapEvec_hpp_
#define Hadrons_MIO_LoadLapEvec_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadLapEvec                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadLapEvecPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadLapEvecPar,
                                    std::string, LapEvecFileName,
                                    int, nVec);
};

template <typename FImpl>
class TLoadLapEvec: public Module<LoadLapEvecPar>
{
public:
    // constructor
    TLoadLapEvec(const std::string name);
    // destructor
    virtual ~TLoadLapEvec(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadLapEvec, TLoadLapEvec<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadLapEvec implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadLapEvec<FImpl>::TLoadLapEvec(const std::string name)
: Module<LoadLapEvecPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadLapEvec<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadLapEvec<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadLapEvec<FImpl>::setup(void)
{
    envCreate(typename DistillationNoise<FImpl>::LapPack, getName(), 1, par().nVec, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadLapEvec<FImpl>::execute(void)
{
    auto & LapEvec4d = envGet(typename DistillationNoise<FImpl>::LapPack, getName() );
    std::string fileName{ par().LapEvecFileName };
    fileName.append( 1, '.' );
    fileName.append( std::to_string( vm().getTrajectory() ) );
    LapEvec4d.read(fileName.c_str(),false);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadLapEvec_hpp_
