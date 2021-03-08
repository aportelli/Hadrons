/*
 * LoadPerambulator.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
 * Author: ferben <ferben@localhost.localdomain>
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

#ifndef Hadrons_MIO_LoadPerambulator_hpp_
#define Hadrons_MIO_LoadPerambulator_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                         LoadPerambulator                                 *
 ******************************************************************************/

class LoadPerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorPar,
                                        std::string, PerambFileName,
                                        std::string, distilNoise,
					int, inversions);
};

template <typename FImpl>
class TLoadPerambulator: public Module<LoadPerambulatorPar>
{
public:
    // constructor
    TLoadPerambulator(const std::string name);
    // destructor
    virtual ~TLoadPerambulator(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadPerambulator, TLoadPerambulator<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadPerambulator<FImpl>::TLoadPerambulator(const std::string name) : Module<LoadPerambulatorPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getInput(void)
{
    return {par().distilNoise};
}

template <typename FImpl>
std::vector<std::string> TLoadPerambulator<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::setup(void)
{
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();	
    int nVec = dilNoise.getNl();	
    int LI = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);	
    int SI = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);	
    const int  Nt{env().getDim(Tdir)};
  
    envCreate(MDistil::PerambTensor, getName(), 1, Nt, nVec, LI, nNoise, par().inversions, SI);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::execute(void)
{
  auto &perambulator = envGet(MDistil::PerambTensor, getName());
  std::string sPerambName{ par().PerambFileName };
  sPerambName.append( 1, '.' );
  sPerambName.append( std::to_string( vm().getTrajectory() ) );
  perambulator.read(sPerambName.c_str());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
