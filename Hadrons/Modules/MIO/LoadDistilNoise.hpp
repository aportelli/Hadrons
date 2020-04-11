/*
 * LoadDistilNoise.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_MIO_LoadDistilNoise_hpp_
#define Hadrons_MIO_LoadDistilNoise_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                         LoadDistilNoise                                 *
 ******************************************************************************/

class LoadDistilNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadDistilNoisePar,
                                        std::string, NoiseFileName,
                                        std::string, DistilParams);
};

template <typename FImpl>
class TLoadDistilNoise: public Module<LoadDistilNoisePar>
{
public:
    // constructor
    TLoadDistilNoise(const std::string name);
    // destructor
    virtual ~TLoadDistilNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadDistilNoise, TLoadDistilNoise<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadDistilNoise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadDistilNoise<FImpl>::TLoadDistilNoise(const std::string name) : Module<LoadDistilNoisePar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadDistilNoise<FImpl>::getInput(void)
{
    return {par().DistilParams};
}

template <typename FImpl>
std::vector<std::string> TLoadDistilNoise<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistilNoise<FImpl>::setup(void)
{
    const MDistil::DistilParameters &dp{envGet(MDistil::DistilParameters,  par().DistilParams)};
    const int Nt{env().getDim(Tdir)};
    envCreate(MDistil::NoiseTensor, getName(), 1, dp.nnoise, Nt, dp.nvec, Ns);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistilNoise<FImpl>::execute(void)
{
  auto &noises = envGet(MDistil::NoiseTensor, getName());
  std::string sNoiseName{ par().NoiseFileName };
  sNoiseName.append( 1, '.' );
  sNoiseName.append( std::to_string( vm().getTrajectory() ) );
  noises.read(sNoiseName.c_str());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
