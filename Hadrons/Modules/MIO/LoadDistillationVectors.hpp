/*
 * LoadDistillationVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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
#ifndef Hadrons_MIO_LoadDistillationVectors_hpp_
#define Hadrons_MIO_LoadDistillationVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/DistillationVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadDistillationVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadDistillationVectorsPar,
                                    std::string,  fileStem,
                                    std::string, distilNoise,
                                    std::string, timeSources);
};

template <typename FImpl>
class TLoadDistillationVectors: public Module<LoadDistillationVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadDistillationVectors(const std::string name);
    // destructor
    virtual ~TLoadDistillationVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int nSourceT_;
};

MODULE_REGISTER_TMP(LoadDistillationVectors, TLoadDistillationVectors<FIMPL>, MIO);

/******************************************************************************
 *                      TLoadDistillationVectors implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadDistillationVectors<FImpl>::TLoadDistillationVectors(const std::string name)
: Module<LoadDistillationVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadDistillationVectors<FImpl>::getInput(void)
{
    std::vector<std::string> in={par().distilNoise};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadDistillationVectors<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistillationVectors<FImpl>::setup(void)
{
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    std::string sourceT = par().timeSources;
    if(par().timeSources.empty())
    {
        nSourceT_=nDT;
    }
    else
    {
        // check whether input is legal, i.e. a number of integers between 0 and (nDT-1)
        std::regex rex("[0-9 ]+");
        std::smatch sm;
        std::regex_match(sourceT, sm, rex);
        if (!sm[0].matched)
        {
            HADRONS_ERROR(Range, "sourceTimes must be list of non-negative integers");
        }
        std::istringstream is(sourceT);
        std::vector<int> iT ( ( std::istream_iterator<int>( is )  ), (std::istream_iterator<int>() ) );
        nSourceT_ = iT.size();
        for (int ii = 0; ii < nSourceT_; ii++)
        {
            if (iT[ii] >= nDT)
            {
                HADRONS_ERROR(Range, "elements of sourceTimes must lie between 0 and nDT");
            }
        }
    }
        
    envCreate(std::vector<FermionField>, getName(), 1, nNoise*nDL*nDS*nSourceT_,
                  envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadDistillationVectors<FImpl>::execute(void)
{
    auto &vec = envGet(std::vector<FermionField>, getName());
    
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        //full time dilution size
    std::vector<unsigned int> dt_list = strToVec<unsigned int>(par().timeSources);
    if(par().timeSources.empty())
    {
        dt_list.resize(nDT);
        std::iota(dt_list.begin(), dt_list.end(), 0);
    }

    LOG(Message) << "Loading time sources (dt) : " <<  dt_list << std::endl;

    unsigned int iD=0;
    for (unsigned int D = 0; D < dilNoise.dilutionSize(); ++D)
    {
        unsigned int dt = dilNoise.dilutionCoordinates(D)[DistillationNoise<FImpl>::Index::t];
        if( std::count(dt_list.begin(), dt_list.end(), dt)!=0 ) // dt is in the list of input time sources
        {
            DistillationVectorsIo::readComponent(vec[iD], par().fileStem, nNoise, nDL, nDS, nDT, D, vm().getTrajectory());
            iD++;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadDistillationVectors_hpp_
