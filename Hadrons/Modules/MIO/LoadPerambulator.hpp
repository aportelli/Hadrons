/*
 * LoadPerambulator.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
 * Author: ferben <ferben@localhost.localdomain>
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

#ifndef Hadrons_MIO_LoadPerambulator_hpp_
#define Hadrons_MIO_LoadPerambulator_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                         LoadPerambulator                                 *
 ******************************************************************************/

class LoadPerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorPar,
                                        std::string, perambFileName,
                                        std::string, distilNoise,
                                        std::string, timeSources);
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
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    const int  Nt{env().getDim(Tdir)};
    int nSourceT;
    std::string sourceT = par().timeSources;
    nSourceT = MDistil::verifyTimeSourcesInput(sourceT,nDT);

    envCreate(MDistil::PerambTensor, getName(), 1, Nt, nVec, nDL, nNoise, nSourceT, nDS);
    envTmp(MDistil::PerambIndexTensor, "PerambTmp", 1, Nt, nVec, nDL, nNoise, nDS);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulator<FImpl>::execute(void)
{
    auto &perambulator = envGet(MDistil::PerambTensor, getName());
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nVec = dilNoise.getNl();        
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    const int  Nt{env().getDim(Tdir)};

    std::string sourceT = par().timeSources;
    std::vector<int> invT;   

    perambulator.MetaData.timeSources = invT;
    envGetTmp(MDistil::PerambIndexTensor, PerambTmp);
    for (int dt = 0; dt < Nt; dt++)
    {
        std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
        //skip dilution indices which are not in invT
        if(it == std::end(invT))
        {
            continue;
        }
        LOG(Message) <<  "reading perambulator dt= " << dt << std::endl;
        int idt=it - std::begin(invT);
        std::string sPerambName {par().perambFileName};
        sPerambName.append(".");
        sPerambName.append(std::to_string(vm().getTrajectory()));
        sPerambName.append("/iDT_");
        sPerambName.append(std::to_string(dt));
        sPerambName.append(".");
        sPerambName.append(std::to_string(vm().getTrajectory()));
        PerambTmp.read(sPerambName.c_str());
        for (int t = 0; t < Nt; t++)
        for (int ivec = 0; ivec < nVec; ivec++)
        for (int idl = 0; idl < nDL; idl++)
        for (int in = 0; in < nNoise; in++)
        for (int ids = 0; ids < nDS; ids++)
        {
            perambulator.tensor(t,ivec,idl,in,idt,ids) = PerambTmp.tensor(t,ivec,idl,in,ids);
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
