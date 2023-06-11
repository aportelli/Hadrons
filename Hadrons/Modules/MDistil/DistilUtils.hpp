/*
 * DistilUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Felix Erben <ferben@ed.ac.uk>
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

#ifndef Hadrons_MDistil_DistilUtils_hpp_
#define Hadrons_MDistil_DistilUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

inline int verifyTimeSourcesInput(std::string sourceT, int nDT)
{
    int nSourceT=0;
    if(sourceT.empty())
    {
        nSourceT=nDT;
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
        nSourceT = iT.size();
        for (int ii = 0; ii < nSourceT; ii++)
        {
            if (iT[ii] >= nDT)
            {
                HADRONS_ERROR(Range, "elements of sourceTimes must lie between 0 and nDT");
            }
        }
        // also check whether input is sorted
        if(!std::is_sorted(iT.begin(), iT.end()))
        {
                HADRONS_ERROR(Range, "elements must be in ascending order");

        }
    }
    return nSourceT;
}

template <typename FImpl>
inline int getSourceTimesFromInput(std::string & sourceT, 
                                   int nDT, 
                                   DistillationNoise<FImpl> & dilNoise,
                                   std::vector<int> & invT)
{
        typedef typename DistillationNoise<FImpl>::Index Index;
    int nSourceT=0;
    std::vector<std::vector<unsigned int>> sourceTimes;
    if(sourceT.empty())
    {
        // create sourceTimes all time-dilution indices
        nSourceT=nDT;
        for (int dt = 0; dt < nDT; dt++)
        {
            std::vector<unsigned int> sT = dilNoise.dilutionPartition(Index::t, dt);
            sourceTimes.push_back(sT);
            invT.push_back(dt);
        }
        LOG(Message) << "Perambulator for all " << nDT << " time-dilution vectors" << std::endl;
    }
    else
    {
        std::istringstream is(sourceT);
        std::vector<int> iT ( ( std::istream_iterator<int>( is )  ), (std::istream_iterator<int>() ) );
        nSourceT = iT.size();
        // create sourceTimes from the chosen subset of time-dilution indices
        for (int dt = 0; dt < nSourceT; dt++)
        {
            std::vector<unsigned int> sT = dilNoise.dilutionPartition(Index::t, iT[dt]);
            sourceTimes.push_back(sT);
            invT.push_back(iT[dt]);
        }
        LOG(Message) << "Perambulator for a subset of " << nSourceT << " time-dilution vectors" << std::endl;
    }
    LOG(Message) << "Source times" << sourceTimes << std::endl;
    return nSourceT;
}

// auxiliar function to print time source logs
inline std::string timeslicesDump(const std::vector<unsigned int> ts)
{
    std::stringstream ss;
    for (auto& t : ts)
        ss << t << " ";
    return ss.str();
};

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
