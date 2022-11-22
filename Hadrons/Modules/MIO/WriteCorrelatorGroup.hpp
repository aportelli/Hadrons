/*
 * WriteCorrelators.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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

#ifndef Hadrons_MIO_WriteCorrelators_hpp_
#define Hadrons_MIO_WriteCorrelators_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to load a single field from disk                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class WriteCorrelatorGroupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WriteCorrelatorGroupPar,
                                    std::vector<std::string>, contractions,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TWriteCorrelatorGroup: public Module<WriteCorrelatorGroupPar>
{
public:
    // constructor
    TWriteCorrelatorGroup(const std::string name);
    // destructor
    virtual ~TWriteCorrelatorGroup(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WriteCorrelatorGroup, ARG(TWriteCorrelatorGroup<FIMPL, FIMPL>), MIO);

/******************************************************************************
 *                TWriteCorrelators implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TWriteCorrelatorGroup<FImpl1, FImpl2>::TWriteCorrelatorGroup(const std::string name)
: Module<WriteCorrelatorGroupPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TWriteCorrelatorGroup<FImpl1, FImpl2>::getInput(void)
{
    return par().contractions;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TWriteCorrelatorGroup<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TWriteCorrelatorGroup<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().output)};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TWriteCorrelatorGroup<FImpl1, FImpl2>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TWriteCorrelatorGroup<FImpl1, FImpl2>::execute(void)
{
    auto &contractionList = par().contractions; // Switch to std::vector with <elem></elem> tags

    HadronsSerializableGroup result(contractionList.size());
    for (const auto &contractionModuleName : contractionList)
    {
        auto &moduleResults = envGet(HadronsSerializable, contractionModuleName);
        result.append(contractionModuleName, moduleResults);
    }

    // Pass a blank string to dump the unpacked group the file, rather than
    // adding a forced and useless outer group around the result
    saveResult(par().output, "", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_WriteCorrelators_hpp_
