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
 *   Module to write an arbitrary collection of Grid::Serializables to disk   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

// Placeholder template argument required by MODULE_REGISTER_TMP
class WriteCorrelatorGroupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WriteCorrelatorGroupPar,
                                    std::vector<std::string>, contractions,
                                    std::string, output);
};

template <typename Placeholder>
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

MODULE_REGISTER_TMP(WriteCorrelatorGroup, TWriteCorrelatorGroup<FIMPL>, MIO);

/******************************************************************************
 *                TWriteCorrelators implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Placeholder>
TWriteCorrelatorGroup<Placeholder>::TWriteCorrelatorGroup(const std::string name)
: Module<WriteCorrelatorGroupPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Placeholder>
std::vector<std::string> TWriteCorrelatorGroup<Placeholder>::getInput(void)
{
    return par().contractions;
}

template <typename Placeholder>
std::vector<std::string> TWriteCorrelatorGroup<Placeholder>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename Placeholder>
std::vector<std::string> TWriteCorrelatorGroup<Placeholder>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().output)};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Placeholder>
void TWriteCorrelatorGroup<Placeholder>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Placeholder>
void TWriteCorrelatorGroup<Placeholder>::execute(void)
{
    LOG(Message) << "Preparing to write collated results into file '" << par().output << "'..." << std::endl;
    auto &contractionList = par().contractions; // Switch to std::vector with <elem></elem> tags

    HadronsSerializableGroup result(contractionList.size());
    for (const auto &contractionModuleName : contractionList)
    {
        auto &moduleResults = envGet(HadronsSerializable, contractionModuleName);
        result.append(contractionModuleName, moduleResults);
        LOG(Message) << "Bundled group '" << contractionModuleName << "' into output." << std::endl;
    }

    // Pass a blank string to dump the unpacked group the file, rather than
    // adding a forced and useless outer group around the result
    saveResult(par().output, "", result);
    LOG(Message) << "Finshed writing collated results into file '" << par().output << "'." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_WriteCorrelators_hpp_
