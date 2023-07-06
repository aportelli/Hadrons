/*
 * WriteResultGroup.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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

#ifndef Hadrons_MIO_WriteResultGroup_hpp_
#define Hadrons_MIO_WriteResultGroup_hpp_

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
class WriteResultGroupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WriteResultGroupPar,
                                    std::vector<std::string>, results,
                                    std::string, output);
};

template <typename Placeholder>
class TWriteResultGroup: public Module<WriteResultGroupPar>
{
public:
    // constructor
    TWriteResultGroup(const std::string name);
    // destructor
    virtual ~TWriteResultGroup(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WriteResultGroup, TWriteResultGroup<void>, MIO);

/******************************************************************************
 *                TWriteResultGroup implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Placeholder>
TWriteResultGroup<Placeholder>::TWriteResultGroup(const std::string name)
: Module<WriteResultGroupPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Placeholder>
std::vector<std::string> TWriteResultGroup<Placeholder>::getInput(void)
{
    return par().results;
}

template <typename Placeholder>
std::vector<std::string> TWriteResultGroup<Placeholder>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename Placeholder>
std::vector<std::string> TWriteResultGroup<Placeholder>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().output)};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Placeholder>
void TWriteResultGroup<Placeholder>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Placeholder>
void TWriteResultGroup<Placeholder>::execute(void)
{
    LOG(Message) << "Preparing to write collated results into file '" << par().output << "'..." << std::endl;
    auto &resultList = par().results; // Switch to std::vector with <elem></elem> tags

    HadronsSerializableGroup resultGroup(resultList.size());
    for (const auto &resultName : resultList)
    {
        auto &result = envGet(HadronsSerializable, resultName);
        resultGroup.append(resultName, result);
        LOG(Message) << "Bundled group '" << resultName << "' into output." << std::endl;
    }

    // Pass a blank string to dump the unpacked group the file, rather than
    // adding a forced and useless outer group around the result
    saveResult(par().output, "", resultGroup);
    LOG(Message) << "Finished writing collated results into file '" << par().output << "'." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_WriteResultGroup_hpp_
