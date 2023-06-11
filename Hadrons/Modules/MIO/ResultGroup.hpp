/*
 * ResultGroup.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MIO_ResultGroup_hpp_
#define Hadrons_MIO_ResultGroup_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to load a single field from disk                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class ResultGroupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ResultGroupPar,
                                    std::vector<std::string>, results);
};

// Placeholder template argument required by MODULE_REGISTER_TMP
template<typename Placeholder>
class TResultGroup: public Module<ResultGroupPar>
{
public:
    // constructor
    TResultGroup(const std::string name);
    // destructor
    virtual ~TResultGroup(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ResultGroup, TResultGroup<void>, MIO);

/******************************************************************************
 *                     TResultGroup implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Placeholder>
TResultGroup<Placeholder>::TResultGroup(const std::string name)
: Module<ResultGroupPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template<typename Placeholder>
std::vector<std::string> TResultGroup<Placeholder>::getInput(void)
{
    return par().results;
}

template<typename Placeholder>
std::vector<std::string> TResultGroup<Placeholder>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

template<typename Placeholder>
std::vector<std::string> TResultGroup<Placeholder>::getOutputFiles(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template<typename Placeholder>
void TResultGroup<Placeholder>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template<typename Placeholder>
void TResultGroup<Placeholder>::execute(void)
{
    LOG(Message) << "Starting result collation into Result Group '" << getName() << "'..." << std::endl;
    auto &resultList  = par().results;
    auto &out         = envGet(HadronsSerializable, getName());
    auto &resultGroup = out.template hold<HadronsSerializableGroup>(resultList.size());

    for (const auto &resultName : resultList)
    {
        auto &result = envGet(HadronsSerializable, resultName);
        resultGroup.append(resultName, result);
        LOG(Message) << "Bundled '" << resultName << "' into Result Group '" << getName() << "'." << std::endl;
    }
    LOG(Message) << "Finished collating results into Result Group '" << getName() << "'." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_ResultGroup_hpp_
