/*
 * SaveFieldSharded.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <Ryan.Hill@ed.ac.uk>
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
#ifndef Hadrons_MIO_SaveFieldSharded_hpp_
#define Hadrons_MIO_SaveFieldSharded_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/FieldIoSharded.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SaveFieldSharded                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class SaveFieldShardedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveFieldShardedPar,
                                    std::string, name,
                                    std::string, fileStem);
};

template <typename Field>
class TSaveFieldSharded: public Module<SaveFieldShardedPar>
{
public:
    // constructor
    TSaveFieldSharded(const std::string name);
    // destructor
    virtual ~TSaveFieldSharded(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SavePropagatorSharded, TSaveFieldSharded<FIMPL::PropagatorField>, MIO);
MODULE_REGISTER_TMP(SaveComplexSharded, TSaveFieldSharded<FIMPL::ComplexField>, MIO);

/******************************************************************************
 *                 TSaveFieldSharded implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSaveFieldSharded<Field>::TSaveFieldSharded(const std::string name)
: Module<SaveFieldShardedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSaveFieldSharded<Field>::getInput(void)
{
    std::vector<std::string> in = {par().name};
    
    return in;
}

template <typename Field>
std::vector<std::string> TSaveFieldSharded<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Field>
std::vector<std::string> TSaveFieldSharded<Field>::getOutputFiles(void)
{
    auto grid = env().getGrid();
    std::vector<std::string> filenames;

    filenames.resize(grid->ProcessorCount());
    for (int rank = 0; rank < grid->ProcessorCount(); ++rank)
    {
        filenames[rank] = resultFilename(par().fileStem + ".rank" + std::to_string(rank));
    }
    return filenames;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveFieldSharded<Field>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveFieldSharded<Field>::execute(void)
{
    auto         &field = envGet(Field, par().name);
    GridBase     *grid  = field.Grid();

    LOG(Message) << "Saving sharded field '" << par().name << "' using stem '" << par().fileStem << "'" << std::endl;
    LOG(Message) << "Field type: " << typeName<Field>() << std::endl;
    std::string filepath = resultFilename(par().fileStem + ".rank" + std::to_string(grid->ThisRank()), "shrd");
    writeShardedFile(filepath, field);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveFieldSharded_hpp_
