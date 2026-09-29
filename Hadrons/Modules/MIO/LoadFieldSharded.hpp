/*
 * LoadFieldSharded.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadFieldSharded_hpp_
#define Hadrons_MIO_LoadFieldSharded_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/FieldIoSharded.hpp>


BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadFieldSharded                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadFieldShardedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadFieldShardedPar,
                                    std::string, name,
                                    unsigned int, Ls,
                                    std::string, fileStem);
};

template <typename FImpl>
class TLoadFieldSharded: public Module<LoadFieldShardedPar>
{
public:
    // constructor
    TLoadFieldSharded(const std::string name);
    // destructor
    virtual ~TLoadFieldSharded(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadPropagatorSharded, TLoadFieldSharded<FIMPL::PropagatorField>, MIO);
MODULE_REGISTER_TMP(LoadComplexSharded, TLoadFieldSharded<FIMPL::ComplexField>, MIO);

/******************************************************************************
 *                 TLoadFieldSharded implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TLoadFieldSharded<Field>::TLoadFieldSharded(const std::string name)
: Module<LoadFieldShardedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TLoadFieldSharded<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TLoadFieldSharded<Field>::getOutput(void)
{
    std::vector<std::string> out = {par().name};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TLoadFieldSharded<Field>::setup(void)
{
    GridBase *grid = nullptr;

    if (par().Ls > 1)
    {
        grid = envGetGrid(Field, par().Ls);
        envCreateLat(Field, par().name, par().Ls);
    }
    else
    {
        grid = envGetGrid(Field);
        envCreateLat(Field, par().name);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TLoadFieldSharded<Field>::execute(void)
{
    auto         &field = envGet(Field, par().name);
    GridBase     *grid  = field.Grid();

    LOG(Message) << "Loading sharded field '" << par().name << "' using stem '" << par().fileStem << "'" << std::endl;
    LOG(Message) << "Field type: " << typeName<Field>() << std::endl;
    std::string filepath = resultFilename(par().fileStem + ".rank" + std::to_string(grid->ThisRank()), "shrd");
    readShardedFile(filepath, field);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadFieldSharded_hpp_
