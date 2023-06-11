/*
 * VectorPackRef.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_VectorPackRef_hpp_
#define Hadrons_MUtilities_VectorPackRef_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Pack fields as a vector of pointers                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class VectorPackRefPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorPackRefPar,
                                    std::vector<std::string>, fields);
};

template <typename Field>
class TVectorPackRef: public Module<VectorPackRefPar>
{
public:
    // constructor
    TVectorPackRef(const std::string name);
    // destructor
    virtual ~TVectorPackRef(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PropagatorVectorPackRef, TVectorPackRef<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                 TVectorPackRef implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TVectorPackRef<Field>::TVectorPackRef(const std::string name)
: Module<VectorPackRefPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TVectorPackRef<Field>::getInput(void)
{
    std::vector<std::string> in = par().fields;
    
    return in;
}

template <typename Field>
std::vector<std::string> TVectorPackRef<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Field>
DependencyMap TVectorPackRef<Field>::getObjectDependencies(void)
{
    DependencyMap dep;

    for (auto &n: par().fields)
    {
        dep.insert({n, getName()});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRef<Field>::setup(void)
{
    envCreate(std::vector<Field *>, getName(), 1, par().fields.size(), nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRef<Field>::execute(void)
{
    LOG(Message) << "Packing " << par().fields.size() << " field pointers" << std::endl;
    LOG(Message) << "fields: " << par().fields << std::endl;
    auto &vec = envGet(std::vector<Field *>, getName());

    for (unsigned int i = 0; i < par().fields.size(); ++i)
    {
        vec[i] = &envGet(Field, par().fields[i]);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_VectorPackRef_hpp_
