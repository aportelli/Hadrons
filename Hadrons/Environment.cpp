/*
 * Environment.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Environment.hpp>
#include <Hadrons/Graph.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

using namespace Grid;
using namespace Hadrons;

#define ERROR_NO_ADDRESS(address)\
HADRONS_ERROR_REF(ObjectDefinition, "no object with address " + std::to_string(address), address);

/******************************************************************************
 *                       Environment implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Environment::Environment(void)
{
    dim_ = GridDefaultLatt().toVector();
    nd_  = dim_.size();
    vol_ = 1.;
    for (auto d: dim_)
    {
        vol_ *= d;
    }
}

// grids ///////////////////////////////////////////////////////////////////////
unsigned int Environment::getNd(void) const
{
    return nd_;
}

std::vector<int> Environment::getDim(void) const
{
    return dim_;
}

int Environment::getDim(const unsigned int mu) const
{
    return dim_[mu];
}

double Environment::getVolume(void) const
{
    return vol_;
}

// random number generator /////////////////////////////////////////////////////
GridParallelRNG * Environment::get4dRng(void)
{
    if (rng4d_ == nullptr)
    {
        rng4d_.reset(new GridParallelRNG(getGrid()));
    }

    return rng4d_.get();
}

GridSerialRNG * Environment::getSerialRng(void)
{
    if (rngSerial_ == nullptr)
    {
        rngSerial_.reset(new GridSerialRNG());
    }

    return rngSerial_.get();
}

// general memory management ///////////////////////////////////////////////////
void Environment::addObject(const std::string name, const int moduleAddress)
{
    if (!hasObject(name))
    {
        ObjInfo info;
        
        info.name   = name;
        info.module = moduleAddress;
        info.data   = nullptr;
        object_.push_back(std::move(info));
        objectAddress_[name] = static_cast<unsigned int>(object_.size() - 1);
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "object '" + name + "' already exists",
                          getObjectAddress(name));
    }
}

void Environment::setObjectStorage(const unsigned int objAddress, 
                                   const Environment::Storage storage)
{
    if (hasObject(objAddress))
    {
        object_[objAddress].storage = storage;
    }
    else
    {
        ERROR_NO_ADDRESS(objAddress);
    }
}

void Environment::setObjectModule(const unsigned int objAddress,
                                  const int modAddress)
{
    if (hasObject(objAddress))
    {
        object_[objAddress].module = modAddress;
    }
    else
    {
        ERROR_NO_ADDRESS(objAddress);
    }
}

void Environment::addObjectDependency(const unsigned int objAddress,
                                      const unsigned int depAddress)
{
    if (objAddress == depAddress)
    {
        HADRONS_ERROR(Logic, "cannot make an object depending on itself");
    }
    if (hasObject(objAddress) and hasObject(depAddress))
    {
        object_[objAddress].dependency.insert(depAddress);
    }
    else if (!hasObject(objAddress))
    {
        ERROR_NO_ADDRESS(objAddress);
    }
    else
    {
        ERROR_NO_ADDRESS(depAddress);
    }
}

void Environment::removeObjectDependency(const unsigned int objAddress,
                                         const unsigned int depAddress)
{
    if (hasObject(objAddress))
    {
        auto &dep = object_[objAddress].dependency;
        auto it   = dep.find(depAddress);

        if (it != dep.end())
        {
            dep.erase(it);
        }
        else
        {
            HADRONS_ERROR(Definition, "object with address "
                          + std::to_string(objAddress) + " does not depend "
                          + "on object with address " + std::to_string(depAddress));
        }
    }
}

const std::set<unsigned int> &Environment::getObjectDependencies(const unsigned int objAddress) const
{
    if (hasObject(objAddress))
    {
        return object_[objAddress].dependency;
    }
    else
    {
        ERROR_NO_ADDRESS(objAddress);
    }
}

bool Environment::hasDependency(const unsigned int objAddress,
                                const unsigned int depAddress) const
{
    if (hasObject(objAddress))
    {
        auto &dep = object_[objAddress].dependency;
        auto it   = dep.find(depAddress);

        return (it != dep.end());
    }
    else
    {
        ERROR_NO_ADDRESS(objAddress);
    }
}

unsigned int Environment::getMaxAddress(void) const
{
    return object_.size();
}

unsigned int Environment::getObjectAddress(const std::string name) const
{
    if (hasObject(name))
    {
        return objectAddress_.at(name);
    }
    else
    {
        HADRONS_ERROR(Definition, "no object with name '" + name + "'");
    }
}

std::string Environment::getObjectName(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].name;
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

std::string Environment::getObjectType(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (object_[address].type)
        {
            return typeName(object_[address].type);
        }
        else
        {
            return "<no type>";
        }
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

std::string Environment::getObjectType(const std::string name) const
{
    return getObjectType(getObjectAddress(name));
}

std::string Environment::getObjectDerivedType(const unsigned int address) const
{
    if (hasObject(address))
    {
        if (object_[address].derivedType)
        {
            return typeName(object_[address].derivedType);
        }
        else
        {
            return "<no type>";
        }
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

std::string Environment::getObjectDerivedType(const std::string name) const
{
    return getObjectDerivedType(getObjectAddress(name));
}

Environment::Size Environment::getObjectSize(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].size;
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

Environment::Size Environment::getObjectSize(const std::string name) const
{
    return getObjectSize(getObjectAddress(name));
}

Environment::Storage Environment::getObjectStorage(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].storage;
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

Environment::Storage Environment::getObjectStorage(const std::string name) const
{
    return getObjectStorage(getObjectAddress(name));
}

int Environment::getObjectModule(const unsigned int address) const
{
    if (hasObject(address))
    {
        return object_[address].module;
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

int Environment::getObjectModule(const std::string name) const
{
    return getObjectModule(getObjectAddress(name));
}

unsigned int Environment::getObjectLs(const unsigned int address) const
{
    if (hasCreatedObject(address))
    {
        return object_[address].Ls;
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

unsigned int Environment::getObjectLs(const std::string name) const
{
    return getObjectLs(getObjectAddress(name));
}

bool Environment::hasObject(const unsigned int address) const
{
    return (address < object_.size());
}

bool Environment::hasObject(const std::string name) const
{
    auto it = objectAddress_.find(name);
    
    return ((it != objectAddress_.end()) and hasObject(it->second));
}

bool Environment::hasCreatedObject(const unsigned int address) const
{
    if (hasObject(address))
    {
        return (object_[address].data != nullptr);
    }
    else
    {
        return false;
    }
}

bool Environment::hasCreatedObject(const std::string name) const
{
    if (hasObject(name))
    {
        return hasCreatedObject(getObjectAddress(name));
    }
    else
    {
        return false;
    }
}

bool Environment::isObject5d(const unsigned int address) const
{
    return (getObjectLs(address) > 1);
}

bool Environment::isObject5d(const std::string name) const
{
    return (getObjectLs(name) > 1);
}

Environment::Size Environment::getTotalSize(void) const
{
    Environment::Size size = 0;
    
    for (auto &o: object_)
    {
        size += o.size;
    }
    
    return size;
}

void Environment::freeObject(const unsigned int address, const bool recursive)
{
    if (hasObject(address))
    {
        if (hasCreatedObject(address))
        {
            LOG(Message) << "Destroying object '" << object_[address].name
                         << "'" << std::endl;
            for (auto &a: object_[address].dependency)
            {
                if (hasCreatedObject(a))
                {
                    if (recursive)
                    {
                        freeObject(a, recursive);
                    }
                    else
                    {
                        HADRONS_ERROR(Logic, "Object '" + object_[address].name + "' cannot be destroyed " +
                                      "because '" + object_[a].name + "' depends on it");
                    }
                }
            }
        }
        object_[address].size = 0;
        object_[address].data.reset(nullptr);
    }
    else
    {
        ERROR_NO_ADDRESS(address);
    }
}

void Environment::freeObject(const std::string name, const bool recursive)
{
    freeObject(getObjectAddress(name), recursive);
}

void Environment::freeSet(const std::set<unsigned int> &objects)
{
    // build object dependency graph from the set
    Graph<unsigned int> graph;

    for (auto &a: objects)
    {
        graph.addVertex(a);
    }
    for (auto &a: objects)
    for (auto &d: getObjectDependencies(a))
    {
        graph.addEdge(d, a);
    }

    // free object in topological order to avoid triggering dependency errors
    std::vector<unsigned int> sorted = graph.topoSort();

    for (auto &a: sorted)
    {
        freeObject(a);
    }
}

void Environment::freeAll(void)
{
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        freeObject(i, true);
    }
}

void Environment::protectObjects(const bool protect)
{
    protect_ = protect;
}

bool Environment::objectsProtected(void) const
{
    return protect_;
}

// print environment content ///////////////////////////////////////////////////
void Environment::printContent(void) const
{
    LOG(Debug) << "Objects: " << std::endl;
    for (unsigned int i = 0; i < object_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getObjectName(i) << " ("
                   << sizeString(getObjectSize(i)) << ")" << std::endl;
    }
}
