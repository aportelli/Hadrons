/*
 * Exceptions.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Exceptions.hpp>
#include <Hadrons/VirtualMachine.hpp>
#include <Hadrons/Module.hpp>

#ifndef ERR_SUFF
#define ERR_SUFF " (" + loc + ")"
#endif

#define CTOR_EXC(name, init) \
Exceptions::name::name(std::string msg, std::string loc)\
:init\
{}

#define CTOR_EXC_REF(name, init) \
Exceptions::name::name(std::string msg, std::string loc, const unsigned int address)\
:init\
{}

using namespace Grid;
using namespace Hadrons;
using namespace Exceptions;

// backtrace cache
std::vector<std::string> HADRONS_NAMESPACE::Exceptions::backtraceStr;

// logic errors
CTOR_EXC(Logic, logic_error(msg + ERR_SUFF))
CTOR_EXC(Definition, Logic("definition error: " + msg, loc))
CTOR_EXC(Implementation, Logic("implementation error: " + msg, loc))
CTOR_EXC(Range, Logic("range error: " + msg, loc))
CTOR_EXC(Size, Logic("size error: " + msg, loc))

// runtime errors
CTOR_EXC(Runtime, runtime_error(msg + ERR_SUFF))
CTOR_EXC(Argument, Runtime("argument error: " + msg, loc))
CTOR_EXC(Database, Runtime("database error: " + msg, loc))
CTOR_EXC(Io, Runtime("IO error: " + msg, loc))
CTOR_EXC(Memory, Runtime("memory error: " + msg, loc))
CTOR_EXC(Parsing, Runtime("parsing error: " + msg, loc))
CTOR_EXC(Program, Runtime("program error: " + msg, loc))
CTOR_EXC(System, Runtime("system error: " + msg, loc))

// virtual machine errors
CTOR_EXC_REF(ObjectDefinition, RuntimeRef("object definition error: " + msg, loc, address));
CTOR_EXC_REF(ObjectType, RuntimeRef("object type error: " + msg, loc, address));

// abort functions
void HADRONS_NAMESPACE::Exceptions::abort(const std::exception& e)
{
    int mod;
    std::string modName;

    // suppress potential extra expections
    try
    {
        auto &vm = VirtualMachine::getInstance();
        mod = vm.getCurrentModule();
        if (mod >= 0)
        {
            modName = vm.getModuleName(mod);
        }
    }
    catch(...)
    {
        mod = -1;
        modName = "";
    }
    LOG(Error) << "FATAL ERROR -- Exception " << typeName(&typeid(e)) 
               << std::endl;
    if (mod >= 0)
    {
        LOG(Error) << "During execution of module '"
                    << modName << "' (address " << mod << ")"
                    << std::endl;
    }
    LOG(Error) << e.what() << std::endl;
    if (!backtraceStr.empty())
    {
        LOG(Error) << "-- BACKTRACE --------------" << std::endl;
        for (auto &s: backtraceStr)
        {
            LOG(Error) << s << std::endl;
        }
        LOG(Error) << "---------------------------" << std::endl;
    }
    LOG(Error) << "Aborting program" << std::endl;
#if defined (GRID_COMMS_MPI) || defined (GRID_COMMS_MPI3) || defined (GRID_COMMS_MPIT)
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(EXIT_FAILURE);
#endif
}
