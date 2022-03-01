/*
 * HDF5Dump.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MIO_HDF5Dump_hpp_
#define Hadrons_MIO_HDF5Dump_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/Modules/MContraction/Meson.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to load a single field from disk                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class HDF5DumpPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HDF5DumpPar,
                                    std::string, name,
                                    std::string, contractions,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class THDF5Dump: public Module<HDF5DumpPar>
{
private:
    typedef typename MContraction::TMeson<FImpl1, FImpl2>::Result MESON_RESULT_T;
public:
    // constructor
    THDF5Dump(const std::string name);
    // destructor
    virtual ~THDF5Dump(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(HDF5Dump, ARG(THDF5Dump<FIMPL, FIMPL>), MIO);

/******************************************************************************
 *                 THDF5Dump implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
THDF5Dump<FImpl1, FImpl2>::THDF5Dump(const std::string name)
: Module<HDF5DumpPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> THDF5Dump<FImpl1, FImpl2>::getInput(void)
{
    auto contractionList = strToVec<std::string>(par().contractions);
    std::vector<std::string> in = contractionList;
    
    return in;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> THDF5Dump<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> THDF5Dump<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().output)};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void THDF5Dump<FImpl1, FImpl2>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void THDF5Dump<FImpl1, FImpl2>::execute(void)
{
    auto contractionList = strToVec<std::string>(par().contractions);

    std::vector<MESON_RESULT_T> result;

    result.clear();
    for (const auto& contractionModuleName : contractionList)
    {
        auto &moduleResults = envGet(std::vector<MESON_RESULT_T>, contractionModuleName);
        for (const auto& resultInstance : moduleResults)
            result.push_back(resultInstance);
    }

    saveResult(par().output, "test", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_HDF5Dump_hpp_
