/*
 * Module.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                       ModuleBase implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
ModuleBase::ModuleBase(const std::string name)
: name_(name)
{}

// access //////////////////////////////////////////////////////////////////////
std::string ModuleBase::getName(void) const
{
    return name_;
}

// random seed override for testing ////////////////////////////////////////////
void ModuleBase::seedOverride(const std::string seed)
{
    if (!seed.empty())
    {
        doSeedOverride_ = true;
        seedOverride_ = seed;
    }
    else
    {
        HADRONS_ERROR(Definition, "provided seed is empty");
    }
}

// get factory registration name if available //////////////////////////////////
std::string ModuleBase::getRegisteredName(void)
{
    HADRONS_ERROR(Definition, "module '" + getName() + "' has no registered type"
                 + " in the factory");
}

// result filename generation //////////////////////////////////////////////////
std::string ModuleBase::resultFilename(const std::string stem, 
                                       const unsigned int traj, 
                                       const std::string ext)
{
    return stem + "." + std::to_string(traj) + "." + ext;
}

std::string ModuleBase::resultFilename(const std::string stem, const std::string ext) const
{
    return resultFilename(stem, vm().getTrajectory(), ext);
}

// result database /////////////////////////////////////////////////////////////
void ModuleBase::generateResultDb(void)
{
    if (db_ and db_->isConnected())
    {
        entryHeader_->traj = vm().getTrajectory();
        for (auto filename: getOutputFiles())
        {
            entryHeader_->filename = filename;
            db_->insert(dbTable_, *entry_, true);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
void ModuleBase::operator()(void)
{
    resetTimers();
    startTimer("_total");
    startTimer("_setup");
    setup();
    stopTimer("_setup");
    startTimer("_execute");
    execute();
    stopAllTimers();
    generateResultDb();
}

// get seed ////////////////////////////////////////////////////////////////////
std::string ModuleBase::getSeed(void)
{
    if (doSeedOverride_)
    {
        if (seedOverride_.empty())
        {
            HADRONS_ERROR(Definition, "seed is empty");
        }
        LOG(Warning) << "seed override: DO NOT USE THIS IN PRODUCTION" << std::endl;

        return seedOverride_;
    }
    else
    {
        return makeSeedString();
    }
}

// make module unique string ///////////////////////////////////////////////////
std::string ModuleBase::makeSeedString(void)
{
    std::string seed;

    if (!vm().getRunId().empty())
    {
        seed += vm().getRunId() + "-";
    }
    seed += getName() + "-" + std::to_string(vm().getTrajectory());

    return seed;
}

// get RNGs seeded from module string //////////////////////////////////////////
GridParallelRNG & ModuleBase::rng4d(void)
{
    auto &r = *env().get4dRng();    
    const std::string seed = getSeed();

    if (seed != seed_)
    {
        seed_ = seed;
        LOG(Message) << "Seeding 4D RNG " << &r << " with string '" 
                     << seed_ << "'" << std::endl;
        r.SeedUniqueString(seed_);
    }

    return r;
}

GridSerialRNG & ModuleBase::rngSerial(void)
{
    auto &r = *env().getSerialRng();
    const std::string seed = getSeed();

    if (seed != seed_)
    {
        seed_ = seed;
        LOG(Message) << "Seeding Serial RNG " << &r << " with string '" 
                     << seed_ << "'" << std::endl;
        r.SeedUniqueString(seed_);
    }

    return r;
}
