/*
 * StatLogger.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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

#ifndef Hadrons_StatLogger_hpp_
#define Hadrons_StatLogger_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>

#if defined GRID_CUDA and !defined GRID_UVM
#define GRID_CUDA_NOUVM 1
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Real-time statistic logger class                     *
 ******************************************************************************/
class StatLogger
{
public:
    struct MemoryEntry: public SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlNotNull<GridTime::rep>, time,
                           SqlNotNull<size_t>, totalCurrent,
                           SqlNotNull<size_t>, envCurrent,
                           SqlNotNull<size_t>, gridCurrent,
                           SqlNotNull<size_t>, commsCurrent,
                           SqlNotNull<size_t>, totalPeak);
    };
    struct DeviceMemoryEntry: public SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlNotNull<GridTime::rep>, time,
                           SqlNotNull<size_t>, totalCurrent,
                           SqlNotNull<size_t>, envCurrent,
                           SqlNotNull<size_t>, gridCurrent,
                           SqlNotNull<size_t>, gridCacheCurrent,
                           SqlNotNull<size_t>, gridTotalCurrent,
                           SqlNotNull<size_t>, evictableCurrent,
                           SqlNotNull<size_t>, hostToDevice,
                           SqlNotNull<size_t>, hostToDeviceTransfers,
                           SqlNotNull<size_t>, deviceToHost,
                           SqlNotNull<size_t>, deviceToHostTransfers);
    };
public:
    // constructor
    StatLogger(void) = default;
    StatLogger(Database &db, const unsigned int periodMs);
    // destructor
    virtual ~StatLogger(void);
    // set and initialise DB
    void setDatabase(Database &db);
    void setPeriod(const unsigned int periodMs);
    unsigned int getPeriod(void) const;
    // logger control
    void start(void);
    void stop(void);
    bool isRunning(void) const;
private:
    // log memory usage
    void logMemory(const GridTime::rep time);
    void logDeviceMemory(const GridTime::rep time);
private:
    Database          *db_{nullptr};
    std::atomic<bool> isRunning_{false};
    std::thread       thread_;
    unsigned int      periodMs_{0};
};

/******************************************************************************
 *                   Utils to query resident memory from OS                   *
 ******************************************************************************/
namespace MemoryUtils
{
    size_t getHostCurrent(void);
    size_t getHostPeak(void);
#if defined GRID_CUDA
    // size_t getDeviceCurrent(void);

#endif
    void   printMemory(void);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_StatLogger_hpp_
