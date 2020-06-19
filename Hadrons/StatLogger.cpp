/*
 * StatLogger.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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

#include <Hadrons/StatLogger.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                         StatLogger implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
StatLogger::StatLogger(Database &db)
{
    setDatabase(db);
}

StatLogger::~StatLogger(void)
{
    if (isRunning())
    {
        stop();
    }
}

void StatLogger::setDatabase(Database &db)
{
    db_ = &db;
    if (db_->isConnected())
    {
        db_->createTable<MemoryEntry>("memory");
        db_->execute(
            "CREATE VIEW IF NOT EXISTS vMemory AS                                                      "
            "SELECT memory.time*1.0e-6 AS timeSec,                                                     "
            "       memory.totalCurrent*0.000000953674316 AS totalCurrentMB,                           "
            "       memory.envCurrent*0.000000953674316 AS envCurrentMB,                               "
            "       memory.gridCurrent*0.000000953674316 AS gridCurrentMB,                             "
            "       memory.commsCurrent*0.000000953674316 AS commsCurrentMB,                           "
            "       (memory.totalCurrent - memory.commsCurrent)*0.000000953674316 AS nocommsCurrentMB, "
            "       memory.totalPeak*0.000000953674316 AS totalPeakMB                                  "
            "FROM memory                                                                               "
            "ORDER BY timeSec;                                                                         "
        );
    }
}

void StatLogger::start(const unsigned int period)
{
    if (isRunning())
    {
        stop();
    }
    if (db_ and db_->isConnected())
    {
        isRunning_.store(true, std::memory_order_release);
        thread_ = std::thread([this, period](void)
        {
            while (isRunning_.load(std::memory_order_acquire))
            {
                auto watch = *GridLogMessage.StopWatch;

                watch.Stop();
                auto time = watch.Elapsed().count();
                watch.Start();
                logMemory(time);
                std::this_thread::sleep_for(std::chrono::milliseconds(period));
            }
        });
    }
}

void StatLogger::stop(void)
{
    if (isRunning())
    {
        isRunning_.store(false, std::memory_order_release);
        thread_.join();
    }
}

bool StatLogger::isRunning(void) const
{
    return (isRunning_.load(std::memory_order_acquire) and thread_.joinable());
}

void StatLogger::logMemory(const GridTime::rep time)
{
    MemoryEntry e;

    e.time         = time;
    e.totalCurrent = MemoryUtils::getCurrentRSS();
    e.envCurrent   = Environment::getInstance().getTotalSize();
    if (Grid::MemoryProfiler::stats)
    {
        e.gridCurrent = Grid::MemoryProfiler::stats->currentlyAllocated;
    }
    else
    {
        e.gridCurrent = 0;
    }
    e.commsCurrent = Grid::GlobalSharedMemory::MAX_MPI_SHM_BYTES;
    e.totalPeak    = MemoryUtils::getPeakRSS();
    if (db_ and db_->isConnected())
    {
        db_->insert("memory", e);
    }
}

void MemoryUtils::printMemory(void)
{
    size_t total, env, comms, peak;

    total = getCurrentRSS();
    env   = Environment::getInstance().getTotalSize();
    comms = Grid::GlobalSharedMemory::MAX_MPI_SHM_BYTES;
    peak  = getPeakRSS();
    LOG(Message) << "Memory: current total " << sizeString(total)
                 << " / environment "        << sizeString(env)
                 << " / comms "              << sizeString(comms);
    if (Grid::MemoryProfiler::stats)
    {
        std::cout << " / grid " 
                  << sizeString(Grid::MemoryProfiler::stats->currentlyAllocated);
    }
    std::cout << " / peak total " << sizeString(peak) << std::endl;
}

/*
 * ALL THE CODE BELOW IS UNDER THE FOLLOWING LICENSE
 * 
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t MemoryUtils::getPeakRSS(void)
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t MemoryUtils::getCurrentRSS(void)
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}
