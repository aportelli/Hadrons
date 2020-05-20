/*
 * MemoryLogger.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MemoryLogger_hpp_
#define Hadrons_MemoryLogger_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>

BEGIN_HADRONS_NAMESPACE

class MemoryLogger
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
public:
    MemoryLogger(void) = default;
    MemoryLogger(Database &db);
    virtual ~MemoryLogger(void);
    static size_t getCurrentRSS(void);
    static size_t getPeakRSS(void);
    static void   print(void);
    void setDatabase(Database &db);
    void start(const unsigned int period);
    void stop(void);
    bool isRunning(void) const;
private:
    Database          *db_{nullptr};
    std::atomic<bool> isRunning_{false};
    std::thread       thread_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_MemoryLogger_hpp_