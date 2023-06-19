/*
 * Test_database_concurrency.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Database.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace Hadrons;

struct TestEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(int, point, double, x, double, f);
};

void test(GridBase *grid, const std::string mode)
{
    std::string filename;

    if (!mode.empty())
    {
        LOG(Message) << "Testing DB concurrency using '" << mode << "' journal mode..." << std::endl;
        filename = "test-concurrency-" + mode + ".db";
    }
    else
    {
        LOG(Message) << "Testing DB concurrency using default SQLite journal mode..." << std::endl;
        filename = "test-concurrency.db";
    }
    if (grid->IsBoss())
    {
        remove(filename.c_str());
    }
    grid->Barrier();

    Database  db(filename, nullptr, mode);
    int       rank  = grid->ThisRank(), npoint = 1000;
    double    dx = 0.01, xi = rank*npoint*dx;
    TestEntry e;

    if (grid->IsBoss())
    {
        db.createTable<TestEntry>("function");
    }
    grid->Barrier();
    for (int i = 0; i < npoint; ++i)
    {
        e.point = i + rank*npoint;
        e.x     = xi + i*dx;
        e.f     = cos(e.x);
        db.insert("function", e);
        LOG(Debug) << i << std::endl;
    }
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    initLogger();

    auto                     &env  = Environment::getInstance();
    auto                     *grid = env.getGrid();
    std::vector<std::string> modes = {"WAL"};

    LOG(Message) << "-- Testing different journal modes..." << std::endl;
    for (auto &m: modes)
    {
        try
        {
            test(grid, m);
        }
        catch(const Hadrons::Exceptions::Database &e)
        {
            LOG(Warning) << "DB exception caught!" << std::endl;
            LOG(Warning) << e.what() << std::endl;
        }
    }
    LOG(Message) << "-- Testing Hadrons default journal mode (will crash if anything goes wrong)..." << std::endl;
    test(grid, HADRONS_SQLITE_DEFAULT_JOURNAL_MODE);

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
