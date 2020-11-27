/*
 * Database.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Database.hpp>

#ifndef HADRONS_SQLITE_MAX_RETRY
#define HADRONS_SQLITE_MAX_RETRY 10
#endif

#ifndef HADRONS_SQLITE_BUSY_TIMEOUT
#define HADRONS_SQLITE_BUSY_TIMEOUT 1000
#endif

using namespace Grid;
using namespace Hadrons;

/******************************************************************************
 *                         QueryResult implementation                         *
 ******************************************************************************/
// row access //////////////////////////////////////////////////////////////////
const std::vector<std::string> & QueryResult::operator[](const unsigned int i) const
{
    return table_.at(i);
}

// column header access ////////////////////////////////////////////////////////
const std::string & QueryResult::colName(const unsigned int j) const
{
    return colName_.at(j);
}

// number of rows and columns //////////////////////////////////////////////////
size_t QueryResult::rows(void) const
{
    return table_.size();
}

size_t QueryResult::cols(void) const
{
    return colName_.size();
}

// broadcast data from boss MPI process ////////////////////////////////////////
void QueryResult::broadcastFromBoss(GridBase *grid)
{
    assert(grid != nullptr);

    int                                 root = grid->BossRank();
    size_t                              r, c, s, k, pos;
    std::vector<std::string::size_type> vs;
    std::string                         cat;

    // broadcast number of rows and columns
    r = rows();
    c = cols();
    grid->Broadcast(root, r);
    grid->Broadcast(root, c);
    if (!grid->IsBoss())
    {
        colName_.resize(c);
        table_.resize(r, std::vector<std::string>(c));
    }

    // pack data
    for (unsigned int j = 0; j < c; ++j)
    {
        cat += colName_[j];
        vs.push_back(colName_[j].size());
    }
    for (unsigned int i = 0; i < r; ++i)
    for (unsigned int j = 0; j < c; ++j)
    {
        cat += table_[i][j];
        vs.push_back(table_[i][j].size());
    }

    // broadcast packed data
    s = cat.size();
    grid->Broadcast(root, s);
    if (!grid->IsBoss())
    {
        cat.resize(s);
    }
    grid->Broadcast(root, &cat[0], s*sizeof(char));
    grid->Broadcast(root, &vs[0], vs.size()*sizeof(std::string::size_type));

    // unpack
    if (!grid->IsBoss())
    {
        pos = 0;
        k   = 0;
        for (unsigned int j = 0; j < c; ++j)
        {
            colName_[j] = cat.substr(pos, vs[k]);
            pos += vs[k];
            k++;
        }
        for (unsigned int i = 0; i < r; ++i)
        for (unsigned int j = 0; j < c; ++j)
        {
            table_[i][j] = cat.substr(pos, vs[k]);
            pos += vs[k];
            k++;
        }
    }
    grid->Barrier();
}

/******************************************************************************
 *                           Database implementation                          *
 ******************************************************************************/
#define BOSS_ONLY if (((grid_ != nullptr) and (grid_->IsBoss())) or (grid_ == nullptr))

// constructor /////////////////////////////////////////////////////////////////
Database::Database(const std::string filename, GridBase *grid, const std::string mode)
{
    if (isConnected() and ((filename != filename_) or (grid != grid_)))
    {
        disconnect();
    }
    setFilename(filename, grid, mode);
}

// destructor //////////////////////////////////////////////////////////////////
Database::~Database(void)
{
    if (isConnected())
    {
        disconnect();
    }
}

// set/get DB filename /////////////////////////////////////////////////////////
void Database::setFilename(const std::string filename, GridBase *grid, const std::string mode)
{
    grid_     = grid;
    filename_ = filename;
    connect();
    if (!mode.empty())
    {
        randomWait(100, grid_);
        execute("PRAGMA journal_mode=" + mode + ";");
#ifndef HADRONS_SQLITE_SYNC
        execute("PRAGMA synchronous = 0;");
#endif
    }
}

std::string Database::getFilename(void) const
{
    return filename_;
}

// test if DB connected ////////////////////////////////////////////////////////
bool Database::isConnected(void) const
{
    return isConnected_;
}

// execute arbitrary SQL statement /////////////////////////////////////////////
#define RETRY_STATUSES ((status == SQLITE_BUSY) or (status == SQLITE_IOERR))

QueryResult Database::execute(const std::string query)
{
    QueryResult result;
    
    BOSS_ONLY
    {
        if (!isConnected())
        {
            HADRONS_ERROR(Database, "no database connected");
        }

        auto callback = [](void *v, int nCol, char **colStr, char **colName)
        {
            std::vector<std::string> line;
            QueryResult              &result = *(static_cast<QueryResult *>(v));

            if (result.colName_.empty())
            {
                for (unsigned int i = 0; i < nCol; ++i)
                {
                    result.colName_.push_back(colName[i]);
                }
            }
            for (unsigned int i = 0; i < nCol; ++i)
            {
                if (colStr[i])
                {
                    line.push_back(colStr[i]);
                }
                else
                {
                    line.push_back("");
                }
            }
            result.table_.push_back(line);

            return SQLITE_OK;
        };

        char *errBuf;
        int  status, attempt = HADRONS_SQLITE_MAX_RETRY;

        do
        {
            status = sqlite3_exec(db_, query.c_str(), callback, &result, &errBuf);
            if ((errBuf != nullptr) and !RETRY_STATUSES)
            {
                std::string errMsg = errBuf;

                sqlite3_free(errBuf);
                HADRONS_ERROR(Database, "error executing query '" + query 
                            + "' in database '" + filename_ + "' (SQLite status " 
                            + std::to_string(status) + ", error '" + errMsg + "')");
                break;
            }
            attempt--;
            if (RETRY_STATUSES)
            {
                LOG(Warning) << "Database '" << filename_ << "' cannot be accessed (SQLite status " 
                             << status << "), randomly retrying in less than 100 ms" << std::endl;
                LOG(Debug) << "Query: '" << query << "'" << std::endl;
                randomWait(100, grid_);
            }
        } while (RETRY_STATUSES and (attempt > 0));
        if (errBuf != nullptr)
        {
            std::string errMsg = errBuf;

            sqlite3_free(errBuf);
            HADRONS_ERROR(Database, "error executing query '" + query 
                        + "' in database '" + filename_ + "' (SQLite status " 
                        + std::to_string(status) + ", error '" + errMsg + "')");
        }
    }
    if (grid_ != nullptr)
    {
        result.broadcastFromBoss(grid_);
    }

    return result;
}

// test if table exists ////////////////////////////////////////////////////////
bool Database::tableExists(const std::string tableName)
{
    QueryResult r = execute("SELECT name FROM sqlite_master WHERE "
                            "type='table' AND name='" + tableName + "';");

    return (r.rows() > 0);
}

// test if table is empty //////////////////////////////////////////////////////
bool Database::tableEmpty(const std::string tableName)
{
    QueryResult r = execute("SELECT COUNT(*) FROM " + tableName + ";");

    return (std::stoi(r[0][0]) == 0);
}

// general tables interface ////////////////////////////////////////////////////
QueryResult Database::getTable(const std::string tableName, const std::string extra)
{
    return execute("SELECT * FROM " + tableName + " " + extra + ";");
}

void Database::insert(const std::string tableName, const SqlEntry &entry, const bool replace)
{
    std::string query;

    query += (replace ? "REPLACE" : "INSERT");
    query += " INTO \"" + tableName + "\" VALUES(";
    query += entry.sqlInsert() + ");";
    execute(query);
}

// key-value tables interface //////////////////////////////////////////////////
void Database::createKeyValueTable(const std::string tableName)
{
    execute("CREATE TABLE " + tableName + " (key TEXT PRIMARY KEY NOT NULL UNIQUE, value BLOB);");
}

std::map<std::string, std::string> Database::getKeyValueTable(const std::string tableName)
{
    std::map<std::string, std::string> map;
    auto                               vec = getTable<KeyValueEntry>(tableName);

    for (auto &e: vec)
    {
        map[e.key] = e.value;
    }

    return map;
}

// private connect/disconnect functions ////////////////////////////////////////
void Database::connect(void)
{
    BOSS_ONLY
    {
        if (!isConnected())
        {
            int status;

            status = sqlite3_open(filename_.c_str(), &db_);
            if (status != SQLITE_OK)
            {
                std::string msg = sqlite3_errmsg(db_);

                HADRONS_ERROR(Io, "cannot connect database in file '" + filename_
                              + "' (SQLite error '" + msg + "')")
            }
            sqlite3_busy_timeout(db_, HADRONS_SQLITE_BUSY_TIMEOUT);
        }
        else
        {
            HADRONS_ERROR(Database, "database already connected");
        }
    }
    isConnected_ = true;
}

void Database::disconnect(void)
{
    BOSS_ONLY
    {
        if (isConnected())
        {
            int status;

            status = sqlite3_close(db_);
            if (status != SQLITE_OK)
            {
                std::string msg = sqlite3_errmsg(db_);

                HADRONS_ERROR(Io, "cannot disconnect database in file '" + filename_
                              + "' (SQLite error '" + msg + "')")
            }
            db_ = nullptr;
        }
        else
        {
            HADRONS_ERROR(Database, "no database connected");
        }
    }
    isConnected_ = false;
}
