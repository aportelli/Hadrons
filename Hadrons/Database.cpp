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

using namespace Grid;
using namespace Hadrons;

#define BOSS_ONLY if (((grid_ != nullptr) and (grid_->IsBoss())) or (grid_ == nullptr))

const std::vector<std::string> & QueryResult::operator[](const unsigned int i) const
{
    return table_.at(i);
}

const std::string & QueryResult::colName(const unsigned int j) const
{
    return colName_.at(j);
}

size_t QueryResult::rows(void) const
{
    return table_.size();
}

size_t QueryResult::cols(void) const
{
    return colName_.size();
}

Database::Database(const std::string filename, GridBase *grid)
: grid_(grid)
{
    if (!filename.empty())
    {
        connect(filename);
    }   
}

Database::~Database(void)
{
    if (db_ != nullptr)
    {
        disconnect();
    }
}

void Database::connect(const std::string filename)
{
    filename_ = filename;

    BOSS_ONLY
    {
        if (db_ == nullptr)
        {
            int status;

            status = sqlite3_open(filename_.c_str(), &db_);
            if (status != SQLITE_OK)
            {
                std::string msg = sqlite3_errmsg(db_);

                HADRONS_ERROR(Io, "cannot connect database in file '" + filename_
                              + "' (SQLite error '" + msg + "')")
            }
        }
        else
        {
            HADRONS_ERROR(Database, "database already connected");
        }
    }
}

void Database::disconnect(void)
{
    BOSS_ONLY
    {
        if (db_ != nullptr)
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
}

QueryResult Database::execute(const std::string query)
{
    QueryResult result;
    
    BOSS_ONLY
    {
        if (db_ == nullptr)
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
                line.push_back(colStr[i]);
            }
            result.table_.push_back(line);

            return SQLITE_OK;
        };

        char *errBuf;

        sqlite3_exec(db_, query.c_str(), callback, &result, &errBuf);
        if (errBuf != nullptr)
        {
            std::string errMsg = errBuf;

            sqlite3_free(errBuf);
            HADRONS_ERROR(Database, "error executing query '" + query 
                          + "' (SQLite error '" + errMsg + "')");
        }
    }

    return result;
}
