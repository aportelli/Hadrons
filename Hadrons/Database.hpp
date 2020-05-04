/*
 * Database.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_Database_hpp_
#define Hadrons_Database_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/SqlEntry.hpp>
#include <Hadrons/sqlite/sqlite3.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 String table class for SQL query results                   *
 ******************************************************************************/
class QueryResult
{
    friend class Database;
public:
    QueryResult(void) = default;
    virtual ~QueryResult(void) = default;
    const std::vector<std::string> & operator[](const unsigned int i) const;
    const std::string & colName(const unsigned int j) const;
    size_t rows(void) const;
    size_t cols(void) const;
private:
    std::vector<std::string>              colName_;
    std::vector<std::vector<std::string>> table_;
};

/******************************************************************************
 *                          Main database class                               *
 ******************************************************************************/
class Database
{
public:
    Database(const std::string filename, GridBase *grid = nullptr);
    virtual ~Database(void);
    QueryResult execute(const std::string query);
    bool tableExists(const std::string tableName);
    template <typename EntryType>
    void createTable(const std::string tableName, const std::string extra = "");
    void insert(const std::string tableName, const SqlEntry &entry, const bool replace = false);
private:
    void connect(const std::string filename);
    void disconnect(void);
private:
    std::string filename_;
    GridBase    *grid_;
    sqlite3     *db_{nullptr};
};

template <typename EntryType>
void Database::createTable(const std::string tableName, const std::string extra)
{
    std::string query;

    query += "CREATE TABLE \"" + tableName + "\" (" + EntryType::sqlSchema();
    query += (extra.empty() ? "" : "," + extra) + ");";
    execute(query);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Database_hpp_
