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
    void broadcastFromBoss(GridBase *grid);
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
    Database(void) = default;
    Database(const std::string filename, GridBase *grid = nullptr);
    virtual ~Database(void);
    void setFilename(const std::string filename, GridBase *grid = nullptr);
    bool isConnected(void) const;
    QueryResult execute(const std::string query);
    bool tableExists(const std::string tableName);
    template <typename EntryType>
    void createTable(const std::string tableName, const std::string extra = "");
    template <typename EntryType>
    std::vector<EntryType> getTable(const std::string tableName, const std::string extra = "");
    template <typename ColType>
    std::vector<ColType> getTableColumn(const std::string tableName, const std::string columnName, const std::string extra = "");
    void insert(const std::string tableName, const SqlEntry &entry, const bool replace = false);
private:
    void connect(void);
    void disconnect(void);
private:
    std::string filename_;
    GridBase    *grid_{nullptr};
    sqlite3     *db_{nullptr};
    bool        isConnected_{false};
};

template <typename EntryType>
void Database::createTable(const std::string tableName, const std::string extra)
{
    std::string query;

    query += "CREATE TABLE " + tableName + " (" + EntryType::sqlSchema();
    query += (extra.empty() ? "" : "," + extra) + ");";
    execute(query);
}

template <typename EntryType>
std::vector<EntryType> Database::getTable(const std::string tableName, 
                                          const std::string extra)
{
    std::vector<EntryType> table;
    EntryType              buf;
    QueryResult            r = execute("SELECT * FROM " + tableName + 
                                       + " " + extra + ";");

    for (unsigned int i = 0; i < r.rows(); ++i)
    {
        buf.deserializeRow(r[i]);
        table.push_back(buf);
    }

    return table;
}

template <typename ColType>
std::vector<ColType> Database::getTableColumn(const std::string tableName, 
                                              const std::string columnName, 
                                              const std::string extra)
{
    std::vector<ColType> column;
    ColType              buf;
    QueryResult          r = execute("SELECT " + columnName + " FROM " + tableName + 
                                     + " " + extra + ";");

    for (unsigned int i = 0; i < r.rows(); ++i)
    {
        column.push_back(SqlEntry::sqlStrTo<ColType>(r[i][0]));
    }

    return column;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Database_hpp_
