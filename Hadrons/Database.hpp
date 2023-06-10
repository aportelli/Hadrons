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

#ifndef HADRONS_SQLITE_DEFAULT_JOURNAL_MODE
#define HADRONS_SQLITE_DEFAULT_JOURNAL_MODE "WAL"
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 String table class for SQL query results                   *
 ******************************************************************************/
class QueryResult
{
    friend class Database;
public:
    // constructor/destructor
    QueryResult(void) = default;
    virtual ~QueryResult(void) = default;
    // row access
    const std::vector<std::string> & operator[](const unsigned int i) const;
    // column header access
    const std::string & colName(const unsigned int j) const;
    // number of rows and columns
    size_t rows(void) const;
    size_t cols(void) const;
private:
    // broadcast data from boss MPI process
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
    struct KeyValueEntry: public SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<std::string>>, key,
                           std::string                       , value);
    };
public:
    // constructors
    Database(void) = default;
    Database(const std::string filename, GridBase *grid = nullptr,
             const std::string mode = HADRONS_SQLITE_DEFAULT_JOURNAL_MODE);
    // destructor
    virtual ~Database(void);
    // set/get DB filename
    void        setFilename(const std::string filename, GridBase *grid = nullptr, 
                            const std::string mode = HADRONS_SQLITE_DEFAULT_JOURNAL_MODE);
    std::string getFilename(void) const;
    // test if DB connected
    bool isConnected(void) const;
    // execute arbitrary SQL statement
    QueryResult execute(const std::string query);
    // test if table exists
    bool tableExists(const std::string tableName);
    // test if table is empty
    bool tableEmpty(const std::string tableName);
    // general tables interface
    template <typename EntryType>
    void createTable(const std::string tableName, const std::string extra = "");
    QueryResult getTable(const std::string tableName, const std::string extra = "");
    template <typename EntryType>
    std::vector<EntryType> getTable(const std::string tableName, const std::string extra = "");
    void insert(const std::string tableName, const SqlEntry &entry, const bool replace = false);
    void insert(const std::string tableName, const std::vector<const SqlEntry *> &entryPtVec, const bool replace = false);
    // key-value tables interface
    void createKeyValueTable(const std::string tableName);
    std::map<std::string, std::string> getKeyValueTable(const std::string tableName);
    template <typename T>
    T getValue(const std::string keyValueTableName, const std::string key);
    template <typename T>
    void insertValue(const std::string keyValueTableName, const std::string key, const T &value, const bool replace = false);
    // get a single column from a table
    template <typename ColType>
    std::vector<ColType> getTableColumn(const std::string tableName, const std::string columnName, const std::string extra = "");
private:
    // private connect/disconnect functions
    void connect(void);
    void disconnect(void);
private:
    std::string filename_;
    GridBase    *grid_{nullptr};
    sqlite3     *db_{nullptr};
    bool        isConnected_{false};
};

/******************************************************************************
 *                  Database class template implementation                    *
 ******************************************************************************/
// general tables interface ///////////////////////////////////////////////////
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
    QueryResult            r = getTable(tableName, extra);

    for (unsigned int i = 0; i < r.rows(); ++i)
    {
        buf.deserializeRow(r[i]);
        table.push_back(buf);
    }

    return table;
}


// key-value tables interface //////////////////////////////////////////////////
template <typename T>
T Database::getValue(const std::string keyValueTableName, const std::string key)
{
    auto vec = getTableColumn<T>(keyValueTableName, "value", "WHERE key = '" + key + "'");

    if (!vec.empty())
    {
        return vec[0];
    }
    else
    {
        HADRONS_ERROR(Database, "no value with key '" + key 
                      + "' in key-value table '" + keyValueTableName 
                      + "' in database '" + filename_ + "'");
    }
}

template <typename T>
void Database::insertValue(const std::string keyValueTableName, const std::string key, 
                           const T &value, const bool replace)
{
    KeyValueEntry e;

    e.key   = key;
    e.value = SqlEntry::sqlStrFrom(value);
    insert(keyValueTableName, e, replace);
}

// get a single column from a table ////////////////////////////////////////////
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
