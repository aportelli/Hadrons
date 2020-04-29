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
#include <Hadrons/sqlite/sqlite3.h>

BEGIN_HADRONS_NAMESPACE

class SqlEntry
{
public:
    template <typename T>
    static std::string strFrom(const T &x);
    template <typename T>
    static std::string sqlType(const T &x);
    template <typename T>
    static typename std::enable_if<std::is_floating_point<T>::value, std::string>::type
    sqlType(void);
    template <typename T>
    static typename std::enable_if<std::is_integral<T>::value, std::string>::type
    sqlType(void);
    template <typename T>
    static typename std::enable_if<!std::is_floating_point<T>::value 
                                   and !std::is_integral<T>::value, std::string>::type
    sqlType(void);
    virtual std::string sqlInsert(void) const = 0;
};

template <typename T>
std::string SqlEntry::strFrom(const T &x)
{
    std::ostringstream stream;
    
    stream << x;

    return stream.str();
}

template <typename T>
std::string SqlEntry::sqlType(const T &x)
{
    return sqlType<T>();
}

template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, std::string>::type
SqlEntry::sqlType(void)
{
    return "REAL";
}

template <typename T>
typename std::enable_if<std::is_integral<T>::value, std::string>::type
SqlEntry::sqlType(void)
{
    return "INTEGER";
}

template <typename T>
typename std::enable_if<!std::is_floating_point<T>::value 
                        and !std::is_integral<T>::value, std::string>::type
SqlEntry::sqlType(void)
{
    return "TEXT";
}

#define HADRONS_SQL_SCHEMA(A, B) schema += "\"" + std::string(#B) + "\" " + sqlType<A>() + " NOT NULL,";
#define HADRONS_SQL_INSERT(A, B)\
if (sqlType<A>() == "TEXT") list += "\"";\
list += strFrom(B);\
if (sqlType<A>() == "TEXT") list += "\"";\
list += ",";

#define HADRONS_SQL_FIELDS(...)\
GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_MEMBER, __VA_ARGS__))\
static std::string sqlSchema(void)\
{\
    std::string schema;\
    \
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_SCHEMA, __VA_ARGS__))\
    schema.pop_back();\
    \
    return schema;\
}\
\
virtual std::string sqlInsert(void) const\
{\
    std::string list;\
    \
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_INSERT, __VA_ARGS__))\
    list.pop_back();\
    \
    return list;\
}

template <typename... Ts>
class MergedSqlEntry: SqlEntry
{
public:
    MergedSqlEntry(const Ts &... entries): pt_{&entries...} {}

    static std::string sqlSchema(void)
    {
        std::array<std::string, sizeof...(Ts)> vec = {Ts::sqlSchema()...};
        std::string                            schema;

        for (auto &s: vec)
        {
            schema += s + ",";
        }
        schema.pop_back();

        return schema;
    }

    virtual std::string sqlInsert(void) const
    {
        std::string list;

        for (auto &e: pt_)
        {
            list += e->sqlInsert() + ",";
        }
        list.pop_back();

        return list;
    }

private:
    std::array<const SqlEntry *, sizeof...(Ts)> pt_;
};

template <typename... Ts>
MergedSqlEntry<Ts...> mergeSqlEntries(const Ts &... entries)
{
    return MergedSqlEntry<Ts...>(entries...);
}

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

class Database
{
public:
    Database(const std::string filename, GridBase *grid = nullptr);
    virtual ~Database(void);
    QueryResult execute(const std::string query);
private:
    void connect(const std::string filename);
    void disconnect(void);
private:
    std::string filename_;
    GridBase    *grid_;
    sqlite3     *db_{nullptr};
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Database_hpp_