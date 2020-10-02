/*
 * SqlEntry.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_SqlEntry_hpp_
#define Hadrons_SqlEntry_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Type traits for SQL column type options                  *
 ******************************************************************************/
template <typename T>
struct SqlColumnOption {};

template <typename T>
struct SqlUnique: SqlColumnOption<SqlUnique<T>>
{
    typedef T type;
    static const std::string option;
};

template <typename T>
const std::string SqlUnique<T>::option = "UNIQUE";

template <typename T>
struct SqlNotNull: SqlColumnOption<SqlNotNull<T>>
{
    typedef T type;
    static const std::string option;
};

template <typename T>
const std::string SqlNotNull<T>::option = "NOT NULL";

template <typename T, typename Enable = void>
struct CppType
{
    typedef T type;
};

// deduce C++ type when decorated with SQL options
template <typename T>
struct CppType<T, typename std::enable_if<std::is_base_of<SqlColumnOption<T>, T>::value>::type>
{
    typedef typename CppType<typename T::type>::type type;
};

/******************************************************************************
 *                          Base class for SQL entries                        *
 ******************************************************************************/
// shortcuts for cumbersome enable_if
#define SER(T, RT)\
typename std::enable_if<std::is_base_of<Serializable, T>::value, RT>::type

#define NOT_SER(T, RT)\
typename std::enable_if<!std::is_base_of<Serializable, T>::value, RT>::type

#define SER_AND_ENUM(T, RT)\
typename std::enable_if<std::is_base_of<Serializable, T>::value and T::isEnum, RT>::type

#define SER_AND_NOT_ENUM(T, RT)\
typename std::enable_if<std::is_base_of<Serializable, T>::value and !T::isEnum, RT>::type

#define NOT_SER_AND_STR(T, RT)\
typename std::enable_if<!std::is_base_of<Serializable, T>::value and std::is_assignable<std::string, T>::value, RT>::type

#define NOT_SER_AND_NOT_STR(T, RT)\
typename std::enable_if<!std::is_base_of<Serializable, T>::value and !std::is_assignable<std::string, T>::value, RT>::type

// base class for SQL rows
class SqlEntry
{
public:
    SqlEntry(void) = default;
    virtual ~SqlEntry(void) = default;
    // string from an arbitrary type
    template <typename T>
    static std::string strFrom(const T &x);
    // XML string from an arbitrary type
    template <typename T>
    static SER(T, std::string) xmlStrFrom(const T &x);
    template <typename T>
    static NOT_SER(T, std::string) xmlStrFrom(const T &x);
    // SQL string from an arbitrary type
    template <typename T>
    static SER_AND_ENUM(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static SER_AND_NOT_ENUM(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static NOT_SER_AND_STR(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static NOT_SER_AND_NOT_STR(T, std::string) sqlStrFrom(const T &x);
    // parse string to an arbitrary type
    template <typename T>
    static T strTo(const std::string str);
    // parse XML string to an arbitrary type
    template <typename T>
    static SER(T, T) xmlStrTo(const std::string str);
    template <typename T>
    static NOT_SER(T, T) xmlStrTo(const std::string str);
    // parse SQL string to an arbitrary type
    template <typename T>
    static SER_AND_ENUM(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static SER_AND_NOT_ENUM(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static NOT_SER_AND_STR(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static NOT_SER_AND_NOT_STR(T, T) sqlStrTo(const std::string str);
    // SQL type (REAL, INTEGER or TEXT) from an arbitrary type
    template <typename T>
    static typename std::enable_if<std::is_floating_point<T>::value, std::string>::type
    sqlType(void);
    template <typename T>
    static typename std::enable_if<std::is_integral<T>::value, std::string>::type
    sqlType(void);
    template <typename T>
    static typename std::enable_if<std::is_base_of<SqlColumnOption<T>, T>::value, std::string>::type
    sqlType(void);
    template <typename T>
    static typename std::enable_if<!std::is_floating_point<T>::value 
                                   and !std::is_integral<T>::value
                                   and !std::is_base_of<SqlColumnOption<T>, T>::value, std::string>::type
    sqlType(void);
    // abstract interface
    virtual std::string sqlInsert(void) const = 0;
    virtual void deserializeRow(const std::vector<std::string> &row) = 0;
    virtual unsigned int cols(void) const = 0;
};

// print SQL entry as CSV
inline std::ostream & operator<<(std::ostream &out, const SqlEntry &e)
{
    out << e.sqlInsert();

    return out;
}

// string from an arbitrary type
template <typename T>
std::string SqlEntry::strFrom(const T &x)
{
    std::ostringstream stream;
    
    stream << x;

    return stream.str();
}

// XML string from an arbitrary type
template <typename T>
SER(T, std::string) SqlEntry::xmlStrFrom(const T &x)
{
    XmlWriter writer("", "");

    write(writer, x.SerialisableClassName(), x);

    return writer.string();
}

template <typename T>
NOT_SER(T, std::string) SqlEntry::xmlStrFrom(const T &x)
{
    XmlWriter writer("", "");

    write(writer, "object", x);

    return writer.string();
}

// SQL string from an arbitrary type
// NB: numerical types are assignable to std::string through char
// implicit conversion! So PODs will use NOT_SER_AND_STR
// that's weird sementically but in principle guaranteed by the Standard
// and avoids an extra POD case below. AP.
template <typename T>
SER_AND_ENUM(T, std::string) SqlEntry::sqlStrFrom(const T &x)
{
    return strFrom(x);
}

template <typename T>
SER_AND_NOT_ENUM(T, std::string) SqlEntry::sqlStrFrom(const T &x)
{
    return xmlStrFrom(x);
}

template <typename T>
NOT_SER_AND_STR(T, std::string) SqlEntry::sqlStrFrom(const T &x)
{
    return strFrom(x);
}

template <typename T>
NOT_SER_AND_NOT_STR(T, std::string) SqlEntry::sqlStrFrom(const T &x)
{
    return xmlStrFrom(x);
}

// parse string to an arbitrary type
template <typename T>
T SqlEntry::strTo(const std::string str)
{
    T                  buf;
    std::istringstream stream(str);
    
    stream >> buf;
    
    return buf;
}

template <>
inline std::string SqlEntry::strTo(const std::string str)
{
    return str;
}

// parse XML string to an arbitrary type
template <typename T>
SER(T, T) SqlEntry::xmlStrTo(const std::string str)
{
    T         buf;
    XmlReader reader(str, true, "");

    read(reader, buf.SerialisableClassName(), buf);

    return buf;
}

template <typename T>
NOT_SER(T, T) SqlEntry::xmlStrTo(const std::string str)
{
    T         buf;
    XmlReader reader(str, true, "");

    read(reader, "object", buf);

    return buf;
}

// parse SQL string to an arbitrary type
template <typename T>
SER_AND_ENUM(T, T) SqlEntry::sqlStrTo(const std::string str)
{
    return strTo<T>(str);
}

template <typename T>
SER_AND_NOT_ENUM(T, T) SqlEntry::sqlStrTo(const std::string str)
{
    return xmlStrTo<T>(str);
}

template <typename T>
NOT_SER_AND_STR(T, T) SqlEntry::sqlStrTo(const std::string str)
{
    return strTo<T>(str);
}

template <typename T>
NOT_SER_AND_NOT_STR(T, T) SqlEntry::sqlStrTo(const std::string str)
{
    return xmlStrTo<T>(str);
}

// SQL type (REAL, INTEGER or TEXT) from an arbitrary type
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
typename std::enable_if<std::is_base_of<SqlColumnOption<T>, T>::value, std::string>::type
SqlEntry::sqlType(void)
{
    return sqlType<typename T::type>() + " " + T::option;
}

template <typename T>
typename std::enable_if<!std::is_floating_point<T>::value 
                        and !std::is_integral<T>::value
                        and !std::is_base_of<SqlColumnOption<T>, T>::value, std::string>::type
SqlEntry::sqlType(void)
{
    return "TEXT";
}   

/******************************************************************************
 *                 "Macro Magic" for SQL entry class declarations             *
 ******************************************************************************/
#define HADRONS_SQL_MEMBER(A, B)      CppType<A>::type B;
#define HADRONS_SQL_BOOL_MEMBER(A, B) bool B{false};
#define HADRONS_SQL_SCHEMA(A, B)      schema += std::string(#B) + " " + sqlType<A>() + ",";
#define HADRONS_SQL_INSERT(A, B)\
if (nullify.B)\
{\
    list += "NULL";\
}\
else\
{\
    if (sqlType<CppType<A>::type>() == "TEXT")\
    {\
        std::string s;\
        s = sqlStrFrom(B);\
        if (!s.empty())\
        {\
            list += "'" + s + "'";\
        }\
        else\
        {\
            list += "NULL";\
        }\
    }\
    else\
    {\
        list += sqlStrFrom(B);\
    }\
}\
list += ",";
#define HADRONS_SQL_DESERIALIZE(A, B) B = sqlStrTo<CppType<A>::type>(*it); it++;
#define HADRONS_SQL_COUNT(A, B) c++;

#define HADRONS_SQL_FIELDS(...)\
GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_MEMBER, __VA_ARGS__))\
struct Nullify\
{\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_BOOL_MEMBER, __VA_ARGS__))\
};\
Nullify nullify;\
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
}\
virtual void deserializeRow(const std::vector<std::string> &row)\
{\
    auto it = row.begin();\
    \
    if (row.size() != cols())\
    {\
        HADRONS_ERROR(Database, "cannot deserialize row (has " + std::to_string(row.size())\
                                + " columns, expected " + std::to_string(cols()) + ")");\
    }\
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_DESERIALIZE, __VA_ARGS__))\
}\
virtual unsigned int cols(void) const\
{\
    unsigned int c = 0;\
    \
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_COUNT, __VA_ARGS__))\
    \
    return c;\
}

/******************************************************************************
 *                      Utility class to merge SQL entries                    *
 ******************************************************************************/
// helper to store an array of pointers on tuple elements
template <unsigned int i, typename... Ts>
struct StorePt
{
    static void apply(std::array<SqlEntry *, sizeof...(Ts)> &array, 
                      std::tuple<Ts...> &t)
    {
        array[i] = &std::get<i>(t);
        StorePt<i - 1, Ts...>::apply(array, t);
    }
};

template <typename... Ts>
struct StorePt<0, Ts...>
{
    static void apply(std::array<SqlEntry *, sizeof...(Ts)> &array,
                      std::tuple<Ts...> &t)
    {
        array[0] = &std::get<0>(t);
    }
};

// merged SQL entries class, derives from SqlEntry
template <typename... Ts>
class MergedSqlEntry: public SqlEntry
{
public:
    typedef std::tuple<typename std::decay<Ts>::type...> Tuple;
public:
    MergedSqlEntry(void) 
    {
        storePt();
    }

    MergedSqlEntry(const MergedSqlEntry<Ts...> &e)
    {
        *this = e;
    }

    MergedSqlEntry(Ts &... entries): data_{entries...} 
    {
        storePt();
    }

    MergedSqlEntry<Ts...> & operator=(const MergedSqlEntry<Ts...> &e)
    {
        if (this != &e)
        {
            data_ = e.data_;
            storePt();
        }

        return *this;
    }

    template <unsigned int i>
    typename std::tuple_element<i, Tuple>::type & getEntry(void)
    {
        return std::get<i>(data_);
    }

    SqlEntry * getEntry(const unsigned int i)
    {
        return pt_[i];
    }

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

        for (auto e: pt_)
        {
            list += e->sqlInsert() + ",";
        }
        list.pop_back();

        return list;
    }

    virtual void deserializeRow(const std::vector<std::string> &row)
    {
        std::vector<std::string> buf;
        auto                     it = row.begin();

        if (row.size() != cols())
        {
            HADRONS_ERROR(Database, "cannot deserialize row (has " + std::to_string(row.size())
                                    + " columns, expected " + std::to_string(cols()) + ")");
        }
        for (unsigned int i = 0; i < pt_.size(); ++i)
        {
            buf.clear();
            for (unsigned int j = 0; j < pt_[i]->cols(); ++j)
            {
                buf.push_back(*it);
                it++;
            }
            pt_[i]->deserializeRow(buf);
        }
    }

    virtual unsigned int cols(void) const
    {
        unsigned int c = 0;

        for (unsigned int i = 0; i < pt_.size(); ++i)
        {
            c += pt_[i]->cols();
        }

        return c;
    }
private:
    void storePt(void)
    {
        StorePt<sizeof...(Ts) - 1, typename std::decay<Ts>::type...>::apply(pt_, data_);
    }
private:
    std::array<SqlEntry *, sizeof...(Ts)> pt_;
    Tuple                                 data_;
};

// function merging SQL entries
template <typename... Ts>
MergedSqlEntry<Ts...> mergeSqlEntries(Ts &... entries)
{
    return MergedSqlEntry<Ts...>(entries...);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_SqlEntry_hpp_
