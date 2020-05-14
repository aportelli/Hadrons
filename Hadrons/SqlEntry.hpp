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

class SqlEntry
{
public:
    template <typename T>
    static std::string strFrom(const T &x);
    template <typename T>
    static SER(T, std::string) xmlStrFrom(const T &x);
    template <typename T>
    static NOT_SER(T, std::string) xmlStrFrom(const T &x);
    template <typename T>
    static SER_AND_ENUM(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static SER_AND_NOT_ENUM(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static NOT_SER_AND_STR(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static NOT_SER_AND_NOT_STR(T, std::string) sqlStrFrom(const T &x);
    template <typename T>
    static T strTo(const std::string str);
    template <typename T>
    static SER(T, T) xmlStrTo(const std::string str);
    template <typename T>
    static NOT_SER(T, T) xmlStrTo(const std::string str);
    template <typename T>
    static SER_AND_ENUM(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static SER_AND_NOT_ENUM(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static NOT_SER_AND_STR(T, T) sqlStrTo(const std::string str);
    template <typename T>
    static NOT_SER_AND_NOT_STR(T, T) sqlStrTo(const std::string str);
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
    virtual std::string sqlInsert(void) const = 0;
    virtual void deserializeRow(const std::vector<std::string> &row) = 0;
    virtual unsigned int cols(void) const = 0;
};

template <typename T>
std::string SqlEntry::strFrom(const T &x)
{
    std::ostringstream stream;
    
    stream << x;

    return stream.str();
}

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

template <typename T>
T SqlEntry::strTo(const std::string str)
{
    T                  buf;
    std::istringstream stream(str);
    
    stream >> buf;
    
    return buf;
}

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
#define HADRONS_SQL_MEMBER(A, B) CppType<A>::type B;
#define HADRONS_SQL_SCHEMA(A, B) schema += std::string(#B) + " " + sqlType<A>() + ",";
#define HADRONS_SQL_INSERT(A, B)\
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
list += ",";
#define HADRONS_SQL_DESERIALIZE(A, B) B = sqlStrTo<CppType<A>::type>(*it); it++;
#define HADRONS_SQL_COUNT(A, B) c++;

#define HADRONS_SQL_FIELDS(...)\
GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_MEMBER, __VA_ARGS__))\
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
template <typename... Ts>
class MergedSqlEntry: public SqlEntry
{
public:
    typedef std::tuple<Ts...> Tuple;
private:
    template <unsigned int i>
    class StorePt
    {
    public:
        static void apply(std::array<SqlEntry *, sizeof...(Ts)> &array, Tuple &t)
        {
            array[i] = &std::get<i>(t);
            StorePt<i-1>::apply(array, t);
        }
    };

    template <>
    class StorePt<0>
    {
    public:
        static void apply(std::array<SqlEntry *, sizeof...(Ts)> &array, Tuple &t)
        {
            array[0] = &std::get<0>(t);
        }
    };
public:
    MergedSqlEntry(void) 
    {
        storePt();
    }

    MergedSqlEntry(const MergedSqlEntry<Ts...> &e)
    {
        data_ = e.data_;
        storePt();
    }

    MergedSqlEntry(Ts &... entries): data_{entries...} 
    {
        storePt();
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
        StorePt<sizeof...(Ts)-1>::apply(pt_, data_);
    }
private:
    std::array<SqlEntry *, sizeof...(Ts)> pt_;
    std::tuple<Ts...>                     data_;
};

template <typename... Ts>
MergedSqlEntry<Ts...> mergeSqlEntries(Ts &... entries)
{
    return MergedSqlEntry<Ts...>(entries...);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_SqlEntry_hpp_
