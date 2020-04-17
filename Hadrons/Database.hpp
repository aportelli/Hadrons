#ifndef Hadrons_Database_hpp_
#define Hadrons_Database_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/sqlite/sqlite3ext.h>

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

#define HADRONS_SQL_SCHEMA(A, B) schema += "\"" + std::string(#B) + "\" " + sqlType<A>() + ",";
#define HADRONS_SQL_INSERT(A, B)\
if (sqlType<A>() == "TEXT") list += "\"";\
list += strFrom(B);\
if (sqlType<A>() == "TEXT") list += "\"";\
list += ",";

#define HADRONS_SQL_FIELDS(cname, ...)\
GRID_MACRO_EVAL(GRID_MACRO_MAP(GRID_MACRO_MEMBER,__VA_ARGS__))\
static inline std::string sqlSchema(void)\
{\
    std::string schema;\
    \
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_SCHEMA, __VA_ARGS__))\
    schema.pop_back();\
    \
    return schema;\
}\
\
inline std::string sqlInsert(void)\
{\
    std::string list;\
    \
    GRID_MACRO_EVAL(GRID_MACRO_MAP(HADRONS_SQL_INSERT, __VA_ARGS__))\
    list.pop_back();\
    \
    return list;\
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Database_hpp_