/*
 * Module.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_Module_hpp_
#define Hadrons_Module_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/VirtualMachine.hpp>

BEGIN_HADRONS_NAMESPACE

// module registration macros
#define MODULE_REGISTER(mod, base, ns)\
class mod: public base\
{\
public:\
    typedef base Base;\
    using Base::Base;\
    virtual std::string getRegisteredName(void)\
    {\
        return std::string(#ns "::" #mod);\
    }\
};\
class ns##mod##ModuleRegistrar\
{\
public:\
    ns##mod##ModuleRegistrar(void)\
    {\
        ModuleFactory &modFac = ModuleFactory::getInstance();\
        modFac.registerBuilder(#ns "::" #mod, [&](const std::string name)\
                              {\
                                  return std::unique_ptr<ns::mod>(new ns::mod(name));\
                              });\
    }\
};\
static ns##mod##ModuleRegistrar ns##mod##ModuleRegistrarInstance;

#define MODULE_REGISTER_TMP(mod, base, ns)\
extern template class base;\
MODULE_REGISTER(mod, ARG(base), ns);

#define HADRONS_MACRO_REDIRECT_12(arg1, arg2, macro, ...) macro
#define HADRONS_MACRO_REDIRECT_23(arg1, arg2, arg3, macro, ...) macro

#define envGetGrid4(latticeType)\
env().template getGrid<typename latticeType::vector_type>()

#define envGetGrid5(latticeType, Ls)\
env().template getGrid<typename latticeType::vector_type>(Ls)

#define envGetGrid(...)\
HADRONS_MACRO_REDIRECT_12(__VA_ARGS__, envGetGrid5, envGetGrid4)(__VA_ARGS__)

#define envGetCoarseGrid4(latticeType, blockSize)\
env().template getCoarseGrid<typename latticeType::vector_type>(blockSize)

#define envGetCoarseGrid5(latticeType, blockSize, Ls)\
env().template getCoarseGrid<typename latticeType::vector_type>(blockSize, Ls)

#define envGetCoarseGrid(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envGetCoarseGrid5, envGetCoarseGrid4)(__VA_ARGS__)

#define envGetRbGrid4(latticeType)\
env().template getRbGrid<typename latticeType::vector_type>()

#define envGetRbGrid5(latticeType, Ls)\
env().template getRbGrid<typename latticeType::vector_type>(Ls)

#define envGetRbGrid(...)\
HADRONS_MACRO_REDIRECT_12(__VA_ARGS__, envGetRbGrid5, envGetRbGrid4)(__VA_ARGS__)

#define envGet(type, name)\
*env().template getObject<type>(name)

#define envGetDerived(base, type, name)\
*env().template getDerivedObject<base, type>(name)

#define envGetTmp(type, var)\
type &var = *env().template getObject<type>(getName() + "_tmp_" + #var)

#define envHasType(type, name)\
env().template isObjectOfType<type>(name)

#define envCreate(type, name, Ls, ...)\
env().template createObject<type>(name, Environment::Storage::standard, Ls, __VA_ARGS__)

#define envCreateDerived(base, type, name, Ls, ...)\
env().template createDerivedObject<base, type>(name, Environment::Storage::standard, Ls, __VA_ARGS__)

#define envCreateLat4(type, name)\
envCreate(type, name, 1, envGetGrid(type))

#define envCreateLat5(type, name, Ls)\
envCreate(type, name, Ls, envGetGrid(type, Ls))

#define envCreateLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envCreateLat5, envCreateLat4)(__VA_ARGS__)

#define envCache(type, name, Ls, ...)\
env().template createObject<type>(name, Environment::Storage::cache, Ls, __VA_ARGS__)

#define envCacheLat4(type, name)\
envCache(type, name, 1, envGetGrid(type))

#define envCacheLat5(type, name, Ls)\
envCache(type, name, Ls, envGetGrid(type, Ls))

#define envCacheLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envCacheLat5, envCacheLat4)(__VA_ARGS__)

#define envTmp(type, name, Ls, ...)\
env().template createObject<type>(getName() + "_tmp_" + name,         \
                                  Environment::Storage::temporary, Ls, __VA_ARGS__)

#define envTmpLat4(type, name)\
envTmp(type, name, 1, envGetGrid(type))

#define envTmpLat5(type, name, Ls)\
envTmp(type, name, Ls, envGetGrid(type, Ls))

#define envTmpLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envTmpLat5, envTmpLat4)(__VA_ARGS__)

/******************************************************************************
 *                            Module class                                    *
 ******************************************************************************/
// base class
class ModuleBase: public TimerArray
{
public:
    struct ResultEntryHeader: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlNotNull<unsigned int>          , traj,
                           SqlUnique<SqlNotNull<std::string>>, filename);
    };
public:
    // constructor
    ModuleBase(const std::string name);
    // destructor
    virtual ~ModuleBase(void) = default;
    // access
    std::string getName(void) const;
    // get factory registration name if available
    virtual std::string getRegisteredName(void);
    // dependencies/products
    virtual std::vector<std::string> getInput(void) = 0;
    virtual std::vector<std::string> getReference(void)
    {
        return std::vector<std::string>(0);
    };
    virtual std::vector<std::string> getOutput(void) = 0;
    virtual std::vector<std::string> getOutputFiles(void)
    {
        return std::vector<std::string>(0);
    };
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name) = 0;
    virtual void saveParameters(XmlWriter &writer, const std::string name) = 0;
    // parameter string
    virtual std::string parString(void) const = 0;
    virtual std::string parClassName(void) const = 0;
    // result filename generation
    static std::string resultFilename(const std::string stem, const unsigned int traj, 
                                      const std::string ext = resultFileExt);
    // result database
    template <typename EntryType>
    void setResultDbEntry(Database &db, const std::string tableName, EntryType &entry);
    void generateResultDb(void);
    // setup
    virtual void setup(void) {};
    // execution
    virtual void execute(void) = 0;
    void operator()(void);
protected:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // virtual machine shortcut
    DEFINE_VM_ALIAS;
    // shortcuts for grid pointers
    template <typename Field>
    GridBase *getGrid4d(const bool redBlack = false);
    template <typename Field>
    GridBase *getGrid5d(const bool redBlack = false, const unsigned int Ls = 1);
    template <typename Field>
    GridBase *getGrid(const bool redBlack = false, const unsigned int Ls = 1);
    template <typename Field>
    GridBase *getGrid(const unsigned int Ls);
    // get RNGs seeded from module string
    GridParallelRNG &rng4d(void);
    GridSerialRNG &rngSerial(void);
    // result file utilities
    std::string resultFilename(const std::string stem, const std::string ext = resultFileExt) const;
    template <typename T>
    void saveResult(const std::string stem, const std::string name, const T &result) const;
private:
    // make module unique string
    std::string makeSeedString(void);
private:
    std::string                             name_, currentTimer_, seed_, dbTable_;
    std::map<std::string, GridStopWatch>    timer_;
    Database                                *db_{nullptr};
    std::unique_ptr<SqlEntry>               entry_{nullptr};
    ResultEntryHeader                       *entryHeader_{nullptr};
};

// derived class, templating the parameter class
template <typename P>
class Module: public ModuleBase
{
public:
    typedef P Par;
public:
    // constructor
    Module(const std::string name);
    // destructor
    virtual ~Module(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name);
    virtual void saveParameters(XmlWriter &writer, const std::string name);
    // parameter strings
    virtual std::string parString(void) const;
    virtual std::string parClassName(void) const;
    // parameter access
    const P &   par(void) const;
    void        setPar(const P &par);
private:
    P par_;
};

// no parameter type
class NoPar {};

template <>
class Module<NoPar>: public ModuleBase
{
public:
    // constructor
    Module(const std::string name): ModuleBase(name) {};
    // destructor
    virtual ~Module(void) = default;
    // parse parameters (do nothing)
    virtual void parseParameters(XmlReader &reader, const std::string name) {};
    virtual void saveParameters(XmlWriter &writer, const std::string name)
    {
        push(writer, "options");
        pop(writer);
    };
    // parameter strings (empty)
    virtual std::string parString(void) const {return "";};
    virtual std::string parClassName(void) const {return "";};
};

/******************************************************************************
 *                           Template implementation                          *
 ******************************************************************************/
template <typename EntryType>
void ModuleBase::setResultDbEntry(Database &db, const std::string tableName, EntryType &entry)
{
    typedef MergedSqlEntry<ResultEntryHeader, EntryType> ResultEntry;

    ResultEntryHeader entryHead;

    if (!db.isConnected())
    {
        HADRONS_ERROR(Database, "result database not connected");
    }
    db_      = &db;
    dbTable_ = tableName;
    entry_.reset(new ResultEntry(mergeSqlEntries(entryHead, entry)));
    if (!db_->tableExists(dbTable_))
    {
        db_->createTable<ResultEntry>(dbTable_);
    }
    ResultEntry *e = dynamic_cast<ResultEntry *>(entry_.get());
    entryHeader_ = &(e->template getEntry<0>());
}

template <typename T>
void ModuleBase::saveResult(const std::string stem, const std::string name, const T &result) const
{
    if (env().getGrid()->IsBoss() and !stem.empty())
    {
        makeFileDir(stem, env().getGrid());
        {
            ResultWriter writer(resultFilename(stem));
            write(writer, name, result);
        }
    }
}

template <typename Field>
GridBase *ModuleBase::getGrid4d(const bool redBlack)
{
    GridBase *g;
    
    if (redBlack)
    {
        g = env().template getRbGrid<typename Field::vector_type>();
    }
    else
    {
        g = env().template getGrid<typename Field::vector_type>();
    }

    return g;
}

template <typename Field>
GridBase *ModuleBase::getGrid5d(const bool redBlack, const unsigned int Ls)
{
    GridBase *grid;

    if (redBlack)
    {
        grid = env().template getRbGrid<typename Field::vector_type>(Ls);
    }
    else
    {
        grid = env().template getGrid<typename Field::vector_type>(Ls);
    }

    return grid;
}

template <typename Field>
GridBase *ModuleBase::getGrid(const bool redBlack, const unsigned int Ls)
{
    if (Ls == 1)
    {
        return getGrid4d<Field>(redBlack);
    }
    else
    {
        return getGrid5d<Field>(redBlack, Ls);
    }
}

template <typename Field>
GridBase *ModuleBase::getGrid(const unsigned int Ls)
{
    return getGrid<Field>(false, Ls);
}

template <typename P>
Module<P>::Module(const std::string name)
: ModuleBase(name)
{}

template <typename P>
void Module<P>::parseParameters(XmlReader &reader, const std::string name)
{
    read(reader, name, par_);
}

template <typename P>
void Module<P>::saveParameters(XmlWriter &writer, const std::string name)
{
    write(writer, name, par_);
}

template <typename P>
std::string Module<P>::parString(void) const
{
    XmlWriter writer("", "");

    write(writer, par_.SerialisableClassName(), par_);

    return writer.string();
}

template <typename P>
std::string Module<P>::parClassName(void) const
{
    return par_.SerialisableClassName();
}

template <typename P>
const P & Module<P>::par(void) const
{
    return par_;
}

template <typename P>
void Module<P>::setPar(const P &par)
{
    par_ = par;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Module_hpp_
