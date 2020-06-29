/*
 * VirtualMachine.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_VirtualMachine_hpp_
#define Hadrons_VirtualMachine_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>
#include <Hadrons/Graph.hpp>
#include <Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

#define DEFINE_VM_ALIAS \
inline VirtualMachine & vm(void) const\
{\
    return VirtualMachine::getInstance();\
}

/******************************************************************************
 *                   Virtual machine for module execution                     *
 ******************************************************************************/
// forward declaration of Module
class ModuleBase;

class VirtualMachine
{
    SINGLETON_DEFCTOR(VirtualMachine);
public:
    typedef SITE_SIZE_TYPE                      Size;
    typedef std::unique_ptr<ModuleBase>         ModPt;
    typedef std::vector<std::set<unsigned int>> GarbageSchedule;
    typedef std::vector<unsigned int>           Program;
    struct MemoryPrint
    {
        Size                 size;
        Environment::Storage storage;
        int                  module;
    };
    struct MemoryProfile
    {
        std::vector<std::map<unsigned int, Size>> module;
        std::vector<MemoryPrint>                  object;
    };
    class GeneticPar: Serializable
    {
    public:
        GeneticPar(void):
            popSize{20}, maxGen{1000}, maxCstGen{100}, mutationRate{.1} {};
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GeneticPar,
                                        unsigned int, popSize,
                                        unsigned int, maxGen,
                                        unsigned int, maxCstGen,
                                        double      , mutationRate);
    };

    // serializable classes for database entries
    struct GlobalEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<std::string>>, name,
                           SqlNotNull<std::string>           , value);
    };

    struct ModuleEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, moduleId,
                           SqlUnique<SqlNotNull<std::string>> , name,
                           SqlNotNull<unsigned int>           , moduleTypeId,
                           std::string                        , parameters);
    };

    struct ModuleTypeEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<unsigned int>           , moduleTypeId,
                           SqlUnique<SqlNotNull<std::string>>, type);
                           
    };

    struct ObjectEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, objectId,
                           SqlUnique<SqlNotNull<std::string>> , name,
                           SqlNotNull<unsigned int>           , objectTypeId,
                           SqlNotNull<SITE_SIZE_TYPE>         , size,
                           SqlNotNull<Environment::Storage>   , storageType,
                           SqlNotNull<unsigned int>           , moduleId);
    };

    struct ObjectTypeEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<unsigned int>           , objectTypeId,
                           SqlUnique<SqlNotNull<std::string>>, type,
                           SqlNotNull<std::string>, baseType);
    };

    struct ScheduleEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, step,
                           SqlUnique<SqlNotNull<unsigned int>>, moduleId);
    };
private:
    struct ModuleInfo
    {
        const std::type_info      *type{nullptr};
        std::string               name;
        ModPt                     data{nullptr};
        std::vector<unsigned int> input, output;
        size_t                    maxAllocated;
    };
public:
    // trajectory counter
    void                setTrajectory(const unsigned int traj);
    unsigned int        getTrajectory(void) const;
    // run tag
    void                setRunId(const std::string id);
    std::string         getRunId(void) const;
    // database
    void                setDatabase(Database &db);
    void                dbRestoreMemoryProfile(void);
    void                dbRestoreModules(void);
    Program             dbRestoreSchedule(void);
    // module management
    void                pushModule(ModPt &pt);
    template <typename M>
    void                createModule(const std::string name);
    template <typename M>
    void                createModule(const std::string name,
                                         const typename M::Par &par);
    void                createModule(const std::string name,
                                     const std::string type,
                                     XmlReader &reader,
                                     const std::string blockName = "options");
    unsigned int        getNModule(void) const;
    ModuleBase *        getModule(const unsigned int address) const;
    ModuleBase *        getModule(const std::string name) const;
    template <typename M>
    M *                 getModule(const unsigned int address) const;
    template <typename M>
    M *                 getModule(const std::string name) const;
    unsigned int        getModuleAddress(const std::string name) const;
    std::string         getModuleName(const unsigned int address) const;
    std::string         getModuleType(const unsigned int address) const;
    std::string         getModuleType(const std::string name) const;
    std::string         getModuleNamespace(const unsigned int address) const;
    std::string         getModuleNamespace(const std::string name) const;
    int                 getCurrentModule(void) const;
    bool                hasModule(const unsigned int address) const;
    bool                hasModule(const std::string name) const;
    // print VM content
    void                printContent(void) const;
    // module graph (could be a const reference if topoSort was const)
    Graph<unsigned int> getModuleGraph(void);
    // dump GraphViz graph
    void                dumpModuleGraph(std::ostream &out);
    void                dumpModuleGraph(void);
    void                dumpModuleGraph(const std::string filename);
    // memory profile
    const MemoryProfile &getMemoryProfile(void);
    void                printMemoryProfile(void) const;
    // garbage collector
    GarbageSchedule     makeGarbageSchedule(const Program &p) const;
    // high-water memory function
    Size                memoryNeeded(const Program &p);
    // genetic scheduler
    Program             schedule(const GeneticPar &par);
    // general execution
    void                executeProgram(const Program &p);
    void                executeProgram(const std::vector<std::string> &p);
    // generate result DB
    void                generateResultDb(void);
private:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // module graph
    void makeModuleGraph(void);
    // memory profile
    void makeMemoryProfile(void);
    void resetProfile(void);
    void resizeProfile(void);
    void updateProfile(const unsigned int address);
    void cleanEnvironment(void);
    void memoryProfile(const std::string name);
    void memoryProfile(const unsigned int address);
    // database handling
    bool         hasDatabase(void) const;
    void         initDatabase(void);
    unsigned int dbInsertModuleType(const std::string type);
    unsigned int dbInsertObjectType(const std::string type, const std::string baseType);
private:
    // general
    std::string                         runId_;
    unsigned int                        traj_;
    // database
    Database                            *db_{nullptr};
    bool                                makeModuleDb_{true}, makeObjectDb_{true}, makeScheduleDb_{true};
    // module and related maps
    std::vector<ModuleInfo>             module_;
    std::map<std::string, unsigned int> moduleAddress_;
    int                                 currentModule_{-1};
    // module graph
    bool                                graphOutdated_{true};
    Graph<unsigned int>                 graph_;
    // memory profile
    bool                                memoryProfileOutdated_{true};
    MemoryProfile                       profile_;     
    // time profile
    GridTime                            totalTime_;
    std::map<std::string, GridTime>     timeProfile_;               
};

/******************************************************************************
 *                   VirtualMachine template implementation                   *
 ******************************************************************************/
// module management ///////////////////////////////////////////////////////////
template <typename M>
void VirtualMachine::createModule(const std::string name)
{
    ModPt pt(new M(name));
    
    pushModule(pt);
}

template <typename M>
void VirtualMachine::createModule(const std::string name,
                               const typename M::Par &par)
{
    ModPt pt(new M(name));
    
    static_cast<M *>(pt.get())->setPar(par);
    pushModule(pt);
}

template <typename M>
M * VirtualMachine::getModule(const unsigned int address) const
{
    if (auto *pt = dynamic_cast<M *>(getModule(address)))
    {
        return pt;
    }
    else
    {
        HADRONS_ERROR(Definition, "module '" + module_[address].name
                     + "' does not have type " + typeid(M).name()
                     + "(has type: " + getModuleType(address) + ")");
    }
}

template <typename M>
M * VirtualMachine::getModule(const std::string name) const
{
    return getModule<M>(getModuleAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_VirtualMachine_hpp_
