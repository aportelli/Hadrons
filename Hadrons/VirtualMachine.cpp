/*
 * VirtualMachine.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/VirtualMachine.hpp>
#include <Hadrons/GeneticScheduler.hpp>
#include <Hadrons/StatLogger.hpp>
#include <Hadrons/ModuleFactory.hpp>

using namespace Grid;
 
using namespace Hadrons;

/******************************************************************************
 *                      VirtualMachine implementation                         *
 ******************************************************************************/
// trajectory counter //////////////////////////////////////////////////////////
void VirtualMachine::setTrajectory(const unsigned int traj)
{
    traj_ = traj;
}

unsigned int VirtualMachine::getTrajectory(void) const
{
    return traj_;
}

// run tag /////////////////////////////////////////////////////////////////////
void VirtualMachine::setRunId(const std::string id)
{
    runId_ = id;
}

std::string VirtualMachine::getRunId(void) const
{
    return runId_;
}

// database ////////////////////////////////////////////////////////////////////
void VirtualMachine::setDatabase(Database &db)
{
    db_ = &db;
    initDatabase();
}

void VirtualMachine::dbRestoreMemoryProfile(void)
{
    if (hasDatabase())
    {
        if (db_->tableExists("objects"))
        {
            auto         table = db_->getTable<ObjectEntry>("objects", "ORDER BY objectId");
            unsigned int nMod;

            auto comp  = [](const ObjectEntry &a, const ObjectEntry &b)
            {
                return (a.moduleId < b.moduleId);
            };

            if (env().getMaxAddress() > 0)
            {
                HADRONS_ERROR(Database, "environment is not empty");
            }
            if (table.size() == 0)
            {
                HADRONS_ERROR(Database, "object table is empty");
            }
            resetProfile();
            nMod = std::max_element(table.begin(), table.end(), comp)->moduleId + 1;
            profile_.module.resize(nMod);
            profile_.object.resize(table.size());
            for (auto &e: table)
            {
                profile_.object[e.objectId].module      = e.moduleId;
                profile_.object[e.objectId].size        = e.size;
                profile_.object[e.objectId].storage     = e.storageType;
                profile_.module[e.moduleId][e.objectId] = e.size;
                env().addObject(e.name, e.moduleId);
                assert(env().getObjectAddress(e.name) == e.objectId);
                env().setObjectStorage(e.objectId, e.storageType);
            }
            memoryProfileOutdated_ = false;
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }
}

void VirtualMachine::dbRestoreModules(void)
{
    if (hasDatabase())
    {
        if (db_->tableExists("modules"))
        {
            std::string prefix    = "Grid::Hadrons::";
            auto        modTable  = db_->getTable<ModuleEntry>("modules", "ORDER BY moduleId");
           
            if (getNModule() > 0)
            {
                HADRONS_ERROR(Database, "module graph is not empty");
            }
            if (modTable.size() == 0)
            {
                HADRONS_ERROR(Database, "module table is empty");
            }
            for (auto &e: modTable)
            {
                auto typeTable = db_->getTable<ModuleTypeEntry>("moduleTypes", 
                    "WHERE moduleTypeId = " + std::to_string(e.moduleTypeId));
                std::string type = typeTable.front().type;

                if ((type.size() > prefix.size()) 
                    and (type.substr(0, prefix.size()) == prefix))
                {
                    type = type.substr(prefix.size(), type.size() - prefix.size());
                }

                std::string par = "<top>" + e.parameters + "</top>";
                XmlReader   reader(par, true, "top");
                auto        &factory = ModuleFactory::getInstance();
                auto        pt       = factory.create(type, e.name);

                pt->parseParameters(reader, pt->parClassName());
                pushModule(pt);
            }
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }
}

VirtualMachine::Program VirtualMachine::dbRestoreSchedule(void)
{
    Program program;

    if (hasDatabase())
    {
        if (db_->tableExists("schedule"))
        {
            auto table = db_->getTable<ScheduleEntry>("schedule");

            if (table.size() == 0)
            {
                HADRONS_ERROR(Database, "schedule table is empty");
            }
            program.resize(table.size());
            for (auto &e: table)
            {
                program[e.step] = e.moduleId;
            }
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }

    return program;
}

bool VirtualMachine::hasDatabase(void) const
{
    return ((db_ != nullptr) and db_->isConnected());
}

void VirtualMachine::initDatabase(void)
{
    db_->execute("PRAGMA foreign_keys = ON;");
    if (!db_->tableExists("global"))
    {
        db_->createTable<GlobalEntry>("global");
    }
    if (!db_->tableExists("moduleTypes"))
    {
        db_->createTable<ModuleTypeEntry>("moduleTypes", "PRIMARY KEY(moduleTypeId)");
    }
    if (!db_->tableExists("modules"))
    {
        db_->createTable<ModuleEntry>("modules", "PRIMARY KEY(moduleId)"
            "FOREIGN KEY(moduleTypeId) REFERENCES moduleTypes(moduleTypeId)");
    }
    else if (!db_->tableEmpty("modules"))
    {
        LOG(Message) << "The module table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeModuleDb_ = false;
    }
    if (!db_->tableExists("objectTypes"))
    {
        db_->createTable<ObjectTypeEntry>("objectTypes", "PRIMARY KEY(objectTypeId)");
    }
    if (!db_->tableExists("objects"))
    {
        db_->createTable<ObjectEntry>("objects", "PRIMARY KEY(objectId)," 
            "FOREIGN KEY(moduleId) REFERENCES modules(moduleId),"
            "FOREIGN KEY(objectTypeId) REFERENCES objectTypes(objectTypeId)");
    }
    else if (!db_->tableEmpty("objects"))
    {
        LOG(Message) << "The object table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeObjectDb_ = false;
    }
    if (!db_->tableExists("schedule"))
    {
        db_->createTable<ScheduleEntry>("schedule", "PRIMARY KEY(step)," 
            "FOREIGN KEY(moduleId) REFERENCES modules(moduleId)");
    }
    else if (!db_->tableEmpty("schedule"))
    {
        LOG(Message) << "The schedule table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeScheduleDb_ = false;
    }
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vModules AS                                     "
        "SELECT moduleId,                                                          "
        "       modules.name,                                                      "
        "       moduleTypes.type AS type,                                          "
        "       modules.parameters                                                 "
        "FROM modules                                                              "
        "INNER JOIN moduleTypes ON modules.moduleTypeId = moduleTypes.moduleTypeId "
        "ORDER BY moduleId;                                                        "
    );
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vObjects AS                                     "
        "SELECT objectId,                                                          "
        "       objects.name,                                                      "
        "       objectTypes.type AS type,                                          "
        "       objectTypes.baseType AS baseType,                                  "
        "       objects.size*1.0/1024/1024 AS sizeMB,                              "
        "       objects.storageType,                                               "
        "       modules.name AS module                                             "
        "FROM objects                                                              "
        "INNER JOIN objectTypes ON objects.objectTypeId = objectTypes.objectTypeId "
        "INNER JOIN modules     ON objects.moduleId = modules.moduleId             "
        "ORDER BY objectId;                                                        "
    );
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vSchedule AS                                    "
        "SELECT step,                                                              "
        "       modules.name AS module                                             "
        "FROM schedule                                                             "
        "INNER JOIN modules ON schedule.moduleId = modules.moduleId                "
        "ORDER BY step;                                                            "
    );
}

unsigned int VirtualMachine::dbInsertModuleType(const std::string type)
{
    QueryResult r = db_->execute("SELECT moduleTypeId FROM moduleTypes "
                                 "WHERE type = '" + type + "';");

    if (r.rows() == 0)
    {
        ModuleTypeEntry e;

        r = db_->execute("SELECT COUNT(*) FROM moduleTypes;");
        e.moduleTypeId = std::stoi(r[0][0]);
        e.type         = type;
        db_->insert("moduleTypes", e);

        return e.moduleTypeId;
    }
    else
    {
        return std::stoi(r[0][0]);
    }
}

unsigned int VirtualMachine::dbInsertObjectType(const std::string type, 
                                             const std::string baseType)
{
    QueryResult r = db_->execute("SELECT objectTypeId FROM objectTypes "
                                 "WHERE type = '" + type + "' "
                                 "AND baseType = '" + baseType + "';");

    if (r.rows() == 0)
    {
        ObjectTypeEntry e;

        r = db_->execute("SELECT COUNT(*) FROM objectTypes;");
        e.objectTypeId = std::stoi(r[0][0]);
        e.type         = type;
        e.baseType     = baseType;
        db_->insert("objectTypes", e);

        return e.objectTypeId;
    }
    else
    {
        return std::stoi(r[0][0]);
    }
}

// module management ///////////////////////////////////////////////////////////
void VirtualMachine::pushModule(VirtualMachine::ModPt &pt)
{
    std::string name = pt->getName();
    
    if (!hasModule(name))
    {
        std::vector<unsigned int> inputAddress;
        unsigned int              address;
        ModuleInfo                m;
        
        // module registration -------------------------------------------------
        m.data = std::move(pt);
        m.type = typeIdPt(*m.data.get());
        m.name = name;
        // input dependencies
        for (auto &in: m.data->getInput())
        {
            if (!env().hasObject(in))
            {
                // if object does not exist, add it with no creator module
                env().addObject(in , -1);
                memoryProfileOutdated_ = true;
            }
            m.input.push_back(env().getObjectAddress(in));
        }
        // reference dependencies
        for (auto &ref: m.data->getReference())
        {
            if (!env().hasObject(ref))
            {
                // if object does not exist, add it with no creator module
                env().addObject(ref , -1);
                memoryProfileOutdated_ = true;
            }
            m.input.push_back(env().getObjectAddress(ref));
        }
        auto inCopy = m.input;
        // if module has inputs with references, they need to be added as
        // an input
        for (auto &in: inCopy)
        {
            int inm = env().getObjectModule(in);

            if (inm > 0)
            {
                if (getModule(inm)->getReference().size() > 0)
                {
                    for (auto &rin: getModule(inm)->getReference())
                    {
                        m.input.push_back(env().getObjectAddress(rin));
                    }
                }
            }
        }
        module_.push_back(std::move(m));
        address              = static_cast<unsigned int>(module_.size() - 1);
        moduleAddress_[name] = address;
        // connecting outputs to potential inputs ------------------------------
        for (auto &out: getModule(address)->getOutput())
        {
            if (!env().hasObject(out))
            {
                // output does not exists, add it
                env().addObject(out, address);
                memoryProfileOutdated_ = true;
                module_[address].output.push_back(env().getObjectAddress(out));
            }
            else
            {
                if (env().getObjectModule(env().getObjectAddress(out)) < 0)
                {
                    // output exists but without creator, correct it
                    env().setObjectModule(env().getObjectAddress(out), address);
                }
                else if (env().getObjectModule(env().getObjectAddress(out)) != address)
                {
                    // output already produced by another module, error
                    HADRONS_ERROR_REF(ObjectDefinition, "object '" + out
                                 + "' is already produced by module '"
                                 + module_[env().getObjectModule(out)].name
                                 + "' (while pushing module '" + name + "')",
                                 env().getObjectAddress(out));
                }
                if (getModule(address)->getReference().size() > 0)
                {
                    // module has references, dependency should be propagated
                    // to children modules; find module with `out` as an input
                    // and add references to their input
                    auto pred = [this, out](const ModuleInfo &n)
                    {
                        auto &in = n.input;
                        auto it  = std::find(in.begin(), in.end(), 
                                             env().getObjectAddress(out));
                        
                        return (it != in.end());
                    };
                    auto it = std::find_if(module_.begin(), module_.end(), pred);
                    while (it != module_.end())
                    {
                        for (auto &ref: getModule(address)->getReference())
                        {
                            it->input.push_back(env().getObjectAddress(ref));
                        }
                        it = std::find_if(++it, module_.end(), pred);
                    }   
                }
            }
        }
        // creating entry in database ------------------------------------------
        if (hasDatabase() and makeModuleDb_)
        {
            ModuleEntry e;

            e.moduleId     = address;
            e.name         = name;
            e.moduleTypeId = dbInsertModuleType(getModuleType(address));
            e.parameters   = getModule(address)->parString();
            db_->insert("modules", e);
        }
        graphOutdated_         = true;
    }
    else
    {
        HADRONS_ERROR(Definition, "module '" + name + "' already exists");
    }
}

unsigned int VirtualMachine::getNModule(void) const
{
    return module_.size();
}

void VirtualMachine::createModule(const std::string name, const std::string type,
                                  XmlReader &reader, const std::string blockName)
{
    auto &factory = ModuleFactory::getInstance();
    auto pt       = factory.create(type, name);
    
    pt->parseParameters(reader, blockName);
    pushModule(pt);
}

ModuleBase * VirtualMachine::getModule(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].data.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

ModuleBase * VirtualMachine::getModule(const std::string name) const
{
    return getModule(getModuleAddress(name));
}

unsigned int VirtualMachine::getModuleAddress(const std::string name) const
{
    if (hasModule(name))
    {
        return moduleAddress_.at(name);
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with name '" + name + "'");
    }
}

std::string VirtualMachine::getModuleName(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].name;
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const unsigned int address) const
{
    if (hasModule(address))
    {
        return typeName(module_[address].type);
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const std::string name) const
{
    return getModuleType(getModuleAddress(name));
}

std::string VirtualMachine::getModuleNamespace(const unsigned int address) const
{
    std::string type = getModuleType(address), ns;
    
    auto pos2 = type.rfind("::");
    auto pos1 = type.rfind("::", pos2 - 2);
    
    return type.substr(pos1 + 2, pos2 - pos1 - 2);
}

std::string VirtualMachine::getModuleNamespace(const std::string name) const
{
    return getModuleNamespace(getModuleAddress(name));
}

int VirtualMachine::getCurrentModule(void) const
{
    return currentModule_;
}

bool VirtualMachine::hasModule(const unsigned int address) const
{
    return (address < module_.size());
}

bool VirtualMachine::hasModule(const std::string name) const
{
    return (moduleAddress_.find(name) != moduleAddress_.end());
}

// print VM content ////////////////////////////////////////////////////////////
void VirtualMachine::printContent(void) const
{
    LOG(Debug) << "Modules: " << std::endl;
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getModuleName(i) << std::endl;
    }
}

// module graph ////////////////////////////////////////////////////////////////
Graph<unsigned int> VirtualMachine::getModuleGraph(void)
{
    if (graphOutdated_)
    {
        makeModuleGraph();
        graphOutdated_ = false;
    }

    return graph_;
}

void VirtualMachine::makeModuleGraph(void)
{
    Graph<unsigned int> graph;
    
    // create vertices
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        graph.addVertex(m);
    }
    // create edges
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            int min = env().getObjectModule(in);

            if (min < 0)
            {
                HADRONS_ERROR_REF(ObjectDefinition, "dependency '" 
                             + env().getObjectName(in) + "' (address " 
                             + std::to_string(in)
                             + ") is not produced by any module", in);
            }
            else
            {
                graph.addEdge(min, m);
            }
        }
    }
    graph_ = graph;
}

// dump GraphViz graph /////////////////////////////////////////////////////////
void VirtualMachine::dumpModuleGraph(std::ostream &out)
{
    makeModuleGraph();
    out << "digraph hadrons {" << std::endl;
    out << "node [shape=record, fontname=\"Courier\", fontsize=\"11\"];" << std::endl;
    out << "graph [fontname = \"Courier\", fontsize=\"11\"];" << std::endl;
    out << "edge [fontname = \"Courier\", fontsize=\"11\"];"<< std::endl;
    for (unsigned int m = 0; m < module_.size(); ++m)
    {

    }
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            int min = env().getObjectModule(in);

            out << min << " -> " << m << " [ label = \""
                << env().getObjectName(in) << "\" ];" << std::endl;
        }
    }
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        out <<  m << " [ label = \"{<f0> " << getModule(m)->getRegisteredName()
            << " |<f1> " << getModuleName(m) << "}\" ];" << std::endl;
    }
    out << "}\n" << std::endl;
}

void VirtualMachine::dumpModuleGraph(void)
{
    dumpModuleGraph(std::cout);
}

void VirtualMachine::dumpModuleGraph(const std::string filename)
{
    std::ofstream f(filename);

    dumpModuleGraph(f);
}

// memory profile //////////////////////////////////////////////////////////////
const VirtualMachine::MemoryProfile & VirtualMachine::getMemoryProfile(void)
{
    if (memoryProfileOutdated_)
    {
        makeMemoryProfile();
        memoryProfileOutdated_ = false;
    }

    return profile_;
}

void VirtualMachine::makeMemoryProfile(void)
{
    bool protect = env().objectsProtected();
    bool hmsg    = HadronsLogMessage.isActive();
    bool gmsg    = GridLogMessage.isActive();
    bool err     = HadronsLogError.isActive();
    auto program = getModuleGraph().topoSort();

    resetProfile();
    profile_.module.resize(getNModule());
    env().protectObjects(false);
    GridLogMessage.Active(false);
    HadronsLogMessage.Active(false);
    for (auto it = program.rbegin(); it != program.rend(); ++it) 
    {
        auto a = *it;

        if (profile_.module[a].empty())
        {
            LOG(Debug) << "Profiling memory for module '" << module_[a].name
                       << "' (" << a << ")" << std::endl;
            memoryProfile(a);
            env().freeAll();
        }
    }
    env().protectObjects(protect);
    GridLogMessage.Active(gmsg);
    HadronsLogMessage.Active(hmsg);
    if (hasDatabase() and makeObjectDb_)
    {
        for (unsigned int i = 0; i < profile_.object.size(); ++i)
        {
            ObjectEntry o;

            o.objectId     = i;
            o.name         = env().getObjectName(i);
            o.objectTypeId = dbInsertObjectType(env().getObjectDerivedType(i),
                                                env().getObjectType(i));
            o.size         = profile_.object[i].size;
            o.moduleId     = profile_.object[i].module;
            o.storageType  = profile_.object[i].storage;
            db_->insert("objects", o);
        }
    }
}

void VirtualMachine::printMemoryProfile(void) const
{
    LOG(Debug) << "Memory profile:" << std::endl;
    LOG(Debug) << "----------------" << std::endl;
    for (unsigned int a = 0; a < profile_.module.size(); ++a)
    {
        LOG(Debug) << getModuleName(a) << " (" << a << ")" << std::endl;
        for (auto &o: profile_.module[a])
        {
            LOG(Debug) << "|__ " << env().getObjectName(o.first) << " ("
                       << profile_.object[o.first].storage << " "
                       << sizeString(o.second) << ")" << std::endl;
        }
        LOG(Debug) << std::endl;
    }
    LOG(Debug) << "----------------" << std::endl;
}

void VirtualMachine::resetProfile(void)
{
    profile_.module.clear();
    profile_.object.clear();
}

void VirtualMachine::resizeProfile(void)
{
    if (env().getMaxAddress() > profile_.object.size())
    {
        MemoryPrint empty;

        empty.size   = 0;
        empty.module = -1;
        profile_.object.resize(env().getMaxAddress(), empty);
    }
}

void VirtualMachine::updateProfile(const unsigned int address)
{
    resizeProfile();
    for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
    {
        int envMod = env().getObjectModule(a);

        if ((env().hasCreatedObject(a) or (envMod == address)) and (profile_.object[a].module == -1))
        {
            profile_.object[a].module   = address;
            profile_.object[a].size     = env().getObjectSize(a);
            profile_.object[a].storage  = env().getObjectStorage(a);
            profile_.module[address][a] = profile_.object[a].size;
            if (envMod < 0)
            {
                env().setObjectModule(a, address);
            }
        }
    }
}

void VirtualMachine::cleanEnvironment(void)
{
    resizeProfile();
    for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
    {
        if (env().hasCreatedObject(a) and (profile_.object[a].module == -1))
        {
            env().freeObject(a);
        }
    }
}

void VirtualMachine::memoryProfile(const unsigned int address)
{
    auto m = getModule(address);

    LOG(Debug) << "Setting up module '" << m->getName() 
               << "' (" << address << ")" << std::endl;
    try
    {
        currentModule_ = address;
        m->setup();
        currentModule_ = -1;
        updateProfile(address);
    }
    catch (Exceptions::ObjectDefinition &exc)
    {
        cleanEnvironment();
        if (!env().hasCreatedObject(exc.getAddress()))
        {
            LOG(Debug) << "Object '" << env().getObjectName(exc.getAddress())
                       << "' missing for setup of '" << m->getName() 
                       << "' (" << address << ")" << std::endl;
            memoryProfile(env().getObjectModule(exc.getAddress()));
        }
        memoryProfile(address);
    }
}

void VirtualMachine::memoryProfile(const std::string name)
{
    memoryProfile(getModuleAddress(name));
}

// garbage collector ///////////////////////////////////////////////////////////
VirtualMachine::GarbageSchedule 
VirtualMachine::makeGarbageSchedule(const Program &p) const
{
    GarbageSchedule freeProg;
    
    freeProg.resize(p.size());
    for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
    {
        if (env().getObjectStorage(a) == Environment::Storage::temporary)
        {
            auto it = std::find(p.begin(), p.end(), env().getObjectModule(a));

            if (it != p.end())
            {
                freeProg[std::distance(p.begin(), it)].insert(a);
            }
        }
        else if (env().getObjectStorage(a) == Environment::Storage::standard)
        {
            auto pred = [a, this](const unsigned int b)
            {
                auto &in = module_[b].input;
                auto it  = std::find(in.begin(), in.end(), a);
                
                return (it != in.end()) or (b == env().getObjectModule(a));
            };
            auto it = std::find_if(p.rbegin(), p.rend(), pred);
            if (it != p.rend())
            {
                freeProg[std::distance(it, p.rend()) - 1].insert(a);
            }
        }
    }

    return freeProg;
}

// high-water memory function //////////////////////////////////////////////////
VirtualMachine::Size VirtualMachine::memoryNeeded(const Program &p)
{
    const MemoryProfile &profile = getMemoryProfile();
    GarbageSchedule     freep    = makeGarbageSchedule(p);
    Size                current = 0, max = 0;

    for (unsigned int i = 0; i < p.size(); ++i)
    {
        for (auto &o: profile.module[p[i]])
        {
            current += o.second;
        }
        max = std::max(current, max);
        for (auto &o: freep[i])
        {
            current -= profile.object[o].size;
        }
    }

    return max;
}

// genetic scheduler ///////////////////////////////////////////////////////////
VirtualMachine::Program VirtualMachine::schedule(const GeneticPar &par)
{
    typedef GeneticScheduler<Size, unsigned int> Scheduler;

    auto graph = getModuleGraph();

    //constrained topological sort using a genetic algorithm
    LOG(Message) << "Scheduling computation..." << std::endl;
    LOG(Message) << "               #module= " << graph.size() << std::endl;
    LOG(Message) << "       population size= " << par.popSize << std::endl;
    LOG(Message) << "       max. generation= " << par.maxGen << std::endl;
    LOG(Message) << "  max. cst. generation= " << par.maxCstGen << std::endl;
    LOG(Message) << "         mutation rate= " << par.mutationRate << std::endl;
    
    unsigned int          gen, prevPeak, nCstPeak = 0;
    std::random_device    rd;
    Scheduler::Parameters gpar;
    
    gpar.popSize      = par.popSize;
    gpar.mutationRate = par.mutationRate;
    gpar.seed         = rd();
    CartesianCommunicator::BroadcastWorld(0, &(gpar.seed), sizeof(gpar.seed));
    Scheduler::ObjFunc memPeak = [this](const Program &p)->Size
    {
        return memoryNeeded(p);
    };
    Scheduler scheduler(graph, memPeak, gpar);
    gen = 0;
    scheduler.initPopulation();
    LOG(Message) << "Start: " << sizeString(scheduler.getMinValue()) 
                 << std::endl;
    do
    {
        scheduler.nextGeneration();
        if (gen != 0)
        {
            if (prevPeak == scheduler.getMinValue())
            {
                nCstPeak++;
            }
            else
            {
                nCstPeak = 0;
            }
        }
        
        prevPeak = scheduler.getMinValue();
        if (gen % 10 == 0)
        {
            LOG(Message) << "Generation " << gen << ": "
                         << sizeString(scheduler.getMinValue()) << std::endl;
        }
        
        gen++;
    } while ((gen < par.maxGen) and (nCstPeak < par.maxCstGen));
    if (hasDatabase() and makeScheduleDb_)
    {
        Program p = scheduler.getMinSchedule();

        for (unsigned int i = 0; i < p.size(); ++i)
        {
            ScheduleEntry s;

            s.step     = i;
            s.moduleId = p[i];
            db_->insert("schedule", s);
        }
    }
    
    return scheduler.getMinSchedule();
}

// general execution ///////////////////////////////////////////////////////////
#define BIG_SEP   "================"
#define SEP       "----------------"
#define SMALL_SEP "................"

void VirtualMachine::executeProgram(const Program &p)
{
    Size            memPeak = 0, sizeBefore, sizeAfter;
    GarbageSchedule freeProg;
    
    // build garbage collection schedule
    LOG(Debug) << "Building garbage collection schedule..." << std::endl;
    freeProg = makeGarbageSchedule(p);
    for (unsigned int i = 0; i < freeProg.size(); ++i)
    {
        std::string msg = "";

        for (auto &a: freeProg[i])
        {
            msg += env().getObjectName(a) + " ";
        }
        msg += "]";
        LOG(Debug) << std::setw(4) << i + 1 << ": [" << msg << std::endl;
    }

    // program execution
    LOG(Debug) << "Executing program..." << std::endl;
    totalTime_ = GridTime::zero();
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        // execute module
        LOG(Message) << SEP << " Measurement step " << i + 1 << "/"
                     << p.size() << " (module '" << module_[p[i]].name
                     << "') " << SEP << std::endl;
        LOG(Message) << SMALL_SEP << " Module execution" << std::endl;
        currentModule_ = p[i];
        (*module_[p[i]].data)();
        currentModule_ = -1;
        sizeBefore = env().getTotalSize();
        // print time profile after execution
        LOG(Message) << SMALL_SEP << " Timings" << std::endl;

        std::map<std::string, GridTime> ctiming, gtiming;
        GridTime                        total;

        ctiming  = module_[p[i]].data->getTimings();
        total    = ctiming.at("_total");
        gtiming["total"]     = ctiming["_total"];   ctiming.erase("_total");
        gtiming["setup"]     = ctiming["_setup"];   ctiming.erase("_setup");
        gtiming["execution"] = ctiming["_execute"]; ctiming.erase("_execute");
        LOG(Message) << "* GLOBAL TIMERS" << std::endl;
        printTimeProfile(gtiming, total);
        if (!ctiming.empty())
        {
            LOG(Message) << "* CUSTOM TIMERS" << std::endl;
            printTimeProfile(ctiming, total);
        }
        timeProfile_[module_[p[i]].name] = total;
        totalTime_ += total;
        // print used memory after execution
        LOG(Message) << SMALL_SEP << " Memory management" << std::endl;
        MemoryUtils::printMemory();
        if (sizeBefore > memPeak)
        {
            memPeak = sizeBefore;
        }
        // garbage collection for step i
        LOG(Message) << "Garbage collection..." << std::endl;
        for (auto &j: freeProg[i])
        {
            env().freeObject(j);
        }
        // print used memory after garbage collection if necessary
        sizeAfter = env().getTotalSize();
        if (sizeBefore != sizeAfter)
        {
            MemoryUtils::printMemory();
        }
        else
        {
            LOG(Message) << "Nothing to free" << std::endl;
        }
    }
    // print total time profile
     LOG(Message) << SEP << " Measurement time profile" << SEP << std::endl;
     LOG(Message) << "Total measurement time: " << totalTime_ << " us" << std::endl;
     LOG(Message) << SMALL_SEP << " Module breakdown" << std::endl;
     printTimeProfile(timeProfile_, totalTime_);
}

void VirtualMachine::executeProgram(const std::vector<std::string> &p)
{
    Program pAddress;
    
    for (auto &n: p)
    {
        pAddress.push_back(getModuleAddress(n));
    }
    executeProgram(pAddress);
}

// generate result DB //////////////////////////////////////////////////////////
void VirtualMachine::generateResultDb(void)
{
    for (auto &m: module_)
    {
        m.data->generateResultDb();
    }
}
