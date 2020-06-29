/*
 * Application.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_Application_hpp_
#define Hadrons_Application_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/VirtualMachine.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Main program manager                               *
 ******************************************************************************/
class Application
{
public:
    // serializable classes for parameter input
    struct TrajRange: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };

    struct DatabasePar: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(DatabasePar,
                                        std::string, applicationDb,
                                        std::string, resultDb,
                                        bool,        restoreModules,
                                        bool,        restoreMemoryProfile,
                                        bool,        restoreSchedule,
                                        bool,        makeStatDb);
        DatabasePar(void): 
        restoreModules{false}, restoreMemoryProfile{false},
        restoreSchedule{false}, makeStatDb{false} {}
    };

    struct GlobalPar: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange,                  trajCounter,
                                        DatabasePar,                database,
                                        VirtualMachine::GeneticPar, genetic,
                                        std::string,                runId,
                                        std::string,                graphFile,
                                        std::string,                scheduleFile,
                                        bool,                       saveSchedule,
                                        int,                        parallelWriteMaxRetry);
        GlobalPar(void): parallelWriteMaxRetry{-1}, saveSchedule{false} {}
    };

    struct ObjectId: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(ObjectId,
                                        std::string, name,
                                        std::string, type);
    };

public:
    // constructors
    Application(void);
    Application(const GlobalPar &par);
    Application(const std::string parameterFileName);
    // destructor
    virtual ~Application(void);
    // access
    void              setPar(const GlobalPar &par);
    const GlobalPar & getPar(void);
    // module creation
    template <typename M>
    void createModule(const std::string name);
    template <typename M>
    void createModule(const std::string name, const typename M::Par &par);
    void createModule(const std::string name, const std::string type, XmlReader &reader);
    // module DB entry for result files
    template <typename EntryType>
    void setResultMetadata(const std::string moduleName, const std::string tableName, const EntryType &entry);
    // generate result DB
    void generateResultDb(void);
    // execute
    void run(void);
    // XML parameter file I/O
    void parseParameterFile(const std::string parameterFileName);
    void saveParameterFile(const std::string parameterFileName, unsigned int prec = 15);
    // schedule computation
    void schedule(void);
    void saveSchedule(const std::string filename);
    void loadSchedule(const std::string filename);
    void printSchedule(void);
    // loop on configurations
    void configLoop(void);
private:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // virtual machine shortcut
    DEFINE_VM_ALIAS;
private:
    long unsigned int       locVol_;
    std::string             parameterFileName_{""};
    GlobalPar               par_;
    VirtualMachine::Program program_;
    Database                db_, resultDb_;
    Grid::MemoryStats       memStats_;
    bool                    scheduled_{false}, loadedSchedule_{false};
};

/******************************************************************************
 *                     Application template implementation                    *
 ******************************************************************************/
// module creation /////////////////////////////////////////////////////////////
template <typename M>
void Application::createModule(const std::string name)
{
    vm().createModule<M>(name);
    scheduled_ = false;
}

template <typename M>
void Application::createModule(const std::string name,
                               const typename M::Par &par)
{
    vm().createModule<M>(name, par);
    scheduled_ = false;
}

// module DB entry /////////////////////////////////////////////////////////////
template <typename EntryType>
void Application::setResultMetadata(const std::string moduleName,
                                   const std::string tableName, 
                                   const EntryType &entry)
{
    auto m = vm().getModule(moduleName);

    m->setResultDbEntry(resultDb_, tableName, entry);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Application_hpp_
