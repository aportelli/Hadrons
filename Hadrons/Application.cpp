/*
 * Application.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#include <Hadrons/Application.hpp>
#include <Hadrons/GeneticScheduler.hpp>
#include <Hadrons/StatLogger.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

#define BIG_SEP "================"
#define SEP     "----------------"

/******************************************************************************
 *                       Application implementation                           *
 ******************************************************************************/
// constructors ////////////////////////////////////////////////////////////////
#define MACOUT(macro)    macro              << " (" << #macro << ")"
#define MACOUTS(macro) HADRONS_STR(macro) << " (" << #macro << ")"

Application::Application(void)
{
    initLogger();
    Grid::MemoryProfiler::stats = &memStats_;
    auto dim = GridDefaultLatt(), mpi = GridDefaultMpi(), loc(dim);

    if (dim.size())
    {
        locVol_ = 1;
        for (unsigned int d = 0; d < dim.size(); ++d)
        {
            loc[d]  /= mpi[d];
            locVol_ *= loc[d];
        }
        LOG(Message) << "====== HADRONS APPLICATION INITIALISATION ======" << std::endl;
        LOG(Message) << "** Dimensions" << std::endl;
        LOG(Message) << "Global lattice: " << dim << std::endl;
        LOG(Message) << "MPI partition : " << mpi << std::endl;
        LOG(Message) << "Local lattice : " << loc << std::endl;
        LOG(Message) << std::endl;
        LOG(Message) << "** Default parameters (and associated C macros)" << std::endl;
        LOG(Message) << "ASCII output precision  : " << MACOUT(DEFAULT_ASCII_PREC) << std::endl;
        LOG(Message) << "Fermion implementation  : " << MACOUTS(FIMPLBASE) << std::endl;
        LOG(Message) << "z-Fermion implementation: " << MACOUTS(ZFIMPLBASE) << std::endl;
        LOG(Message) << "Scalar implementation   : " << MACOUTS(SIMPLBASE) << std::endl;
        LOG(Message) << "Gauge implementation    : " << MACOUTS(GIMPLBASE) << std::endl;
        LOG(Message) << "Eigenvector base size   : " 
                     << MACOUT(HADRONS_DEFAULT_LANCZOS_NBASIS) << std::endl;
        LOG(Message) << "Schur decomposition     : " << MACOUTS(HADRONS_DEFAULT_SCHUR) << std::endl;
        LOG(Message) << std::endl;
    }
}

Application::~Application(void)
{
    Grid::MemoryProfiler::stats = nullptr;
}

Application::Application(const Application::GlobalPar &par)
: Application()
{
    setPar(par);
}

Application::Application(const std::string parameterFileName)
: Application()
{
    parameterFileName_ = parameterFileName;
}

// access //////////////////////////////////////////////////////////////////////
void Application::setPar(const Application::GlobalPar &par)
{
    par_ = par;
    if (!getPar().database.applicationDb.empty())
    {
        LOG(Message) << "Connecting to application database in file '" 
                     << getPar().database.applicationDb << "'..." << std::endl;
        db_.setFilename(getPar().database.applicationDb, isGridInit() ? env().getGrid() : nullptr);
        vm().setDatabase(db_);
        if (getPar().database.restoreMemoryProfile)
        {
            vm().dbRestoreMemoryProfile();
            LOG(Message) << "Memory profile restored from application database" << std::endl;
        }
        if (getPar().database.restoreModules)
        {
            vm().dbRestoreModules();
            LOG(Message) << "Modules restored from application database" << std::endl;
        }
        if (getPar().database.restoreSchedule)
        {
            program_ = vm().dbRestoreSchedule();
            loadedSchedule_ = true;
            scheduled_      = true;
            LOG(Message) << "Schedule restored from application database" << std::endl;
        }
    }
    if (!getPar().database.resultDb.empty())
    {
        LOG(Message) << "Connecting to result database in file '" 
                     << getPar().database.resultDb << "'..." << std::endl;
        resultDb_.setFilename(getPar().database.resultDb, isGridInit() ? env().getGrid() : nullptr);
    }
}

const Application::GlobalPar & Application::getPar(void)
{
    return par_;
}

// module creation /////////////////////////////////////////////////////////////
void Application::createModule(const std::string name, const std::string type, 
                               XmlReader &reader)
{
    vm().createModule(name, type, reader);
}

// generate result DB //////////////////////////////////////////////////////////
void Application::generateResultDb(void)
{
    auto range = par_.trajCounter;
    
    for (unsigned int t = range.start; t < range.end; t += range.step)
    {
        vm().setTrajectory(t);
        vm().generateResultDb();
    }
}

// execute /////////////////////////////////////////////////////////////////////
void Application::run(void)
{
    Database   statDb;
    StatLogger statLogger;

    LOG(Message) << "====== HADRONS APPLICATION START ======" << std::endl;
    if (!parameterFileName_.empty() and (vm().getNModule() == 0))
    {
        parseParameterFile(parameterFileName_);
    }
    if (getPar().runId.empty())
    {
        HADRONS_ERROR(Definition, "run id is empty");
    }
    LOG(Message) << "RUN ID '" << getPar().runId << "'" << std::endl;
    BinaryIO::latticeWriteMaxRetry = getPar().parallelWriteMaxRetry;
    LOG(Message) << "Attempt(s) for resilient parallel I/O: " 
                 << BinaryIO::latticeWriteMaxRetry << std::endl;
    vm().setRunId(getPar().runId);
    if (!getPar().database.statDbBase.empty())
    {
        std::string        statDbFilename;
        std::ostringstream oss;
        auto now      = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        auto nowLocal = *std::localtime(&now);

        oss << std::put_time(&nowLocal, "%Y%m%d-%H%M%S");
        statDbFilename = getPar().database.statDbBase + getPar().runId + "-stat-" + oss.str() + ".db";
        LOG(Message) << "Logging run statistics in '" << statDbFilename << "'" << std::endl;
        if (env().getGrid()->IsBoss())
        {
            statDb.setFilename(statDbFilename);
            statLogger.setDatabase(statDb);
            statLogger.start(500);
        }
    }
    if (getPar().saveSchedule or getPar().scheduleFile.empty())
    {
        schedule();
        if (getPar().saveSchedule)
        {
            std::string filename;

            filename = (getPar().scheduleFile.empty()) ? 
                         "hadrons.sched" : getPar().scheduleFile;
            saveSchedule(filename);
        }
    }
    else
    {
        loadSchedule(getPar().scheduleFile);
    }
    printSchedule();
    vm().printMemoryProfile();
    if (!getPar().graphFile.empty())
    {
        makeFileDir(getPar().graphFile, env().getGrid());
        vm().dumpModuleGraph(getPar().graphFile);
    }
    configLoop();
    if (!getPar().database.statDbBase.empty() and env().getGrid()->IsBoss())
    {
        statLogger.stop();
    }
}

// parse parameter file ////////////////////////////////////////////////////////
void Application::parseParameterFile(const std::string parameterFileName)
{
    XmlReader reader(parameterFileName, false, HADRONS_XML_TOPLEV);
    GlobalPar par;
    ObjectId  id;
    
    LOG(Message) << "Building application from '" << parameterFileName << "'..." << std::endl;
    read(reader, "parameters", par);
    setPar(par);
    if (!par.database.restoreModules)
    {
        if (!push(reader, "modules"))
        {
            HADRONS_ERROR(Parsing, "Cannot open node 'modules' in parameter file '" 
                                + parameterFileName + "'");
        }
        if (!push(reader, "module"))
        {
            HADRONS_ERROR(Parsing, "Cannot open node 'modules/module' in parameter file '" 
                                + parameterFileName + "'");
        }
        do
        {
            read(reader, "id", id);
            createModule(id.name, id.type, reader);
        } while (reader.nextElement("module"));
        pop(reader);
        pop(reader);
    }
    else
    {
        LOG(Message) << "XML module list ignored (restored from database '"
                     << par.database.applicationDb << "')" << std::endl;
    }
}

void Application::saveParameterFile(const std::string parameterFileName, unsigned int prec)
{
    LOG(Message) << "Saving application to '" << parameterFileName << "'..." << std::endl;
    if (env().getGrid()->IsBoss())
    {
        XmlWriter          writer(parameterFileName, HADRONS_XML_TOPLEV);
        ObjectId           id;
        const unsigned int nMod = vm().getNModule();

        writer.setPrecision(prec);
        write(writer, "parameters", getPar());
        push(writer, "modules");
        for (unsigned int i = 0; i < nMod; ++i)
        {
            push(writer, "module");
            id.name = vm().getModuleName(i);
            id.type = vm().getModule(i)->getRegisteredName();
            write(writer, "id", id);
            vm().getModule(i)->saveParameters(writer, "options");
            pop(writer);
        }
        pop(writer);
        pop(writer);
    }
}

// schedule computation ////////////////////////////////////////////////////////
void Application::schedule(void)
{
    if (!scheduled_ and !loadedSchedule_)
    {
        program_   = vm().schedule(par_.genetic);
        scheduled_ = true;
    }
}

void Application::saveSchedule(const std::string filename)
{
    LOG(Message) << "Saving current schedule to '" << filename << "'..."
                 << std::endl;
    if (env().getGrid()->IsBoss())
    {
        TextWriter               writer(filename);
        std::vector<std::string> program;
        
        if (!scheduled_)
        {
            HADRONS_ERROR(Definition, "Computation not scheduled");
        }

        for (auto address: program_)
        {
            program.push_back(vm().getModuleName(address));
        }
        write(writer, "schedule", program);
    }
}

void Application::loadSchedule(const std::string filename)
{
    TextReader               reader(filename);
    std::vector<std::string> program;
    
    LOG(Message) << "Loading schedule from '" << filename << "'..."
                 << std::endl;
    read(reader, "schedule", program);
    program_.clear();
    for (auto &name: program)
    {
        program_.push_back(vm().getModuleAddress(name));
    }
    loadedSchedule_ = true;
    scheduled_      = true;
}

void Application::printSchedule(void)
{
    if (!scheduled_ and !loadedSchedule_)
    {
        HADRONS_ERROR(Definition, "Computation not scheduled");
    }
    auto peak = vm().memoryNeeded(program_);
    LOG(Message) << "Schedule (memory needed: " << sizeString(peak) << "):"
                 << std::endl;
    for (unsigned int i = 0; i < program_.size(); ++i)
    {
        LOG(Message) << std::setw(4) << i + 1 << ": "
                     << vm().getModuleName(program_[i]) << std::endl;
    }
}

// loop on configurations //////////////////////////////////////////////////////
void Application::configLoop(void)
{
    auto range = par_.trajCounter;
    
    for (unsigned int t = range.start; t < range.end; t += range.step)
    {
        LOG(Message) << BIG_SEP << " Starting measurement for trajectory " << t
                     << " " << BIG_SEP << std::endl;
        vm().setTrajectory(t);
        vm().executeProgram(program_);
    }
    LOG(Message) << BIG_SEP << " End of measurement " << BIG_SEP << std::endl;
    env().freeAll();
}
