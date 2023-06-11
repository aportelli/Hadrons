/*
 * Perambulator.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: felixerben <46817371+felixerben@users.noreply.github.com>
 * Author: felixerben <felix.erben@ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
* Author: Antonin Portelli <antonin.portelli@me.com>
* Author: Felix Erben <felix.erben@ed.ac.uk>
* Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
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

#ifndef Hadrons_MDistil_Perambulator_hpp_
#define Hadrons_MDistil_Perambulator_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/DistillationVectors.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE

#define START_P_TIMER(name) if (this) this->startTimer(name)
#define STOP_P_TIMER(name)  if (this) this->stopTimer(name)
// #define GET_TIMER(name)   ((this != nullptr) ? this->getDTimer(name) : 0.)

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                             Perambulator                                    *
 ******************************************************************************/

GRID_SERIALIZABLE_ENUM(pMode, undef, perambOnly, 0, inputSolve, 1, outputSolve, 2, saveSolveOnly, 3);

class PerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambulatorPar,
                                    std::string, lapEigenPack,
                                    std::string, solver,
                                    int,         sourceBatchSize,
                                    std::string, perambOutFileName,
                                    std::string, unsmSolveOutFileName,
                                    std::string, unsmSolve,
                                    std::string, distilNoise,
                                    std::string, timeSources,
                                    pMode, perambMode,
                                    std::string, nVec);
};

template <typename FImpl>
class TPerambulator: public Module<PerambulatorPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    // constructor
    TPerambulator(const std::string name);
    // destructor
    virtual ~TPerambulator(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(Perambulator, TPerambulator<FIMPL>, MDistil);
MODULE_REGISTER_TMP(ZPerambulator, TPerambulator<ZFIMPL>, MDistil);

/******************************************************************************
 *                 TPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPerambulator<FImpl>::TPerambulator(const std::string name) : Module<PerambulatorPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getInput(void)
{
    std::vector<std::string> in={par().lapEigenPack, par().distilNoise};
    pMode perambMode{par().perambMode};
    if(perambMode == pMode::inputSolve)
    {
        in.push_back(par().unsmSolve);
    }
    else
    {
        in.push_back(par().solver);
    }
    return in;
}


template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getOutput(void)
{
    std::vector<std::string> out{ getName() };
    pMode perambMode{par().perambMode};
    if(perambMode == pMode::outputSolve)
    {
        out.push_back( getName()+"_unsm_solve" );
    }
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::setup(void)
{
    GridCartesian * grid4d = envGetGrid(FermionField);
    GridCartesian * grid3d = envGetSliceGrid(FermionField,grid4d->Nd() -1);
    const int  Nt{env().getDim(Tdir)};
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);
    pMode perambMode{par().perambMode};
    // get nVec from DilutedNoise class, unless specified here. This is useful (and allowed) only in the inputSolve mode, 
    // where an already computed unsm solve can be recycled to compute a new perambulator with a smaller nVec.   
    int nVec=0;
    if(par().nVec.empty())
    {
        nVec = dilNoise.getNl();
    }
    else
    {
        nVec = std::stoi(par().nVec);
    }
    if(nVec > dilNoise.getNl())
    {
        HADRONS_ERROR(Argument, "nVec cannot be larger than the one specified in DilutedNoise");
    }
    if(perambMode != pMode::inputSolve && nVec < dilNoise.getNl())
    {
        HADRONS_ERROR(Argument, "only perambMode = inputSolve supports a different nVec to the one specified in DilutedNoise");
    }
    // If we run with reduced nVec we still need the same DilutedNoise object, we have to keep track of two different values of nDL:
    // the one used for the original unsm solve, and the one used in this module execution 
    int nDL_reduced=nDL;
    if(nDL>nVec)
    {
        nDL_reduced=nVec;
    }
    // Read in and verify the format of the vector of time sources used in this module execution. 
    // It must be a subset of available time dilution indices.
    int nSourceT;
    std::string sourceT = par().timeSources;
    nSourceT = verifyTimeSourcesInput(sourceT, nDT);

    int sourceBatchSize = par().sourceBatchSize;
    int nSources = nNoise * nDL_reduced * nDS * nSourceT;
    if(nSources % sourceBatchSize != 0)
    {
        HADRONS_ERROR(Argument, "sourceBatchSize must divide nNoise * nD");
    }
    //need a vector here
    envTmp(std::vector<FermionField>,  "dist_source_vec",  1, sourceBatchSize, envGetGrid(FermionField));
    envTmp(std::vector<FermionField>,  "fermion4dtmp_vec", 1, sourceBatchSize, envGetGrid(FermionField));

    // Perambulator dimensions need to use the reduced value for nDL     
    envCreate(PerambTensor, getName(), 1, Nt, nVec, nDL_reduced, nNoise, nSourceT, nDS);
    envTmp(PerambIndexTensor, "PerambDT",1,Nt,nVec,nDL_reduced,nNoise,nDS);
    if(perambMode == pMode::outputSolve)
    {
        LOG(Message)<< "setting up output field for unsm solves to " << getName()+"_unsm_solve" << std::endl;
        if(!par().unsmSolveOutFileName.empty())
        {
            LOG(Message)<< "also saving solves to stem '" << par().unsmSolveOutFileName << "'" << std::endl;
        }
        envCreate(std::vector<FermionField>, getName()+"_unsm_solve", 1, nNoise*nDL*nDS*nSourceT,
        envGetGrid(FermionField));
    }
    else if(perambMode == pMode::inputSolve)
    {
        LOG(Message)<< "setting up input solve from diluted noise '" << par().distilNoise << "' and reduced nVec=" << par().nVec <<  std::endl;
    }

    envTmp(FermionField,         "fermion3dtmp", 1, grid3d);
    envTmpLat(ColourVectorField, "cv4dtmp");
    envTmp(ColourVectorField,    "cv3dtmp", 1, grid3d);
    envTmp(ColourVectorField,    "evec3d",  1, grid3d);
    const int Ntlocal{grid4d->LocalDimensions()[3]};
    envTmp(std::vector<ColourVectorField>,    "evec3d_tmp",  1, Ntlocal * nVec,  grid3d);

    // No solver needed if an already existing solve is recycled    
    if(perambMode != pMode::inputSolve)
    {
        Ls_ = env().getObjectLs(par().solver);
        envTmp(std::vector<FermionField>, "v5dtmp_vec", Ls_, sourceBatchSize, envGetGrid(FermionField, Ls_));
        envTmp(std::vector<FermionField>, "v5dtmp_sol_vec", Ls_,  sourceBatchSize, envGetGrid(FermionField, Ls_));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::execute(void)
{
    const int Nt{env().getDim(Tdir)};
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nVec=0;
    if(par().nVec.empty())
    {
        nVec = dilNoise.getNl();
    }
    else
    {
        nVec = std::stoi(par().nVec);
    }
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    int nD = nDL * nDS * nDT;
    // If we run with reduced nVec we still need the same DilutedNoise object, but a smaller laplacian dilution dimension
    int nDL_reduced=nDL;
    if(nDL>nVec)
    {
        nDL_reduced=nVec;
    }
    auto &perambulator = envGet(PerambTensor, getName());
    auto &epack = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEigenPack);

    pMode perambMode{par().perambMode};
    LOG(Message)<< "Mode " << perambMode << std::endl;
    LOG(Message)<< "Source batch size = " << par().sourceBatchSize << std::endl;

    envGetTmp(FermionField,      fermion3dtmp);
    envGetTmp(ColourVectorField, cv4dtmp);
    envGetTmp(ColourVectorField, cv3dtmp);
    envGetTmp(ColourVectorField, evec3d);
    GridCartesian * grid4d = envGetGrid(FermionField);
    GridCartesian * grid3d = envGetSliceGrid(FermionField,grid4d->Nd() -1);
    const int Ntlocal{grid4d->LocalDimensions()[3]};
    const int Ntfirst{grid4d->LocalStarts()[3]};

    envGetTmp(std::vector<ColourVectorField>, evec3d_tmp);
    for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
    {
        for (int ivec = 0; ivec < nVec; ivec++)
        {
            int jvec= ivec + nVec * (t-Ntfirst);
            ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
            evec3d_tmp[jvec] = evec3d;
        }
    }


    std::string sourceT = par().timeSources;
    std::vector<int> invT;
    getSourceTimesFromInput(sourceT,nDT,dilNoise,invT);
    perambulator.MetaData.timeSources = invT;
    
    int sourceBatchSize = par().sourceBatchSize;
    envGetTmp(std::vector<FermionField>, dist_source_vec);
    envGetTmp(std::vector<FermionField>, fermion4dtmp_vec);

    int idt,dt,dk,ds,dIndexSolve = 0; 
    std::array<unsigned int, 3> index;
    int iSource=0;
    std::vector<int> sourceIndices(sourceBatchSize);
    for (int inoise = 0; inoise < nNoise; inoise++)
    {
        for (int d = 0; d < nD; d++)
        {
            // create batched sources
            // save all the inoise, d indices looped through to unroll them later
            index = dilNoise.dilutionCoordinates(d);
            dt = index[DistillationNoise<FImpl>::Index::t];
            dk = index[DistillationNoise<FImpl>::Index::l];
            // Skip laplacian dilution indices which are larger than (reduced) number of eigenvectors used for the perambulator 
            if(dk>=nDL_reduced)
            {
                continue;
            }
            ds = index[DistillationNoise<FImpl>::Index::s];
            std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
            // Skip dilution indices which are not in invT
            if(it == std::end(invT))
            {
                continue;
            }
            idt=it - std::begin(invT);
            // there is no reason to use this batched, but it works just as well
            if(perambMode == pMode::inputSolve)
            {
                START_P_TIMER("input solve");
                auto &solveIn = envGet(std::vector<FermionField>, par().unsmSolve);
                // Index of the solve just has the reduced time dimension & uses nDL from solveIn
                dIndexSolve = ds + nDS * dk + nDL * nDS * idt;
                fermion4dtmp_vec[iSource] = solveIn[inoise+nNoise*dIndexSolve];
                STOP_P_TIMER("input solve");
                LOG(Message) << "re-using source vector: noise " << inoise << " dilution (d_t,d_k,d_alpha) : (" << dt << ","<< dk << "," << ds << ")" << std::endl;
            } 
            else 
            {
                // Fill batched vector of distillation sources
                // also set solution batch vector to zero
                dist_source_vec[iSource] = dilNoise.makeSource(d,inoise);
                fermion4dtmp_vec[iSource]=0;
            }
            sourceIndices[iSource]=inoise+nNoise*d;
            iSource++;
            // if batch size not yet reached, continue batching
            if(iSource < sourceBatchSize)
            {
                continue;
            }
            LOG(Message) <<  "batch built from indices: " << sourceIndices << std::endl;
            if(perambMode != pMode::inputSolve)
            {
                // Solve on batched vector
                auto &solver=envGet(Solver, par().solver);
                auto &mat = solver.getFMat();
                if (Ls_ == 1)
                {
                    START_P_TIMER("solver");
                    solver(fermion4dtmp_vec, dist_source_vec);
                    STOP_P_TIMER("solver");
                }
                else
                {
                    START_P_TIMER("solver handling");
                    envGetTmp(std::vector<FermionField>,      v5dtmp_vec);
                    envGetTmp(std::vector<FermionField>,      v5dtmp_sol_vec);
                    for (iSource = 0; iSource < sourceBatchSize; iSource ++)
                    {  
                        mat.ImportPhysicalFermionSource(dist_source_vec[iSource], v5dtmp_vec[iSource]);
                    }
                    STOP_P_TIMER("solver handling");
                    START_P_TIMER("solver");
                    solver(v5dtmp_sol_vec, v5dtmp_vec);
                    STOP_P_TIMER("solver");
                    START_P_TIMER("solver handling");
                    for (iSource = 0; iSource < sourceBatchSize; iSource ++)
                    {  
                        mat.ExportPhysicalFermionSolution(v5dtmp_sol_vec[iSource], fermion4dtmp_vec[iSource]);
                    }
                    STOP_P_TIMER("solver handling");
                }
                if(perambMode == pMode::outputSolve)
                {
                    for (iSource = 0; iSource < sourceBatchSize; iSource ++)
                    {
                        START_P_TIMER("output solve");
                        int in = sourceIndices[iSource] % nNoise;
                        int id = sourceIndices[iSource] / nNoise;
                        index = dilNoise.dilutionCoordinates(id);
                        dt = index[DistillationNoise<FImpl>::Index::t];
                        dk = index[DistillationNoise<FImpl>::Index::l];
                        ds = index[DistillationNoise<FImpl>::Index::s];
                        std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
                        idt=it - std::begin(invT);
                        // Index of the solve just has the reduced time dimension                         
                        dIndexSolve = ds + nDS * dk + nDL * nDS * idt;
                        auto &solveOut = envGet(std::vector<FermionField>, getName()+"_unsm_solve");
                        solveOut[in+nNoise*dIndexSolve] = fermion4dtmp_vec[iSource];                     
                        STOP_P_TIMER("output solve");
                    }
                }
                if(perambMode == pMode::saveSolveOnly or (perambMode == pMode::outputSolve and !par().unsmSolveOutFileName.empty()))
                {
                    for (iSource = 0; iSource < sourceBatchSize; iSource ++)
                    {
                        START_P_TIMER("save solve");
                        int in = sourceIndices[iSource] % nNoise;
                        int id = sourceIndices[iSource] / nNoise;
                        index = dilNoise.dilutionCoordinates(id);
                        dt = index[DistillationNoise<FImpl>::Index::t];
                        dk = index[DistillationNoise<FImpl>::Index::l];
                        ds = index[DistillationNoise<FImpl>::Index::s];
                        std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
                        idt=it - std::begin(invT);
                        // Index of the solve has the full time dimension
                        dIndexSolve = dilNoise.dilutionIndex(dt,dk,ds);
                        std::string sFileName(par().unsmSolveOutFileName);
                        sFileName.append("_noise");
                        sFileName.append(std::to_string(in));
                        DistillationVectorsIo::writeComponent(sFileName, fermion4dtmp_vec[iSource], "unsmSolve", nNoise, nDL, nDS, nDT, invT, in+nNoise*dIndexSolve, vm().getTrajectory());
                        STOP_P_TIMER("save solve");
                    }
                }
            }
            START_P_TIMER("perambulator computation");
            for (iSource = 0; iSource < sourceBatchSize; iSource ++)
            {
                int in = sourceIndices[iSource] % nNoise;
                int id = sourceIndices[iSource] / nNoise;
                index = dilNoise.dilutionCoordinates(id);
                dt = index[DistillationNoise<FImpl>::Index::t];
                dk = index[DistillationNoise<FImpl>::Index::l];
                ds = index[DistillationNoise<FImpl>::Index::s];
                std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
                idt=it - std::begin(invT);
                for (int is = 0; is < Ns; is++)
                {
                    cv4dtmp = peekSpin(fermion4dtmp_vec[iSource],is);
                    for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
                    {
                        ExtractSliceLocal(cv3dtmp,cv4dtmp,0,t-Ntfirst,Tdir); 
                        for (int ivec = 0; ivec < nVec; ivec++)
                        {
                            int jvec= ivec + nVec * (t-Ntfirst);
                            evec3d = evec3d_tmp[jvec];
                            pokeSpin(perambulator.tensor(t, ivec, dk, in,idt,ds),static_cast<Complex>(innerProduct(evec3d, cv3dtmp)),is);
                        }
                    }
                }
            }
            STOP_P_TIMER("perambulator computation");
            iSource=0;
        }
    }

    // Now share my timeslice data with other members of the grid
    const int NumSlices{grid4d->_processors[Tdir] / grid3d->_processors[Tdir]};
    if (NumSlices > 1)
    {
        START_P_TIMER("perambulator sharing");
        LOG(Debug) <<  "Sharing perambulator data with other nodes" << std::endl;
        const int MySlice {grid4d->_processor_coor[Tdir]};
        const int TensorSize {static_cast<int>(perambulator.tensor.size() * PerambTensor::Traits::count)};
        if (TensorSize != perambulator.tensor.size() * PerambTensor::Traits::count)
        {
            HADRONS_ERROR(Range, "peramb size overflow");
        }
        const int SliceCount {TensorSize/NumSlices};
        using InnerScalar = typename PerambTensor::Traits::scalar_type;
        InnerScalar * const PerambData {EigenIO::getFirstScalar( perambulator.tensor )};
        // Zero data for all timeslices other than the slice computed by 3d boss nodes
        for (int Slice = 0 ; Slice < NumSlices ; ++Slice)
        {
            if (!grid3d->IsBoss() || Slice != MySlice)
            {
                InnerScalar * const SliceData {PerambData + Slice * SliceCount};
                for (int j = 0 ; j < SliceCount ; ++j)
                {
                    SliceData[j] = 0.;
                }
            }
        }
        grid4d->GlobalSumVector(PerambData, TensorSize);
        STOP_P_TIMER("perambulator sharing");
    }

    // Save the perambulator to disk from the boss node
    if (grid4d->IsBoss() && !par().perambOutFileName.empty())
    {
        START_P_TIMER("perambulator io");
        envGetTmp(PerambIndexTensor, PerambDT);
        std::vector<std::string> nHash = dilNoise.generateHash();
        PerambDT.MetaData.noiseHashes = nHash;
        PerambDT.MetaData.Version = "0.1";
        for (int dt = 0; dt < Nt; dt++)
        {
            std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
            // Skip dilution indices which are not in invT
            if(it == std::end(invT))
            {
                continue;
            }
            LOG(Message) <<  "saving perambulator dt= " << dt << std::endl;
            idt=it - std::begin(invT);
            std::string sPerambName {par().perambOutFileName};
            sPerambName.append(".");
            sPerambName.append(std::to_string(vm().getTrajectory()));
            sPerambName.append("/iDT_");
            sPerambName.append(std::to_string(dt));
            sPerambName.append(".");
            sPerambName.append(std::to_string(vm().getTrajectory()));
            makeFileDir(sPerambName, grid4d);
            for (int t = 0; t < Nt; t++)
            for (int ivec = 0; ivec < nVec; ivec++)
            for (int idl = 0; idl < nDL_reduced; idl++)
            for (int in = 0; in < nNoise; in++)
            for (int ids = 0; ids < nDS; ids++)
            {
                PerambDT.tensor(t,ivec,idl,in,ids) = perambulator.tensor(t,ivec,idl,in,idt,ids);
            }
            PerambDT.MetaData.timeDilutionIndex = dt;
            PerambDT.write(sPerambName.c_str());
        }
        STOP_P_TIMER("perambulator io");
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Perambulator_hpp_
