/*
 * Perambulator.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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

#include <Hadrons/Modules/MDistil/Distil.hpp>
#include <Hadrons/DilutedNoise.hpp>


BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                             Perambulator                                    *
 ******************************************************************************/



class PerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambulatorPar,
                                    std::string, lapEigenPack,
                                    std::string, solver,
                                    std::string, perambFileName,
                                    std::string, unsmearedSolveFileName,
                                    std::string, unsmearedSolve,
                                    std::string, distilNoise,
                                    std::string, timeSources,
                                    pMode, perambMode,
                                    int, nVec);
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
    std::vector<std::string> out={par().lapEigenPack, par().solver, par().distilNoise};
    pMode perambMode{par().perambMode};
    if(perambMode == pMode::inputSolve)
    {
        out.push_back(par().unsmearedSolve);
    }
    return out;
}


template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getOutput(void)
{
    // Always return perambulator with name of module
    std::string objName{ getName() };
    std::vector<std::string> output{ objName };
    pMode perambMode{par().perambMode};
    if(perambMode == pMode::outputSolve)
    {
        LOG(Message)<< "unsmeared solves are an output" << std::endl;
        objName.append( "_unsmeared_solve" );
        output.push_back( objName );
    }
    return output;
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
    int nSourceT;
    std::string sourceT = par().timeSources;
    if(par().timeSources.empty())
    {
	nSourceT=nDT;
    }
    else
    {
	// check whether input is legal, i.e. a number of integers between 0 and (nDT-1)
	std::regex rex("[0-9 ]+");
	std::smatch sm;
	std::regex_match(sourceT, sm, rex);
	if (!sm[0].matched)
	{
	    HADRONS_ERROR(Range, "sourceTimes must be list of non-negative integers");
	}
	std::istringstream is(sourceT);
	std::vector<int> iT ( ( std::istream_iterator<int>( is )  ), (std::istream_iterator<int>() ) );
	nSourceT = iT.size();
        for (int ii = 0; ii < nSourceT; ii++)
	{
	    if (iT[ii] >= nDT)
	    {
		HADRONS_ERROR(Range, "elements of sourceTimes must lie between 0 and nDT");
	    }
	}
    }
	    
    std::string objName{ getName() };
    envCreate(PerambTensor, objName, 1, Nt, par().nVec, nDL, nNoise, nSourceT, nDS);
    pMode perambMode{par().perambMode};
    if(perambMode == pMode::outputSolve)
    {
        LOG(Message)<< "setting up output field for unsmeared solves" << std::endl;
        objName.append( "_unsmeared_solve" );
        envCreate(std::vector<FermionField>, objName, 1, nNoise*nDL*nDS*nSourceT,
                  envGetGrid(FermionField));
    }
    
    envTmpLat(FermionField,      "dist_source");
    envTmpLat(FermionField,      "fermion4dtmp");
    envTmp(FermionField,         "fermion3dtmp", 1, grid3d);
    envTmpLat(ColourVectorField, "cv4dtmp");
    envTmp(ColourVectorField,    "cv3dtmp", 1, grid3d);
    envTmp(ColourVectorField,    "evec3d",  1, grid3d);
    
    Ls_ = env().getObjectLs(par().solver);
    envTmpLat(FermionField, "v5dtmp", Ls_);
    envTmpLat(FermionField, "v5dtmp_sol", Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::execute(void)
{
    const int Nt{env().getDim(Tdir)};
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();	
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);	
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);	
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);	
    int nD = nDL * nDS * nDT;

    auto &solver=envGet(Solver, par().solver);
    auto &mat = solver.getFMat();
    envGetTmp(FermionField, v5dtmp);
    envGetTmp(FermionField, v5dtmp_sol);
    std::string objName{ getName() };
    auto &perambulator = envGet(PerambTensor, objName);
    auto &epack = envGet(LapEvecs, par().lapEigenPack);
    
    objName.append( "_unsmeared_solve" );
    envGetTmp(FermionField,      dist_source);
    envGetTmp(FermionField,      fermion4dtmp);
    envGetTmp(FermionField,      fermion3dtmp);
    envGetTmp(ColourVectorField, cv4dtmp);
    envGetTmp(ColourVectorField, cv3dtmp);
    envGetTmp(ColourVectorField, evec3d);
    GridCartesian * grid4d = envGetGrid(FermionField);
    GridCartesian * grid3d = envGetSliceGrid(FermionField,grid4d->Nd() -1);
    const int Ntlocal{grid4d->LocalDimensions()[3]};
    const int Ntfirst{grid4d->LocalStarts()[3]};

    pMode perambMode{par().perambMode};
    LOG(Message)<< "Mode " << perambMode << std::endl;

    std::vector<FermionField> solveIn;
    if(perambMode == pMode::inputSolve)
    {
        LOG(Message) << "unsmeared solves are an input" << std::endl;
        solveIn         = envGet(std::vector<FermionField>, par().unsmearedSolve);
    }
    
    std::string sourceT = par().timeSources;
    int nSourceT;
    std::vector<int> invT;
    if(par().timeSources.empty())
    {
	// create sourceTimes all time-dilution indices
	nSourceT=nDT;
        for (int dt = 0; dt < nDT; dt++)
	{
	    std::vector<unsigned int> sT = dilNoise.timeSlices(dt);
	    perambulator.MetaData.sourceTimes.push_back(sT);
	    invT.push_back(dt);
	}
        LOG(Message) << "Computing inversions on all " << nDT << " time-dilution vectors" << std::endl;
    }
    else
    {
	std::istringstream is(sourceT);
	std::vector<int> iT ( ( std::istream_iterator<int>( is )  ), (std::istream_iterator<int>() ) );
	nSourceT = iT.size();
	// create sourceTimes from the chosen subset of time-dilution indices
        for (int dt = 0; dt < nSourceT; dt++)
	{
	    std::vector<unsigned int> sT = dilNoise.timeSlices(iT[dt]);
	    perambulator.MetaData.sourceTimes.push_back(sT);
	    invT.push_back(iT[dt]);
	}
        LOG(Message) << "Computing inversions on a subset of " << nSourceT << " time-dilution vectors" << std::endl;
    }
    LOG(Message) << "Source times" << perambulator.MetaData.sourceTimes << std::endl;

    int nVecFullSolve = 0;
    if(perambMode == pMode::inputSolve)
    {
        nVecFullSolve   = solveIn.size()/nNoise/nDS/nSourceT;
        LOG(Message) << "Using solve originally computed with Nvec=" << nVecFullSolve << std::endl;
    }

    for (int inoise = 0; inoise < nNoise; inoise++)
    {
        for (int d = 0; d < nD; d++)
        {
            std::array<unsigned int, 3> index = dilNoise.dilutionCoordinates(d);
   	    int dt = index[0];
   	    int dk = index[1];
   	    int ds = index[2];
   	    int idt = 0; 
            bool dtExists = std::find(std::begin(invT), std::end(invT), dt) != std::end(invT);
   	    if(!dtExists)
   	    {
   	        //skip dilution indices which are not in invT
   	        continue;
   	    }
   	    idt+=1;
   	    // index of the solve just has the reduced time dimension 
   	    int dIndexSolve = ds + nDS * dk + nVecFullSolve * nDL * idt;
            if(perambMode == pMode::inputSolve)
   	    {
                fermion4dtmp = solveIn[inoise+nNoise*dIndexSolve];
   	    } 
   	    else 
   	    {
                LOG(Message) <<  "LapH source vector from noise " << inoise << " and dilution component (d_t,d_k,d_alpha) : (" << dt << ","<< dk << "," << ds << ")" << std::endl;
                dist_source = dilNoise.makeSource(d,inoise);
                fermion4dtmp=0;
                if (Ls_ == 1)
		{
                    solver(fermion4dtmp, dist_source);
		}
		else
                {
                    mat.ImportPhysicalFermionSource(dist_source, v5dtmp);
                    solver(v5dtmp_sol, v5dtmp);
                    mat.ExportPhysicalFermionSolution(v5dtmp_sol, fermion4dtmp);
                }
                if(perambMode == pMode::outputSolve)
                {
                    auto &solveOut = envGet(std::vector<FermionField>, objName);
                    solveOut[inoise+nNoise*dIndexSolve] = fermion4dtmp;
                }
   	    }
            for (int is = 0; is < Ns; is++)
            {
                cv4dtmp = peekSpin(fermion4dtmp,is);
                for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
                {
                    ExtractSliceLocal(cv3dtmp,cv4dtmp,0,t-Ntfirst,Tdir); 
   	            for (int ivec = 0; ivec < par().nVec; ivec++)
                    {
                        ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
                        pokeSpin(perambulator.tensor(t, ivec, dk, inoise,idt,ds),static_cast<Complex>(innerProduct(evec3d, cv3dtmp)),is);
                    }
                }
            }
        }
    }
    // Now share my timeslice data with other members of the grid
    const int NumSlices{grid4d->_processors[Tdir] / grid3d->_processors[Tdir]};
    if (NumSlices > 1)
    {
        LOG(Debug) <<  "Sharing perambulator data with other nodes" << std::endl;
        const int MySlice {grid4d->_processor_coor[Tdir]};
        const int TensorSize {static_cast<int>(perambulator.tensor.size() * PerambTensor::Traits::count)};
        if (TensorSize != perambulator.tensor.size() * PerambTensor::Traits::count){
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
    }
    
    // Save the perambulator to disk from the boss node
    if (grid4d->IsBoss())
    {
        std::string sPerambName {par().perambFileName};
        sPerambName.append(".");
        sPerambName.append(std::to_string(vm().getTrajectory()));
        perambulator.write(sPerambName.c_str());
    }

    // Also save the unsmeared sink if specified
    // Alternative: Make this an A2AVector and use A2AVectorIO::read
    std::string sFileName(par().unsmearedSolveFileName);
    if(perambMode == pMode::outputSolve && !sFileName.empty())
    {
        sFileName.append(1, '.');
        const std::string sTrajNum{std::to_string(vm().getTrajectory())};
        sFileName.append(sTrajNum);
        auto &solveOut = envGet(std::vector<FermionField>, objName);
	// TODO: Which writer do we want here??
	/*#ifdef HAVE_HDF5
            using Default_Writer = Grid::Hdf5Writer;
        #else
            using Default_Writer = Grid::BinaryWriter;
        #endif
	Default_Writer w( sFileName );*/
	Grid::BinaryWriter w( sFileName );
        write( w, "unsmearedSolve", solveOut );
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Perambulator_hpp_
