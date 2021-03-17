/*
 * LapEvec.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
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

#ifndef Hadrons_MDistil_LapEvec_hpp_
#define Hadrons_MDistil_LapEvec_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 
 Laplacian eigenvectors - parameters

 Computes the eigenvectors of the 3D-Laplacian. The gauge field provided to 
 this module should be built from stout-smeared gauge links, where the smearing
 is only applied to the spatial components of the gauge field,
 i.e. rho_{4i} = rho_{i4} = rho_{44} = 0. 

 Chebyshev-preconditioning is needed for convergence of the nvec lowest 
 eigenvectors.
 
 ******************************************************************************/

struct ChebyshevParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ChebyshevParameters,
                                    int, PolyOrder,
                                    double, alpha,
                                    double, beta)
    ChebyshevParameters() = default;
    template <class ReaderClass> ChebyshevParameters(Reader<ReaderClass>& Reader){read(Reader,"Chebyshev",*this);}
};

struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
                                    int, Nvec,
                                    int, Nk,
                                    int, Np,
                                    int, MaxIt,
                                    double, resid,
                                    int, IRLLog)
    LanczosParameters() = default;
    template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

// These are the actual parameters passed to the module during construction

struct LapEvecPar: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar
                                    ,std::string,         gauge
                                    ,ChebyshevParameters, Cheby
                                    ,LanczosParameters,   Lanczos
                                    ,std::string,         FileName)
};

/******************************************************************************
 
 Laplacian eigenvectors - Module (class) definition
 
 ******************************************************************************/

template <typename FImpl>
class TLapEvec: public Module<LapEvecPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TLapEvec(const std::string name);
    // destructor
    virtual ~TLapEvec(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<FIMPL>, MDistil);

/******************************************************************************
 TLapEvec implementation
 ******************************************************************************/

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLapEvec<FImpl>::TLapEvec(const std::string name) : Module<LapEvecPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLapEvec<FImpl>::getInput(void)
{
    return std::vector<std::string>{par().gauge};
}

template <typename FImpl>
std::vector<std::string> TLapEvec<FImpl>::getOutput(void)
{
    return {getName()}; // This is the higher dimensional eigenpack
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLapEvec<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    // Temporaries
    envTmpLat(GaugeField, "Umu_smear");
    envTmp(LatticeGaugeField, "UmuNoTime", 1, gridLD);
    envTmp(ColourVectorField,  "src",1,gridLD);
    envTmp(std::vector<LapEvecs>,  "eig", 1, Ntlocal);
    // Output objects
    envCreate(LapEvecs, getName(), 1, par().Lanczos.Nvec, gridHD);
}

/*************************************************************************************
 
 -Grad^2 (Peardon, 2009, pg 2, equation 3, https://arxiv.org/abs/0905.2160)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 *************************************************************************************/

//template<typename FImpl> //would this be desired? 
template<typename Field, typename GaugeField>
class Laplacian3D : public LinearOperatorBase<Field>, public LinearFunction<Field> {
    typedef typename GaugeField::vector_type vCoeff_t;
public:
    int          nd; // number of spatial dimensions
    std::vector< LapEvec::ColourMatrixField > U;
    // Construct this operator given a gauge field and the number of dimensions it should act on
    Laplacian3D( GaugeField& gf, int dimSpatial = Tdir ) : nd{dimSpatial}
    {
        if (dimSpatial<1)
        {
            HADRONS_ERROR(Range,"Must be at least one spatial dimension");
        }
        for (int mu = 0 ; mu < nd ; mu++)
            U.push_back(PeekIndex<LorentzIndex>(gf,mu));
    }
    
    // Apply this operator to "in", return result in "out"
    void operator()(const Field& in, Field& out) {
        if (nd > in.Grid()->Nd())
        {
            HADRONS_ERROR(Range,"nd too large");
        }
        conformable( in, out );
        out = ( ( Real ) ( 2 * nd ) ) * in;
        Field tmp_(in.Grid());
        //typedef typename GaugeField::vector_type vCoeff_t;
        for (int mu = 0 ; mu < nd ; mu++)
        {
            out -= U[mu] * Cshift( in, mu, 1);
            tmp_ = adj( U[mu] ) * in;
            out -= Cshift(tmp_,mu,-1);
        }
    }
    
    void OpDiag (const Field &in, Field &out) { HADRONS_ERROR(Definition, "OpDiag() undefined"); };
    void OpDir  (const Field &in, Field &out,int dir,int disp) { HADRONS_ERROR(Definition, "OpDir() undefined"); };
    void OpDirAll  (const Field &in, std::vector<Field> &out){ HADRONS_ERROR(Definition, "OpDirAll() undefined"); };
    void Op     (const Field &in, Field &out) { HADRONS_ERROR(Definition, "Op() undefined"); };
    void AdjOp  (const Field &in, Field &out) { HADRONS_ERROR(Definition, "AdjOp() undefined"); };
    void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2) { HADRONS_ERROR(Definition, "HermOpAndNorm() undefined"); };
    void HermOp(const Field &in, Field &out) { operator()(in,out); };
};

template<typename Field>
class Laplacian3DHerm : public LinearFunction<Field> {
public:
    OperatorFunction<Field>   & poly_;
    LinearOperatorBase<Field> &Linop_;
    Laplacian3DHerm(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop)
    : poly_{poly}, Linop_{linop} {}
    void operator()(const Field& in, Field& out)
    {
        poly_(Linop_,in,out);
    }
};

/******************************************************************************
 Calculate low-mode eigenvalues of the Laplacian
 ******************************************************************************/

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLapEvec<FImpl>::execute(void)
{
    const ChebyshevParameters &ChebPar{par().Cheby};
    const LanczosParameters   &LPar{par().Lanczos};
    
    // Disable IRL logging if requested
    LOG(Message) << "IRLLog=" << LPar.IRLLog << std::endl;
    const int PreviousIRLLogState{GridLogIRL.isActive()};
    GridLogIRL.Active( LPar.IRLLog == 0 ? 0 : 1 );
    
    // Stout smeared gauge field
    envGetTmp(GaugeField, Umu_smear);
    Umu_smear = envGet(GaugeField, par().gauge); 
    
    ////////////////////////////////////////////////////////////////////////
    // Invert nabla operator separately on each time-slice
    ////////////////////////////////////////////////////////////////////////
    
    auto & eig4d = envGet(LapEvecs, getName() );
    envGetTmp(std::vector<LapEvecs>, eig);   // Eigenpack for each timeslice
    envGetTmp(GaugeField, UmuNoTime); // Gauge field without time dimension
    envGetTmp(ColourVectorField, src);
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    uint32_t ConvergenceErrors{0};
    const int NtFull{env().getDim(Tdir)};
    TimesliceEvals Evals{ NtFull, LPar.Nvec };
    for (int t = 0; t < NtFull; t++)
    {
        for (int v = 0; v < LPar.Nvec; v++)
	{
            Evals.tensor( t, v ) = 0;
	}
    }
    for (int t = 0; t < Ntlocal; t++ )
    {
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        LOG(Message) << " Compute eigenpack, local timeslice = " << t << " / " << Ntlocal << std::endl;
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        eig[t].resize(LPar.Nk+LPar.Np,gridLD);
        
        // Construct smearing operator
        ExtractSliceLocal(UmuNoTime,Umu_smear,0,t,Tdir); // switch to 3d/4d objects
        Laplacian3D<ColourVectorField,GaugeField> Nabla(UmuNoTime);
        LOG(Message) << "Chebyshev preconditioning to order " << ChebPar.PolyOrder
                     << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
        Chebyshev<ColourVectorField> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.PolyOrder);
        
        // Construct source vector according to Test_dwf_compressed_lanczos.cc
        src = 11.0; // NB: This is a dummy parameter and just needs to be non-zero
        RealD nn = norm2(src);
        nn = Grid::sqrt(nn);
        src = src * (1.0/nn);
        
        Laplacian3DHerm<ColourVectorField> NablaCheby(Cheb,Nabla);
        ImplicitlyRestartedLanczos<ColourVectorField>
        IRL(NablaCheby,Nabla,LPar.Nvec,LPar.Nk,LPar.Nk+LPar.Np,LPar.resid,LPar.MaxIt);
        int Nconv = 0;
        IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
        if (Nconv < LPar.Nvec)
        {
            // NB: Can't assert here since we are processing local slices - i.e. not all nodes would assert
            ConvergenceErrors = 1;
            LOG(Error) << "MDistil::LapEvec : Not enough eigenvectors converged. If this occurs in practice, we should modify the eigensolver to iterate once more to ensure the second convergence test does not take us below the requested number of eigenvectors" << std::endl;
        }
        if( Nconv != LPar.Nvec )
	{
            eig[t].resize(LPar.Nvec, gridLD);
	}
	RotateEigen( eig[t].evec ); // Rotate the eigenvectors into our phase convention
        
        for (int i=0;i<LPar.Nvec;i++)
	{
            InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,Tdir);
            if(t==0 && Ntfirst==0)
	    {
                eig4d.eval[i] = eig[t].eval[i]; // TODO: Discuss: is this needed? Is there a better way?
	    }
            if(gridLD->IsBoss()) // Only do this on one node per timeslice, so a global sum will work
	    {
                Evals.tensor(t + Ntfirst,i) = eig[t].eval[i];
	    }
        }
    }
    GridLogIRL.Active( PreviousIRLLogState );
    gridHD->GlobalSum(ConvergenceErrors);
    if(ConvergenceErrors!=0)
    {
        HADRONS_ERROR(Program,"The eingensolver failed to find enough eigenvectors on at least one node");
    }
    // Now write out the 4d eigenvectors
    std::string sEigenPackName(par().FileName);
    if( !sEigenPackName.empty() )
    {
        eig4d.record.solverXml = parString();
        ModuleBase * b{vm().getModule(par().gauge)};
        std::string sOperatorXml{ "<module><id><type>" };
        sOperatorXml.append( b->getRegisteredName() );
        sOperatorXml.append( "</type></id><options>" );
        sOperatorXml.append( b->parString() );
        sOperatorXml.append( "</options></module>" );
        eig4d.record.operatorXml = sOperatorXml;
        sEigenPackName.append(1, '.');
        std::size_t NameLen{ sEigenPackName.length() };
        const std::string sTrajNum{std::to_string(vm().getTrajectory())};
        sEigenPackName.append(sTrajNum);
        eig4d.write(sEigenPackName,false);
        // Communicate eig[t].evec to boss-node, save into new object evecs
        gridHD->GlobalSumVector(EigenIO::getFirstScalar(Evals.tensor),
                                static_cast<int>(EigenIO::getScalarCount(Evals.tensor)));
        if(gridHD->IsBoss())
        {
            sEigenPackName.resize(NameLen);
            sEigenPackName.append("evals.");
            sEigenPackName.append(sTrajNum);
            Evals.write( sEigenPackName );
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
