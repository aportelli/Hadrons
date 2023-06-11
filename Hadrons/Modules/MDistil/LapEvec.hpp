/*
 * LapEvec.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Author Name <43034299+mmphys@users.noreply.github.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
 * Author: ferben <ferben@c180030.wlan.net.ed.ac.uk>
 * Author: ferben <ferben@c183011.wlan.net.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: ferben <ferben@localhost.localdomain>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/DilutedNoise.hpp>

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
                                    int, polyOrder,
                                    double, alpha,
                                    double, beta)
    ChebyshevParameters() = default;
    template <class ReaderClass> ChebyshevParameters(Reader<ReaderClass>& Reader){read(Reader,"Chebyshev",*this);}
};

struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
                                    int, nVec,
                                    int, nK,
                                    int, nP,
                                    int, maxIt,
                                    double, resid,
                                    int, irlLog)
    LanczosParameters() = default;
    template <class ReaderClass> LanczosParameters(Reader<ReaderClass>& Reader){read(Reader,"Lanczos",*this);}
};

// These are the actual parameters passed to the module during construction

struct LapEvecPar: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LapEvecPar
                                    ,std::string,         gauge
                                    ,ChebyshevParameters, cheby
                                    ,LanczosParameters,   lanczos
                                    ,std::string,         fileName)
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
    virtual void RotateEigen(std::vector<ColourVectorField> & evec);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LapEvec, TLapEvec<FIMPL>, MDistil);

/*************************************************************************************
 Rotate eigenvectors into our phase convention
 First component of first eigenvector is real and positive
 *************************************************************************************/

template <typename FImpl>
inline void TLapEvec<FImpl>::RotateEigen(std::vector<ColourVectorField> & evec)
{
    ColourVector cv0;
    auto grid = evec[0].Grid();
    Coordinate siteFirst(grid->Nd(),0);
    for( int k = 0 ; k < evec.size() ; k++ )
    {
        peekSite(cv0, evec[k], siteFirst);
        const std::complex<Real> cplx0{cv0()()(0).real(), cv0()()(0).imag()};
        const Real cplx0_mag{ std::abs(cplx0) };
        const std::complex<Real> std_phase{std::conj(cplx0/cplx0_mag)};
        LOG(Message) << "RotateEigen() : Vector " << k <<  " Site 0 : |" << cplx0 << "|=" << cplx0_mag
                 << " => phase=" << (std::arg(std_phase) / M_PI) << " pi" << std::endl;
        {
            const Grid::Complex phase{std_phase.real(),std_phase.imag()};
            evec[k] *= phase;
            // Get rid of the rounding error in imaginary phase on the very first site
            peekSite(cv0, evec[k], siteFirst);
            cv0()()(0).imag(0); // this should be zero after the phase multiply - force it to be so
            pokeSite(cv0, evec[k], siteFirst);
        }
    }
}

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
    envTmp(std::vector<typename DistillationNoise<FImpl>::LapPack>,  "eig", 1, Ntlocal);
    // Output objects
    envCreate(typename DistillationNoise<FImpl>::LapPack, getName(), 1, par().lanczos.nVec, gridHD);
}

/*************************************************************************************
 
 -Grad^2 (Peardon, 2009, pg 2, equation 3, https://arxiv.org/abs/0905.2160)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 *************************************************************************************/

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
    using LinearFunction<Field>::operator();
    void operator()(const Field& in, Field& out) {
        if (nd > in.Grid()->Nd())
        {
            HADRONS_ERROR(Range,"nd too large");
        }
        conformable( in, out );
        out = ( ( Real ) ( 2 * nd ) ) * in;
        Field tmp_(in.Grid());
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
    using LinearFunction<Field>::operator();
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
    const ChebyshevParameters &ChebPar{par().cheby};
    const LanczosParameters   &LPar{par().lanczos};
    
    // Disable IRL logging if requested
    LOG(Message) << "irlLog=" << LPar.irlLog << std::endl;
    const int PreviousIRLLogState{GridLogIRL.isActive()};
    GridLogIRL.Active( LPar.irlLog == 0 ? 0 : 1 );
    
    // Stout smeared gauge field
    envGetTmp(GaugeField, Umu_smear);
    Umu_smear = envGet(GaugeField, par().gauge); 
    
    ////////////////////////////////////////////////////////////////////////
    // Invert nabla operator separately on each time-slice
    ////////////////////////////////////////////////////////////////////////
    
    auto & eig4d = envGet(typename DistillationNoise<FImpl>::LapPack, getName() );
    envGetTmp(std::vector<typename DistillationNoise<FImpl>::LapPack>, eig);   // Eigenpack for each timeslice
    envGetTmp(GaugeField, UmuNoTime); // Gauge field without time dimension
    envGetTmp(ColourVectorField, src);
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    uint32_t ConvergenceErrors{0};
    const int NtFull{env().getDim(Tdir)};
    TimesliceEvals Evals{ NtFull, LPar.nVec };
    for (int t = 0; t < NtFull; t++)
    {
        for (int v = 0; v < LPar.nVec; v++)
        {
            Evals.tensor( t, v ) = 0;
        }
    }
    for (int t = 0; t < Ntlocal; t++ )
    {
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        LOG(Message) << " Compute eigenpack, local timeslice = " << t << " / " << Ntlocal << std::endl;
        LOG(Message) << " Lanczos residual = " << LPar.resid << std::endl;
        LOG(Message) << " Number of Lap eigenvectors (nvec) = " << LPar.nVec << std::endl;
        LOG(Message) << "------------------------------------------------------------" << std::endl;
        eig[t].resize(LPar.nK+LPar.nP,gridLD);
        
        // Construct smearing operator
        ExtractSliceLocal(UmuNoTime,Umu_smear,0,t,Tdir); // switch to 3d/4d objects
        Laplacian3D<ColourVectorField,GaugeField> Nabla(UmuNoTime);
        LOG(Message) << "Chebyshev preconditioning to order " << ChebPar.polyOrder
                     << " with parameters (alpha,beta) = (" << ChebPar.alpha << "," << ChebPar.beta << ")" << std::endl;
        Chebyshev<ColourVectorField> Cheb(ChebPar.alpha,ChebPar.beta,ChebPar.polyOrder);
        
        // Construct source vector according to Test_dwf_compressed_lanczos.cc
        src = 11.0; // NB: This is a dummy parameter and just needs to be non-zero
        RealD nn = norm2(src);
        nn = Grid::sqrt(nn);
        src = src * (1.0/nn);
        
        Laplacian3DHerm<ColourVectorField> NablaCheby(Cheb,Nabla);
        ImplicitlyRestartedLanczos<ColourVectorField>
        IRL(NablaCheby,Nabla,LPar.nVec,LPar.nK,LPar.nK+LPar.nP,LPar.resid,LPar.maxIt);
        int Nconv = 0;
        IRL.calc(eig[t].eval,eig[t].evec,src,Nconv);
        if (Nconv < LPar.nVec)
        {
            // NB: Can't assert here since we are processing local slices - i.e. not all nodes would assert
            ConvergenceErrors = 1;
            LOG(Error) << "MDistil::LapEvec : Not enough eigenvectors converged. If this occurs in practice, we should modify the eigensolver to iterate once more to ensure the second convergence test does not take us below the requested number of eigenvectors" << std::endl;
        }
        if( Nconv != LPar.nVec )
        {
            eig[t].resize(LPar.nVec, gridLD);
        }
        RotateEigen( eig[t].evec ); // Rotate the eigenvectors into our phase convention
        
        for (int i=0;i<LPar.nVec;i++)
        {
            InsertSliceLocal(eig[t].evec[i],eig4d.evec[i],0,t,Tdir);
            if(t==0 && Ntfirst==0)
            {
                eig4d.eval[i] = eig[t].eval[i]; 
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
    std::string sEigenPackName(par().fileName);
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
            Evals.MetaData.Version="0.1";
            Evals.write( sEigenPackName );
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_LapEvec_hpp_
