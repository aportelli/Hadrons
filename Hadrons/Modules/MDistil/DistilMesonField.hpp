#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Modules/MDistil/DistilMatrix.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>
#include <Hadrons/Modules/MNoise/ExactDistillation.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilMesonField                                 *
 * Eliminates DistilVectors module. Receives LapH eigenvectors and 
 * perambulator/noise (as left/right objs). Computes MesonFields by 
 * block (and chunking it) and save them to H5 files.
 * 
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

class DistilMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPar,
                                    std::string,                outPath,
                                    std::string,                mesonFieldType,
                                    std::string,                lapEigenPack,
                                    std::string,                leftNoise,
                                    std::string,                rightNoise,
                                    std::vector<std::string>,   noisePairs,
                                    std::string,                leftTimeSources,
                                    std::string,                rightTimeSources,
                                    std::string,                leftPeramb,
                                    std::string,                rightPeramb,
                                    unsigned int,               blockSize,
                                    unsigned int,               cacheSize,
                                    std::string,                gamma,
                                    std::vector<std::string>,   momenta)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DmfComputation<FImpl, Complex, HADRONS_DISTIL_IO_TYPE>    Computation;
public:
    typedef typename Computation::Index Index;
    typedef typename Computation::DistilVector DistilVector;
    typedef typename Computation::DistillationNoise DistillationNoise;
public:
    // constructor
    TDistilMesonField(const std::string name);
    // destructor
    virtual ~TDistilMesonField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string                         outputMFPath_;
    bool                                hasPhase_{false};
    std::map<Side, unsigned int>        dilutionSize_ls_;
    std::vector<std::vector<RealF>>     momenta_;
    std::vector<Gamma::Algebra>         gamma_;  
    bool                                isExact_=false;
    std::map<Side,std::string>          dmf_type_;
};

MODULE_REGISTER_TMP(DistilMesonField, TDistilMesonField<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilMesonField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMesonField<FImpl>::TDistilMesonField(const std::string name)
: Module<DistilMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getInput(void)
{   
    std::vector<std::string> in = {par().lapEigenPack, par().leftNoise, par().rightNoise};
    std::string c = par().mesonFieldType;
    // check mesonfield type
    if(!(c=="phi-phi" || c=="phi-rho" || c=="rho-phi" || c=="rho-rho"))
    {
        HADRONS_ERROR(Argument,"Bad meson field type");
    }

    dmf_type_.emplace(Side::left  , c.substr(0,3));
    dmf_type_.emplace(Side::right , c.substr(4,7));

    //require peramb dependency if phi case
    for(Side s : sides)
    {
        if(dmf_type_.at(s)=="phi")
        {
            in.push_back( s==Side::left ? par().leftPeramb : par().rightPeramb);
        }
    }

    if( ( vm().getModuleType(par().leftNoise) =="Grid::Hadrons::MNoise::ExactDistillation" ) &&
        ( vm().getModuleType(par().rightNoise)=="Grid::Hadrons::MNoise::ExactDistillation" ) )
    {
        isExact_=true;
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::setup(void)
{
    GridCartesian *g    = envGetGrid(FermionField);
    GridCartesian *g3d  = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    outputMFPath_       = par().outPath;
    dilutionSize_ls_    = { {Side::left,noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s)},
                            {Side::right,noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s)} };
    
    // time source input validation
    MDistil::verifyTimeSourcesInput(par().leftTimeSources,noisel.dilutionSize(Index::t));
    MDistil::verifyTimeSourcesInput(par().rightTimeSources,noiser.dilutionSize(Index::t));

    // parse momenta
    momenta_.clear();
    for(auto &p_string : par().momenta)
    {
        auto p = strToVec<RealF>(p_string);
        if (p.size() != g->Nd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(g->Nd() - 1));
        }
        momenta_.push_back(p);
    }

    //parse gamma
    gamma_.clear();
    if (par().gamma == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gamma);
    }

    envTmpLat(ComplexField,             "coor");
    envCache(std::vector<ComplexField>, "phasename",    1, momenta_.size(), g );
    envTmp(DistilVector,                "dvl",          1, noisel.dilutionSize() , g);
    envTmp(DistilVector,                "dvr",          1, noiser.dilutionSize() , g);
    envTmp(Computation,                 "computation",  1, dmf_type_, g, g3d, noisel, noiser, par().blockSize, 
                par().cacheSize, env().getDim(g->Nd() - 1), momenta_.size(), gamma_.size(), isExact_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    // temps
    envGetTmp(DistilVector,     dvl);
    envGetTmp(DistilVector,     dvr);
    envGetTmp(Computation,      computation);

    //start
    GridCartesian *g        = envGetGrid(FermionField);
    auto &epack             = envGet(typename DistillationNoise::LapPack, par().lapEigenPack);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");
    const unsigned int nVec = epack.evec.size();
    const unsigned int nd   = g->Nd();
    const unsigned int nt   = env().getDim(nd - 1);

    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int)>  MetadataFn;
    std::map<Side, DistilVector & > dist_vecs = {{Side::left,dvl}  ,{Side::right,dvr}};
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    std::vector<std::vector<int>>       noise_pairs;

    // nvec check (assume nvec cant be different on different sides)
    std::map<Side, DistillationNoise & > noises = {{Side::left,noisel},{Side::right,noiser}};
    if((noisel.getNl() != nVec) || (noiser.getNl() != nVec))
    {
        HADRONS_ERROR(Size, "Incompatibility between number of Laplacian eigenvectors and Laplacian subspace size in noises.");
    }

    // fetch time sources 
    std::vector<unsigned int> tsourcel = strToVec<unsigned int>(par().leftTimeSources);
    std::vector<unsigned int> tsourcer = strToVec<unsigned int>(par().rightTimeSources);
    std::map<Side, std::vector<unsigned int>> timeDilSource = {{Side::left,tsourcel},{Side::right,tsourcer}};
    for(Side s : sides)     //in case it's empty, include all possible time sources
    {
        if(timeDilSource.at(s).empty()){
            timeDilSource.at(s).resize(noises.at(s).dilutionSize(Index::t));
            std::iota( timeDilSource.at(s).begin() , timeDilSource.at(s).end() , 0);    //creates sequence from 0 to TI-1
        }
    }

    // parse timeSource from perambulators and validate input against that
    std::map<Side, std::string> peramb_input = {{Side::left,par().leftPeramb},{Side::right,par().rightPeramb}};
    std::map<Side, std::vector<int>> ts_peramb;
    for(Side s : sides)
    {
        if(computation.isPhi(s))
        {
            auto & inPeramb = envGet(PerambTensor , peramb_input.at(s));
            ts_peramb.emplace(s , inPeramb.MetaData.timeSources);
            if( !std::includes(ts_peramb.at(s).begin(), ts_peramb.at(s).end(),
                             timeDilSource.at(s).begin(), timeDilSource.at(s).end()) )  //check if input time source is compatible with peramb's (subset of it)
            {
                std::string errside = (s==Side::left) ? "left" : "right";
                HADRONS_ERROR(Argument,"Time sources are not available on " + errside + " perambulator");
            }
        }
    }

    std::string filepath = par().outPath + "/" + par().mesonFieldType + "/";
    Hadrons::mkdir(filepath);

    //auxiliar lambda functions for names and metadata
    auto filenameDmfFn = [this, filepath](const unsigned int m, const unsigned int o, const int nl, const int nr)
    {
        std::stringstream ss;
        ss << gamma_[o] << "_p";
        for (unsigned int mu = 0; mu < momenta_[m].size(); ++mu)
            ss << momenta_[m][mu] << ((mu == momenta_[m].size() - 1) ? "" : "_");

        std::string filename = filepath + ss.str(); 
        if(!isExact_)
        {
            filename += "_n" + std::to_string(nl) + "_" + std::to_string(nr) ;
        }
        filename += "." + std::to_string(vm().getTrajectory()) + ".h5";
        return filename;
    };

    auto metadataDmfFn = [this, &nt, &nVec, &noisel, &noiser](const unsigned int m, const unsigned int o, const int nl, const int nr)
    {
        DistilMesonFieldMetadata<FImpl> md;
        for (auto pmu: momenta_[m])
        {
            md.Momentum.push_back(pmu);
        }
        md.Operator         = gamma_[o];
        md.Nt               = nt;   
        md.Nvec             = nVec;     //nvec is the same for both sides
        md.NoisePair        = {nl,nr};
        md.MesonFieldType   = par().mesonFieldType;
        if(isExact_)
        {
            md.NoiseHashesLeft     = {"0"}; // exact distil convention
            md.NoiseHashesRight    = {"0"};
        }
        else
        {
            md.NoiseHashesLeft   = noisel.generateHash();
            md.NoiseHashesRight  = noiser.generateHash();
        }
        return md;
    };

    // prepare noise pairs for execution
    if(isExact_)
    {
        noise_pairs.push_back({0,0});
    }
    else
    {
        for(auto &npair : par().noisePairs)
        {
            noise_pairs.push_back(strToVec<int>(npair));
        }
    }

    //compute momentum phases
    if (!hasPhase_)
    {
        startTimer("momentum phases");
        envGetTmp(ComplexField, coor);
        Complex           i(0.0,1.0);
        for (unsigned int j = 0; j < momenta_.size(); ++j)
        {
            phase[j] = Zero();
            for(unsigned int mu = 0; mu < momenta_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                phase[j] = phase[j] + (momenta_[j][mu]/env().getDim()[mu])*coor;
            }
            phase[j] = exp((Real)(2*M_PI)*i*phase[j]);
        }
        hasPhase_ = true;
        stopTimer("momentum phases");
    }
    
    LOG(Message) << "Selected time-dilution partitions :"         << std::endl;
    LOG(Message) << " Left : " << timeslicesFn(timeDilSource.at(Side::left)) << std::endl;
    LOG(Message) << " Right : " << timeslicesFn(timeDilSource.at(Side::right)) << std::endl;
    LOG(Message) << "Left/right Laplacian-spin dilution sizes : " 
        << dilutionSize_ls_.at(Side::left) << "/" << dilutionSize_ls_.at(Side::right) << std::endl;
    LOG(Message) << "Meson field type : " << par().mesonFieldType << std::endl;
    LOG(Message) << "Momenta :" << std::endl;
    for (auto &p: momenta_)
    {
        LOG(Message) << " " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears :" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << " " << g << std::endl;
    }
    LOG(Message) << "Block size : " << par().blockSize << std::endl;
    LOG(Message) << "Cache block size : " << par().cacheSize << std::endl;

    //execution
    for(auto &npair : noise_pairs)
    {
        LOG(Message) << "Noise pair : " << npair[0] << " " << npair[1] << std::endl;
        //computation of distillation vectors (phi or rho)
        if(computation.isPhi(Side::left) || computation.isPhi(Side::right)) //if theres at least one phi, populate peramb
        {
            std::map<Side, PerambTensor&> peramb;
            for(Side s : sides)
            {
                if(computation.isPhi(s)){
                    PerambTensor &perambtemp = envGet( PerambTensor , s==Side::left ? par().leftPeramb : par().rightPeramb);
                    peramb.emplace(s , perambtemp);
                }
            }
            computation.makeDistVecs(dist_vecs, npair, epack, timeDilSource, peramb);
        }
        else
        {
            computation.makeDistVecs(dist_vecs, npair, epack, timeDilSource);
        }

        // computing mesonfield blocks and saving to disk
        computation.execute(filenameDmfFn, metadataDmfFn, gamma_, dist_vecs, npair, phase, timeDilSource, this);

        LOG(Message) << "Meson fields saved at " << outputMFPath_ << std::endl;
    }
    LOG(Message) << "A2AUtils::MesonField kernel executed " << computation.global_counter << " times over " << par().cacheSize << "^2 cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) : "          << computation.global_flops/computation.global_counter        << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) : "           << computation.global_bytes/computation.global_counter        << " GB/s/node "    << std::endl;
    LOG(Message) << "Average IO speed (write) : "             << computation.global_iospeed/computation.global_counter      << " MB/s "    << std::endl;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
