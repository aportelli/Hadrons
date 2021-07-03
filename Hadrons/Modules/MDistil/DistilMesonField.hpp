#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/DistilMatrix.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

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
    std::map<Side, unsigned int>        dilutionSize_ls_;
    std::vector<std::vector<RealF>>     momenta_;
    std::vector<Gamma::Algebra>         gamma_;  
    bool                                isExact_=false;
    std::map<Side,std::string>          dmf_type_;
    std::vector<unsigned int>           tsourcel_;
    std::vector<unsigned int>           tsourcer_;
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
    GridCartesian *g3d  = envGetSliceGrid(FermionField, g->Nd() - 1);
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    outputMFPath_       = par().outPath;
    dilutionSize_ls_    = { {Side::left,noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s)},
                            {Side::right,noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s)} };

    if(par().blockSize > dilutionSize_ls_.at(Side::left) or par().blockSize > dilutionSize_ls_.at(Side::right))
    {
         HADRONS_ERROR(Size, "blockSize needs to be <= Laplacian-spin space dimensions.");
    }

    if( ( vm().getModuleType(par().leftNoise) =="Grid::Hadrons::MNoise::ExactDistillation" ) &&
        ( vm().getModuleType(par().rightNoise)=="Grid::Hadrons::MNoise::ExactDistillation" ) )
    {
        isExact_=true;
    }

    if(par().noisePairs.empty() and !isExact_)
    {
        HADRONS_ERROR(Size, "Missing noise pairs input for stochastic distillation.");
    }

    // time source input validation
    MDistil::verifyTimeSourcesInput(par().leftTimeSources,noisel.dilutionSize(Index::t));
    MDistil::verifyTimeSourcesInput(par().rightTimeSources,noiser.dilutionSize(Index::t));
    tsourcel_ = strToVec<unsigned int>(par().leftTimeSources);
    tsourcer_ = strToVec<unsigned int>(par().rightTimeSources);

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
    envTmp(std::vector<ComplexField>,   "phase",        1, momenta_.size(), g );
    envTmp(DistilVector,                "dvl",          1, dilutionSize_ls_.at(Side::left), g);
    envTmp(DistilVector,                "dvr",          1, dilutionSize_ls_.at(Side::right), g);
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
    envGetTmp(std::vector<ComplexField>, phase);

    //start
    GridCartesian *g        = envGetGrid(FermionField);
    auto &epack             = envGet(typename DistillationNoise::LapPack, par().lapEigenPack);
    const unsigned int nVec = epack.evec.size();
    const unsigned int nd   = g->Nd();
    const unsigned int nt   = env().getDim(nd - 1);

    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int)>  MetadataFn;
    std::map<Side, DistilVector & > dist_vecs = {{Side::left,dvl}  ,{Side::right,dvr}};
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    std::vector<std::vector<unsigned int>>       noise_pairs;

    // nvec check (assume nvec cant be different on different sides)
    std::map<Side, DistillationNoise & > noises = {{Side::left,noisel},{Side::right,noiser}};
    if((noisel.getNl() != nVec) || (noiser.getNl() != nVec))
    {
        HADRONS_ERROR(Size, "Incompatibility between number of Laplacian eigenvectors and Laplacian subspace size in noises.");
    }

    // fetch time sources input
    std::map<Side, std::vector<unsigned int>> time_sources = {{Side::left,tsourcel_},{Side::right,tsourcer_}};
    std::map<Side, std::string> peramb_input = {{Side::left,par().leftPeramb},{Side::right,par().rightPeramb}}; // perambulators time sources
    std::map<Side, std::vector<int>> ts_peramb;
    for(Side s : sides)     
    {
        if(computation.isPhi(s))
        {
            auto & inPeramb = envGet(PerambTensor , peramb_input.at(s));
            ts_peramb.emplace(s , inPeramb.MetaData.timeSources);
            if(time_sources.at(s).empty())  //in case it's empty and it's a phi, include all peramb time sources
            {
                for(auto tperamb : ts_peramb.at(s))
                    time_sources.at(s).push_back(static_cast<unsigned int>(tperamb));
            }
            else    // if it's not empty, validate it against peamb time sources (check if is subset of it)
            {
                if( !std::includes(ts_peramb.at(s).begin(), ts_peramb.at(s).end(),
                                time_sources.at(s).begin(), time_sources.at(s).end()) )
                {
                    std::string errside = (s==Side::left) ? "left" : "right";
                    HADRONS_ERROR(Argument,"Time sources are not available on " + errside + " perambulator");
                }
            }
        }
        else
        {
            if(time_sources.at(s).empty())   //in case it's empty and it's a rho, include all time sources
            {
                time_sources.at(s).resize(noises.at(s).dilutionSize(Index::t));
                std::iota( time_sources.at(s).begin() , time_sources.at(s).end() , 0);    //creates sequence from 0 to TI-1
            }
        }
    }

    std::string filepath = par().outPath + "/" + par().mesonFieldType + "/";
    Hadrons::mkdir(filepath);

    //auxiliar lambda functions for names and metadata
    auto filenameDmfFn = [this, filepath](const unsigned int m, const unsigned int o, const unsigned int nl, const unsigned int nr)
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

    auto metadataDmfFn = [this, &nt, &nVec, &noisel, &noiser](const unsigned int m, const unsigned int o, const unsigned int nl, const unsigned int nr)
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
            md.NoiseHashLeft   = noisel.generateHash();
            md.NoiseHashRight  = noiser.generateHash();
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
            noise_pairs.push_back(strToVec<unsigned int>(npair));
        }
    }

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
    stopTimer("momentum phases");
    
    LOG(Message) << "Selected time-dilution partitions :"         << std::endl;
    LOG(Message) << " Left : " << MDistil::timeslicesDump(time_sources.at(Side::left)) << std::endl;
    LOG(Message) << " Right : " << MDistil::timeslicesDump(time_sources.at(Side::right)) << std::endl;
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
            computation.execute(filenameDmfFn, metadataDmfFn, gamma_, dist_vecs, npair, phase, time_sources, epack, this, peramb);
        }
        else
        {
            computation.execute(filenameDmfFn, metadataDmfFn, gamma_, dist_vecs, npair, phase, time_sources, epack, this);
        }
        LOG(Message) << "Meson fields saved to " << outputMFPath_ << std::endl;
    }
    LOG(Message) << "A2AUtils::MesonField kernel executed "   << computation.global_counter << " times over " << 
                    DISTIL_NT_CHUNK_SIZE  << "x" << par().cacheSize << "x" << par().cacheSize << " cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) : "          << computation.global_flops/computation.global_counter    << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) : "           << computation.global_bytes/computation.global_counter    << " GB/s/node "    << std::endl;
    LOG(Message) << "Average IO speed (write) : "             << computation.global_iospeed/computation.global_counter  << " MB/s "    << std::endl;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
