#ifndef Hadrons_MDistil_DistilMesonFieldPinned_hpp_
#define Hadrons_MDistil_DistilMesonFieldPinned_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/DistilMatrix.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/*******************************************************************
 *                         DistilMesonFieldPinned                        *
 * Receives LapH eigenvectors and perambulator/noise.              *
 * Computes DistilMesonFieldPinneds by block (of spin-lap dilution size) *
 * and save them to H5 files.                                      *
 *******************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

class DistilMesonFieldPinnedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPinnedPar,
                                    std::string,                outPath,
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
                                    std::string,                deltaT,
                                    std::string,                pinSide,
                                    std::string,                gamma,
                                    std::vector<std::string>,   momenta)
};

template <typename FImpl>
class TDistilMesonFieldPinned: public Module<DistilMesonFieldPinnedPar>
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
    TDistilMesonFieldPinned(const std::string name);
    // destructor
    virtual ~TDistilMesonFieldPinned(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string                         outputMFPath_;
    std::map<Side, unsigned int>        dilSizeLS_;
    std::vector<std::vector<RealF>>     momenta_;
    std::vector<Gamma::Algebra>         gamma_;  
    bool                                isExact_=false;
    std::map<Side,std::string>          dmfType_;
    std::vector<unsigned int>           tSourceL_;
    std::vector<unsigned int>           tSourceR_;
    Side                                pinned_side_;
    std::vector<unsigned int>           delta_t_list_;
};

MODULE_REGISTER_TMP(DistilMesonFieldPinned, TDistilMesonFieldPinned<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilMesonFieldPinned implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMesonFieldPinned<FImpl>::TDistilMesonFieldPinned(const std::string name)
: Module<DistilMesonFieldPinnedPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonFieldPinned<FImpl>::getInput(void)
{   
    std::vector<std::string> in = {par().lapEigenPack, par().leftNoise, par().rightNoise};

    //define meson field type (if a peramb object, set to phi, rho otherwise)
    dmfType_.emplace(Side::left   , par().leftPeramb.empty()  ? "rho" : "phi");
    dmfType_.emplace(Side::right  , par().rightPeramb.empty() ? "rho" : "phi");

    //require peramb dependency if phi case
    for(Side s : sides)
    {
        if(dmfType_.at(s)=="phi")
        {
            in.push_back( s==Side::left ? par().leftPeramb : par().rightPeramb);
        }
    }
    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilMesonFieldPinned<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonFieldPinned<FImpl>::setup(void)
{
    if( envHasDerivedType(DistillationNoise, ExactDistillationPolicy<FImpl>, par().leftNoise)
        and envHasDerivedType(DistillationNoise, ExactDistillationPolicy<FImpl>, par().rightNoise) )
    {
        isExact_ = true;
    }

    GridCartesian *g            = envGetGrid(FermionField);
    GridCartesian *g3d          = envGetSliceGrid(FermionField, g->Nd() - 1);
    DistillationNoise &noisel   = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser   = envGet( DistillationNoise , par().rightNoise);
    outputMFPath_   = par().outPath;
    dilSizeLS_      = { {Side::left,noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s)},
                        {Side::right,noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s)} };

    if(par().blockSize > dilSizeLS_.at(Side::left) or par().blockSize > dilSizeLS_.at(Side::right))
    {
         HADRONS_ERROR(Size, "blockSize needs to be <= Laplacian-spin space dimensions.");
    }

    if(par().noisePairs.empty() and !isExact_)
    {
        HADRONS_ERROR(Size, "Missing noise pairs input for stochastic distillation.");
    }

    // time source input validation
    MDistil::verifyTimeSourcesInput(par().leftTimeSources,noisel.dilutionSize(Index::t));
    MDistil::verifyTimeSourcesInput(par().rightTimeSources,noiser.dilutionSize(Index::t));
    tSourceL_ = strToVec<unsigned int>(par().leftTimeSources);
    tSourceR_ = strToVec<unsigned int>(par().rightTimeSources);

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

    //parse pinSide
    if(par().pinSide=="right")
    {
        pinned_side_ = Side::right;
    }
    else if(par().pinSide=="left")
    {
        pinned_side_ = Side::left;
    }
    else if( par().pinSide.empty() )    //default option
    {
        if( dmfType_.at(Side::left)=="rho" and dmfType_.at(Side::right)=="phi" )
        {
            pinned_side_ = Side::left;
        }
        else if( dmfType_.at(Side::left)=="phi" and dmfType_.at(Side::right)=="phi" )
        {
            pinned_side_ = Side::right;
        }
        else
        {
            pinned_side_ = Side::right;
        }
    }
    else
    {
        HADRONS_ERROR(Argument, "pinSide input is not valid (right or left).");
    }

    //parse deltaT
    if(dmfType_.at(pinned_side_)=="phi")
    {
        delta_t_list_ = strToVec<unsigned int>(par().deltaT);
    }
    else
    {
        delta_t_list_ = {0};    //pinning a rho field implies in delta_t=0 by definition
    }

    std::map<Side, unsigned int> dilSizeT = { {Side::left, pinned_side_==Side::left ? 1 : DISTILVECTOR_TIME_BATCH_SIZE},
                                                {Side::right, pinned_side_==Side::right ? 1 : DISTILVECTOR_TIME_BATCH_SIZE} };  //the pinned distilvector has always time-dilution dimension 1

    envTmpLat(ComplexField,             "coor");
    envTmp(std::vector<ComplexField>,   "phase",        1, momenta_.size(), g );
    envTmp(DistilVector,                "dvl",          1, dilSizeT.at(Side::left)*dilSizeLS_.at(Side::left), g);
    envTmp(DistilVector,                "dvr",          1, dilSizeT.at(Side::right)*dilSizeLS_.at(Side::right), g);
    envTmp(Computation,                 "computation",  1, dmfType_, g, g3d, noisel, noiser, par().blockSize, 
                par().cacheSize, env().getDim(g->Nd() - 1), momenta_.size(), gamma_.size(), isExact_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonFieldPinned<FImpl>::execute(void)
{
    // temps
    envGetTmp(DistilVector, dvl);
    envGetTmp(DistilVector, dvr);
    envGetTmp(Computation,  computation);
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

    // nvec check against noises (and assuming nvec cannot be different on different sides!)
    std::map<Side, DistillationNoise & > noises = {{Side::left,noisel},{Side::right,noiser}};
    if((noisel.getNl() != nVec) || (noiser.getNl() != nVec))
    {
        HADRONS_ERROR(Size, "Incompatibility between number of Laplacian eigenvectors and size of Laplacian subspace in the noises.");
    }

    // fetch time sources input
    std::map<Side, std::vector<unsigned int>> time_sources = {{Side::left,tSourceL_},{Side::right,tSourceR_}};
    std::map<Side, std::string> peramb_input = {{Side::left,par().leftPeramb},{Side::right,par().rightPeramb}}; // perambulator time sources
    std::map<Side, std::vector<int>> ts_peramb;
    for(Side s : sides)     
    {
        if(computation.isPhi(s))
        {
            auto & inPeramb = envGet(PerambTensor , peramb_input.at(s));
            ts_peramb.emplace(s , inPeramb.MetaData.timeSources);
            if(time_sources.at(s).empty())  //in case it's empty and it's a phi, include all available peramb time sources
            {
                for(auto tperamb : ts_peramb.at(s))
                {
                    time_sources.at(s).push_back(static_cast<unsigned int>(tperamb));
                }
            }
            else    // if it's not empty, validate it against peramb time sources (check if it is subset of that)
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
                std::iota( time_sources.at(s).begin() , time_sources.at(s).end() , 0);
            }
        }
        if(time_sources.at(s).size()%DISTILVECTOR_TIME_BATCH_SIZE != 0){ //distil vector batch size only supports being a divisor of number of time sources
            std::string errside = (s==Side::left) ? "left" : "right";
            HADRONS_ERROR(Range, "Number of time sources (" + errside + ") not divisible by distil vector batch size.");
        }
    }

    std::string filepath = par().outPath + "/" + dmfType_.at(Side::left) + "-" + dmfType_.at(Side::right) + "." + std::to_string(vm().getTrajectory()) + "/";
    Hadrons::mkdir(filepath);

    //auxiliar lambda expressions for names and metadata
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
        filename += ".h5";
        return filename;
    };

    auto metadataDmfFn = [this, &nt, &nVec, &noisel, &noiser](const unsigned int m, const unsigned int o, const unsigned int nl, const unsigned int nr, DilutionMap lmap, DilutionMap rmap)
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
        md.MesonFieldType   = dmfType_.at(Side::left) + "-" + dmfType_.at(Side::right);
        md.PinnedSide       = (pinned_side_==Side::left) ? "left" : "right";
        if(isExact_)
        {
            md.NoiseHashLeft   = noisel.generateHash()[nl];
            md.NoiseHashRight  = noiser.generateHash()[nr];
        }
        md.TimeDilutionLeft  = lmap[Index::t];
        md.TimeDilutionRight = rmap[Index::t];
        md.LapDilutionLeft   = lmap[Index::l];
        md.LapDilutionRight  = rmap[Index::l];
        md.SpinDilutionLeft  = lmap[Index::s];
        md.SpinDilutionRight = rmap[Index::s];
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
    
    if(isExact_)
    {
        LOG(Message) << "Exact distillation" << std::endl;
    }
    LOG(Message) << "Distil vector batch size (time-dilution direction) : " << DISTILVECTOR_TIME_BATCH_SIZE << std::endl;
    std::string pinned_side_str = (pinned_side_==Side::left) ? "left" : "right";
    LOG(Message) << "Pinned side : " << pinned_side_str << std::endl;
    LOG(Message) << "DeltaT list : " << delta_t_list_ << std::endl;
    LOG(Message) << "Selected time-dilution partitions :"   << std::endl;
    LOG(Message) << " Left : " << MDistil::timeslicesDump(time_sources.at(Side::left)) << std::endl;
    LOG(Message) << " Right : " << MDistil::timeslicesDump(time_sources.at(Side::right)) << std::endl;
    LOG(Message) << "Left/right Laplacian-spin dilution sizes : " 
        << dilSizeLS_.at(Side::left) << "/" << dilSizeLS_.at(Side::right) << std::endl;
    LOG(Message) << "Meson field type : " << dmfType_.at(Side::left) + "-" + dmfType_.at(Side::right) << std::endl;
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
        std::map<Side, unsigned int> noise_idx = {{Side::left,npair[0]},{Side::right,npair[1]}};
        LOG(Message) << "Noise pair : " << noise_idx.at(Side::left) << " " << noise_idx.at(Side::right) << std::endl;
        //computation of distillation vectors (phi or rho)
        if(computation.isPhi(Side::left) || computation.isPhi(Side::right)) //if theres at least one phi, populate peramb map
        {
            std::map<Side, PerambTensor&> peramb;
            for(Side s : sides)
            {
                if(computation.isPhi(s)){
                    PerambTensor &perambtemp = envGet( PerambTensor , s==Side::left ? par().leftPeramb : par().rightPeramb);
                    peramb.emplace(s , perambtemp);
                }
            }
            computation.executePinned(filenameDmfFn, metadataDmfFn, gamma_, dist_vecs, noise_idx, phase, time_sources, epack, this, pinned_side_, delta_t_list_, peramb);
        }
        else
        {
            computation.executePinned(filenameDmfFn, metadataDmfFn, gamma_, dist_vecs, noise_idx, phase, time_sources, epack, this, pinned_side_, delta_t_list_);
        }
        LOG(Message) << "Meson fields saved to " << outputMFPath_ << std::endl;
    }
    LOG(Message) << "A2AUtils::MesonField kernel executed " << computation.blockCounter_ << " times over "
        << nt  << "x" << par().cacheSize << "x" << par().cacheSize << " cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) : "    << computation.blockFlops_/computation.blockCounter_    << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) : "     << computation.blockBytes_/computation.blockCounter_    << " GB/s/node "    << std::endl;
    LOG(Message) << "Average IO speed (write) : "       << computation.blockIoSpeed_/computation.blockCounter_  << " MB/s "    << std::endl;
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonFieldPinned_hpp_
