#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Modules/MDistil/DistilMatrix.hpp>

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
                                    std::string,                lapEvec,
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
    typedef DmfHelper<FImpl>         Helper;
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
    unsigned int blockSize_;
    unsigned int cacheSize_;
    std::map<std::string,std::string>   dmf_type_;
    std::vector<std::vector<int>>       noisePairs_;           // read from extermal object (diluted noise class)
    std::string                         outputMFStem_;
    bool                                hasPhase_{false};
    std::map<std::string, unsigned int> dilutionSize_ls_;
    std::map<std::string, std::string>  perambInput_ ;
    std::vector<std::string>            sides_       ;
    std::vector<std::vector<RealF>>     momenta_;
    std::vector<Gamma::Algebra>         gamma_;
    
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
    sides_           = {"left","right"};
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getInput(void)
{   
    std::vector<std::string> in = {par().lapEvec, par().leftNoise, par().rightNoise};

    std::string c = par().mesonFieldType;
    // check mesonfield type
    if(!(c=="phi-phi" || c=="phi-rho" || c=="rho-phi" || c=="rho-rho"))
    {
        HADRONS_ERROR(Argument,"Bad meson field type");
    }

    dmf_type_.emplace("left"  , c.substr(0,3));
    dmf_type_.emplace("right" , c.substr(4,7));

    // should I do a softer check?
    for(auto s : sides_)
    {
        if(dmf_type_.at(s)=="phi")
        {
            in.push_back( s=="left" ? par().leftPeramb : par().rightPeramb);
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
    outputMFStem_       = par().outPath;
    GridCartesian *g    = envGetGrid(FermionField);
    GridCartesian *g3d  = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)

    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    perambInput_    = {{"left",par().leftPeramb},{"right",par().rightPeramb}};

    blockSize_ = par().blockSize;
    cacheSize_ = par().cacheSize;

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

    dilutionSize_ls_ = { {"left",noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s)} , {"right",noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s)} };

    envTmpLat(ComplexField,             "coor");
    envCache(std::vector<ComplexField>, "phasename",    1, momenta_.size(), g );
    envTmp(DistilVector,                "dvl",          1, noisel.dilutionSize() , g);
    envTmp(DistilVector,                "dvr",          1, noiser.dilutionSize() , g);
    envTmp(Computation,                 "computation",  1, dmf_type_, momenta_, gamma_, g, g3d, blockSize_ , cacheSize_, env().getDim(g->Nd() - 1));
    envTmp(Helper,                      "helper"   ,    1, noisel, noiser , dmf_type_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    GridCartesian *g        = envGetGrid(FermionField);
    auto &epack             = envGet(typename DistillationNoise::LapPack, par().lapEvec);
    const unsigned int nVec = epack.evec.size();
    const unsigned int nd   = g->Nd();
    const unsigned int nt   = env().getDim(nd - 1);

    // temps
    envGetTmp(DistilVector,     dvl);
    envGetTmp(DistilVector,     dvr);
    envGetTmp(Computation,      computation);
    envGetTmp(Helper,           helper);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");

    // do not use operator []!! similar but better way to do that?
    std::map<std::string, DistilVector & > distVectors = {{"left",dvl}  ,{"right",dvr}};
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    std::map<std::string, DistillationNoise & >   noises = {{"left",noisel},{"right",noiser}};
    
    if((noisel.getNl() != nVec) || (noiser.getNl() != nVec))
    {
        HADRONS_ERROR(Size, "Incompatibility between number of Laplacian eigenvectors and Laplacian subspace size in noises.");
    }

    //encapsulate this in helper or computation class
    std::vector<unsigned int> lSources = strToVec<unsigned int>(par().leftTimeSources);
    std::vector<unsigned int> rSources = strToVec<unsigned int>(par().rightTimeSources);

    std::map<std::string, std::vector<unsigned int>> timeDilSource = {{"left",lSources},{"right",rSources}};
    //in case it's empty, include all possible time sources
    for(auto s : sides_)
    {
        if(timeDilSource.at(s).empty()){
            timeDilSource.at(s).resize(noises.at(s).dilutionSize(Index::t));
            std::iota( timeDilSource.at(s).begin() , timeDilSource.at(s).end() , 0);    //create sequence from 0 to TI-1
        }
        std::sort(timeDilSource.at(s).begin(), timeDilSource.at(s).end());  //guarantee they are ordered

        if( timeDilSource.at(s).size() > noises.at(s).dilutionSize(Index::t) )
        {
            HADRONS_ERROR(Argument,"Invalid number of time sources.");
        }
        else if( !std::none_of( timeDilSource.at(s).cbegin() , timeDilSource.at(s).cend(), 
                                [s,&noises](int dt){ return dt >= noises.at(s).dilutionSize(Index::t); }) ) //checks if any element is larger than time dilution size
        {
            HADRONS_ERROR(Argument,"Invalid value for one or more time sources.");
        }
    }

    // parse timeSource from Perambs
    std::map<std::string, std::vector<int>> peramb_st;
    for(auto & s : sides_)
    {
        if(dmf_type_.at(s)=="phi")
        {
            auto & inPeramb = envGet(PerambTensor , perambInput_.at(s));
            peramb_st.emplace(s , inPeramb.MetaData.timeSources);
            std::sort( peramb_st.at(s).begin() , peramb_st.at(s).end() );   //guarantee they are ordered
            if( !std::includes(peramb_st.at(s).begin(), peramb_st.at(s).end(),
                             timeDilSource.at(s).begin(), timeDilSource.at(s).end()) )  //check if input time source is compatible with peramb's (subset of it)
            {
                HADRONS_ERROR(Argument,"Time sources are not available on " +s+ " perambulator");
            }
        }
    }

    noisePairs_ = helper.parseNoisePairs(par().noisePairs);

    //compute momentum phase
    if (!hasPhase_)
    {
        startTimer("momentum phases");
        envGetTmp(ComplexField, coor);
        helper.computePhase(momenta_, coor, env().getDim(), phase);
        hasPhase_ = true;
        stopTimer("momentum phases");
    }
    
    LOG(Message) << "Selected time-dilution partitions:"         << std::endl;
    LOG(Message) << "Left:" << timeDilSource.at("left") << std::endl;
    LOG(Message) << "Right:" << timeDilSource.at("right") << std::endl;
    LOG(Message) << "Meson field type: "    << par().mesonFieldType << std::endl;
    LOG(Message) << "Selected block size: " << par().blockSize << std::endl;
    LOG(Message) << "Selected cache size: " << par().cacheSize << std::endl;
    LOG(Message) << "Lap-spin dilution size (left x right): " << dilutionSize_ls_.at("left") << " x " << dilutionSize_ls_.at("right") << std::endl;

    for(auto &inoise : noisePairs_)
    {
        LOG(Message) << "Noise pair: " << inoise << std::endl;
        LOG(Message) << "Gamma:" << gamma_ << std::endl;
        LOG(Message) << "Momenta:" << momenta_ << std::endl;

        //computation of distvectors
        if(dmf_type_.at("left")=="phi" || dmf_type_.at("right")=="phi")
        {
            std::map<std::string, PerambTensor&> peramb; // = { {"left", perambl} , {"right", perambr} };
            for(auto s : sides_)
            {
                if(dmf_type_.at(s)=="phi"){
                    PerambTensor &perambtemp = envGet( PerambTensor , s=="left" ? par().leftPeramb : par().rightPeramb);
                    peramb.emplace(s , perambtemp);
                }
            }
            computation.distVec(distVectors, noises, inoise, epack, timeDilSource, peramb);
        }
        else
        {
            computation.distVec(distVectors, noises, inoise, epack, timeDilSource);
        }

        // computing mesonfield blocks and saving to disk
        LOG(Message) << "Time-dilution blocks computation starting..." << std::endl;
        computation.execute(outputMFStem_, distVectors, noises, inoise, phase, dilutionSize_ls_, timeDilSource, helper.timeSliceMap(noises.at("left")), helper.timeSliceMap(noises.at("right")), this);

        LOG(Message) << "Meson fields saved at " << outputMFStem_ << std::endl;
    }
    LOG(Message) << "A2AUtils::MesonField kernel executed " << computation.global_counter << " times over " << cacheSize_ << "^2 cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) "          << computation.global_flops/computation.global_counter        << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) "           << computation.global_bytes/computation.global_counter        << " GB/s/node "    << std::endl;
    LOG(Message) << "Average IO speed (write) "             << computation.global_iospeed/computation.global_counter      << " MB/s "    << std::endl;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
