#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include "DmfTemp.hpp"

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilMesonField                                 *
 * Eliminates DistilVectors module. Receives LapH eigenvectors and 
 * perambulator/noise (as left/right fields). Computes MesonFields by 
 * block(and chunking it) and save them to H5 file.
 * 
 * For now, do not load anything from disk. Trying phi phi case 
 * with full-dilution.
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

class DistilMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPar,
                                    std::string,                outPath,
                                    std::string,                mesonFieldCase,
                                    std::string,                lapEvec,
                                    std::string,                leftPeramb,
                                    std::string,                rightPeramb,
                                    std::string,                leftNoise,
                                    std::string,                rightNoise,
                                    std::string,                gamma,
                                    int,                        blockSize,
                                    int,                        cacheSize,
                                    std::string,                leftTimeDilSources,
                                    std::string,                rightTimeDilSources,
                                    std::vector<std::string>,   momenta,
                                    std::vector<std::string>,   noisePairs)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DmfComputation<FImpl, FermionField, Complex, ComplexF>    Computation;
    typedef DmfHelper<FImpl, FermionField>         Helper;
public:
    typedef typename Computation::Index Index;
    typedef typename Computation::DistilVector DistilVector;
    typedef typename Computation::DistillationNoise DistillationNoise;
public:
    // constructor
    TDistilMesonField(const std::string name);
    // destructorœ
    virtual ~TDistilMesonField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    int blockSize_;
    int cacheSize_;
    std::map<std::string,std::string>   dmf_case_;
    // std::string                        momphName_;
    std::vector<Gamma::Algebra>         gamma_;
    std::vector<std::vector<RealF>>     momenta_;
    int                                 nExt_;
    int                                 nStr_;
    int                                 eff_nt_;
    std::vector<std::vector<int>>       noisePairs_;           // read from extermal object (diluted noise class)
    TimeSliceMap                        st_;
    std::string                         outputMFStem_;
    bool                                hasPhase_{false};
    std::map<std::string, int>          dilutionSize_ls_;
    std::map<std::string, std::string>  noiseInput_  ;
    std::map<std::string, std::string>  perambInput_ ;
    std::vector<std::string>            sides_       ;
    // DmfHelper<FImpl>                   *helper_;
    
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
    return{par().lapEvec, par().leftPeramb, par().rightPeramb, par().leftNoise, par().rightNoise};
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
    std::string c = par().mesonFieldCase;
    // check mesonfield case
    if(!(c=="phi phi" || c=="phi rho" || c=="rho phi" || c=="rho rho"))
    {
        HADRONS_ERROR(Argument,"Bad meson field case");
    }

    dmf_case_.emplace("left"  , c.substr(0,3));
    dmf_case_.emplace("right" , c.substr(4,7));

    outputMFStem_       = par().outPath;
    GridCartesian *g    = envGetGrid(FermionField);
    GridCartesian *g3d  = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)

    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    noiseInput_     = {{"left",par().leftNoise},{"right",par().rightNoise}}; //apparently not used
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
    nExt_       = momenta_.size(); //noise pairs computed independently, but can optmize embedding it into nExt??

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
    nStr_       = gamma_.size();

    //not taking into account different spin/lap dilution on each side, just different time dilutions
    // if( noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s) != noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s) )
    // {
    //     HADRONS_ERROR(Argument,"Spin-Lap dilution spaces do not have same dimension on both sides.");
    // }
    // dilutionSize_ls_ = noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s);

    dilutionSize_ls_ = { {"left",noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s)} , {"right",noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s)} };

    envTmpLat(ComplexField,             "coor");
    envCache(std::vector<ComplexField>, "phasename",    1, momenta_.size(), g );
    envTmp(DistilVector,                "dvl",          1, noisel.dilutionSize() , g);
    envTmp(DistilVector,                "dvr",          1, noiser.dilutionSize() , g);
    envTmp(Computation,                 "computation",  1, dmf_case_ , g, g3d, blockSize_ , cacheSize_, env().getDim(g->Nd() - 1));
    envTmp(Helper,                      "helper"   ,    1, noisel, noiser , dmf_case_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    GridCartesian *g        = envGetGrid(FermionField);

    auto &epack             = envGet(LapEvecs, par().lapEvec);
    int nVec                = epack.evec.size();
    const unsigned int nd   = g->Nd();
    const int nt            = env().getDim(nd - 1);

    // temps
    envGetTmp(DistilVector,     dvl);
    envGetTmp(DistilVector,     dvr);
    envGetTmp(Computation,      computation);
    envGetTmp(Helper,           helper);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");

    // do not use operator []!! similar but better way to do that?
    std::map<std::string, DistilVector & > distVector = {{"left",dvl}  ,{"right",dvr}};
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    std::map<std::string, DistillationNoise & >   noise = {{"left",noisel},{"right",noiser}};
    
    //encapsulate this in helper or computation class
    std::vector<int> lSources = strToVec<int>(par().leftTimeDilSources);
    std::vector<int> rSources = strToVec<int>(par().rightTimeDilSources);
    if(lSources.empty())
    {
        lSources.resize(noisel.dilutionSize(Index::t));
        std::iota(std::begin(lSources), std::end(lSources), 0);
    }
    if(rSources.empty())
    {
        rSources.resize(noiser.dilutionSize(Index::t));
        std::iota(std::begin(rSources), std::end(rSources), 0);
    }
    std::map<std::string, std::vector<int>> timeDilSource = {{"left",lSources},{"right",rSources}};

    // todo: do the following only in the necessary cases
    PerambTensor &perambl = envGet( PerambTensor , par().leftPeramb);
    PerambTensor &perambr = envGet( PerambTensor , par().rightPeramb);
    std::map<std::string, PerambTensor&> peramb = {{"left",perambl},{"right",perambr}};
    // parse source times
    std::map<std::string, std::vector<int>> peramb_st = {{"left",{}},{"right",{}}}; //source times from the peramb objs
    for(auto & s : sides_)
    {
        if(dmf_case_.at(s)=="phi")
        {
            auto & inPeramb = envGet(PerambTensor , perambInput_.at(s));
            peramb_st.at(s) = inPeramb.MetaData.timeSources;
            //TODO: check timeDilSource are subsets of this
        }
    }

    eff_nt_     = helper.computeEffTimeDimension();
    noisePairs_ = helper.parseNoisePairs(par().noisePairs);

    // ObjArray_LR<DistillationNoise*> n;
    // n = {&noisel, &noiser};
    // std::cout << "dilution size= " << n[0]->dilutionSize() << std::endl;

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
    LOG(Message) << " Left:" << timeDilSource.at("left") << std::endl;
    LOG(Message) << " Right:" << timeDilSource.at("right") << std::endl;
    LOG(Message) << "Meson field case: "    << par().mesonFieldCase << std::endl;
    LOG(Message) << "EffTime dimension = "     << eff_nt_ << std::endl;
    LOG(Message) << "Selected block size: " << par().blockSize << std::endl;
    LOG(Message) << "Selected cache size: " << par().cacheSize << std::endl;
    LOG(Message) << "Lap-spin dilution size (left x right): " << dilutionSize_ls_.at("left") << " x " << dilutionSize_ls_.at("right") << std::endl;

    for(auto &inoise : noisePairs_)
    {
        // set up io object and metadata for all gamma/momenta -> turn into method
        std::vector<A2AMatrixIo<ComplexF>> matrixIoTable;
        DistilMesonFieldMetadata<FImpl> md;
        for(int iExt=0; iExt<nExt_; iExt++)
        for(int iStr=0; iStr<nStr_; iStr++)
        {
            // metadata;
            md.momentum             = momenta_[iExt];
            md.gamma                = gamma_[iStr];
            md.noise_pair           = inoise;
            md.leftTimeDilSources   = timeDilSource.at("left");
            md.rightTimeDilSources  = timeDilSource.at("right");

            std::stringstream ss;
            ss << md.gamma << "_";
            for (unsigned int mu = 0; mu < md.momentum.size(); ++mu)
                ss << md.momentum[mu] << ((mu == md.momentum.size() - 1) ? "" : "_");
            std::string groupName = ss.str();

            // io init
            std::string outStem = outputMFStem_ + "/noise" + std::to_string(inoise[0]) + "_" + std::to_string(inoise[1]) + "/";
            Hadrons::mkdir(outStem);
            std::string mfName = groupName+"_"+dmf_case_.at("left")+"-"+dmf_case_.at("right")+".h5";
            A2AMatrixIo<ComplexF> matrixIo(outStem+mfName, groupName, eff_nt_, dilutionSize_ls_.at("left"), dilutionSize_ls_.at("right"));
            matrixIoTable.push_back(matrixIo);
            //initialize file with no outputName group (containing atributes of momentum and gamma) but no dataset inside
            if(g->IsBoss())
            {
                startTimer("IO: total");
                startTimer("IO: file creation");
                matrixIoTable.back().initFile(md);
                stopTimer("IO: file creation");
                stopTimer("IO: total");
            }
        }

        LOG(Message) << "Noise pair: " << inoise << std::endl;
        LOG(Message) << "Gamma:" << std::endl << gamma_ << std::endl;
        LOG(Message) << "momenta:" << std::endl << momenta_ << std::endl;

        //computation of distvectors
        computation.distVec(distVector, noise, inoise, peramb, epack, timeDilSource);

        // computing mesonfield blocks and saving to disk
        LOG(Message) << "Time-dilution blocks computation starting..." << std::endl;
        computation.execute(matrixIoTable, st_, distVector, noise, gamma_, phase, nExt_, nStr_, dilutionSize_ls_, eff_nt_, timeDilSource, this);

        LOG(Message) << "Meson fields saved at " << outputMFStem_ << std::endl;
    }
    LOG(Message) << "A2AUtils::MesonField kernel executed " << computation.global_counter << " times over " << cacheSize_ << "^2 cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) "          << computation.global_flops/computation.global_counter        << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) "           << computation.global_bytes/computation.global_counter        << " GB/s/node "    << std::endl;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
