#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>

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

class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma,
                                    std::vector<int>,   noise_pair,
                                    )
};

class DistilMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPar,
                                    std::string,    OutputStem,
                                    std::string,    MesonFieldCase,
                                    std::string,    LapEvec,
                                    std::string,    LeftPeramb,
                                    std::string,    RightPeramb,
                                    std::string,    LeftDPar,
                                    std::string,    RightDPar,
                                    std::vector<std::string>, NoisePairs,
                                    std::vector<std::string>, SourceTimesLeft,
                                    std::vector<std::string>, SourceTimesRight,
                                    std::string,    Gamma,
                                    std::vector<std::string>, Momenta,
                                    int,            BlockSize,
                                    int,            CacheSize,
                                    std::string,    LeftNoise,
                                    std::string,    RightNoise,)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DNoise;
    typedef std::vector<std::set<unsigned int>> TDilutionMap;
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
    Vector<ComplexF>                    blockbuf_;
    Vector<Complex>                     cachebuf_;
    int                                 eff_nt_;
    std::vector<std::vector<int>>       noisePairs_;           // read from extermal object (diluted noise class)
    TDilutionMap                        st_;
    std::string                         outputMFStem_;
    bool                                hasPhase_{false};
    int                                 dilutionSize_LS_;
    std::map<std::string, std::string>  noiseInput_  ;
    std::map<std::string, std::string>  perambInput_ ;
    std::vector<std::string>            sides        ;
    
};

MODULE_REGISTER_TMP(DistilMesonField, TDistilMesonField<FIMPL>, MDistil);

// aux class
template <typename FImpl>
class DMesonFieldHelper
{
public:
    typedef typename TDistilMesonField<FImpl>::DNoise DNoise;
    typedef typename TDistilMesonField<FImpl>::TDilutionMap TDilutionMap;
private:
    int nt_;
    int nd_;
    std::string mfCase_;
    TDilutionMap timeMapl_, timeMapr_;
public:
    
    DMesonFieldHelper(DNoise &nl, DNoise &nr, std::string in_case)
    : mfCase_(in_case) , timeMapl_(nl.getMap()[0]) , timeMapr_(nr.getMap()[0])
    {
        nt_ = nr.getNt();
        nd_ = nr.getGrid()->Nd();
        assert( timeMapl_.size() == timeMapr_.size() );  //number of partitions should be the same (?)

        // check mesonfield case
        if(mfCase_=="phi phi" || mfCase_=="phi rho" || mfCase_=="rho phi" || mfCase_=="rho rho")
        {
            LOG(Message) << "Meson field case checked: " << mfCase_ << std::endl;
        }
        else
        {
            HADRONS_ERROR(Argument,"Bad meson field case");
        }
    }

    int computeTimeDimension(TDilutionMap st)
    {
        // compute eff_nt (<=nt_), the number of non-zero timeslices in the final object, when there's at least one rho involved
        int eff_nt = 1;
        if(mfCase_=="rho rho" || mfCase_=="rho phi" || mfCase_=="phi rho")
        {
            for(auto &e : st)
                e.size() > eff_nt ? eff_nt = e.size() : NULL;      //get the highest possible eff_nt from st
        }
        else
        {
            eff_nt = nt_;
        }
        return eff_nt;
    }

    TDilutionMap getSourceTimes()
    {
        // find maximal subset, may be useless in the future
        std::vector<std::set<unsigned int>> st ;
        st.clear();
        for(int p=0 ; p<timeMapl_.size() ; p++)
        {
            std::set<unsigned int> outset;
            outset.clear();
            std::set_intersection(timeMapl_[p].begin(), timeMapl_[p].end(), 
                                    timeMapr_[p].begin(), timeMapr_[p].end(),
                                    std::inserter(outset, outset.begin()));
            st.push_back(outset);
        }
        return st;
    }

    std::vector<std::vector<RealF>> parseMomenta(std::vector<std::string> inputP)
    {
        std::vector<std::vector<RealF>> m;
        m.clear();
        for(auto &p_string : inputP)
        {
            auto p = strToVec<RealF>(p_string);

            if (p.size() != nd_ - 1)
            {
                HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                    + " components instead of " 
                                    + std::to_string(nd_ - 1));
            }
            m.push_back(p);
        }
        return(m);
    }

    std::vector<Gamma::Algebra> parseGamma(std::string inputG)
    {  
        std::vector<Gamma::Algebra> g;
        if (inputG == "all")
        {
            g = {
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
            g = strToVec<Gamma::Algebra>(inputG);
        }
        return(g);
    }

    std::vector<std::vector<int>> parseNoisePairs(std::vector<std::string> inputN , std::map<std::string,std::string> caseMap , std::map<std::string,int> noiseDim )
    {
        const std::vector<std::string> sides = {"left","right"};
        std::vector<std::vector<int>> nPairs;
        nPairs.clear();
        for(auto &npair : inputN)
        {
            nPairs.push_back(strToVec<int>(npair));
            std::map<std::string, int>  noiseMapTemp = { {"left", nPairs.back()[0]} , {"right",nPairs.back()[1]} };
            for(auto &side : sides){
                if(caseMap.at(side)=="phi")    // turn this into macro?
                {
                    if( noiseMapTemp.at(side) >= noiseDim.at(side) )    // verify if input noise number is valid ( < tensor nnoise dimension)
                    {
                        HADRONS_ERROR(Size,"Noise pair element " + std::to_string(noiseMapTemp.at(side)) + "(>=" +std::to_string(noiseDim.at(side)) + ") unavailable in input tensor");
                    }
                }
                else
                {
                    if( noiseMapTemp.at(side) >= noiseDim.at(side) )
                    {
                        HADRONS_ERROR(Size,"Noise pair element " + std::to_string(noiseMapTemp.at(side)) + "(>=" +std::to_string(noiseDim.at(side)) + ") unavailable in input tensor");
                    }
                }
                if( noiseMapTemp.at(side) < 0)
                {
                    HADRONS_ERROR(Size,"Negative noise pair element");
                }
            }
        }
        return(nPairs);
    }

    void computePhase(std::vector<std::vector<RealF>> momenta_, typename FImpl::ComplexField &coor, std::vector<int> dim, std::vector<typename FImpl::ComplexField> &phase)
    {
        Complex           i(0.0,1.0);
        for (unsigned int j = 0; j < momenta_.size(); ++j)
        {
            phase[j] = Zero();
            for(unsigned int mu = 0; mu < momenta_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                phase[j] = phase[j] + (momenta_[j][mu]/dim[mu])*coor;
            }
            phase[j] = exp((Real)(2*M_PI)*i*phase[j]);
        }
    }
};

/******************************************************************************
 *                 TDistilMesonField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMesonField<FImpl>::TDistilMesonField(const std::string name)
: Module<DistilMesonFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getInput(void)
{   
    return{par().LapEvec, par().LeftPeramb, par().RightPeramb, par().LeftDPar, par().RightDPar, par().LeftNoise, par().RightNoise};
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
    std::map<std::string, std::string>  noiseInput_     = {{"left",par().LeftNoise},{"right",par().RightNoise}};
    std::map<std::string, std::string>  perambInput_    = {{"left",par().LeftPeramb},{"right",par().RightPeramb}};
    std::vector<std::string>            sides           = {"left","right"};

    outputMFStem_ = par().OutputStem;
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)
    // hard-coded dilution scheme (assuming dpL=dpR for now)
    const DistilParameters &dpL = envGet(DistilParameters, par().LeftDPar);
    const DistilParameters &dpR = envGet(DistilParameters, par().RightDPar);

    DNoise &noisel = envGet( DNoise , par().LeftNoise);
    DNoise &noiser = envGet( DNoise , par().RightNoise);

    DMesonFieldHelper<FImpl> helper(noisel, noiser, par().MesonFieldCase);

    dmf_case_.emplace("left" , par().MesonFieldCase.substr(0,3));   //left
    dmf_case_.emplace("right" , par().MesonFieldCase.substr(4,7));  //right
    
    LOG(Message) << "Time dimension = " << eff_nt_ << std::endl;
    LOG(Message) << "Selected block size: " << par().BlockSize << std::endl;
    LOG(Message) << "Selected cache size: " << par().CacheSize << std::endl;

    // parse source times
    // outermost dimension is the time-dilution index, innermost one are the non-zero source timeslices
    // in phi phi, save all timeslices, but in the other cases save only the non-zero  ones...
    st_ = helper.getSourceTimes();
    eff_nt_ = helper.computeTimeDimension(st_);
    
    // parse and validate input
    std::map<std::string, int> noiseDimension = {{"left",0},{"right",0}}; ;
    for(auto &side : sides)
    {
        if(dmf_case_.at(side)=="phi")
        {
            auto &inPeramb = envGet(PerambTensor , perambInput_.at(side));
            noiseDimension.at(side) = inPeramb.tensor.dimensions().at(3);
        }
        else
        {
            auto &inNoise = envGet(DNoise , noiseInput_.at(side));
            noiseDimension.at(side) = inNoise.size();
        }
    }
    noisePairs_ = helper.parseNoisePairs(par().NoisePairs , dmf_case_ , noiseDimension);
    
    blockSize_ = {par().BlockSize};
    cacheSize_ ={par().CacheSize};

    // momenta and gamma parse -> turn into method
    momenta_ = helper.parseMomenta(par().Momenta);
    gamma_ = helper.parseGamma(par().Gamma);
    nExt_ = momenta_.size(); //noise pairs computed independently, but can optmize embedding it into nExt??
    nStr_ = gamma_.size();
    
    //populate matrix sets
    blockbuf_.resize(nExt_*nStr_*eff_nt_*par().BlockSize*par().BlockSize);
    cachebuf_.resize(nExt_*nStr_*env().getDim(g->Nd() - 1)*par().CacheSize*par().CacheSize);
    
    assert( noisel.getMap()[1].size()*noisel.getMap()[2].size() == noiser.getMap()[1].size()*noiser.getMap()[2].size() );
    dilutionSize_LS_ = noisel.getMap()[1].size()*noisel.getMap()[2].size();
    
    envTmp(FermionField,                    "fermion3dtmp",         1, g3d);
    envTmp(ColourVectorField,               "fermion3dtmp_nospin",  1, g3d);
    envTmp(ColourVectorField,               "evec3d",               1, g3d);
    envTmp(std::vector<FermionField>,       "left",                 1, st_.size()*dilutionSize_LS_, g); //stL.size()*dpL.LI*dpL.SI, g);
    envTmp(std::vector<FermionField>,       "right",                1, st_.size()*dilutionSize_LS_, g); //stR.size()*dpR.LI*dpR.SI, g);
    envTmpLat(ComplexField, "coor");
    envCache(std::vector<ComplexField>,     "phasename",            1, momenta_.size(), g);
    envTmpLat(FermionField,                 "fermion4dtmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    // temps
    envGetTmp(FermionField,                 fermion3dtmp);
    envGetTmp(ColourVectorField,            fermion3dtmp_nospin);
    envGetTmp(FermionField,                 fermion4dtmp);
    envGetTmp(ColourVectorField,            evec3d);
    envGetTmp(std::vector<FermionField>,    left);
    envGetTmp(std::vector<FermionField>,    right);

    int vol = env().getGrid()->_gsites;
    const unsigned int nd = env().getGrid()->Nd();
    const int nt = env().getDim(nd - 1);
    const int Ntlocal = env().getGrid()->LocalDimensions()[nd - 1];
    const int Ntfirst = env().getGrid()->LocalStarts()[nd - 1];

    // hard-coded dilution schem &  assuming dilution_left == dilution_right; other cases?...
    // replace by noise class
    const DistilParameters &dpL = envGet(DistilParameters, par().LeftDPar);
    const DistilParameters &dpR = envGet(DistilParameters, par().RightDPar);

    DistillationNoise<FImpl> &noisel = envGet( DistillationNoise<FImpl> , par().LeftNoise);
    DistillationNoise<FImpl> &noiser = envGet( DistillationNoise<FImpl> , par().RightNoise);
    DMesonFieldHelper<FImpl> helper(noisel, noiser,  par().MesonFieldCase);
    noisel.dumpDilutionMap();

    int nInversions = MIN(dpL.inversions,dpR.inversions);   //number of inversions on perambulator
    assert(nInversions >= st_.size());   //just check this on the case where there is a perambulator, check somewhere else

    auto &epack = envGet(LapEvecs, par().LapEvec);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");

    //compute momentum phase
    if (!hasPhase_)
    {
        startTimer("momentum phases");
        envGetTmp(ComplexField, coor);
        helper.computePhase(momenta_, coor, env().getDim(), phase);
        hasPhase_ = true;
        stopTimer("momentum phases");
    }

    int nVec        = dpL.nvec;
    int nNoiseLeft  = dpL.nnoise;
    int nNoiseRight = dpR.nnoise;

    // do not use operator []!! similar but better way to do that? maybe map of pointers?
    std::map<std::string, std::vector<FermionField>&>       lrDistVector    = {{"left",left},{"right",right}};
    std::map<std::string, DNoise&>                          lrNoise         = {{"left",noisel},{"right",noiser}};

    long    global_counter = 0;
    double  global_flops = 0.0;
    double  global_bytes = 0.0;

    for(auto &inoise : noisePairs_)
    {
        // set up io object and metadata for all gamma/momenta
        std::vector<A2AMatrixIo<ComplexF>> matrixIoTable;
        DistilMesonFieldMetadata md;
        for(int iExt=0; iExt<nExt_; iExt++)
        for(int iStr=0; iStr<nStr_; iStr++)
        {
            // metadata;
            md.momentum = momenta_[iExt];
            md.gamma = gamma_[iStr];
            md.noise_pair = inoise;
            // md.dilution_time = noisel.getMap()[0];
            // md.dilution_lap = noisel.getMap()[1];
            // md.dilution_spin = noisel.getMap()[2];

            std::stringstream ss;
            ss << md.gamma << "_";
            for (unsigned int mu = 0; mu < md.momentum.size(); ++mu)
                ss << md.momentum[mu] << ((mu == md.momentum.size() - 1) ? "" : "_");
            std::string groupName = ss.str();

            //init file here (do not create dataset yet)
            //IO configuration for fixed test gamma and momentum
            
            std::string outputStem = outputMFStem_ + "/noise" + std::to_string(inoise[0]) + "_" + std::to_string(inoise[1]) + "/";
            Hadrons::mkdir(outputStem);
            std::string mfName = groupName+"_"+dmf_case_.at("left")+"-"+dmf_case_.at("right")+".h5";
            A2AMatrixIo<ComplexF> matrixIo(outputStem+mfName, groupName, eff_nt_, dilutionSize_LS_, dilutionSize_LS_);  // automatise name choice according to momenta_ and gamma_
            matrixIoTable.push_back(matrixIo);

            //initialize file with no outputName group (containing atributes of momentum and gamma) but no dataset inside
            if(env().getGrid()->IsBoss())
            {
                startTimer("IO: total");
                startTimer("IO: file creation");
                matrixIoTable.back().initFile(md);
                stopTimer("IO: file creation");
                stopTimer("IO: total");
            }
        }

        std::map<std::string,int&> iNoise = {{"left",inoise[0]},{"right",inoise[1]}}; //do not use []

        LOG(Message) << "Noise pair: " << inoise << std::endl;

        LOG(Message) << "Gamma:" << std::endl;
        LOG(Message) << gamma_ << std::endl;
        LOG(Message) << "Momenta:" << std::endl;
        LOG(Message) << momenta_ << std::endl;

        for(auto &side : sides)    // computation, still ignoring gamma5 hermiticity
        {
            for(int iD=0 ; iD<lrNoise.at(side).dilutionSize() ; iD++)  // computation of phi or rho
            {
                //each dtL, dtR pair corresponds to a different dataset here
                // for (int dt = 0; dt<st_.size(); dt++)         //loop over time dilution index
                // {
                    std::array<unsigned int,3> c = lrNoise.at(side).dilutionCoordinates(iD);
                    int dt = c[0] , dk = c[1] , ds = c[2];
                    lrDistVector.at(side)[iD] = Zero();
                    if(dmf_case_.at(side)=="phi")
                    {
                        auto &inTensor = envGet(PerambTensor , perambInput_.at(side));
                        for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)   //loop over (local) timeslices
                        {
                            fermion3dtmp = Zero();
                            for (int k = 0; k < nVec; k++)
                            {
                                ExtractSliceLocal(evec3d,epack.evec[k],0,t-Ntfirst,nd - 1);
                                fermion3dtmp += evec3d * inTensor.tensor(t, k, dk, iNoise.at(side), dt, ds);
                            }
                            InsertSliceLocal(fermion3dtmp,lrDistVector.at(side)[iD],0,t-Ntfirst,nd - 1);
                        }
                    }
                    else if(dmf_case_.at(side)=="rho"){
                        lrDistVector.at(side)[iD] = lrNoise.at(side).makeSource(iD, iNoise.at(side));
                    }
                // }
            }
        }

        // computing mesonfield blocks and saving to disk
        for (int dtL = 0; dtL < lrNoise.at("left").getMap()[0].size() ; dtL++)
        for (int dtR = 0; dtR < lrNoise.at("right").getMap()[0].size() ; dtR++)
        {
            if(!(par().MesonFieldCase=="rho rho" && dtL!=dtR))
            {
                std::string datasetName = "dtL"+std::to_string(dtL)+"_dtR"+std::to_string(dtR);
                LOG(Message) << "- Computing dilution dataset " << datasetName << "..." << std::endl;

                // int nblocki = left.size()/blockSize_ + (((left.size() % blockSize_) != 0) ? 1 : 0);
                // int nblockj = right.size()/blockSize_ + (((right.size() % blockSize_) != 0) ? 1 : 0);
                int nblocki = dilutionSize_LS_/blockSize_ + (((dilutionSize_LS_ % blockSize_) != 0) ? 1 : 0);
                int nblockj = dilutionSize_LS_/blockSize_ + (((dilutionSize_LS_ % blockSize_) != 0) ? 1 : 0);


                // loop over blocks in the current time-dilution block
                for(int iblock=0 ; iblock<dilutionSize_LS_ ; iblock+=blockSize_) //set according to memory size
                for(int jblock=0 ; jblock<dilutionSize_LS_ ; jblock+=blockSize_)
                {
                    int iblockSize = MIN(dilutionSize_LS_-iblock,blockSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                    int jblockSize = MIN(dilutionSize_LS_-jblock,blockSize_);
                    A2AMatrixSet<ComplexF> block(blockbuf_.data(), nExt_ , nStr_ , eff_nt_, iblockSize, jblockSize);

                    LOG(Message) << "Distil matrix block " 
                    << jblock/blockSize_ + nblocki*iblock/blockSize_ + 1 
                    << "/" << nblocki*nblockj << " [" << iblock << " .. " 
                    << iblock+iblockSize-1 << ", " << jblock << " .. " << jblock+jblockSize-1 << "]" 
                    << std::endl;

                    LOG(Message) << "Block size = "         << eff_nt_*iblockSize*jblockSize*sizeof(ComplexF) << "MB/momentum/gamma" << std::endl;
                    LOG(Message) << "Cache blocks size = "   << nt*cacheSize_*cacheSize_*sizeof(ComplexD) << "MB/momentum/gamma" << std::endl;  //remember to change this in case I change chunk size from nt to something else

                    double flops       = 0.0;
                    double bytes       = 0.0;
                    double time_kernel = 0.0;
                    double nodes    = env().getGrid()->NodeCount();

                    // loop over cache_ blocks in the current block
                    for(int icache=0 ; icache<iblockSize ; icache+=cacheSize_)   //set according to cache_ size
                    for(int jcache=0 ; jcache<jblockSize ; jcache+=cacheSize_)
                    {
                        int icacheSize = MIN(iblockSize-icache,cacheSize_);      // icacheSize is the size of the current cache_ block (indexed by ii); N_ii-ii is the size of the remainder cache_ block
                        int jcacheSize = MIN(jblockSize-jcache,cacheSize_);
                        A2AMatrixSet<Complex> blockCache(cachebuf_.data(), nExt_, nStr_, nt, icacheSize, jcacheSize);

                        double timer = 0.0;
                        startTimer("kernel");
                        A2Autils<FImpl>::MesonField(blockCache, &left[dtL*dilutionSize_LS_+iblock+icache], &right[dtR*dilutionSize_LS_+jblock+jcache], gamma_, phase, nd - 1, &timer);
                        stopTimer("kernel");

                        time_kernel += timer;
                        
                        // nExt is currently # of momenta , nStr is # of gamma matrices
                        flops += vol*(2*8.0+6.0+8.0*nExt_)*icacheSize*jcacheSize*nStr_;
                        bytes += vol*(12.0*sizeof(ComplexD))*icacheSize*jcacheSize
                              +  vol*(2.0*sizeof(ComplexD)*nExt_)*icacheSize*jcacheSize*nStr_;

                        // std::cout<< "block dimensions " << block.dimensions().at(2) << std::endl;
                        // std::cout<< "blockcache dimensions " << blockCache.dimensions().at(2) << std::endl << std::cin.get();

                        // loop through the cacheblock (inside them) and point blockCache to block
                        startTimer("cache copy");
                        if(par().MesonFieldCase=="phi phi")
                        {
                            thread_for_collapse( 5, iExt ,nExt_,{
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int t=0 ; t<nt ; t++)
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                                block(iExt,iStr,t,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,t,iicache,jjcache);
                            });
                        }
                        else
                        {
                            std::vector<std::vector<unsigned int>> v;
                            v.clear();
                            for(int i ; i<st_.size() ; i++){
                                std::vector<unsigned int> temp(st_[i].begin() , st_[i].end());
                                v.push_back(temp);
                            }

                            // for ( uint64_t iExt=0;iExt<nExt_;iExt++)
                            thread_for_collapse( 5, iExt ,nExt_,{
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int it=0 ; it<eff_nt_ ; it++)  //only wish to copy non-zero timeslices to block
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                                block(iExt,iStr,it,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,v[dtL][it],iicache,jjcache);
                            });
                        }
                        stopTimer("cache copy");
                    }

                    LOG(Message) << "Kernel perf (flops) " << flops/time_kernel/1.0e3/nodes 
                                << " Gflop/s/node " << std::endl;
                    LOG(Message) << "Kernel perf (read) " << bytes/time_kernel*0.000931322574615478515625/nodes //  1.0e6/1024/1024/1024/nodes
                                << " GB/s/node "  << std::endl;
                    global_counter++;
                    global_flops += flops/time_kernel/1.0e3/nodes ;
                    global_bytes += bytes/time_kernel*0.000931322574615478515625/nodes ; // 1.0e6/1024/1024/1024/nodes

                    // saving current block to disk
                    LOG(Message) << "Writing block to disk" << std::endl;
                    startTimer("IO: total");
                    startTimer("IO: write block");
                    double ioTime = -getDTimer("IO: write block");
#ifdef HADRONS_A2AM_PARALLEL_IO
                    //parallel io
                    int inode = env().getGrid()->ThisRank();
                    int nnode = env().getGrid()->RankCount(); 
                    LOG(Message) << "Starting parallel IO. Rank count=" << nnode  << std::endl;
                    env().getGrid()->Barrier();
                    for(int ies=inode ; ies<nExt_*nStr_ ; ies+=nnode){
                        int iExt = ies/nStr_;
                        int iStr = ies%nStr_;
                        if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName, cacheSize_);   //set surface chunk size as cacheSize_ (the chunk itself is 3D)
                        }
                        else{
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName);
                        }
                    }
                    env().getGrid()->Barrier();
#else
                    // serial io, can remove later
                    LOG(Message) << "Starting serial IO" << std::endl;
                    for(int iExt=0; iExt<nExt_; iExt++)
                    for(int iStr=0; iStr<nStr_; iStr++)
                    {
                        if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName, cacheSize_);   //set surface chunk size as cacheSize_ (the chunk itself is 3D)
                        }
                        else{
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName);
                        }
                    }
#endif
                    stopTimer("IO: total");
                    stopTimer("IO: write block");
                    ioTime    += getDTimer("IO: write block");
                    int bytesBlockSize  = static_cast<double>(nExt_*nStr_*eff_nt_*iblockSize*jblockSize*sizeof(ComplexF));
                    LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                    << ioTime  << " us (" 
                                    << bytesBlockSize/ioTime*0.95367431640625 // 1.0e6/1024/1024
                                    << " MB/s)" << std::endl;
                }
            }
        }
        LOG(Message) << "Meson fields saved at " << outputMFStem_ << std::endl;
    }
    LOG(Message) << "MesonField kernel executed " << global_counter << " times on " << cacheSize_ << "^2 cache blocks" << std::endl;
    LOG(Message) << "Average kernel perf (flops) " << global_flops/global_counter << " Gflop/s/node " << std::endl;
    LOG(Message) << "Average kernel perf (read) " << global_bytes/global_counter  << " GB/s/node "  << std::endl;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
