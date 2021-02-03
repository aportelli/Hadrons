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



class DistilMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPar,
                                    std::string,    outputStem,
                                    std::string,    mesonFieldCase,
                                    std::string,    lapEvec,
                                    std::string,    leftPeramb,
                                    std::string,    rightPeramb,
                                    std::string,    leftNoise,
                                    std::string,    rightNoise,
                                    std::vector<std::string>, noisePairs,
                                    std::string,    gamma,
                                    std::vector<std::string>, momenta,
                                    int,            blockSize,
                                    int,            cacheSize,)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef std::vector<std::vector<unsigned int>> TimeSliceMap;
    typedef typename DistillationNoise::Index Index;
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
    TimeSliceMap                        st_;
    std::string                         outputMFStem_;
    bool                                hasPhase_{false};
    int                                 dilutionSize_LS_;
    std::map<std::string, std::string>  noiseInput_  ;
    std::map<std::string, std::string>  perambInput_ ;
    std::vector<std::string>            sides        ;
    
};

MODULE_REGISTER_TMP(DistilMesonField, TDistilMesonField<FIMPL>, MDistil);

// metadata class
class DistilMesonFieldMetadata: Serializable
{
public:
    typedef typename TDistilMesonField<FImpl>::TimeSliceMap TimeSliceMap;
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma,
                                    std::vector<int>,   noise_pair,
                                    TimeSliceMap,       time_dilution,
                                    )
};

// aux class
template <typename FImpl>
class DMesonFieldHelper
{
public:
    typedef typename TDistilMesonField<FImpl>::DistillationNoise DistillationNoise;
    typedef typename TDistilMesonField<FImpl>::TimeSliceMap TimeSliceMap;
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename TDistilMesonField<FImpl>::Index Index;
private:
    int nt_;
    int nd_;
    std::map<std::string,std::string> dmfCase_;
    TimeSliceMap timeMapl_, timeMapr_;
    const std::vector<std::string> sides = {"left","right"};
public:
    
    DMesonFieldHelper(DistillationNoise &nl, DistillationNoise &nr, std::string in_case)
    : timeMapl_( timeSliceMap(nl) ) , timeMapr_( timeSliceMap(nr) )
    {
        nt_ = nr.getNt();
        nd_ = nr.getGrid()->Nd();
        assert( timeMapl_.size() == timeMapr_.size() );  //number of partitions should be the same (?)

        // check mesonfield case
        if(!(in_case=="phi phi" || in_case=="phi rho" || in_case=="rho phi" || in_case=="rho rho"))
        {
            HADRONS_ERROR(Argument,"Bad meson field case");
        }

        dmfCase_.emplace("left"  , in_case.substr(0,3));
        dmfCase_.emplace("right" , in_case.substr(4,7));
    }

    std::map<std::string,std::string> getValidCase(){
        return dmfCase_;
    }

    int computeTimeDimension(TimeSliceMap st)
    {
        // compute eff_nt (<=nt_), the number of non-zero timeslices in the final object, when there's at least one rho involved
        int eff_nt = 1;
        if(dmfCase_.at("left")=="rho" || dmfCase_.at("right")=="rho")
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

    TimeSliceMap getSourceTimes()
    {
        // find intersection, may be useless in the future
        TimeSliceMap st;
        for(int p=0 ; p<timeMapl_.size() ; p++)
        {
            std::vector<unsigned int> temp;
            std::set_intersection(timeMapl_[p].begin(), timeMapl_[p].end(), 
                                    timeMapr_[p].begin(), timeMapr_[p].end(),
                                    std::back_inserter(temp));
            st.push_back(temp);
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
        
        std::vector<std::vector<int>> nPairs;
        nPairs.clear();
        for(auto &npair : inputN)
        {
            nPairs.push_back(strToVec<int>(npair));
            std::map<std::string, int>  noiseMapTemp = { {"left", nPairs.back()[0]} , {"right",nPairs.back()[1]} };
            for(auto &side : sides){
                if( noiseMapTemp.at(side) >= noiseDim.at(side) )    // verify if input noise number is valid ( < tensor nnoise dimension)
                {
                    HADRONS_ERROR(Size,"Noise pair element " + std::to_string(noiseMapTemp.at(side)) + "(>=" +std::to_string(noiseDim.at(side)) + ") unavailable in input tensor");
                }
                if( noiseMapTemp.at(side) < 0)
                {
                    HADRONS_ERROR(Size,"Negative noise pair element");
                }
            }
        }
        return(nPairs);
    }

    void computePhase(std::vector<std::vector<RealF>> momenta, ComplexField &coor, std::vector<int> dim, std::vector<ComplexField> &phase)
    {
        Complex           i(0.0,1.0);
        for (unsigned int j = 0; j < momenta.size(); ++j)
        {
            phase[j] = Zero();
            for(unsigned int mu = 0; mu < momenta[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                phase[j] = phase[j] + (momenta[j][mu]/dim[mu])*coor;
            }
            phase[j] = exp((Real)(2*M_PI)*i*phase[j]);
        }
    }

    TimeSliceMap timeSliceMap(DistillationNoise &n)
    {
        TimeSliceMap m;
        for(unsigned int it=0 ; it<n.dilutionSize(Index::t) ; it++)
        {
            std::vector<unsigned int> temp = n.timeSlices(it);
            m.push_back(temp);
        }
        return m;
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
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    DMesonFieldHelper<FImpl> helper(noisel, noiser, par().mesonFieldCase);

    dmf_case_ = helper.getValidCase();

    std::map<std::string, std::string>  noiseInput_     = {{"left",par().leftNoise},{"right",par().rightNoise}};
    std::map<std::string, std::string>  perambInput_    = {{"left",par().leftPeramb},{"right",par().rightPeramb}};
    std::vector<std::string>            sides           = {"left","right"};

    outputMFStem_ = par().outputStem;
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)

    // parse source times
    // outermost dimension is the time-dilution index, innermost one are the non-zero source timeslices
    // in phi phi, save all timeslices, but in the other cases save only the non-zero  ones...
    st_ = helper.getSourceTimes();
    for(auto i : st_)
        std::cout << i << std::endl;

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
            auto &inNoise = envGet(DistillationNoise , noiseInput_.at(side));
            noiseDimension.at(side) = inNoise.size();
        }
    }
    noisePairs_ = helper.parseNoisePairs(par().noisePairs , dmf_case_ , noiseDimension);
    
    blockSize_ = {par().blockSize};
    cacheSize_ ={par().cacheSize};

    // momenta and gamma parse
    momenta_ = helper.parseMomenta(par().momenta);
    gamma_ = helper.parseGamma(par().gamma);
    nExt_ = momenta_.size(); //noise pairs computed independently, but can optmize embedding it into nExt??
    nStr_ = gamma_.size();
    
    //populate matrix sets
    blockbuf_.resize(nExt_*nStr_*eff_nt_*blockSize_*blockSize_);
    cachebuf_.resize(nExt_*nStr_*env().getDim(g->Nd() - 1)*cacheSize_*cacheSize_);
    
    assert( noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s) == noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s) );
    dilutionSize_LS_ = noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s);
    
    envTmp(FermionField,                    "fermion3dtmp",         1, g3d);
    envTmp(ColourVectorField,               "evec3d",               1, g3d);
    envTmp(std::vector<FermionField>,       "dvl",                  1, st_.size()*dilutionSize_LS_, g);
    envTmp(std::vector<FermionField>,       "dvr",                  1, st_.size()*dilutionSize_LS_, g);
    envTmpLat(ComplexField, "coor");
    envCache(std::vector<ComplexField>,     "phasename",            1, momenta_.size(), g);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    // temps
    envGetTmp(FermionField,                 fermion3dtmp);
    envGetTmp(ColourVectorField,            evec3d);
    envGetTmp(std::vector<FermionField>,    dvl);
    envGetTmp(std::vector<FermionField>,    dvr);
    auto &epack = envGet(LapEvecs, par().lapEvec);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");

    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    DMesonFieldHelper<FImpl> helper(noisel, noiser,  par().mesonFieldCase);
    // do not use operator []!! similar but better way to do that? maybe map of pointers?
    std::map<std::string, std::vector<FermionField>&>       distVector    = {{"left",dvl}  ,{"right",dvr}};
    std::map<std::string, DistillationNoise&>               noise         = {{"left",noisel},{"right",noiser}};

    int nVec = epack.evec.size();
    int vol = env().getGrid()->_gsites;
    const unsigned int nd = env().getGrid()->Nd();
    const int nt = env().getDim(nd - 1);
    const int Ntlocal = env().getGrid()->LocalDimensions()[nd - 1];
    const int Ntfirst = env().getGrid()->LocalStarts()[nd - 1];

    //compute momentum phase
    if (!hasPhase_)
    {
        startTimer("momentum phases");
        envGetTmp(ComplexField, coor);
        helper.computePhase(momenta_, coor, env().getDim(), phase);
        hasPhase_ = true;
        stopTimer("momentum phases");
    }

    long    global_counter = 0;
    double  global_flops = 0.0;
    double  global_bytes = 0.0;

    LOG(Message) << "Meson field case: " << par().mesonFieldCase << std::endl;
    LOG(Message) << "Time dimension = " << eff_nt_ << std::endl;
    LOG(Message) << "Selected block size: " << par().blockSize << std::endl;
    LOG(Message) << "Selected cache size: " << par().cacheSize << std::endl;

    for(auto &inoise : noisePairs_)
    {
        // set up io object and metadata for all gamma/momenta -> turn into method
        std::vector<A2AMatrixIo<ComplexF>> matrixIoTable;
        DistilMesonFieldMetadata md;
        for(int iExt=0; iExt<nExt_; iExt++)
        for(int iStr=0; iStr<nStr_; iStr++)
        {
            // metadata;
            md.momentum = momenta_[iExt];
            md.gamma = gamma_[iStr];
            md.noise_pair = inoise;
            md.dilution_time = st_;

            std::stringstream ss;
            ss << md.gamma << "_";
            for (unsigned int mu = 0; mu < md.momentum.size(); ++mu)
                ss << md.momentum[mu] << ((mu == md.momentum.size() - 1) ? "" : "_");
            std::string groupName = ss.str();

            // io init
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
        LOG(Message) << "momenta:" << std::endl;
        LOG(Message) << momenta_ << std::endl;

        for(auto &side : sides)    // computation
        {
            for(int iD=0 ; iD<noise.at(side).dilutionSize() ; iD++)  // computation of phi or rho
            {
                std::array<unsigned int,3> c = noise.at(side).dilutionCoordinates(iD);
                int dt = c[0] , dk = c[1] , ds = c[2];
                distVector.at(side)[iD] = Zero();
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
                        InsertSliceLocal(fermion3dtmp,distVector.at(side)[iD],0,t-Ntfirst,nd - 1);
                    }
                }
                else if(dmf_case_.at(side)=="rho"){
                    distVector.at(side)[iD] = noise.at(side).makeSource(iD, iNoise.at(side));
                }
            }
        }

        // computing mesonfield blocks and saving to disk
        for (int dtL = 0; dtL < noise.at("left").dilutionSize(Index::t) ; dtL++)
        for (int dtR = 0; dtR < noise.at("right").dilutionSize(Index::t) ; dtR++)
        {
            if(!(par().mesonFieldCase=="rho rho" && dtL!=dtR))
            {
                std::string datasetName = "dtL"+std::to_string(dtL)+"_dtR"+std::to_string(dtR);
                LOG(Message) << "- Computing dilution dataset " << datasetName << "..." << std::endl;

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
                    LOG(Message) << "Cache block size = "   << nt*cacheSize_*cacheSize_*sizeof(ComplexD) << "MB/momentum/gamma" << std::endl;  //remember to change this in case I change chunk size from nt to something else

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
                        // assuming certain indexation here! (dt must be the slowest index for this to work; otherwise will have to compute l/r block at each contraction)
                        A2Autils<FImpl>::MesonField(blockCache, &dvl[dtL*dilutionSize_LS_+iblock+icache], &dvr[dtR*dilutionSize_LS_+jblock+jcache], gamma_, phase, nd - 1, &timer);
                        stopTimer("kernel");
                        time_kernel += timer;

                        // nExt is currently # of momenta , nStr is # of gamma matrices
                        flops += vol*(2*8.0+6.0+8.0*nExt_)*icacheSize*jcacheSize*nStr_;
                        bytes += vol*(12.0*sizeof(ComplexD))*icacheSize*jcacheSize
                              +  vol*(2.0*sizeof(ComplexD)*nExt_)*icacheSize*jcacheSize*nStr_;

                        // loop through the cacheblock (inside them) and point blockCache to block
                        startTimer("cache copy");
                        if(par().mesonFieldCase=="phi phi")
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
                            thread_for_collapse( 5, iExt ,nExt_,{
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int it=0 ; it<eff_nt_ ; it++)  //only wish to copy non-zero timeslices to block
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                                block(iExt,iStr,it,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,st_[dtL][it],iicache,jjcache);
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
