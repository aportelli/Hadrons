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
                                    std::string,    gamma,
                                    int,            blockSize,
                                    int,            cacheSize,
                                    std::vector<std::string>, sourceTimes,
                                    std::vector<std::string>, momenta,
                                    std::vector<std::string>, noisePairs,)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef std::vector<std::vector<FermionField>> DistilVector;
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
    int                                 dilutionSize_ls_;
    std::map<std::string, std::string>  noiseInput_  ;
    std::map<std::string, std::string>  perambInput_ ;
    std::vector<std::string>            sides_       ;
    // DmfHelper<FImpl>                   *helper_;
    
};

MODULE_REGISTER_TMP(DistilMesonField, TDistilMesonField<FIMPL>, MDistil);

// metadata class
template <typename FImpl>
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
class DmfHelper
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
    TimeSliceMap noiseTimeMapl_, noiseTimeMapr_;
    const std::vector<std::string> sides = {"left","right"};
public:
    
    DmfHelper(DistillationNoise & nl, DistillationNoise & nr, std::string in_case)
    : noiseTimeMapl_( timeSliceMap(nl) ) , noiseTimeMapr_( timeSliceMap(nr) )
    {
        nt_ = nr.getNt();
        nd_ = nr.getGrid()->Nd();

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

    int computeEffTimeDimension(TimeSliceMap st)
    {
        // compute eff_nt (<=nt_), the number of non-zero timeslices in the final object, when there's at least one rho involved
        int eff_nt = 1;
        if(dmfCase_.at("left")=="rho" || dmfCase_.at("right")=="rho")
        {
            for(auto &e : st)
                e.size() > eff_nt ? eff_nt = e.size() : 0;      //get the highest possible eff_nt from st
        }
        else
        {
            eff_nt = nt_;
        }
        return eff_nt;
    }

    TimeSliceMap getSourceTimes(std::map<std::string, TimeSliceMap> perambTimeMap , TimeSliceMap st_input)
    {//check if noise_st_i contains peramb_st_i (case==phi), take the intersection between the l/r intersection result, check if input is subset of that
        TimeSliceMap st,st_dependencies;
        std::map<std::string, TimeSliceMap> inter = { {"left",{}},{"right",{}} };
        std::map<std::string , TimeSliceMap> noiseTimeMap = { {"left",noiseTimeMapl_},{"right",noiseTimeMapr_} };
        for(auto &s : sides)
        {
            if(dmfCase_.at(s)=="phi")
            {
                // std::cout << "noise/peramb intersection, " << s << std::endl; 
                inter.at(s) = getIntersectionMap(noiseTimeMap.at(s) , perambTimeMap.at(s));
                
                if(inter.at(s).empty())
                {
                    HADRONS_ERROR(Argument,"Time dilution not compatible between noise and perambulator.");
                }
            }
            else
            {
                // std::cout << "no noise/peramb intersection (rho case), " << s << std::endl; 
                inter.at(s) = noiseTimeMap.at(s);
            }
        }

        // std::cout << "l/r dependencies intersection" << std::endl; 
        st_dependencies  = getIntersectionMap(inter.at("left") , inter.at("right"));

        // std::cout << "dependencies/input intersection" << std::endl; 
        st = getIntersectionMap(st_dependencies , st_input);

        return st;
    }
    
    TimeSliceMap getIntersectionMap(TimeSliceMap m1, TimeSliceMap m2)
    {
        TimeSliceMap inter;
        // std::cout << "m1:" << std::endl;
        // for(auto e : m1)
        //     std::cout << e << std::endl;

        // std::cout << "m2:" << std::endl;
        // for(auto e : m1)
        //     std::cout << e << std::endl;

        // std::cout << "inter:" << std::endl;
        for(unsigned int p=0 ; p<m1.size() ; p++)
        {
            std::vector<unsigned int> temp;
            std::set_intersection(m1[p].begin(), m1[p].end(), 
                                m2[p].begin(), m2[p].end(),
                                std::back_inserter(temp));
            // std::cout << temp << std::endl;
            inter.push_back(temp);
        }
        // std::cin.get();
        return inter;
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

    std::vector<std::vector<int>> parseNoisePairs(std::vector<std::string> inputN)
    {
        
        std::vector<std::vector<int>> nPairs;
        nPairs.clear();
        for(auto &npair : inputN)
        {
            nPairs.push_back(strToVec<int>(npair));
            std::map<std::string, int>  noiseMapTemp = { {"left", nPairs.back()[0]} , {"right",nPairs.back()[1]} };
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

    TimeSliceMap timeSliceMap(DistillationNoise & n)
    {
        TimeSliceMap m;
        for(unsigned int it=0 ; it<n.dilutionSize(Index::t) ; it++)
        {
            std::vector<unsigned int> temp = n.timeSlices(it);
            m.push_back(temp);
        }
        return m;
    }

    void dumpMap(TimeSliceMap m)
    {
        std::string o = "{";
        int i=0;
        for(auto & d : m)
        {
            std::string s = "";
            for (auto e: d)
            {
                s += std::to_string(e) + " "; 
            }
            s.pop_back();
            LOG(Message) << "  " << i << ": {" << s << "}" << std::endl;
            i++;
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
    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    DmfHelper<FImpl> helper(noisel, noiser, par().mesonFieldCase);
    dmf_case_ = helper.getValidCase();

    noiseInput_     = {{"left",par().leftNoise},{"right",par().rightNoise}}; //apparently not used
    perambInput_    = {{"left",par().leftPeramb},{"right",par().rightPeramb}};

    outputMFStem_       = par().outputStem;
    GridCartesian *g    = envGetGrid(FermionField);
    GridCartesian *g3d  = envGetSliceGrid(FermionField, g->Nd() - 1);  // 3d grid (as a 4d one with collapsed time dimension)
    
    blockSize_ = par().blockSize;
    cacheSize_ = par().cacheSize;

    // momenta and gamma parse
    momenta_    = helper.parseMomenta(par().momenta);
    gamma_      = helper.parseGamma(par().gamma);
    nExt_       = momenta_.size(); //noise pairs computed independently, but can optmize embedding it into nExt??
    nStr_       = gamma_.size();
       
    //not taking into account different spin/lap dilution on each side, just different time dilutions
    assert( noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s) == noiser.dilutionSize(Index::l)*noiser.dilutionSize(Index::s) );
    dilutionSize_ls_ = noisel.dilutionSize(Index::l)*noisel.dilutionSize(Index::s);
    envTmp(std::vector<std::vector<FermionField>>,                "dvl",          1, noisel.dilutionSize(Index::t) , std::vector<FermionField>(dilutionSize_ls_, g) ); // temp. setting this to full dilution size
    envTmp(std::vector<std::vector<FermionField>>,                "dvr",          1, noiser.dilutionSize(Index::t) , std::vector<FermionField>(dilutionSize_ls_, g) );
    envTmp(FermionField,                "fermion3dtmp", 1, g3d );
    envTmp(ColourVectorField,           "evec3d",       1, g3d );
    envTmpLat(ComplexField,             "coor");
    envCache(std::vector<ComplexField>, "phasename",    1, momenta_.size(), g );
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    // auto & inPeramb = envGet(PerambTensor , par().leftPeramb);
    // for(auto e : inPeramb.MetaData.sourceTimes)
    //     std::cout << e << std::endl;
    // std::cin.get();

    auto &epack = envGet(LapEvecs, par().lapEvec);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");

    int nVec = epack.evec.size();
    int vol = env().getGrid()->_gsites;
    const unsigned int nd = env().getGrid()->Nd();
    const int nt = env().getDim(nd - 1);
    const int Ntlocal = env().getGrid()->LocalDimensions()[nd - 1];
    const int Ntfirst = env().getGrid()->LocalStarts()[nd - 1];

    // temps
    envGetTmp(std::vector<std::vector<FermionField>>,         dvl);
    envGetTmp(std::vector<std::vector<FermionField>>,         dvr);
    envGetTmp(ColourVectorField,    evec3d);
    envGetTmp(FermionField,         fermion3dtmp);
    // do not use operator []!! similar but better way to do that?
    std::map<std::string, std::vector<std::vector<FermionField>> & > distVector = {{"left",dvl}  ,{"right",dvr}};

    DistillationNoise &noisel = envGet( DistillationNoise , par().leftNoise);
    DistillationNoise &noiser = envGet( DistillationNoise , par().rightNoise);
    std::map<std::string, DistillationNoise & >   noise = {{"left",noisel},{"right",noiser}};

    DmfHelper<FImpl> helper(noisel, noiser,  par().mesonFieldCase);

    // parse source times
    std::map<std::string, TimeSliceMap> peramb_st = {{"left",{}},{"right",{}}} ;
    for(auto & s : sides_)
    {
        if(dmf_case_.at(s)=="phi")
        {
            auto & inPeramb = envGet(PerambTensor , perambInput_.at(s));
            peramb_st.at(s) = inPeramb.MetaData.sourceTimes;
        }
    }
    TimeSliceMap st_input;
    for(auto e : par().sourceTimes)
    {
        std::vector<unsigned int> temp = strToVec<unsigned int>(e);
        st_input.push_back(temp);
    }
    st_         = helper.getSourceTimes(peramb_st , st_input);
    eff_nt_     = helper.computeEffTimeDimension(st_);
    noisePairs_ = helper.parseNoisePairs(par().noisePairs);

    //populate matrix sets
    blockbuf_.resize(nExt_*nStr_*eff_nt_*blockSize_*blockSize_);
    cachebuf_.resize(nExt_*nStr_*env().getDim(env().getGrid()->Nd() - 1)*cacheSize_*cacheSize_);

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
    
    LOG(Message) << "Source times:"         << std::endl;
    helper.dumpMap(st_);
    LOG(Message) << "Meson field case: "    << par().mesonFieldCase << std::endl;
    LOG(Message) << "EffTime dimension = "     << eff_nt_ << std::endl;
    LOG(Message) << "Selected block size: " << par().blockSize << std::endl;
    LOG(Message) << "Selected cache size: " << par().cacheSize << std::endl;

    for(auto &inoise : noisePairs_)
    {
        // set up io object and metadata for all gamma/momenta -> turn into method
        std::vector<A2AMatrixIo<ComplexF>> matrixIoTable;
        DistilMesonFieldMetadata<FImpl> md;
        for(int iExt=0; iExt<nExt_; iExt++)
        for(int iStr=0; iStr<nStr_; iStr++)
        {
            // metadata;
            md.momentum         = momenta_[iExt];
            md.gamma            = gamma_[iStr];
            md.noise_pair       = inoise;
            md.time_dilution    = st_;

            std::stringstream ss;
            ss << md.gamma << "_";
            for (unsigned int mu = 0; mu < md.momentum.size(); ++mu)
                ss << md.momentum[mu] << ((mu == md.momentum.size() - 1) ? "" : "_");
            std::string groupName = ss.str();

            // io init
            std::string outputStem = outputMFStem_ + "/noise" + std::to_string(inoise[0]) + "_" + std::to_string(inoise[1]) + "/";
            Hadrons::mkdir(outputStem);
            std::string mfName = groupName+"_"+dmf_case_.at("left")+"-"+dmf_case_.at("right")+".h5";
            A2AMatrixIo<ComplexF> matrixIo(outputStem+mfName, groupName, eff_nt_, dilutionSize_ls_, dilutionSize_ls_);
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

        for(auto &s : sides_)    // computation
        {
            for(int iD=0 ; iD<noise.at(s).dilutionSize() ; iD++)  // computation of phi or rho
            {
                std::array<unsigned int,3> c = noise.at(s).dilutionCoordinates(iD);
                unsigned int dt = c[0] , dl = c[1] , ds = c[2];
                unsigned int iD_ls = ds + noise.at(s).dilutionSize(Index::s) * dl;
                // std::cout << noise.at(s).dilutionSize(Index::s) << " " << iD << " " << iD_ls << " " << dt << std::endl;
                distVector.at(s)[dt][iD_ls] = Zero();
                if(dmf_case_.at(s)=="phi")
                {
                    auto &inTensor = envGet(PerambTensor , perambInput_.at(s));
                    for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)   //loop over (local) timeslices
                    {
                        fermion3dtmp = Zero();
                        for (int k = 0; k < nVec; k++)
                        {
                            ExtractSliceLocal(evec3d,epack.evec[k],0,t-Ntfirst,nd - 1);
                            fermion3dtmp += evec3d * inTensor.tensor(t, k, dl, iNoise.at(s), dt, ds);
                        }
                        InsertSliceLocal(fermion3dtmp,distVector.at(s)[dt][iD_ls],0,t-Ntfirst,nd - 1);
                    }
                }
                else if(dmf_case_.at(s)=="rho"){
                    distVector.at(s)[dt][iD_ls] = noise.at(s).makeSource(iD, iNoise.at(s));
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

                int nblocki = dilutionSize_ls_/blockSize_ + (((dilutionSize_ls_ % blockSize_) != 0) ? 1 : 0);
                int nblockj = dilutionSize_ls_/blockSize_ + (((dilutionSize_ls_ % blockSize_) != 0) ? 1 : 0);

                // loop over blocls in the current time-dilution block
                for(int iblock=0 ; iblock<dilutionSize_ls_ ; iblock+=blockSize_) //set according to memory size
                for(int jblock=0 ; jblock<dilutionSize_ls_ ; jblock+=blockSize_)
                {
                    int iblockSize = MIN(dilutionSize_ls_-iblock,blockSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                    int jblockSize = MIN(dilutionSize_ls_-jblock,blockSize_);
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
                        A2Autils<FImpl>::MesonField(blockCache, &dvl[dtL][iblock+icache], &dvr[dtR][jblock+jcache], gamma_, phase, nd - 1, &timer);
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
                            // for(int iExt=0 ; iExt<nExt_ ; iExt++)
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int t=0 ; t<nt ; t++)
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                            {
                                block(iExt,iStr,t,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,t,iicache,jjcache);
                            }
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
