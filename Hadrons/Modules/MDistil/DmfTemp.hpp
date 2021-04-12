#ifndef atwork_hpp
#define atwork_hpp

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

// template <typename T>
// using ObjArray_LR = std::array<T, 2>;

using TimeSliceMap = std::vector<std::vector<unsigned int>>; // this is here because TimeSliceMap is a return type in methods below

static void printMap(TimeSliceMap &m)
{
    for(auto& d : m)
    {
        for(auto& t : d)
        {
            std::cout << t << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//computation class declaration
template <typename FImpl, typename Field, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef typename std::vector<Field> DistilVector;
    typedef typename DistillationNoise::Index Index;
public:
    //todo: maybe reposition later
    long    global_counter = 0;
    double  global_flops   = 0.0;
    double  global_bytes   = 0.0;
private:
    const std::vector<std::string> sides_ = {"left","right"};
    std::map<std::string,std::string> dmfCase_;
    GridCartesian *     g_;
    GridCartesian *     g3d_;
    ColourVectorField   evec3d_;
    Field               tmp3d_;
    Vector<Tio>         bBuf_;
    Vector<T>           cBuf_;
    unsigned int        bSize_;
    unsigned int        cSize_;
    unsigned int        nt_;
    unsigned int        nd_;
public:
    DmfComputation( std::map<std::string,std::string> c,
                    GridCartesian * g, GridCartesian * g3d,
                    unsigned int blockSize,
                    unsigned int cacheSize,
                    unsigned int nt);
public:
    void execute(std::vector<A2AMatrixIo<Tio>> io_table,
                TimeSliceMap st,
                std::map<std::string, DistilVector&> dv,
                std::map<std::string, DistillationNoise&> n,
                std::vector<Gamma::Algebra> gamma,
                std::vector<ComplexField> ph,
                const unsigned int n_ext,
                const unsigned int n_str,
                std::map<std::string, int>  dil_size_ls,
                const unsigned int eff_nt,
                std::map<std::string, std::vector<int>>   timeDilSource,
                TimerArray * tarray
                );
    void distVec(std::map<std::string, DistilVector&> & dv,
                std::map<std::string, DistillationNoise&> n,
                const std::vector<int> inoise,
                std::map<std::string, PerambTensor&> & peramb,
                const LapEvecs            & epack,
                std::map<std::string, std::vector<int>>   timeDilSource
                );
};

// aux class declaration
template <typename FImpl, typename Field>
class DmfHelper
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef typename std::vector<Field> DistilVector;
    typedef typename DistillationNoise::Index Index;
private:
    int nt_;
    int nd_;
    std::map<std::string,std::string> dmfCase_;
    TimeSliceMap noiseTimeMapl_, noiseTimeMapr_;
    const std::vector<std::string> sides = {"left","right"};
    unsigned int effTime(TimeSliceMap tmap);
public:
    DmfHelper(DistillationNoise & nl, DistillationNoise & nr, std::map<std::string,std::string> c);
    int computeEffTimeDimension();
    std::vector<std::vector<int>> parseNoisePairs(std::vector<std::string> inputN);
    void computePhase(std::vector<std::vector<RealF>> momenta, ComplexField &coor, std::vector<int> dim, std::vector<ComplexField> &phase);
    TimeSliceMap timeSliceMap(DistillationNoise & n);
    void dumpMap(TimeSliceMap m);
};

// metadata class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma,
                                    std::vector<int>,   noise_pair,
                                    std::vector<int>,   leftTimeDilSources,
                                    std::vector<int>,   rightTimeDilSources,
                                    )
};

//####################################
//# computation class implementation #
//####################################

template <typename FImpl, typename Field, typename T, typename Tio>
void DmfComputation<FImpl,Field,T,Tio>
::distVec(std::map<std::string, DistilVector&>&     dv,
          std::map<std::string, DistillationNoise&> n,
          const std::vector<int>                    inoise,
          std::map<std::string, PerambTensor&>&     peramb,
          const LapEvecs&                           epack,
          std::map<std::string, std::vector<int>>   timeDilSource)
{
    const int nd = g_->Nd();
    const int nVec = epack.evec.size();
    const int Ntfirst = g_->LocalStarts()[nd - 1];
    const int Ntlocal = g_->LocalDimensions()[nd - 1];
    std::map<std::string,int> iNoise = {{"left",inoise[0]},{"right",inoise[1]}};

    // std::cout << "tDilSLeft=" << timeDilSource.at("left") << std::endl;
    // std::cout << "tDilSRight=" << timeDilSource.at("right") << std::endl;
    // std::cin.get();

    for(std::string s : sides_)    // computation
    for(int iD=0 ; iD<n.at(s).dilutionSize() ; iD++)  // computation of phi or rho
    {
        std::array<unsigned int,3> c = n.at(s).dilutionCoordinates(iD);
        unsigned int dt = c[0] , dl = c[1] , ds = c[2];
        dv.at(s)[iD] = Zero();
        std::vector<int> tDilS = timeDilSource.at(s);
        if(std::find(tDilS.begin(), tDilS.end(), dt) != tDilS.end()) // if time source is available, compute that block
        {
            // std::cout << "dt=" << dt << std::endl;
            // std::cin.get();
            if(dmfCase_.at(s)=="phi")
            {
                for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)   //loop over (local) timeslices
                {
                    tmp3d_ = Zero();
                    for (int k = 0; k < nVec; k++)
                    {
                        ExtractSliceLocal(evec3d_,epack.evec[k],0,t-Ntfirst,nd - 1);
                        // std::cout << t<< " " <<  k<< " " <<  dl<< " " <<  iNoise.at(s)<< " " <<  dt<< " " <<  ds << std::endl;
                        tmp3d_ += evec3d_ * peramb.at(s).tensor(t, k, dl, iNoise.at(s), dt, ds);
                    }
                    InsertSliceLocal(tmp3d_,dv.at(s)[iD],0,t-Ntfirst,nd - 1);
                }
            }
            else if(dmfCase_.at(s)=="rho"){
                dv.at(s)[iD] = n.at(s).makeSource(iD, iNoise.at(s));
            }
        }
    }
}

template <typename FImpl, typename Field, typename T, typename Tio>
DmfComputation<FImpl,Field,T,Tio>
::DmfComputation(std::map<std::string,std::string> c,
                 GridCartesian * g, GridCartesian * g3d,
                 unsigned int blockSize,
                 unsigned int cacheSize,
                 unsigned int nt)
: dmfCase_(c), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d)
, nt_(nt), nd_(g->Nd()), bSize_(blockSize) , cSize_(cacheSize)
{
}

template <typename FImpl, typename Field, typename T, typename Tio>
void DmfComputation<FImpl,Field,T,Tio>
::execute(std::vector<A2AMatrixIo<Tio>> io_table,
          TimeSliceMap st,
          std::map<std::string, DistilVector&> dv,
          std::map<std::string, DistillationNoise&> n,
          std::vector<Gamma::Algebra> gamma,
          std::vector<ComplexField> ph,
          const unsigned int n_ext,
          const unsigned int n_str,
          std::map<std::string, int>  dil_size_ls,
          const unsigned int eff_nt,
          std::map<std::string, std::vector<int>>   timeDilSource,
          TimerArray * tarray)
{
    bBuf_.resize(n_ext*n_str*eff_nt*bSize_*bSize_);
    cBuf_.resize(n_ext*n_str*nt_*cSize_*cSize_);
    
    const unsigned int vol = g_->_gsites;
    std::string dmfcase = dmfCase_.at("left") + " " + dmfCase_.at("right");
    // computing mesonfield blocks and saving to disk
    for (auto& dtL : timeDilSource.at("left"))
    for (auto& dtR : timeDilSource.at("right"))
    {
        std::vector<unsigned int> pL;
        std::vector<unsigned int> pR;
        std::vector<unsigned int> stInter;
        if(dmfCase_.at("left")=="phi")
        {
            pL.resize(nt_);
            std::iota(std::begin(pL), std::end(pL), 0); //phi: filling with all time slices <-> intersects with non-empty
        }
        else
        {
            pL =  n.at("left").timeSlices(dtL);
        }
        if(dmfCase_.at("right")=="phi")
        {
            pR.resize(nt_);
            std::iota(std::begin(pR), std::end(pR), 0); //phi: filling with all time slices <-> intersects with non-empty
        }
        else
        {
            pR = n.at("right").timeSlices(dtR);
        }

        std::set_intersection(pL.begin(), pL.end(), 
                        pR.begin(), pR.end(),
                        std::back_inserter(stInter));
        // std::cout << "pR=" << pR << " pL" << pL << " stInter=" << stInter << std::endl;
        // if(!stInter.empty())

        if( !stInter.empty() ) // only execute rho rho case when partitions have at least one time slice in common
        {
            LOG(Message) << "################### dtL_" << dtL << " dtR_" << dtR << " ################### " << std::endl; 
            LOG(Message) << "At least one rho found. Time slices to be saved=" << stInter << "..." << std::endl;
            std::string datasetName = "dtL"+std::to_string(dtL)+"_dtR"+std::to_string(dtR);
            // LOG(Message) << "Computing dilution dataset " << datasetName << "..." << std::endl;

            int nblocki = dil_size_ls.at("left")/bSize_ + (((dil_size_ls.at("left") % bSize_) != 0) ? 1 : 0);
            int nblockj = dil_size_ls.at("right")/bSize_ + (((dil_size_ls.at("right") % bSize_) != 0) ? 1 : 0);

            // loop over blocls in the current time-dilution block
            for(int iblock=0 ; iblock<dil_size_ls.at("left") ; iblock+=bSize_) //set according to memory size
            for(int jblock=0 ; jblock<dil_size_ls.at("right") ; jblock+=bSize_)
            {
                int iblockSize = MIN(dil_size_ls.at("left")-iblock,bSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                int jblockSize = MIN(dil_size_ls.at("right")-jblock,bSize_);
                A2AMatrixSet<Tio> block(bBuf_.data(), n_ext , n_str , eff_nt, iblockSize, jblockSize);

                LOG(Message) << "Distil matrix block " 
                << jblock/bSize_ + nblocki*iblock/bSize_ + 1 
                << "/" << nblocki*nblockj << " [" << iblock << " .. " 
                << iblock+iblockSize-1 << ", " << jblock << " .. " << jblock+jblockSize-1 << "]" 
                << std::endl;

                LOG(Message) << "Block size = "         << eff_nt*iblockSize*jblockSize*sizeof(Tio) << "MB/momentum/gamma" << std::endl;
                LOG(Message) << "Cache block size = "   << nt_*cSize_*cSize_*sizeof(T) << "MB/momentum/gamma" << std::endl;  //remember to change this in case I change chunk size from nt_ to something else

                double flops        = 0.0;
                double bytes        = 0.0;
                double time_kernel  = 0.0;
                double nodes        = g_->NodeCount();

                // loop over cache_ blocks in the current block
                for(int icache=0 ; icache<iblockSize ; icache+=cSize_)   //set according to cache_ size
                for(int jcache=0 ; jcache<jblockSize ; jcache+=cSize_)
                {
                    int icacheSize = MIN(iblockSize-icache,cSize_);      // icacheSize is the size of the current cache_ block (indexed by ii); N_ii-ii is the size of the remainder cache_ block
                    int jcacheSize = MIN(jblockSize-jcache,cSize_);
                    A2AMatrixSet<T> blockCache(cBuf_.data(), n_ext, n_str, nt_, icacheSize, jcacheSize);

                    double timer = 0.0;
                    tarray->startTimer("kernel");
                    // assuming certain indexation here! (dt must be the slowest index for this to work; otherwise will have to compute l/r block at each contraction)
                    unsigned int iDl = n.at("left").dilutionIndex(dtL,0,0) , iDr = n.at("right").dilutionIndex(dtR,0,0);
                    A2Autils<FImpl>::MesonField(blockCache, &dv.at("left")[iDl+iblock+icache], &dv.at("right")[iDr+jblock+jcache], gamma, ph, nd_ - 1, &timer);
                    tarray->stopTimer("kernel");
                    time_kernel += timer;

                    // nExt is currently # of momenta , nStr is # of gamma matrices
                    flops += vol*(2*8.0+6.0+8.0*n_ext)*icacheSize*jcacheSize*n_str;
                    bytes += vol*(12.0*sizeof(T))*icacheSize*jcacheSize
                            +  vol*(2.0*sizeof(T)*n_ext)*icacheSize*jcacheSize*n_str;

                    // loop through the cacheblock (inside them) and point blockCache to block
                    tarray->startTimer("cache copy");
                    if(dmfcase=="phi phi")
                    {
                        thread_for_collapse( 5, iExt ,n_ext,{
                        // for(int iExt=0 ; iExt<n_ext ; iExt++)
                        for(int iStr=0 ;iStr<n_str ; iStr++)
                        for(int t=0 ; t<nt_ ; t++)
                        for(int iicache=0 ; iicache<icacheSize ; iicache++)
                        for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                        {
                            block(iExt,iStr,t,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,t,iicache,jjcache);
                        }
                        });
                    }
                    else
                    {
                        thread_for_collapse( 5, iExt ,n_ext,{
                        for(int iStr=0 ;iStr<n_str ; iStr++)
                        for(int it=0 ; it<stInter.size() ; it++)  //only wish to copy non-zero timeslices to block; guaranteed to fit because of neff = sup (partition sizes)
                        for(int iicache=0 ; iicache<icacheSize ; iicache++)
                        for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                            block(iExt,iStr,it,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,stInter[it],iicache,jjcache);
                        });
                    }
                    tarray->stopTimer("cache copy");
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
                tarray->startTimer("IO: total");
                tarray->startTimer("IO: write block");
                double ioTime = -tarray->getDTimer("IO: write block");
    #ifdef HADRONS_A2AM_PARALLEL_IO
                //parallel io
                int inode = g_->ThisRank();
                int nnode = g_->RankCount(); 
                LOG(Message) << "Starting parallel IO. Rank count=" << nnode  << std::endl;
                g_->Barrier();
                for(int ies=inode ; ies<n_ext*n_str ; ies+=nnode){
                    int iExt = ies/n_str;
                    int iStr = ies%n_str;
                    if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                        io_table[iStr + n_str*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName, cSize_);   //set surface chunk size as cSize_ (the chunk itself is 3D)
                    }
                    else{
                        io_table[iStr + n_str*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName);
                    }
                }
                g_->Barrier();
    #else
                // serial io, can remove later
                LOG(Message) << "Starting serial IO" << std::endl;
                for(int iExt=0; iExt<n_ext; iExt++)
                for(int iStr=0; iStr<n_str; iStr++)
                {
                    if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                        matrixIoTable[iStr + n_str*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName, cSize_);   //set surface chunk size as cSize_ (the chunk itself is 3D)
                    }
                    else{
                        matrixIoTable[iStr + n_str*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName);
                    }
                }
    #endif
                tarray->stopTimer("IO: total");
                tarray->stopTimer("IO: write block");
                ioTime    += tarray->getDTimer("IO: write block");
                int bytesBlockSize  = static_cast<double>(n_ext*n_str*eff_nt*iblockSize*jblockSize*sizeof(Tio));
                LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                << ioTime  << " us (" 
                                << bytesBlockSize/ioTime*0.95367431640625 // 1.0e6/1024/1024
                                << " MB/s)" << std::endl;
            }
        }
    }
}

//############################
//# helper class implementation #
//############################
template <typename FImpl, typename Field>
DmfHelper<FImpl,Field>::DmfHelper(DistillationNoise & nl, DistillationNoise & nr, std::map<std::string,std::string> c)
: noiseTimeMapl_( timeSliceMap(nl) ) , noiseTimeMapr_( timeSliceMap(nr) ) , dmfCase_(c)
{
    nt_ = nr.getNt();
    nd_ = nr.getGrid()->Nd();
}

template <typename FImpl, typename Field>
unsigned int DmfHelper<FImpl,Field>::effTime(TimeSliceMap tmap)
{
    // assuming it's a rho
    // compute eff_nt (1<=eff_nt<=nt_), the time extensin in the final object when there's at least one rho involved
    int eff_nt = 1;
    for(auto &e : tmap)
        e.size() > eff_nt ? eff_nt = e.size() : 0;      //get the highest possible eff_nt from st
    return eff_nt;
}

template <typename FImpl, typename Field>
int DmfHelper<FImpl,Field>::computeEffTimeDimension()
{
    if(dmfCase_.at("left")=="rho" || dmfCase_.at("right")=="rho")
    {
        unsigned int left_neff = (dmfCase_.at("left")=="phi") ? 0  : effTime(noiseTimeMapl_);
        unsigned int right_neff = (dmfCase_.at("right")=="phi") ? 0 : effTime(noiseTimeMapr_);
        std::cout << left_neff << " " << right_neff << std::endl;
        return MAX(left_neff,right_neff);    // if it's from phi (=0) it will be ignored by MAX
    }
    else
    {
        return nt_;
    }
}

template <typename FImpl, typename Field>
std::vector<std::vector<int>> DmfHelper<FImpl,Field>::parseNoisePairs(std::vector<std::string> inputN)
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

template <typename FImpl, typename Field>
void DmfHelper<FImpl,Field>::computePhase(std::vector<std::vector<RealF>> momenta, ComplexField &coor, std::vector<int> dim, std::vector<ComplexField> &phase)
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

template <typename FImpl, typename Field>
TimeSliceMap DmfHelper<FImpl,Field>::timeSliceMap(DistillationNoise & n)
{
    TimeSliceMap m;
    for(unsigned int it=0 ; it<n.dilutionSize(Index::t) ; it++)
    {
        std::vector<unsigned int> temp = n.timeSlices(it);
        m.push_back(temp);
    }
    return m;
}

template <typename FImpl, typename Field>
void DmfHelper<FImpl,Field>::dumpMap(TimeSliceMap m)
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

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif