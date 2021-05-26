#ifndef Distil_matrix_hpp_
#define Distil_matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/TimerArray.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

#define DISTIL_MATRIX_NAME "DistilMesonField"

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

using TimeSliceMap = std::vector<std::vector<unsigned int>>;
using DilutionMap  = std::array<std::vector<std::vector<unsigned int>>,3>;

// metadata serialiser class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    unsigned int,               Nt,
                                    unsigned int,               Nvec,
                                    std::vector<RealF>,         Momentum,
                                    Gamma::Algebra,             Operator,               // just gamma matrices for now, but could turn into more general operators in the future
                                    std::vector<int>,           NoisePair,
                                    std::string,                MesonFieldType,
                                    std::vector<std::vector<std::string>>,   NoiseHashes)
};

auto timeslicesFn = [](const std::vector<unsigned int> ts)
{
    std::stringstream ss;
    for (auto& t : ts)
        ss << t << " ";
    return ss.str();
};

//metadata io class
class DistilMetadataIo
{
private:
    std::string filename_, metadataname_;
public:
    // constructor
    DistilMetadataIo(std::string filename, std::string metadataname):filename_(filename) , metadataname_(metadataname) {}
    //methods
    template <typename T>
    void write2dMetadata(const std::string name, const std::vector<std::vector<T>>& data)
    {
        // auxiliar variable-length struct (see hdf5 variable-length documentation)
        typedef struct  {
            size_t len; /* Length of VL data (in base type units) */      
            void *p;    /* Pointer to VL data */        
        } VlStorage;
#ifdef HAVE_HDF5

        Hdf5Reader  reader(filename_, false);
        push(reader, metadataname_);    //creates main h5 group
        H5NS::Group &subgroup = reader.getGroup();

        H5NS::VarLenType vl_type(Hdf5Type<T>::type());
        std::vector<VlStorage> vl_data(data.size());
        for(unsigned int i=0 ; i<data.size() ; i++)
        {
            vl_data.at(i).len = data.at(i).size();
            vl_data.at(i).p = (void*) &data.at(i).front();
        }

        hsize_t         attrDim = data.size();
        H5NS::DataSpace attrSpace(1, &attrDim);
        H5NS::Attribute attr = subgroup.createAttribute(name, vl_type, attrSpace);
        attr.write(vl_type, &vl_data.front());
#else
        HADRONS_ERROR(Implementation, "distillation I/O needs HDF5 library");
#endif
    }
};

//computation class declaration
template <typename FImpl, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef typename std::vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int)>  MetadataFn;
public:
    //todo: maybe reposition later
    long    global_counter = 0;
    double  global_flops   = 0.0;
    double  global_bytes   = 0.0;
    double  global_iospeed = 0.0;
private:
    const std::vector<std::string>      sides_ = {"left","right"};
    std::map<std::string,std::string>   dmfType_;
    GridCartesian*                      g_;
    GridCartesian*                      g3d_;
    ColourVectorField                   evec3d_;
    FermionField                        tmp3d_;
    Vector<Tio>                         bBuf_;
    Vector<T>                           cBuf_;
    const unsigned int                  bSize_;
    const unsigned int                  cSize_;
    const unsigned int                  nt_;
    const unsigned int                  nd_;
    const unsigned int                  next_;
    const unsigned int                  nstr_;
public:
    DmfComputation(std::map<std::string,std::string>    c,
                   GridCartesian*                       g,
                   GridCartesian*                       g3d,
                   const unsigned int                   blockSize,
                   const unsigned int                   cacheSize,
                   const unsigned int                   nt,
                   const unsigned int                   next,
                   const unsigned int                   nstr);
public:
    void execute(const FilenameFn                                   &filenameDmfFn,
                 const MetadataFn                                   &metadataDmfFn,
                 std::vector<Gamma::Algebra>                        gamma_,
                 std::map<std::string, DistilVector&>               dv,
                 std::map<std::string, DistillationNoise&>          n,
                 std::vector<int>                                   inoise,
                 std::vector<ComplexField>                          ph,
                 std::map<std::string, unsigned int>                dil_size_ls,
                 std::map<std::string, std::vector<unsigned int>>   timeDilSource,
                 DilutionMap                                        leftMap,
                 DilutionMap                                        rightMap,
                 TimerArray*                                        tarray);
    void distVec(std::map<std::string, DistilVector&>               dv,
                 std::map<std::string, DistillationNoise&>          n,
                 std::vector<int>                                   inoise,
                 LapPack&                                           epack,
                 std::map<std::string, std::vector<unsigned int>>   timeDilSource,
                 std::map<std::string, PerambTensor&>               peramb={});
private:
    void makePhiComponent(FermionField&         phiComponent,
                          DistillationNoise&    n,
                          const int             inoise,
                          const unsigned int    iD,
                          PerambTensor&         peramb,
                          LapPack&              epack);
    void makeRhoComponent(FermionField&         rhoComponent,
                          DistillationNoise&    n,
                          const int             inoise,
                          const unsigned int    iD);
};

//####################################
//# computation class implementation #
//####################################
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePhiComponent(FermionField&        phiComponent,
                   DistillationNoise&   n,
                   const int            inoise,
                   const unsigned int   iD,
                   PerambTensor&        peramb,
                   LapPack&             epack)
{
    std::array<unsigned int,3> c = n.dilutionCoordinates(iD);
    unsigned int dt = c[Index::t] , dl = c[Index::l] , ds = c[Index::s];
    std::vector<int> p_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(p_ts.begin(), p_ts.end(), dt);
    unsigned int idt = std::distance(p_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj
    const unsigned int nVec = epack.evec.size();
    const unsigned int Ntfirst = g_->LocalStarts()[nd_ - 1];
    const unsigned int Ntlocal = g_->LocalDimensions()[nd_ - 1];
    for (unsigned int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
    {
        tmp3d_ = Zero();
        for (unsigned int k = 0; k < nVec; k++)
        {
            ExtractSliceLocal(evec3d_,epack.evec[k],0,t-Ntfirst,nd_ - 1);
            tmp3d_ += evec3d_ * peramb.tensor(t, k, dl, inoise, idt, ds);
        }
        InsertSliceLocal(tmp3d_,phiComponent,0,t-Ntfirst,nd_ - 1);
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRhoComponent(FermionField&        rhoComponent,
                   DistillationNoise&   n,
                   const int            inoise,
                   const unsigned int   iD)   
{
    rhoComponent = n.makeSource(iD, inoise);
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::distVec(std::map<std::string, DistilVector&>              dv,
          std::map<std::string, DistillationNoise&>         n,
          std::vector<int>                                  inoise,
          LapPack&                                          epack,
          std::map<std::string, std::vector<unsigned int>>  timeDilSource,
          std::map<std::string, PerambTensor&>              peramb)
{
    std::map<std::string,int> iNoise = {{"left",inoise[0]},{"right",inoise[1]}};
    for(std::string s : sides_)    // computation
    for(unsigned int iD=0 ; iD<n.at(s).dilutionSize() ; iD++)  // computation of phi or rho
    {
        dv.at(s)[iD] = Zero();
        std::vector<unsigned int> tDilS = timeDilSource.at(s);
        unsigned int dt = n.at(s).dilutionCoordinates(iD)[Index::t];
        if(std::find(tDilS.begin(), tDilS.end(), dt) != tDilS.end()) // if time source is available, compute that block
        {
            if(dmfType_.at(s)=="phi")
            {
                makePhiComponent(dv.at(s)[iD] , n.at(s) , iNoise.at(s) , iD , peramb.at(s), epack);
            }
            else if(dmfType_.at(s)=="rho"){
                makeRhoComponent(dv.at(s)[iD] , n.at(s) , iNoise.at(s), iD);
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
DmfComputation<FImpl,T,Tio>
::DmfComputation(std::map<std::string,std::string>  c,
                 GridCartesian*                     g,
                 GridCartesian*                     g3d,
                 const unsigned int                 blockSize,
                 const unsigned int                 cacheSize,
                 const unsigned int                 nt,
                 const unsigned int                 next,
                 const unsigned int                 nstr)
: dmfType_(c), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d)
, nt_(nt) , nd_(g->Nd()), bSize_(blockSize) , cSize_(cacheSize)
, next_(next) , nstr_(nstr)
{
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::execute(const FilenameFn                                  &filenameDmfFn,
          const MetadataFn                                  &metadataDmfFn,
          std::vector<Gamma::Algebra>                       gamma_,
          std::map<std::string, DistilVector&>              dv,
          std::map<std::string, DistillationNoise&>         n,
          std::vector<int>                                  inoise,
          std::vector<ComplexField>                         ph,
          std::map<std::string, unsigned int>               dil_size_ls,
          std::map<std::string, std::vector<unsigned int>>  timeDilSource,
          DilutionMap                                       leftMap,
          DilutionMap                                       rightMap,
          TimerArray*                                       tarray)
{
    std::vector<std::vector<unsigned int>> timeDilutionPairList;
    bool fileIsInit = false;
    cBuf_.resize(next_*nstr_*nt_*cSize_*cSize_);
    
    
    const unsigned int vol = g_->_gsites;
    std::string dmfType = dmfType_.at("left") + " " + dmfType_.at("right");
    // computing mesonfield blocks and saving to disk
    for (unsigned int dtL : timeDilSource.at("left"))
    for (unsigned int dtR : timeDilSource.at("right"))
    {
        std::map<std::string,std::vector<unsigned int>> p = { {"left",{}} , {"right",{}}};
        for(auto s : sides_)
        {
            if(dmfType_.at(s)=="phi")
            {
                p.at(s).resize(nt_);
                std::iota(std::begin(p.at(s)), std::end(p.at(s)), 0); //phi: filling with all time slices <-> intersects with non-empty
            }
            else
            {
                p.at(s) =  n.at(s).timeSlices(s=="left" ? dtL : dtR);
            }
        }

        std::vector<unsigned int> stInter;
        std::set_intersection(p.at("left").begin(), p.at("left").end(), 
                              p.at("right").begin(), p.at("right").end(),
                              std::back_inserter(stInter));

        const int nt_sparse = stInter.size();
        bBuf_.resize(next_*nstr_*nt_sparse*bSize_*bSize_); //does Hadrons environment know about this?

        if( !stInter.empty() ) // only execute case when partitions have at least one time slice in common
        {
            timeDilutionPairList.push_back({dtL,dtR});
            LOG(Message) << "------------------------ " << dtL << "-" << dtR << " ------------------------" << std::endl; 
            LOG(Message) << "Saving time slices : " << timeslicesFn(stInter) << std::endl;
            LOG(Message) << "Time extension in file : " << nt_sparse << std::endl;
            std::string datasetName = std::to_string(dtL)+"-"+std::to_string(dtR);

            unsigned int nblocki = dil_size_ls.at("left")/bSize_ + (((dil_size_ls.at("left") % bSize_) != 0) ? 1 : 0);
            unsigned int nblockj = dil_size_ls.at("right")/bSize_ + (((dil_size_ls.at("right") % bSize_) != 0) ? 1 : 0);

            // loop over blocks in the current time-dilution block
            for(unsigned int iblock=0 ; iblock<dil_size_ls.at("left") ; iblock+=bSize_) //set according to memory size
            for(unsigned int jblock=0 ; jblock<dil_size_ls.at("right") ; jblock+=bSize_)
            {
                unsigned int iblockSize = MIN(dil_size_ls.at("left")-iblock,bSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                unsigned int jblockSize = MIN(dil_size_ls.at("right")-jblock,bSize_);
                A2AMatrixSet<Tio> block(bBuf_.data(), next_ , nstr_ , nt_sparse, iblockSize, jblockSize);

                LOG(Message) << "Distil matrix block " 
                << jblock/bSize_ + nblocki*iblock/bSize_ + 1 
                << "/" << nblocki*nblockj << " [" << iblock << " .. " 
                << iblock+iblockSize-1 << ", " << jblock << " .. " << jblock+jblockSize-1 << "]" 
                << std::endl;

                LOG(Message) << "Block size : "         << nt_sparse*iblockSize*jblockSize*sizeof(Tio)/1024. << "KB/momentum/gamma" << std::endl;
                LOG(Message) << "Cache block size : "   << DISTIL_NT_CHUNK_SIZE*cSize_*cSize_*sizeof(T) << "B/momentum/gamma" << std::endl;  //remember to change this in case I change chunk size from nt_ to something else

                double flops        = 0.0;
                double bytes        = 0.0;
                double time_kernel  = 0.0;
                double nodes        = g_->NodeCount();

                // loop over cache_ blocks in the current block
                for(unsigned int icache=0 ; icache<iblockSize ; icache+=cSize_)   //set according to cache_ size
                for(unsigned int jcache=0 ; jcache<jblockSize ; jcache+=cSize_)
                {
                    unsigned int icacheSize = MIN(iblockSize-icache,cSize_);      // icacheSize is the size of the current cache_ block (indexed by ii); N_ii-ii is the size of the remainder cache_ block
                    unsigned int jcacheSize = MIN(jblockSize-jcache,cSize_);
                    A2AMatrixSet<T> blockCache(cBuf_.data(), next_, nstr_, nt_, icacheSize, jcacheSize);

                    double timer = 0.0;
                    tarray->startTimer("kernel");
                    // assuming certain indexation here! (dt must be the slowest index for this to work; otherwise will have to compute l/r block at each contraction)
                    unsigned int iDl = n.at("left").dilutionIndex(dtL,0,0) , iDr = n.at("right").dilutionIndex(dtR,0,0);
                    A2Autils<FImpl>::MesonField(blockCache, &dv.at("left")[iDl+iblock+icache], &dv.at("right")[iDr+jblock+jcache], gamma_, ph, nd_ - 1, &timer);
                    tarray->stopTimer("kernel");
                    time_kernel += timer;

                    // nExt is currently # of momenta , nStr is # of gamma matrices
                    flops += vol*(2*8.0+6.0+8.0*next_)*icacheSize*jcacheSize*nstr_;
                    bytes += vol*(12.0*sizeof(T))*icacheSize*jcacheSize
                            +  vol*(2.0*sizeof(T)*next_)*icacheSize*jcacheSize*nstr_;

                    // loop through the cacheblock (inside them) and point blockCache to block
                    tarray->startTimer("cache copy");
                    unsigned int stSize = stInter.size();
                    thread_for_collapse(5,iext,next_,{
                    for(unsigned int istr=0;istr<nstr_;istr++)
                    for(unsigned int it=0;it<stSize;it++)
                    for(unsigned int iicache=0;iicache<icacheSize;iicache++)
                    for(unsigned int jjcache=0;jjcache<jcacheSize;jjcache++)
                    {
                        block(iext,istr,it,icache+iicache,jcache+jjcache)=blockCache(iext,istr,stInter[it],iicache,jjcache);
                    }
                    });
                    tarray->stopTimer("cache copy");
                }

                LOG(Message) << "Kernel perf (flops) " << flops/time_kernel/1.0e3/nodes 
                            << " Gflop/s/node " << std::endl;
                LOG(Message) << "Kernel perf (read) " << bytes/time_kernel*0.000931322574615478515625/nodes //  *1.0e6/1024/1024/1024/nodes
                            << " GB/s/node "  << std::endl;
                global_counter++;
                global_flops += flops/time_kernel/1.0e3/nodes ;
                global_bytes += bytes/time_kernel*0.000931322574615478515625/nodes ; // 1.0e6/1024/1024/1024/nodes

                // saving current block to disk
                double ioTime = (tarray!=nullptr) ? -tarray->getDTimer("IO: write block") : 0.0;
                tarray->startTimer("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
                //parallel io
                unsigned int inode = g_->ThisRank();
                unsigned int nnode = g_->RankCount(); 
                LOG(Message) << "Starting parallel IO. Rank count=" << nnode  << std::endl;
                g_->Barrier();
                for(unsigned int ies=inode ; ies<next_*nstr_ ; ies+=nnode){
                    unsigned int iext = ies/nstr_;
                    unsigned int istr = ies%nstr_;
                    
                    // metadata;
                    DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,inoise[0],inoise[1]);                    
                    A2AMatrixIo<HADRONS_DISTIL_IO_TYPE> matrixIo(filenameDmfFn(iext,istr,inoise[0],inoise[1]), DISTIL_MATRIX_NAME, nt_sparse, dil_size_ls.at("left"), dil_size_ls.at("right"));

                    tarray->startTimer("IO: write block");
                    if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                        if( (dtL==timeDilSource.at("left")[0]) && (dtR==timeDilSource.at("right")[0]) )     //execute this once per block
                        {
                            tarray->startTimer("IO: file creation");
                            matrixIo.initFile(md);
                            tarray->stopTimer("IO: file creation");
                        }
                        matrixIo.saveBlock(block, iext , istr , iblock, jblock, datasetName, stInter, cSize_);   //sets 2D chunk size and creates dataset
                    }
                    else{
                        matrixIo.saveBlock(block, iext , istr , iblock, jblock, datasetName);
                    }
                    tarray->stopTimer("IO: write block");
                }
                g_->Barrier();
#endif
                tarray->stopTimer("IO: total");
                ioTime    += tarray->getDTimer("IO: write block");
                unsigned int bytesBlockSize  = static_cast<double>(next_*nstr_*nt_sparse*iblockSize*jblockSize*sizeof(Tio));
                double iospeed = bytesBlockSize/ioTime*0.95367431640625;     // 1.0e6/1024/1024
                LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                << ioTime  << " us (" << iospeed << " MB/s)" << std::endl;
                global_iospeed += iospeed;
            }
        }
    }
    //saving dilution schemes (2d ragged metadata)
    if(g_->IsBoss())
    {
        tarray->startTimer("IO: total");
        for(unsigned int iext=0 ; iext<next_ ; iext++)
        for(unsigned int istr=0 ; istr<nstr_ ; istr++)
        {
            const std::string distilname = DISTIL_MATRIX_NAME;
            DistilMetadataIo mdIo(filenameDmfFn(iext,istr,inoise[0],inoise[1]), distilname+"/Metadata" );
            mdIo.write2dMetadata("TimeDilutionPairs", timeDilutionPairList);
            //schemes
            mdIo.write2dMetadata("TimeDilutionLeft"  , leftMap[Index::t] );
            mdIo.write2dMetadata("TimeDilutionRight" , rightMap[Index::t]);
            mdIo.write2dMetadata("LapDilutionLeft"   , leftMap[Index::l] );
            mdIo.write2dMetadata("LapDilutionRight"  , rightMap[Index::l]);
            mdIo.write2dMetadata("SpinDilutionLeft"  , leftMap[Index::s] );
            mdIo.write2dMetadata("SpinDilutionRight" , rightMap[Index::s]);
        }
        tarray->stopTimer("IO: total");
    }
}

//############################
// aux class declaration
//############################
template <typename FImpl>
class DmfHelper
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef typename std::vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
private:
    unsigned int nd_;
    std::map<std::string,std::string> dmfType_;
    const std::vector<std::string> sides = {"left","right"};
public:
    DmfHelper(const unsigned int nd, std::map<std::string,std::string> c);
    std::vector<std::vector<int>> parseNoisePairs(std::vector<std::string> inputN);
    void computePhase(std::vector<std::vector<RealF>> momenta, ComplexField &coor, std::vector<int> dim, std::vector<ComplexField> &phase);
    DilutionMap getMap(DistillationNoise & n);
};

//############################
//# helper class implementation #
//############################
template <typename FImpl>
DmfHelper<FImpl>::DmfHelper(const unsigned int nd, std::map<std::string,std::string>   c)
: dmfType_(c) , nd_(nd)
{
}

template <typename FImpl>
std::vector<std::vector<int>> DmfHelper<FImpl>::parseNoisePairs(std::vector<std::string> inputN)
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

template <typename FImpl>
void DmfHelper<FImpl>::computePhase(std::vector<std::vector<RealF>> momenta,
                                    ComplexField&                   coor,
                                    std::vector<int>                dim,
                                    std::vector<ComplexField>&      phase)
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

template <typename FImpl>
DilutionMap DmfHelper<FImpl>::getMap(DistillationNoise & n)
{
    DilutionMap m;
    for(auto &id : { Index::t, Index::l, Index::s })
    for(unsigned int it=0 ; it<n.dilutionSize(id) ; it++)
    {
        std::vector<unsigned int> temp = n.getDilutionPartition(id,it);
        m[id].push_back(temp);
    }
    return m;
}


END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Distil_matrix_hpp_
