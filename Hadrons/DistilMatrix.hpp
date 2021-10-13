#ifndef Distil_matrix_hpp_
#define Distil_matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

// make this an input?
#ifndef DISTILVECTOR_TIME_BATCH_SIZE
#define DISTILVECTOR_TIME_BATCH_SIZE 1
#endif

#define DISTIL_MATRIX_NAME      "DistilMesonField"
#define METADATA_NAME           "Metadata"
#define DILUTION_METADATA_NAME  "DilutionSchemes"

#define START_TIMER(name) if (tarray) tarray->startTimer(name)
#define STOP_TIMER(name)  if (tarray) tarray->stopTimer(name)
#define GET_TIMER(name)   ((tarray != nullptr) ? tarray->getDTimer(name) : 0.)

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

// Dimensions:
//   0 - ext - external field (momentum phase, etc...)
//   1 - str - spin-color structure or covariant derivative operators
//   2 - i   - left  A2A mode index
//   3 - j   - right A2A mode index
template <typename T>
using DistilMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 4, Eigen::RowMajor>>;

using DilutionMap  = std::array<std::vector<std::vector<unsigned int>>,3>;
enum Side {left = 0, right = 1};
const std::vector<Side> sides =  {Side::left,Side::right};

// metadata serialiser class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    unsigned int,               Nt,
                                    unsigned int,               Nvec,
                                    std::vector<RealF>,         Momentum,
                                    Gamma::Algebra,             Operator,               // can turn into more general operators in the future
                                    std::vector<unsigned int>,  NoisePair,
                                    std::string,                MesonFieldType,
                                    std::vector<std::string>,   NoiseHashLeft,
                                    std::vector<std::string>,   NoiseHashRight)
};

//metadata io class to deal with 2d ragged arrays (possibly remove after Mike's serialisation revision)
class DistilMetadataIo
{
private:
    std::string fileName_, metadataName_;
public:
    DistilMetadataIo(std::string filename, std::string metadataname):fileName_(filename) , metadataName_(metadataname) {}
    template <typename T>
    void write2dMetadata(const std::string name, const std::vector<std::vector<T>>& data, const std::string newgroupname="")
    {
#ifdef HAVE_HDF5
        // variable-length struct
        typedef struct  {
            size_t len; /* Length of VL data (in base type units) */      
            void *p;    /* Pointer to VL data */        
        } VlStorage;
        
        Hdf5Reader  reader(fileName_, false);
        push(reader, metadataName_);
        H5NS::Group &subgroup = reader.getGroup();

        if(!newgroupname.empty())
        {
            H5NS::Exception::dontPrint();
            try{
                subgroup = subgroup.openGroup(newgroupname);
            } catch (...) {
                subgroup = subgroup.createGroup(newgroupname);
            }
        }

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

/******************************************************************************
 *                  Class to handle Distil matrix block HDF5 I/O                 *
 ******************************************************************************/

// test purposes
template <typename FImpl>
class TestMd: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestMd,
                                    std::string,    name)
};


template <typename T>
class DistilMatrixIo
{
public:
    // constructors
    DistilMatrixIo(void) = default;
    DistilMatrixIo(std::string filename, std::string dataname, const unsigned int nt=0, const unsigned int ni = 0,
                const unsigned int nj = 0);
    // destructor
    ~DistilMatrixIo(void) = default;
    // access
    unsigned int getNi(void) const;
    unsigned int getNj(void) const;
    unsigned int getNt(void) const;
    size_t       getSize(void) const;

    // block I/O
    template <typename MetadataType>
    void initFile(const MetadataType &d);

    void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                   const unsigned int blockSizei, const unsigned int blockSizej, std::string t_name, std::string datasetName, const unsigned int chunkSize);

    void saveBlock(const DistilMatrixSet<T> &m, const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j, std::string datasetName,
                               const unsigned int t, const unsigned int chunkSize);

    // template <template <class> class Vec, typename VecT>
    // void load(Vec<VecT> &v, double *tRead = nullptr, GridBase *grid = nullptr, std::string datasetName="");
private:
    std::string  filename_{""}, dataname_{""};
    unsigned int nt_{0}, ni_{0}, nj_{0};
};

// constructor /////////////////////////////////////////////////////////////////
template <typename T>
DistilMatrixIo<T>::DistilMatrixIo(std::string filename, std::string dataname, 
                            const unsigned int nt, const unsigned int ni,
                            const unsigned int nj)
: filename_(filename), dataname_(dataname), ni_(ni), nj_(nj), nt_(nt)
{}

template <typename T>
template <typename MetadataType>
void DistilMatrixIo<T>::initFile(const MetadataType &d)
{
#ifdef HAVE_HDF5
    Hdf5Writer writer(filename_);
    push(writer, dataname_);    //creates main h5 group
    write(writer, "Metadata", d);
#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void DistilMatrixIo<T>::saveBlock(const T *data, 
                               const unsigned int i, 
                               const unsigned int j,
                               const unsigned int blockSizei,
                               const unsigned int blockSizej,
                               std::string t_name,
                               std::string datasetName,
                               const unsigned int chunkSize)
{
#ifdef HAVE_HDF5
    
    H5NS::H5File file(filename_,H5F_ACC_RDWR);
    H5NS::Group rootgroup = file.openGroup(DISTIL_MATRIX_NAME);

    // create t group if it doesnt exist, or open it
    H5NS::Group tgroup;
    // if ( H5Lexists( file.getId(), t_name, H5P_DEFAULT ) > 0 )
    // {
    //     tgroup = rootgroup.openGroup(t_name);
    //     // grp=h5file.openGroup("A");
    // }
    // else
    // {
    //     tgroup = rootgroup.createGroup(t_name);
    //     // grp=h5file.createGroup("A");
    // }

    H5NS::Exception::dontPrint();
    try{
        tgroup = rootgroup.openGroup(t_name);
    } catch (...) {                         //WHICH EXCEPTION SHOULD IT CATCH??
        tgroup = rootgroup.createGroup(t_name);
    }

    H5NS::DataSet        dataset;
    //create zeroed dataset T1,T2, or open it
    H5NS::Exception::dontPrint();
    try{
        dataset = tgroup.openDataSet(datasetName);
    } catch (...) {                         //WHICH EXCEPTION SHOULD IT CATCH??
        H5NS::DataSpace      dataspace;
        std::vector<hsize_t>    dim = {static_cast<hsize_t>(ni_), 
                                    static_cast<hsize_t>(nj_)},
                                chunk = {static_cast<hsize_t>(chunkSize), 
                                    static_cast<hsize_t>(chunkSize)};
        dataspace.setExtentSimple(dim.size(), dim.data());
        H5NS::DSetCreatPropList     plist;
        plist.setChunk(chunk.size(), chunk.data());
        plist.setFletcher32();
        dataset = tgroup.createDataSet(datasetName, Hdf5Type<T>::type(), dataspace, plist);
    }


    std::vector<hsize_t> count = {blockSizei, blockSizej},
                         offset = {static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1},
                         block  = {1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);



#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void DistilMatrixIo<T>::saveBlock(const DistilMatrixSet<T> &m,
                               const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j, std::string datasetName,
                               const unsigned int t, const unsigned int chunkSize)
{
    unsigned int blockSizei = m.dimension(2);
    unsigned int blockSizej = m.dimension(3);
    unsigned int nstr       = m.dimension(1);
    size_t       offset     = ((ext*nstr + str)*nt_ + t)*blockSizei*blockSizej;

    std::string t_name = std::to_string(t);

    saveBlock(m.data() + offset, i, j, blockSizei, blockSizej, t_name, datasetName, chunkSize);
}

/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

//computation class declaration
template <typename FImpl, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef std::vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int)>  MetadataFn;
public:
    long    blockCounter_ = 0;
    double  blockFlops_ = 0.0, blockBytes_ = 0.0, blockIoSpeed_ = 0.0;
private:
    std::map<Side,std::string>          dmfType_;
    GridCartesian*                      g_;
    GridCartesian*                      g3d_;
    ColourVectorField                   evec3d_;
    FermionField                        tmp3d_;
    FermionField                        tmp4d_;
    std::vector<Tio>                    bBuf_;
    std::vector<Tio>                    bufPinnedT_;
    std::vector<T>                      cBuf_;
    const unsigned int                  blockSize_; //eventually turns into io chunk size
    const unsigned int                  cacheSize_;
    const unsigned int                  nt_;
    const unsigned int                  nd_;
    const unsigned int                  nExt_;
    const unsigned int                  nStr_;
    const bool                          isExact_;
    const bool                          onlyDiag_;
    bool                                isInitFile_=false;
    std::map<Side, unsigned int>        dilSizeLS_;
    std::map<Side, DistillationNoise&>  distilNoise_;
    const unsigned int                  dvBatchSize_ = DISTILVECTOR_TIME_BATCH_SIZE;
public:
    DmfComputation(std::map<Side,std::string>   mf_type,
                   GridCartesian*               g,
                   GridCartesian*               g3d,
                   DistillationNoise&           nl,
                   DistillationNoise&           nr,
                   const unsigned int           block_size,
                   const unsigned int           cache_size,
                   const unsigned int           nt,
                   const unsigned int           n_ext,
                   const unsigned int           n_str,
                   const bool                   is_exact,
                   const bool                   only_diag);
    bool isPhi(Side s);
    bool isRho(Side s);
    DilutionMap getMap(Side s);
private:
    void makePhiComponent(FermionField&         phi_component,
                          DistillationNoise&    n,
                          const unsigned int    n_idx,
                          const unsigned int    D,
                          PerambTensor&         peramb,
                          LapPack&              epack);
    void makeRhoComponent(FermionField&         rho_component,
                          DistillationNoise&    n,
                          const unsigned int    n_idx,
                          const unsigned int    D);
    void makeDvLapSpinBlock(std::map<Side, DistilVector&>               dv,
                                      std::map<Side, unsigned int>      n_idx,
                                      LapPack&                          epack,
                                      Side                              s,
                                      unsigned int                      dt,
                                      unsigned int                      iibatch,
                                      std::map<Side, PerambTensor&>     peramb={});
    void makeDvLapSpinBatch(std::map<Side, DistilVector&>               dv,
                                      std::map<Side, unsigned int>      n_idx,
                                      LapPack&                          epack,
                                      Side                              s,
                                      std::vector<unsigned int>         dt_list,
                                      std::map<Side, PerambTensor&>     peramb);
    std::vector<unsigned int> fetchDvBatchIdxs(unsigned int               ibatch,
                                               std::vector<unsigned int>  time_dil_sources);
    void makePinnedPhiComponent(FermionField&                phi_component,
                            DistillationNoise&               n,
                            const unsigned int               n_idx,
                            const unsigned int               D,
                            const unsigned int               delta_t,
                            PerambTensor&                    peramb,
                            LapPack&                         epack);
    void makePinnedRhoComponent(FermionField&            rho_component,
                            DistillationNoise&           n,
                            const unsigned int           n_idx,
                            const unsigned int           t,
                            const unsigned int           D,
                            LapPack&                     epack);

    void makePinnedDvLapSpinBlock(std::map<Side, DistilVector&>     dv,
                               std::vector<unsigned int>            dt_list,
                               std::map<Side, unsigned int>         n_idx,
                               LapPack&                             epack,
                               Side                                 s,
                               const unsigned int                   delta_t,
                               std::map<Side, PerambTensor&>        peramb);
public:
    void execute(const FilenameFn                               &filenameDmfFn,
                 const MetadataFn                               &metadataDmfFn,
                 std::vector<Gamma::Algebra>                    gamma,
                 std::map<Side, DistilVector&>                  dv,
                 std::map<Side, unsigned int>                   n_idx,
                 std::vector<ComplexField>                      ph,
                 std::map<Side, std::vector<unsigned int>>      time_dil_source,
                 LapPack&                                       epack,
                 TimerArray*                                    tarray,
                 std::map<Side, PerambTensor&>                  peramb={});
};

//####################################
//# computation class implementation #
//####################################

template <typename FImpl, typename T, typename Tio>
DmfComputation<FImpl,T,Tio>
::DmfComputation(std::map<Side,std::string>     mf_type,
                 GridCartesian*                 g,
                 GridCartesian*                 g3d,
                 DistillationNoise&             nl,
                 DistillationNoise&             nr,
                 const unsigned int             block_size,
                 const unsigned int             cache_size,
                 const unsigned int             nt,
                 const unsigned int             n_ext,
                 const unsigned int             n_str,
                 const bool                     is_exact,
                 const bool                     only_diag)
: dmfType_(mf_type), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d), tmp4d_(g)
, nt_(nt) , nd_(g->Nd()), blockSize_(block_size) , cacheSize_(cache_size)
, nExt_(n_ext) , nStr_(n_str) , isExact_(is_exact) , onlyDiag_(only_diag)
{
    cBuf_.resize(nExt_*nStr_*nt_*cacheSize_*cacheSize_);
    bBuf_.resize(nExt_*nStr_*nt_*blockSize_*blockSize_); //maximum size

    const int nt_sparse = 1; //full time dilution , todo: adapt to dilution
    bufPinnedT_.resize(nExt_*nStr_*nt_sparse*blockSize_*blockSize_);

    distilNoise_ = std::map<Side, DistillationNoise&> ({{Side::left,nl},{Side::right,nr}});
    dilSizeLS_ = { {Side::left,nl.dilutionSize(Index::l)*nl.dilutionSize(Index::s)} ,
                             {Side::right,nr.dilutionSize(Index::l)*nr.dilutionSize(Index::s)} };

}

template <typename FImpl, typename T, typename Tio>
DilutionMap DmfComputation<FImpl,T,Tio>::getMap(Side s)
{
    DilutionMap m;
    for(auto dil_idx : { Index::t, Index::l, Index::s })
    for(unsigned int it=0 ; it<distilNoise_.at(s).dilutionSize(dil_idx) ; it++)
    {
        std::vector<unsigned int> temp = distilNoise_.at(s).dilutionPartition(dil_idx,it);
        m[dil_idx].push_back(temp);
    }
    return m;
}

template <typename FImpl, typename T, typename Tio>
bool DmfComputation<FImpl,T,Tio>::isPhi(Side s)
{
    return (dmfType_.at(s)=="phi" ? true : false);
}

template <typename FImpl, typename T, typename Tio>
bool DmfComputation<FImpl,T,Tio>::isRho(Side s)
{
    return (dmfType_.at(s)=="rho" ? true : false);
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePhiComponent(FermionField&            phi_component,
                   DistillationNoise&       n,
                   const unsigned int       n_idx,
                   const unsigned int       D,
                   PerambTensor&            peramb,
                   LapPack&                 epack)
{
    std::array<unsigned int,3> d_coor = n.dilutionCoordinates(D);
    unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
    std::vector<int> peramb_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(peramb_ts.begin(), peramb_ts.end(), dt);
    unsigned int idt = std::distance(peramb_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj 
    const unsigned int nVec = epack.evec.size();
    const unsigned int Nt_first = g_->LocalStarts()[nd_ - 1];
    const unsigned int Nt_local = g_->LocalDimensions()[nd_ - 1];
    for (unsigned int t = Nt_first; t < Nt_first + Nt_local; t++)
    {
        tmp3d_ = Zero();
        for (unsigned int k = 0; k < nVec; k++)
        {
            ExtractSliceLocal(evec3d_,epack.evec[k],0,t-Nt_first,nd_ - 1);
            tmp3d_ += evec3d_ * peramb.tensor(t, k, dk, n_idx, idt, ds);
        }
        InsertSliceLocal(tmp3d_,phi_component,0,t-Nt_first,nd_ - 1);
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRhoComponent(FermionField&        rho_component,
                   DistillationNoise&   n,
                   const unsigned int   n_idx,
                   const unsigned int   D)   
{
    rho_component = n.makeSource(D, n_idx);
}

// lap-spin blocks have fixed dimensions of (lap-spin dilution size left)x(lap-spin dilution size right)
// and each is identified by the starting posision in time dilution space, (T1,dtR)
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinBlock(std::map<Side, DistilVector&>                  dv,
                               std::map<Side, unsigned int>         n_idx,
                               LapPack&                             epack,
                               Side                                 s,
                               unsigned int                         dt,
                               unsigned int                         iibatch,
                               std::map<Side, PerambTensor&>        peramb)
{
    unsigned int D_offset = distilNoise_.at(s).dilutionIndex(dt,0,0);    // t is the slowest index
    unsigned int iD_offset = iibatch*dilSizeLS_.at(s);
    for(unsigned int iD=iD_offset ; iD<iD_offset+dilSizeLS_.at(s) ; iD++)
    {
        unsigned int D = iD + D_offset - iD_offset;
        if(isPhi(s))
        {
            makePhiComponent(dv.at(s)[iD] , distilNoise_.at(s) , n_idx.at(s) , D , peramb.at(s), epack);
        }
        else if(isRho(s))
        {
            makeRhoComponent(dv.at(s)[iD] , distilNoise_.at(s) , n_idx.at(s) , D);
        }
    }
}

// a batch is composed of several lap-spin blocks
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinBatch(std::map<Side, DistilVector&>              dv,
                               std::map<Side, unsigned int>     n_idx,
                               LapPack&                         epack,
                               Side                             s,
                               std::vector<unsigned int>        dt_list,
                               std::map<Side, PerambTensor&>    peramb)
{
    for(unsigned int idt=0 ; idt<dt_list.size() ; idt++)
    {
        makeDvLapSpinBlock(dv,n_idx,epack,s,dt_list[idt],idt,peramb);
    }
}

// fetch time dilution indices (sources) in dv batch ibatch
template <typename FImpl, typename T, typename Tio>
std::vector<unsigned int> DmfComputation<FImpl,T,Tio>
::fetchDvBatchIdxs(unsigned int ibatch, std::vector<unsigned int> time_dil_sources)
{
    std::vector<unsigned int> batch_dt;
    for(unsigned int dt=ibatch*dvBatchSize_; dt<(ibatch+1)*dvBatchSize_; dt++){
        batch_dt.push_back(time_dil_sources[dt]);
    }
    return batch_dt;
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePinnedPhiComponent(FermionField&              phi_component,
                   DistillationNoise&               n,
                   const unsigned int               n_idx,
                   const unsigned int               D,
                   const unsigned int               delta_t,
                   PerambTensor&                    peramb,
                   LapPack&                         epack)
{
    std::array<unsigned int,3> d_coor = n.dilutionCoordinates(D);
    unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
    
    //compute at pin_time, insert at t
    const unsigned int pin_time = (dt + delta_t)%nt_;     //TODO: generalise to dilution
    const unsigned int t        = pin_time;               //TODO: generalise to dilution

    std::vector<int> peramb_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(peramb_ts.begin(), peramb_ts.end(), dt);
    unsigned int idt_peramb = std::distance(peramb_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj 
    const unsigned int nVec = epack.evec.size();
    const unsigned int Nt_first = g_->LocalStarts()[nd_ - 1];
    const unsigned int Nt_local = g_->LocalDimensions()[nd_ - 1];
    if( (pin_time>=Nt_first) && (pin_time<Nt_first+Nt_local) )
    {
        tmp3d_ = Zero();
        for (unsigned int k = 0; k < nVec; k++)
        {
            ExtractSliceLocal(evec3d_,epack.evec[k],0,pin_time-Nt_first,nd_ - 1);
            tmp3d_ += evec3d_ * peramb.tensor(pin_time, k, dk, n_idx, idt_peramb, ds);
        }
        InsertSliceLocal(tmp3d_,phi_component,0,t-Nt_first,nd_ - 1);        //shift this byt delta_t
        // std::cout << phi_component << std::endl;
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePinnedRhoComponent(FermionField&          rho_component,
                   DistillationNoise&           n,
                   const unsigned int           n_idx,
                   const unsigned int           D,
                   const unsigned int           delta_t,
                   LapPack&                     epack)
{
    // abstract this to makePinnedSource() kind of method in DilutedNoise?
    std::array<unsigned int,3> d_coor = n.dilutionCoordinates(D);
    unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];

    typename DistillationNoise::NoiseType   noise(n.getNoise()[n_idx].data(), nt_, epack.eval.size(), Ns);

    //compute at pin_time, insert at t
    const unsigned int pin_time = (dt + delta_t)%nt_;     //TODO: generalise to dilution
    const unsigned int t        = dt;               //TODO: generalise to dilution

    // adapt to dilution! (add an it loop and untrivialise ik and is loops)
    
    const unsigned int Nt_first = g_->LocalStarts()[nd_ - 1];
    const unsigned int Nt_local = g_->LocalDimensions()[nd_ - 1];
    if( (pin_time>=Nt_first) && (pin_time<Nt_first+Nt_local) )
    {
        for (unsigned int ik = dk; ik < dk+1; ik++)
        {
            for (unsigned int is = ds; is < ds+1; is++)
            {
                ExtractSliceLocal(evec3d_, epack.evec[ik], 0, pin_time - Nt_first, nd_-1);
                evec3d_ = evec3d_*noise(dt, ik, is);
                tmp3d_  = Zero();
                pokeSpin(tmp3d_, evec3d_, is);
                tmp4d_ = Zero();
                InsertSliceLocal(tmp3d_, tmp4d_, 0, t - Nt_first, nd_-1);
                rho_component += tmp4d_;
            }
        }
    }
    // std::cout << rho_component << std::endl;
    // std::cin.get();
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePinnedDvLapSpinBlock(std::map<Side, DistilVector&>                dv,
                               std::vector<unsigned int>                dt_list,
                               std::map<Side, unsigned int>             n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               const unsigned int                       delta_t,
                               std::map<Side, PerambTensor&>            peramb)
{
    for(unsigned int D=0 ; D<dilSizeLS_.at(s) ; D++)
        dv.at(s)[D] = Zero();

    for(unsigned int D=0 ; D<distilNoise_.at(s).dilutionSize() ; D++)
    {
        std::array<unsigned int,3> d_coor = distilNoise_.at(s).dilutionCoordinates(D);
        unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
        if( std::count(dt_list.begin(), dt_list.end(), dt) )
        {
            const unsigned int Dpinned = distilNoise_.at(s).dilutionIndex(0,dk,ds);
            std::vector<unsigned int>::iterator itr_dt = std::find(dt_list.begin(), dt_list.end(), dt);
            unsigned int idt = std::distance(dt_list.begin(), itr_dt);
            if(isPhi(s))
            {
                makePinnedPhiComponent(dv.at(s)[Dpinned] , distilNoise_.at(s) , n_idx.at(s) , D , delta_t , peramb.at(s), epack);
            }
            else if(isRho(s))
            {
                makePinnedRhoComponent(dv.at(s)[Dpinned] , distilNoise_.at(s) , n_idx.at(s), D, delta_t , epack);
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::execute(const FilenameFn                              &filenameDmfFn,
          const MetadataFn                              &metadataDmfFn,
          std::vector<Gamma::Algebra>                   gamma,
          std::map<Side, DistilVector&>                 dv,
          std::map<Side, unsigned int>                  n_idx,
          std::vector<ComplexField>                     ph,
          std::map<Side, std::vector<unsigned int>>     time_dil_source,
          LapPack&                                      epack,
          TimerArray*                                   tarray,
          std::map<Side, PerambTensor&>                 peramb)
{
    std::vector<std::vector<unsigned int>> time_dil_pair_list;
    const unsigned int vol = g_->_gsites;

    //make these input parameters
    std::vector<unsigned int> delta_t_list={0,1,2,3,4,5,6,7};

    Side pinned_side, fixed_side;
    pinned_side=Side::left;
    // hard-coding this for now to not complicate inputs
    // if(isRho(Side::left) && isPhi(Side::right))
    // {
    //     pinned_side = Side::left;
    // }
    // else if(isPhi(Side::left) && isPhi(Side::right))
    // {
    //     pinned_side = Side::right;
    // }
    // else
    // {
    //     pinned_side = Side::right;
    // }
    fixed_side = (pinned_side==Side::right ? Side::left : Side::right);

    if(isRho(pinned_side))  //if the pinned side has a rho field, force delta_t to be 0, as the others should yield trivial zeros
    {
        delta_t_list={0};
    }

    for(auto delta_t : delta_t_list)
    {        
        START_TIMER("distil vectors");
        makePinnedDvLapSpinBlock(dv, time_dil_source.at(pinned_side), n_idx, epack, pinned_side, delta_t, peramb);
        STOP_TIMER("distil vectors");

        for (unsigned int ibatchFixed=0 ; ibatchFixed<time_dil_source.at(fixed_side).size()/dvBatchSize_ ; ibatchFixed++)   //loop over left dv batches
        {
            START_TIMER("distil vectors");
            std::vector<unsigned int> batch_dtFixed;
            batch_dtFixed = fetchDvBatchIdxs(ibatchFixed,time_dil_source.at(fixed_side));
            makeDvLapSpinBatch(dv, n_idx, epack, fixed_side, batch_dtFixed, peramb);
            STOP_TIMER("distil vectors");
            for (unsigned int idtFixed=0 ; idtFixed<batch_dtFixed.size() ; idtFixed++)
            {
                unsigned int Tfixed = batch_dtFixed[idtFixed];

                std::vector<unsigned int> pinned_times={};
                if( (isRho(Side::left) && isRho(Side::right)) 
                    || onlyDiag_)
                {
                    pinned_times.clear();
                    pinned_times.push_back(Tfixed); //only diagonal rhorho blocks don't vanish
                }
                else
                {
                    for(auto t : time_dil_source.at(pinned_side))
                    {
                        pinned_times.push_back(t+delta_t);
                    }
                }

                if(pinned_side==Side::right)
                {
                    LOG(Message) << "------------------------ " << Tfixed << " X ( " << MDistil::timeslicesDump(time_dil_source.at(pinned_side)) << ") ------------------------" << std::endl; 
                }
                else
                {
                    LOG(Message) << "------------------------ ( " << MDistil::timeslicesDump(time_dil_source.at(pinned_side)) << ") X " << Tfixed << " ------------------------" << std::endl; 
                }
                LOG(Message) << "Time shift : " << delta_t << std::endl; 

                unsigned int nblocki = dilSizeLS_.at(Side::left)/blockSize_ + (((dilSizeLS_.at(Side::left) % blockSize_) != 0) ? 1 : 0);
                unsigned int nblockj = dilSizeLS_.at(Side::right)/blockSize_ + (((dilSizeLS_.at(Side::right) % blockSize_) != 0) ? 1 : 0);

                // loop over blocks within the current time-dilution block
                for(unsigned int i=0 ; i<dilSizeLS_.at(Side::left) ; i+=blockSize_) //set according to memory size
                for(unsigned int j=0 ; j<dilSizeLS_.at(Side::right) ; j+=blockSize_)
                {
                    double flops=0.0, bytes=0.0, time_kernel=0.0, nodes=g_->NodeCount();
                    // iblock_size is the size of the current block (indexed by i); N_i-i is the size of the possible remainder block
                    unsigned int iblock_size = MIN(dilSizeLS_.at(Side::left)-i,blockSize_);
                    unsigned int jblock_size = MIN(dilSizeLS_.at(Side::right)-j,blockSize_);
                    A2AMatrixSet<Tio> block(bBuf_.data(), nExt_ , nStr_ , nt_, iblock_size, jblock_size);

                    LOG(Message) << "Distil matrix block " 
                    << j/blockSize_ + nblocki*i/blockSize_ + 1 
                    << "/" << nblocki*nblockj << " [" << i << " .. " 
                    << i+iblock_size-1 << ", " << j << " .. " << j+jblock_size-1 << "]" 
                    << std::endl;

                    // loop over cache blocks within the current block
                    for(unsigned int ii=0 ; ii<iblock_size ; ii+=cacheSize_)
                    for(unsigned int jj=0 ; jj<jblock_size ; jj+=cacheSize_)
                    {
                        unsigned int icache_size = MIN(iblock_size-ii,cacheSize_);      
                        unsigned int jcache_size = MIN(jblock_size-jj,cacheSize_);
                        A2AMatrixSet<T> cache(cBuf_.data(), nExt_, nStr_, nt_, icache_size, jcache_size);

                        double timer = 0.0;
                        START_TIMER("kernel");
                        unsigned int dv_idxLeftOffset = (pinned_side==Side::right) ? idtFixed*dilSizeLS_.at(Side::left) : 0;
                        unsigned int dv_idxRightOffset = (pinned_side==Side::left) ? idtFixed*dilSizeLS_.at(Side::right) : 0;
                        A2Autils<FImpl>::MesonField(cache, &dv.at(Side::left)[dv_idxLeftOffset +i+ii], &dv.at(Side::right)[dv_idxRightOffset +j+jj], gamma, ph, nd_ - 1, &timer);
                        STOP_TIMER("kernel");
                        time_kernel += timer;

                        flops += vol*(2*8.0+6.0+8.0*nExt_)*icache_size*jcache_size*nStr_;
                        bytes += vol*(12.0*sizeof(T))*icache_size*jcache_size
                                +  vol*(2.0*sizeof(T)*nExt_)*icache_size*jcache_size*nStr_;

                        // copy from cache
                        START_TIMER("cache copy");
                        thread_for_collapse(5,iext,nExt_,{
                        // for(unsigned int iext=0;iext<nExt_;iext++)
                        for(unsigned int istr=0;istr<nStr_;istr++)
                        for(unsigned int t=0;t<nt_;t++)     // TODO: could loop just through the times I know are non zero now...
                        for(unsigned int iii=0;iii<icache_size;iii++)
                        for(unsigned int jjj=0;jjj<jcache_size;jjj++)
                        {
                            block(iext,istr,t,ii+iii,jj+jjj) = cache(iext,istr,t,iii,jjj);
                        }
                        });
                        STOP_TIMER("cache copy");
                    }

                    LOG(Message) << "Kernel perf (flops) " << flops/time_kernel/1.0e3/nodes 
                                << " Gflop/s/node " << std::endl;
                    LOG(Message) << "Kernel perf (read) " << bytes/time_kernel*1.0e6/1024/1024/1024/nodes
                                << " GB/s/node "  << std::endl;
                    blockCounter_++;
                    blockFlops_ += flops/time_kernel/1.0e3/nodes ;
                    blockBytes_ += bytes/time_kernel*1.0e6/1024/1024/1024/nodes;

                    // saving current block to disk
#ifdef HADRONS_A2AM_PARALLEL_IO     //only one implemented, todo: add exception
                    //parallel io
                    unsigned int inode = g_->ThisRank();
                    unsigned int nnode = g_->RankCount(); 
                    for(unsigned int it=0 ; it<time_dil_source.at(pinned_side).size() ; it++)   // TODO: generalise to both sides; include this in the node loop?
                    {
                        
                        // TODO: generalise to dilution
                        unsigned int t = (time_dil_source.at(pinned_side)[it] + delta_t)%nt_;
                        const unsigned int Tpinned = time_dil_source.at(pinned_side)[it];

                        DistilMatrixSet<Tio> block_pinned(bBuf_.data(), nExt_ , nStr_ , iblock_size, jblock_size);
                        // START_TIMER("cache copy");
                        // thread_for_collapse(4,iext,nExt_,{
                        // // for(unsigned int iext=0;iext<nExt_;iext++)
                        // for(unsigned int istr=0;istr<nStr_;istr++)
                        // for(unsigned int ii=0;ii<iblock_size;ii++)
                        // for(unsigned int jj=0;jj<jblock_size;jj++)
                        // {
                        //     block_pinned(iext,istr,ii,jj) = block(iext,istr,t,ii,jj);
                        // }
                        // });
                        // STOP_TIMER("cache copy");

                        std::string dataset_name;
                        dataset_name = std::to_string( (pinned_side==Side::right) ? Tfixed : Tpinned ) 
                            + "-" + std::to_string( (pinned_side==Side::right) ? Tpinned : Tfixed );   


                        LOG(Message) << "Starting parallel IO. Rank count=" << nnode << ". Sparse " << dataset_name << " block" << std::endl;
                        double ioTime = -GET_TIMER("IO: write block");
                        START_TIMER("IO: total");
                        g_->Barrier();
                        for(unsigned int k=inode ; k<nExt_*nStr_ ; k+=nnode){
                            unsigned int iext = k/nStr_;
                            unsigned int istr = k%nStr_;
                            // metadata;
                            DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right));

                            DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(filenameDmfFn(iext, istr, n_idx.at(Side::left), n_idx.at(Side::right)),
                                    DISTIL_MATRIX_NAME, nt_, dilSizeLS_.at(Side::left), dilSizeLS_.at(Side::right));

                            if( ( Tfixed==time_dil_source.at(fixed_side).front() ) &&
                                ( Tpinned==time_dil_source.at(pinned_side).front() ) &&
                                ( t==(time_dil_source.at(pinned_side).front() + delta_t_list.front())%nt_ ) && 
                                (i==0) && (j==0) )  //executes only once per file
                            {
                                START_TIMER("IO: file creation");
                                matrix_io.initFile(md);
                                STOP_TIMER("IO: file creation");
                            }
                            START_TIMER("IO: write block");
                            matrix_io.saveBlock(block_pinned, iext, istr, i, j, dataset_name, t, blockSize_);
                            STOP_TIMER("IO: write block");
                        }
                        g_->Barrier();
                        STOP_TIMER("IO: total");
                        ioTime    += GET_TIMER("IO: write block");
                        unsigned int bytesBlockSize  = static_cast<double>(nExt_*nStr_*iblock_size*jblock_size*sizeof(Tio));
                        double iospeed = bytesBlockSize/ioTime*1.0e6/1024/1024;
                        LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                        << ioTime  << " us (" << iospeed << " MB/s) (chunking "
                                        << blockSize_ << "x" << blockSize_ << ")" << std::endl;
                        blockIoSpeed_ += iospeed;
                    }
                }
#endif              
            }
        }
    }
    
    // //write dilution schemes and time source pairs
    // if(g_->IsBoss())
    // {
    //     START_TIMER("IO: total");
    //     for(unsigned int iext=0 ; iext<nExt_ ; iext++)
    //     for(unsigned int istr=0 ; istr<nStr_ ; istr++)
    //     {
    //         DistilMetadataIo mdIo(filenameDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right)),
    //                 std::string(DISTIL_MATRIX_NAME) + "/" + std::string(METADATA_NAME) );
    //         mdIo.write2dMetadata("TimeSourcePairs", time_dil_pair_list);
    //         //  dilution schemes (2d ragged metadata) - replace by grid serialisation code when support to ragged vectors is done
    //         DilutionMap lmap = getMap(Side::left);
    //         DilutionMap rmap = getMap(Side::right);
    //         mdIo.write2dMetadata("TimeDilutionLeft" , lmap[Index::t], DILUTION_METADATA_NAME);
    //         mdIo.write2dMetadata("TimeDilutionRight", rmap[Index::t], DILUTION_METADATA_NAME);
    //         mdIo.write2dMetadata("LapDilutionLeft"  , lmap[Index::l], DILUTION_METADATA_NAME);
    //         mdIo.write2dMetadata("LapDilutionRight" , rmap[Index::l], DILUTION_METADATA_NAME);
    //         mdIo.write2dMetadata("SpinDilutionLeft" , lmap[Index::s], DILUTION_METADATA_NAME);
    //         mdIo.write2dMetadata("SpinDilutionRight", rmap[Index::s], DILUTION_METADATA_NAME);
    //     }
    //     STOP_TIMER("IO: total");
    // }

    

    // TestMd<FImpl> md;
    // md.name = "heyo";


    // bBuf_.resize(1*1*12*12);
    // DistilMatrixSet<Tio> test_block(bBuf_.data(), 1 , 1 , 6, 6);
    // DistilMatrixSet<Tio> test_block_2(bBuf_.data(), 1 , 1 , 6, 6);

    // thread_for_collapse(4,iext,1,{
    // // for(unsigned int iext=0;iext<nExt_;iext++)
    // for(unsigned int istr=0;istr<1;istr++)
    // for(unsigned int ii=0;ii<6;ii++)
    // for(unsigned int jj=0;jj<6;jj++)
    // {
    //     test_block(iext,istr,ii,jj) = 0;
    // }
    // });

    // // test_block(0,0,0,0) = 66.6;
    // test_block(0,0,3,3) = 33.3;

    // DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io("test.h5", DISTIL_MATRIX_NAME, 12, 12);
    // matrix_io.initFile(md);
    // matrix_io.saveBlock(test_block, 0, 0, 5, 5, "0-0", "t666" , 6);

    // test_block_2(0,0,2,2) = 66.6;
    // matrix_io.saveBlock(test_block_2, 0, 0, 5, 5, "0-0", "t555" , 6);



}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Distil_matrix_hpp_
