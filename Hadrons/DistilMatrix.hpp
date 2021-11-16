#ifndef Distil_matrix_hpp_
#define Distil_matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>
#include <Hadrons/DistillationVectors.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

// number of distil vector time-dilution components held in memory at the same time
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

// Dimensions of distil matrix set:
//   0 : ext - external field (momentum phase, etc...)
//   1 : str - spin-color structure or covariant derivative operators
//   2 : i   - left  spin-Laplacian mode index
//   3 : j   - right spin-Laplacian mode index
template <typename T>
using DistilMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 4, Eigen::RowMajor>>;
using namespace MDistil;

// Dimensions of distil matrix:
//   0 : i   - left  spin-Laplacian mode index
//   1 : j   - right spin-Laplacian mode index
template <typename T>
using DistilMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

using DilutionMap  = std::array<std::vector<std::vector<unsigned int>>,3>;

enum Side {left = 0, right = 1};
const std::vector<Side> sides =  {Side::left,Side::right};  //for easy looping over sides

// metadata serialiser class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    unsigned int,                               Nt,
                                    unsigned int,                               Nvec,
                                    std::vector<RealF>,                         Momentum,
                                    Gamma::Algebra,                             Operator,               // potentially more general operators in the future
                                    std::vector<unsigned int>,                  NoisePair,
                                    std::string,                                MesonFieldType,
                                    std::string,                                NoiseHashLeft,
                                    std::string,                                NoiseHashRight,
                                    std::vector<std::vector<unsigned int>>,     TimeDilutionLeft,
                                    std::vector<std::vector<unsigned int>>,     TimeDilutionRight,
                                    std::vector<std::vector<unsigned int>>,     LapDilutionLeft,
                                    std::vector<std::vector<unsigned int>>,     LapDilutionRight,
                                    std::vector<std::vector<unsigned int>>,     SpinDilutionLeft,
                                    std::vector<std::vector<unsigned int>>,     SpinDilutionRight,
                                    std::string,                                RelativeSide,
                                    )
};

/******************************************************************************
 *                  Class to handle Distil matrix block HDF5 I/O              *
 ******************************************************************************/
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
    size_t       getSize(void) const;

    // block I/O
    template <typename MetadataType>
    void initFile(const MetadataType &d);

    void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                   const unsigned int blockSizei, const unsigned int blockSizej, std::string t_name, std::string datasetName, const unsigned int chunkSize);

    void saveBlock(const DistilMatrixSet<T> &m, const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j, std::string datasetName,
                               const unsigned int t, const unsigned int chunkSize);

    template <template <class> class Vec, typename VecT>
    void load(Vec<VecT> &v, const unsigned int t, const std::string dataset_name, double *tRead = nullptr, GridBase *grid = nullptr);
private:
    std::string  filename_{""}, dataname_{""};
    unsigned int nt_{0}, ni_{0}, nj_{0};
};

// implementation /////////////////////////////////////////////////////////////////
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
size_t DistilMatrixIo<T>::getSize(void) const
{
    return nt_*ni_*nj_*sizeof(T);
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

template <typename T>
template <template <class> class Vec, typename VecT>
void DistilMatrixIo<T>::load(Vec<VecT> &v, const unsigned int t, const std::string dataset_name, double *tRead, GridBase *grid)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;

    if (!(grid) or grid->IsBoss())
    {
        Hdf5Reader reader(filename_);
        push(reader, dataname_);
        auto &root_group = reader.getGroup();
        std::string t_name = std::to_string(t);
        auto t_group = root_group.openGroup(t_name);
        dataset = t_group.openDataSet(dataset_name);
        datatype = dataset.getCompType();
        dataspace = dataset.getSpace();
        hdim.resize(dataspace.getSimpleExtentNdims());
        dataspace.getSimpleExtentDims(hdim.data());
        if ((ni_ * nj_ != 0) and
            ((hdim[0] != ni_) or (hdim[0] != nj_)))
        {
            HADRONS_ERROR(Size, "distil matrix size mismatch (got "
                + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + ", expected "
                + std::to_string(ni_) + "x" + std::to_string(nj_));
        }
        ni_ = hdim[0];
        nj_ = hdim[1];
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), &ni_, sizeof(unsigned int));
        grid->Broadcast(grid->BossRank(), &nj_, sizeof(unsigned int));
    }

    DistilMatrix<T>         buf(ni_, nj_);
    int broadcastSize =  sizeof(T) * buf.size();
    std::vector<hsize_t> count    = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)},
                         stride   = {1, 1},
                         block    = {1, 1},
                         memCount = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)};
    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    std::cout << "Loading timeslice";
    std::cout.flush();
    *tRead = 0.;

    std::vector<hsize_t> offset = {0, 0};
    
    if (!(grid) or grid->IsBoss())
    {
        dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                    stride.data(), block.data());
    }
    if (tRead) *tRead -= usecond();
    if (!(grid) or grid->IsBoss())
    {
        dataset.read(buf.data(), datatype, memspace, dataspace);
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), buf.data(), broadcastSize);
    }
    if (tRead) *tRead += usecond();
    v[0] = buf.template cast<VecT>();
    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

//####################################
//# computation class declaration    #
//####################################

template <typename FImpl, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef std::vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  TrajNumberFn;
    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int, DilutionMap, DilutionMap)>  MetadataFn;
public:
    long    blockCounter_ = 0;
    double  blockFlops_ = 0.0, blockBytes_ = 0.0, blockIoSpeed_ = 0.0;
private:
    std::map<Side,std::string>          dmfType_;
    GridCartesian*                      g_;
    GridCartesian*                      g3d_;
    ColourVectorField                   evec3d_;
    FermionField                        tmp3d_, tmp4d_;
    std::vector<Tio>                    bBuf_; //potentially large objects
    std::vector<T>                      cBuf_;
    const unsigned int                  blockSize_; //eventually turns into io chunk size
    const unsigned int                  cacheSize_;
    const unsigned int                  traj_;
    const unsigned int                  nt_, nd_, nExt_, nStr_;
    const bool                          isExact_;
    bool                                isInitFile_=false;
    std::map<Side, unsigned int>        dilSizeLS_;
    std::map<Side, DistillationNoise&>  distilNoise_;
    const unsigned int                  dvBatchSize_ = DISTILVECTOR_TIME_BATCH_SIZE;
    std::map<Side, std::string>         vectorStem_={{Side::left,""},{Side::right,""}};
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
                   const unsigned int           traj,
                   const std::string            left_vector_stem="",
                   const std::string            right_vector_stem="");
    bool isPhi(Side s);
    bool isRho(Side s);
    DilutionMap fetchDilutionMap(Side s);
private:
    void makePhiComponent(FermionField&         phi_component,
                          DistillationNoise&    n,
                          const unsigned int    n_idx,
                          const unsigned int    D,
                          PerambTensor&         peramb,
                          LapPack&              epack);
    void loadPhiComponent(FermionField&            phi_component,
                          DistillationNoise&       n,
                          const unsigned int       D,
                          std::string              vector_stem,
                          PerambTensor&            peramb,
                          LapPack&                 epack);
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
                                               std::vector<unsigned int>  time_dil_sources,
                                               unsigned int                    shift=0);
    void makeRelativePhiComponent(FermionField&                phi_component,
                            DistillationNoise&               n,
                            const unsigned int               n_idx,
                            const unsigned int               D,
                            const unsigned int               delta_t,
                            PerambTensor&                    peramb,
                            LapPack&                         epack);
    void makeRelativeRhoComponent(FermionField&            rho_component,
                            DistillationNoise&           n,
                            const unsigned int           n_idx,
                            const unsigned int           t,
                            const unsigned int           D,
                            LapPack&                     epack);
    void makeRelativeDvLapSpinBlock(std::map<Side, DistilVector&>     dv,
                               std::vector<unsigned int>            dt_list,
                               std::map<Side, unsigned int>         n_idx,
                               LapPack&                             epack,
                               Side                                 s,
                               const unsigned int                   delta_t,
                               std::map<Side, PerambTensor&>        peramb);
public:
    void executeRelative(const FilenameFn                               &filenameDmfFn,
                       const MetadataFn                               &metadataDmfFn,
                       std::vector<Gamma::Algebra>                    gamma,
                       std::map<Side, DistilVector&>                  dv,
                       std::map<Side, unsigned int>                   n_idx,
                       std::vector<ComplexField>                      ph,
                       std::map<Side, std::vector<unsigned int>>      time_dil_source,
                       LapPack&                                       epack,
                       TimerArray*                                    tarray,
                       Side                                           relative_side,
                       std::vector<unsigned int>                      delta_t_list,
                       std::map<Side, PerambTensor&>                  peramb={});
    void executeFixed(const FilenameFn                               &filenameDmfFn,
                 const MetadataFn                               &metadataDmfFn,
                 std::vector<Gamma::Algebra>                    gamma,
                 std::map<Side, DistilVector&>                  dv,
                 std::map<Side, unsigned int>                   n_idx,
                 std::vector<ComplexField>                      ph,
                 std::map<Side, std::vector<unsigned int>>      time_dil_source,
                 LapPack&                                       epack,
                 TimerArray*                                    tarray,
                 bool                                           only_diag,
                 const unsigned int                             diag_shift,
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
                 const unsigned int             traj,
                 const std::string              left_vector_stem,
                 const std::string              right_vector_stem)
: dmfType_(mf_type), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d), tmp4d_(g)
, nt_(nt) , nd_(g->Nd()), blockSize_(block_size) , cacheSize_(cache_size)
, nExt_(n_ext) , nStr_(n_str) , isExact_(is_exact), traj_(traj)
{
    cBuf_.resize(nExt_*nStr_*nt_*cacheSize_*cacheSize_);
    bBuf_.resize(nExt_*nStr_*nt_*blockSize_*blockSize_); //maximum size

    distilNoise_ = std::map<Side, DistillationNoise&> ({{Side::left,nl},{Side::right,nr}});
    dilSizeLS_ = { {Side::left,nl.dilutionSize(Index::l)*nl.dilutionSize(Index::s)} ,
                             {Side::right,nr.dilutionSize(Index::l)*nr.dilutionSize(Index::s)} };

    vectorStem_ = { {Side::left,left_vector_stem} , {Side::right,right_vector_stem}};
}

template <typename FImpl, typename T, typename Tio>
DilutionMap DmfComputation<FImpl,T,Tio>::fetchDilutionMap(Side s)
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
::loadPhiComponent(FermionField&            phi_component,
                   DistillationNoise&       n,
                   const unsigned int       D,
                   std::string              vector_stem,
                   PerambTensor&            peramb,
                   LapPack&                 epack)
{
    const unsigned int nVec = epack.evec.size();
    DistillationVectorsIo::readComponent(phi_component, vector_stem, n.size(),
				    n.dilutionSize(Index::l), n.dilutionSize(Index::s), n.dilutionSize(Index::t), peramb.MetaData.timeSources, D, traj_);
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
// and each is identified by the starting posision in time dilution space, (dtL,dtR)
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
            if(vectorStem_.at(s).empty())
            {
                makePhiComponent(dv.at(s)[iD] , distilNoise_.at(s) , n_idx.at(s) , D , peramb.at(s), epack);
            }
            else
            {
                loadPhiComponent(dv.at(s)[iD] , distilNoise_.at(s) , D , vectorStem_.at(s), peramb.at(s), epack);
            }
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
::fetchDvBatchIdxs(unsigned int ibatch, std::vector<unsigned int> time_dil_sources, const unsigned int shift)
{
    std::vector<unsigned int> batch_dt;
    for(unsigned int dt=ibatch*dvBatchSize_; dt<(ibatch+1)*dvBatchSize_; dt++){
        batch_dt.push_back( (time_dil_sources[dt]+shift)%distilNoise_.at(Side::right).dilutionSize(Index::t) );
    }
    return batch_dt;
}

// makeRelative methods build distil vectors with reorganised time slices
// (in order to compute multiple time-dilution blocks at different t with a single call of MesonField kernel)
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativePhiComponent(FermionField&              phi_component,
                   DistillationNoise&               n,
                   const unsigned int               n_idx,
                   const unsigned int               D,
                   const unsigned int               delta_t,
                   PerambTensor&                    peramb,
                   LapPack&                         epack)
{
    std::array<unsigned int,3> d_coor = n.dilutionCoordinates(D);
    unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
    
    //compute at relative_time, insert at t=relative_time
    const unsigned int relative_time = (dt + delta_t)%nt_;               //TODO: generalise to dilution
    const unsigned int t        = relative_time;                         //TODO: generalise to dilution 

    std::vector<int> peramb_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(peramb_ts.begin(), peramb_ts.end(), dt);
    unsigned int idt_peramb = std::distance(peramb_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj 
    const unsigned int nVec = epack.evec.size();
    const unsigned int Nt_first = g_->LocalStarts()[nd_ - 1];
    const unsigned int Nt_local = g_->LocalDimensions()[nd_ - 1];
    if( (relative_time>=Nt_first) and (relative_time<Nt_first+Nt_local) )
    {
        tmp3d_ = Zero();
        for (unsigned int k = 0; k < nVec; k++)
        {
            ExtractSliceLocal(evec3d_,epack.evec[k],0,relative_time-Nt_first,nd_ - 1);
            tmp3d_ += evec3d_ * peramb.tensor(relative_time, k, dk, n_idx, idt_peramb, ds);
        }
        InsertSliceLocal(tmp3d_,phi_component,0,t-Nt_first,nd_ - 1);
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativeRhoComponent(FermionField&          rho_component,
                   DistillationNoise&           n,
                   const unsigned int           n_idx,
                   const unsigned int           D,
                   const unsigned int           delta_t,
                   LapPack&                     epack)
{
    // abstract this to makeRelativeSource() kind of method in DilutedNoise?
    std::array<unsigned int,3> d_coor = n.dilutionCoordinates(D);
    unsigned int dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];

    typename DistillationNoise::NoiseType   noise(n.getNoise()[n_idx].data(), nt_, epack.eval.size(), Ns);

    //compute at relative_time, insert at t
    const unsigned int relative_time = (dt + delta_t)%nt_;     //TODO: generalise to dilution
    const unsigned int t        = relative_time;               //TODO: generalise to dilution

    // adapt to dilution! (add an it loop and untrivialise ik and is loops)
    
    const unsigned int Nt_first = g_->LocalStarts()[nd_ - 1];
    const unsigned int Nt_local = g_->LocalDimensions()[nd_ - 1];
    if( (relative_time>=Nt_first) and (relative_time<Nt_first+Nt_local) )
    {
        for (unsigned int ik = dk; ik < dk+1; ik++)
        {
            for (unsigned int is = ds; is < ds+1; is++)
            {
                ExtractSliceLocal(evec3d_, epack.evec[ik], 0, relative_time - Nt_first, nd_-1);
                evec3d_ = evec3d_*noise(dt, ik, is);
                tmp3d_  = Zero();
                pokeSpin(tmp3d_, evec3d_, is);
                tmp4d_ = Zero();
                InsertSliceLocal(tmp3d_, tmp4d_, 0, t - Nt_first, nd_-1);
                rho_component += tmp4d_;
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativeDvLapSpinBlock(std::map<Side, DistilVector&>                dv,
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
        if( std::count(dt_list.begin(), dt_list.end(), dt)!=0 )    //if dt is in dt_list
        {
            const unsigned int Drelative = distilNoise_.at(s).dilutionIndex(0,dk,ds);
            std::vector<unsigned int>::iterator itr_dt = std::find(dt_list.begin(), dt_list.end(), dt);
            unsigned int idt = std::distance(dt_list.begin(), itr_dt);
            if(isPhi(s))
            {
                makeRelativePhiComponent(dv.at(s)[Drelative] , distilNoise_.at(s) , n_idx.at(s) , D , delta_t , peramb.at(s), epack);
            }
            else if(isRho(s))
            {
                makeRelativeRhoComponent(dv.at(s)[Drelative] , distilNoise_.at(s) , n_idx.at(s), D, delta_t , epack);
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::executeRelative(const FilenameFn                              &filenameDmfFn,
                const MetadataFn                              &metadataDmfFn,
                std::vector<Gamma::Algebra>                   gamma,
                std::map<Side, DistilVector&>                 dv,
                std::map<Side, unsigned int>                  n_idx,
                std::vector<ComplexField>                     ph,
                std::map<Side, std::vector<unsigned int>>     time_dil_source,
                LapPack&                                      epack,
                TimerArray*                                   tarray,
                Side                                          relative_side,
                std::vector<unsigned int>                     delta_t_list,
                std::map<Side, PerambTensor&>                 peramb)
{
    const unsigned int vol = g_->_gsites;
    Side anchored_side = (relative_side==Side::right ? Side::left : Side::right);

    for(auto delta_t : delta_t_list)
    {        
        START_TIMER("distil vectors");
        makeRelativeDvLapSpinBlock(dv, time_dil_source.at(relative_side), n_idx, epack, relative_side, delta_t, peramb);
        STOP_TIMER("distil vectors");

        //loop over left dv batches
        for (unsigned int ibatchAnchored=0 ; ibatchAnchored<time_dil_source.at(anchored_side).size()/dvBatchSize_ ; ibatchAnchored++)
        {
            START_TIMER("distil vectors");
            std::vector<unsigned int> batch_dtAnchored;
            batch_dtAnchored = fetchDvBatchIdxs(ibatchAnchored,time_dil_source.at(anchored_side));
            makeDvLapSpinBatch(dv, n_idx, epack, anchored_side, batch_dtAnchored, peramb);
            STOP_TIMER("distil vectors");
            for (unsigned int idtAnchored=0 ; idtAnchored<batch_dtAnchored.size() ; idtAnchored++)
            {
                unsigned int Tanchored = batch_dtAnchored[idtAnchored];

                if(relative_side==Side::right)
                {
                    LOG(Message) << "------------ " << Tanchored << " X ( " << MDistil::timeslicesDump(time_dil_source.at(relative_side)) << ") ------------" << std::endl; 
                }
                else
                {
                    LOG(Message) << "------------ ( " << MDistil::timeslicesDump(time_dil_source.at(relative_side)) << ") X " << Tanchored << " ------------" << std::endl; 
                }
                LOG(Message) << "Time shift (deltaT) : " << delta_t << std::endl; 

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
                        unsigned int dv_idxLeftOffset = (relative_side==Side::right) ? idtAnchored*dilSizeLS_.at(Side::left) : 0;
                        unsigned int dv_idxRightOffset = (relative_side==Side::left) ? idtAnchored*dilSizeLS_.at(Side::right) : 0;
                        A2Autils<FImpl>::MesonField(cache, &dv.at(Side::left)[dv_idxLeftOffset +i+ii], &dv.at(Side::right)[dv_idxRightOffset +j+jj], gamma, ph, nd_ - 1, &timer);
                        STOP_TIMER("kernel");
                        time_kernel += timer;

                        flops += vol*(2*8.0+6.0+8.0*nExt_)*icache_size*jcache_size*nStr_;
                        bytes += vol*(12.0*sizeof(T))*icache_size*jcache_size
                                +  vol*(2.0*sizeof(T)*nExt_)*icache_size*jcache_size*nStr_;

                        // copy from cache
                        START_TIMER("cache copy");
                        thread_for_collapse(5,iext,nExt_,{
                        for(unsigned int istr=0;istr<nStr_;istr++)
                        for(unsigned int t=0;t<nt_;t++)     // TODO(low priority): could loop just over the times I know are non zero now...
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

                    // io section
                    unsigned int inode = g_->ThisRank();
                    unsigned int nnode = g_->RankCount(); 
                    LOG(Message) << "Starting parallel IO. Rank count=" << nnode << std::endl;
                    for(unsigned int it=0 ; it<time_dil_source.at(relative_side).size() ; it++)
                    {
                        // TODO: generalise to dilution
                        unsigned int t = (time_dil_source.at(relative_side)[it] + delta_t)%nt_;
                        const unsigned int Trelative = time_dil_source.at(relative_side)[it];

                        DistilMatrixSet<Tio> block_relative(bBuf_.data(), nExt_ , nStr_ , iblock_size, jblock_size);

                        std::string dataset_name;
                        dataset_name = std::to_string( (relative_side==Side::right) ? Tanchored : Trelative ) 
                            + "-" + std::to_string( (relative_side==Side::right) ? Trelative : Tanchored );

                        std::vector<unsigned int> relative_partition = distilNoise_.at(relative_side).dilutionPartition(Index::t,Trelative);
                        std::vector<unsigned int> anchored_partition = distilNoise_.at(anchored_side).dilutionPartition(Index::t,Tanchored);

                        if( !( isRho(relative_side) and std::count(relative_partition.begin(), relative_partition.end(), t)==0  ) 
                            and !(isRho(anchored_side) and std::count(anchored_partition.begin(), anchored_partition.end(), t)==0) )
                        {
                            LOG(Message)    << "Saving block " << dataset_name << " , t=" << t << std::endl;

                            double ioTime = -GET_TIMER("IO: write block");
                            START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
                            g_->Barrier();
                            for(unsigned int k=inode ; k<nExt_*nStr_ ; k+=nnode)
                            {
                                unsigned int iext = k/nStr_;
                                unsigned int istr = k%nStr_;

                                // metadata;
                                DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(filenameDmfFn(iext, istr, n_idx.at(Side::left), n_idx.at(Side::right)),
                                        DISTIL_MATRIX_NAME, nt_, dilSizeLS_.at(Side::left), dilSizeLS_.at(Side::right));

                                //executes once per file
                                if( ( Tanchored==time_dil_source.at(anchored_side).front() ) and      //first time-dilution idx at one side
                                    ( Trelative==time_dil_source.at(relative_side).front() ) and    // same as above for the other side
                                    ( t==(time_dil_source.at(relative_side).front() + delta_t_list.front())%nt_ ) and    //first time slice
                                    (i==0) and (j==0) )  //first IO block
                                {
                                    //fetch metadata
                                    DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right),fetchDilutionMap(Side::left),fetchDilutionMap(Side::right));
                                    //init file and write metadata
                                    START_TIMER("IO: file creation");
                                    matrix_io.initFile(md);
                                    STOP_TIMER("IO: file creation");
                                }
                                START_TIMER("IO: write block");
                                matrix_io.saveBlock(block_relative, iext, istr, i, j, dataset_name, t, blockSize_);
                                STOP_TIMER("IO: write block");
                            }
                            g_->Barrier();
#else
    HADRONS_ERROR(Implementation, "DistilMesonField serial IO not implemented.");
#endif              
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
                }
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::executeFixed(const FilenameFn                         &filenameDmfFn,
          const MetadataFn                              &metadataDmfFn,
          std::vector<Gamma::Algebra>                   gamma,
          std::map<Side, DistilVector&>                 dv,
          std::map<Side, unsigned int>                  n_idx,
          std::vector<ComplexField>                     ph,
          std::map<Side, std::vector<unsigned int>>     time_dil_source,
          LapPack&                                      epack,
          TimerArray*                                   tarray,
          bool                                          only_diag,
          const unsigned int                            diag_shift,
          std::map<Side, PerambTensor&>                 peramb)
{
    const unsigned int vol = g_->_gsites;

    //loop over left dv batches
    for (unsigned int ibatchL=0 ; ibatchL<time_dil_source.at(Side::left).size()/dvBatchSize_ ; ibatchL++)   //loop over left dv batches
    {
        std::vector<unsigned int> batch_dtL = fetchDvBatchIdxs(ibatchL,time_dil_source.at(Side::left));
        for (unsigned int idtL=0 ; idtL<batch_dtL.size() ; idtL++)
        {
            unsigned int dtL = batch_dtL[idtL];
            for (unsigned int ibatchR=0 ; ibatchR<time_dil_source.at(Side::right).size()/dvBatchSize_ ; ibatchR++)  //loop over right dv batches
            {
                std::vector<unsigned int> batch_dtR = fetchDvBatchIdxs(ibatchR,time_dil_source.at(Side::right), diag_shift);
                for (unsigned int idtR=0 ; idtR<batch_dtR.size() ; idtR++)
                {
                    unsigned int dtR = batch_dtR[idtR];
                    // fetch necessary time slices for this time-dilution block
                    std::map<Side,std::vector<unsigned int>> time_partition = { {Side::left,{}} , {Side::right,{}}};
                    for(auto s : sides)
                    {
                        if(isPhi(s))
                        {
                            time_partition.at(s).resize(nt_);
                            //phi: fill with all time slices so the intersection with a vector X is equal to X
                            std::iota(std::begin(time_partition.at(s)), std::end(time_partition.at(s)), 0); 
                        }
                        else
                        {
                            time_partition.at(s) = distilNoise_.at(s).dilutionPartition(Index::t, s==Side::left ? dtL : dtR);
                        }
                    }
                    std::vector<unsigned int> ts_intersection;
                    std::set_intersection(time_partition.at(Side::left).begin(), time_partition.at(Side::left).end(), 
                                        time_partition.at(Side::right).begin(), time_partition.at(Side::right).end(),
                                        std::back_inserter(ts_intersection));

                    // only execute when partitions have at least one time slice in common; only computes diagonal (up to diag_shift) when onlydiag true
                    if( !ts_intersection.empty() and 
                        ( !only_diag or ((dtL+diag_shift)%distilNoise_.at(Side::right).dilutionSize(Index::t)==dtR) ) )
                    {
                        LOG(Message) << "------------------------ " << dtL << "-" << dtR << " ------------------------" << std::endl; 
                        LOG(Message) << "Saving time slices : " << MDistil::timeslicesDump(ts_intersection) << std::endl;

                        START_TIMER("distil vectors");
                        makeDvLapSpinBatch(dv, n_idx, epack, Side::left, batch_dtL, peramb);
                        makeDvLapSpinBatch(dv, n_idx, epack, Side::right, batch_dtR, peramb);
                        STOP_TIMER("distil vectors");

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
                                unsigned int dv_idxL = idtL*dilSizeLS_.at(Side::left) +i+ii;
                                unsigned int dv_idxR = idtR*dilSizeLS_.at(Side::right) +j+jj;
                                A2Autils<FImpl>::MesonField(cache, &dv.at(Side::left)[dv_idxL], &dv.at(Side::right)[dv_idxR], gamma, ph, nd_ - 1, &timer);
                                STOP_TIMER("kernel");
                                time_kernel += timer;

                                flops += vol*(2*8.0+6.0+8.0*nExt_)*icache_size*jcache_size*nStr_;
                                bytes += vol*(12.0*sizeof(T))*icache_size*jcache_size
                                        +  vol*(2.0*sizeof(T)*nExt_)*icache_size*jcache_size*nStr_;

                                // copy cache to block
                                START_TIMER("cache copy");
                                unsigned int ts_size = ts_intersection.size();
                                thread_for_collapse(5,iext,nExt_,{
                                for(unsigned int istr=0;istr<nStr_;istr++)
                                for(unsigned int it=0;it<ts_size;it++)
                                for(unsigned int iii=0;iii<icache_size;iii++)
                                for(unsigned int jjj=0;jjj<jcache_size;jjj++)
                                {
                                    block(iext,istr,ts_intersection[it],ii+iii,jj+jjj)=cache(iext,istr,ts_intersection[it],iii,jjj);
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

                            // io section
                            unsigned int inode = g_->ThisRank();
                            unsigned int nnode = g_->RankCount(); 
                            LOG(Message) << "Starting parallel IO. Rank count=" << nnode << std::endl;
                            for(unsigned int t=0 ; t<nt_ ; t++)
                            {
                                // TODO: generalise to dilution
                                // unsigned int t = ts_intersection[it];
                                std::vector<unsigned int> left_partition = distilNoise_.at(Side::left).dilutionPartition(Index::t,dtL);
                                std::vector<unsigned int> right_partition = distilNoise_.at(Side::right).dilutionPartition(Index::t,dtR);

                                if( !( isRho(Side::left) and std::count(left_partition.begin(), left_partition.end(), t)==0  ) 
                                    and !(isRho(Side::right) and std::count(right_partition.begin(), right_partition.end(), t)==0) )
                                {
                                    DistilMatrixSet<Tio> block_relative(bBuf_.data(), nExt_ , nStr_ , iblock_size, jblock_size);
                                    std::string dataset_name = std::to_string(dtL)+"-"+std::to_string(dtR);
                                    LOG(Message)    << "Saving block block " << dataset_name << " , t=" << t << std::endl;

                                    double ioTime = -GET_TIMER("IO: write block");
                                    START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
                                    g_->Barrier();
                                    for(unsigned int k=inode ; k<nExt_*nStr_ ; k+=nnode)
                                    {
                                        unsigned int iext = k/nStr_;
                                        unsigned int istr = k%nStr_;

                                        // metadata;
                                        DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(filenameDmfFn(iext, istr, n_idx.at(Side::left), n_idx.at(Side::right)),
                                                DISTIL_MATRIX_NAME, nt_, dilSizeLS_.at(Side::left), dilSizeLS_.at(Side::right));

                                        //executes once per file
                                        if( ( dtL==time_dil_source.at(Side::left).front() ) and      //first time-dilution idx at one side
                                            ( dtR==((time_dil_source.at(Side::right).front()+diag_shift)%distilNoise_.at(Side::right).dilutionSize(Index::t)) ) and    // same as above for the other side
                                            ( t==ts_intersection.front() ) and    //first time slice
                                            (i==0) and (j==0) )  //first IO block
                                        {
                                            //fetch metadata
                                            DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right),
                                                                    fetchDilutionMap(Side::left),fetchDilutionMap(Side::right));
                                            //init file and write metadata
                                            START_TIMER("IO: file creation");
                                            matrix_io.initFile(md);
                                            STOP_TIMER("IO: file creation");
                                        }
                                        START_TIMER("IO: write block");
                                        matrix_io.saveBlock(block_relative, iext, istr, i, j, dataset_name, t, blockSize_);
                                        STOP_TIMER("IO: write block");
                                    }
                                    g_->Barrier();
#else
    HADRONS_ERROR(Implementation, "DistilMesonField serial IO not implemented.");
#endif              
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
                        }
                    }
                }
            }
        }
    }
}

// END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Distil_matrix_hpp_
