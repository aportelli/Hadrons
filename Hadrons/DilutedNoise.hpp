/*
 * DilutedNoise.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_DilutedNoise_hpp_
#define Hadrons_DilutedNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Grid/util/Sha.h>

BEGIN_HADRONS_NAMESPACE

class DilutedNoise
{
public:
    DilutedNoise(void) = default;
    virtual ~DilutedNoise(void) = default;
    virtual void resize(const int nNoise) = 0;
    virtual int  size(void) const = 0;
    virtual int  dilutionSize(void) const = 0;
};

/******************************************************************************
 *                 Abstract container for distillation noise                  *
 ******************************************************************************/
template <typename FImpl>
class DistillationNoise: public DilutedNoise
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef EigenPack<ColourVectorField> LapPack;
    typedef typename ComplexField::scalar_type Type;
    typedef std::array<std::vector<std::set<unsigned int>>, 3> DilutionMap;
    typedef Eigen::TensorMap<Eigen::Tensor<Type, 3, Eigen::RowMajor>> NoiseType;
    enum Index {t = 0, l = 1, s = 2};
public:
    DistillationNoise(GridCartesian *g, GridCartesian *g3d, const LapPack &pack, 
                      const unsigned int nNoise = 1);
    virtual ~DistillationNoise(void) = default;
    std::vector<Vector<Type>> & getNoise(void);
    const std::vector<Vector<Type>> & getNoise(void) const;
    unsigned int dilutionIndex(const unsigned t, const unsigned l, const unsigned s) const;
    std::array<unsigned int, 3> dilutionCoordinates(const unsigned int d) const;
    std::vector<unsigned int> dilutionPartition(const Index p, const unsigned i);
    const FermionField & makeSource(const unsigned int d, const unsigned int i);
    // access
    virtual void resize(const int nNoise);
    virtual int  size(void) const;
    virtual int  dilutionSize(void) const;
    virtual int  dilutionSize(const Index ind) const = 0;
    GridCartesian *getGrid(void) const;
    unsigned int getNt(void) const;
    unsigned int getNl(void) const;
    unsigned int getNs(void) const;
    // dump dilution map
    void dumpDilutionMap(void);
    // generate noise
    void generateNoise(GridSerialRNG &rng);
    virtual std::vector<std::string> generateHash(void);
    void save(const std::string filename, const std::string distilname, const unsigned int traj);
protected:
    virtual void buildMap(void) = 0;
    DilutionMap  &getMap(const bool createIfEmpty = true);
    bool         mapEmpty(void) const;
    void load(const std::string filestem, const std::string distilname);
protected:
    DilutionMap                    map_;
    GridCartesian                  *grid_, *grid3d_;
    size_t                         noiseSize_;
private:
    const LapPack                  &pack_;
    FermionField                   src_, tmp3d_, tmp4d_;
    ColourVectorField              evec3d_;
    std::vector<Vector<Type>>      noise_;
};

template <typename FImpl>
DistillationNoise<FImpl>::DistillationNoise(GridCartesian *g, 
                                            GridCartesian *g3d,
                                            const LapPack &pack, 
                                            const unsigned int nNoise)
: DilutedNoise()
, grid_(g), grid3d_(g3d), pack_(pack)
, src_(grid_), tmp3d_(grid3d_), tmp4d_(grid_), evec3d_(grid3d_)
{
    noiseSize_ = getNt()*getNl()*getNs();
    resize(nNoise);
}

template <typename FImpl>
std::vector<Vector<typename DistillationNoise<FImpl>::Type>> &
DistillationNoise<FImpl>::getNoise()
{
    return noise_;
}

template <typename FImpl>
const std::vector<Vector<typename DistillationNoise<FImpl>::Type>> &
DistillationNoise<FImpl>::getNoise() const
{
    return noise_;
}

template <typename FImpl>
unsigned int DistillationNoise<FImpl>::dilutionIndex(const unsigned t,
                                                     const unsigned l,
                                                     const unsigned s) const
{
    unsigned int nl = dilutionSize(Index::l), ns = dilutionSize(Index::s);

    return s + ns*(l + nl*t);
}

template <typename FImpl>
std::array<unsigned int, 3> DistillationNoise<FImpl>::dilutionCoordinates(const unsigned int d) const
{
    std::array<unsigned int, 3> c;
    unsigned int                nl = dilutionSize(Index::l), ns = dilutionSize(Index::s), buf;

    buf         = d;
    c[Index::s] = buf%ns;
    buf         = (buf - c[Index::s])/ns;
    c[Index::l] = buf%nl;
    c[Index::t] = (buf - c[Index::l])/nl;

    return c;
}

template <typename FImpl>
typename DistillationNoise<FImpl>::DilutionMap & DistillationNoise<FImpl>::getMap(const bool createIfEmpty)
{
    if (mapEmpty() and createIfEmpty)
    {
        buildMap();
    }

    return map_;
}

template <typename FImpl>
std::vector<unsigned int> DistillationNoise<FImpl>::dilutionPartition(const Index p, const unsigned ip)
{
    std::vector<unsigned int> s;
    DilutionMap               &map = getMap();
    for (auto i: map[p][ip])
    {
        s.push_back(i);
    }

    return s;
}

template <typename FImpl>
const typename DistillationNoise<FImpl>::FermionField & 
DistillationNoise<FImpl>::makeSource(const unsigned int d, const unsigned int i)
{
    const int   tDir    = grid_->Nd() - 1;
    const int   nt      = grid_->GlobalDimensions()[tDir];
    const int   ntLocal = grid_->LocalDimensions()[tDir];
    const int   tFirst  = grid_->LocalStarts()[tDir];
    NoiseType   noise(noise_[i].data(), nt, pack_.eval.size(), Ns);
    auto        c = dilutionCoordinates(d);
    std::string cstr;
    DilutionMap &map = getMap();

    for (auto j: c)
    {
        cstr += std::to_string(j) + " ";
    }
    cstr.pop_back();
    // LOG(Message) << "Making distillation source for dilution index " << d    //do we want to keep these messages?
                //  << " ~ (" << cstr << ")" << std::endl;
    src_ = Zero();
    for (int it: map[Index::t][c[Index::t]])
    { 
        if(it >= tFirst && it < tFirst + ntLocal)
        {
            for (int ik: map[Index::l][c[Index::l]])
            {
                for (int is: map[Index::s][c[Index::s]])
                {
                    // LOG(Message) << it << " " << ik << " " << is << std::endl;
                    ExtractSliceLocal(evec3d_, pack_.evec[ik], 0, it - tFirst, tDir);
                    // LOG(Message) << noise(it, ik, is) << std::endl;
                    evec3d_ = evec3d_*noise(it, ik, is);
                    tmp3d_  = Zero();
                    pokeSpin(tmp3d_, evec3d_, is);
                    tmp4d_ = Zero();
                    InsertSliceLocal(tmp3d_, tmp4d_, 0, it - tFirst, tDir);
                    src_ += tmp4d_;
                }
            }
        }
    }

    return src_;
}

template <typename FImpl>
void DistillationNoise<FImpl>::resize(const int nNoise)
{
    noise_.resize(nNoise, Vector<Type>(noiseSize_));
}

template <typename FImpl>
int DistillationNoise<FImpl>::size(void) const
{
    return noise_.size();
}

template <typename FImpl>
int DistillationNoise<FImpl>::dilutionSize(void) const
{
    return dilutionSize(Index::t)*dilutionSize(Index::l)*dilutionSize(Index::s);
}

template <typename FImpl>
GridCartesian * DistillationNoise<FImpl>::getGrid(void) const
{
    return grid_;
}

template <typename FImpl>
unsigned int DistillationNoise<FImpl>::getNt(void) const
{
    return grid_->GlobalDimensions()[grid_->Nd() - 1];
}

template <typename FImpl>
unsigned int DistillationNoise<FImpl>::getNl(void) const
{
    return pack_.eval.size();
}

template <typename FImpl>
unsigned int DistillationNoise<FImpl>::getNs(void) const
{
    return Ns;
}

template <typename FImpl>
void DistillationNoise<FImpl>::dumpDilutionMap(void)
{
    DilutionMap &map = getMap();

    auto dump = [this,&map](const unsigned int index)
    {
        for (unsigned int i = 0; i < map[index].size(); ++i)
        {
            std::string s;

            for (auto t: map[index][i])
            {
                s += std::to_string(t) + " "; 
            }
            s.pop_back();
            LOG(Message) << "  " << i << ": {" << s << "}" << std::endl;
        }
    };

    LOG(Message) << "Time index sets:" << std::endl;
    dump(Index::t);
    LOG(Message) << "Laplacian index sets:" << std::endl;
    dump(Index::l);
    LOG(Message) << "Spin index sets:" << std::endl;
    dump(Index::s);
}

template <typename FImpl>
void DistillationNoise<FImpl>::generateNoise(GridSerialRNG &rng)
{
    const Type shift(1., 1.);
    const double invSqrt2 = 0.7071067812;
    Type      eta;

    for (auto &n: noise_)
    for (unsigned int i = 0; i < n.size(); ++i)
    {
        bernoulli(rng, eta);
        n[i] = (2.*eta - shift)*invSqrt2;
    }
}

template <typename FImpl>
bool DistillationNoise<FImpl>::mapEmpty(void) const
{
    bool empty = false;

    for (auto &m: map_)
    {
        empty = empty or m.empty();
    }

    return empty;
}

template <typename FImpl>
std::vector<std::string> DistillationNoise<FImpl>::generateHash(void)
{
    if (mapEmpty())
    {
        buildMap();
    }

    // noise hash by hit
    std::vector<std::string> shash(noise_.size());
    for(unsigned int i=0 ; i<noise_.size() ; i++)
    {
        auto hn = GridChecksum::sha256( &noise_[i].front() , sizeof(Type)*noiseSize_ );
        shash[i] = GridChecksum::sha256_string(hn);
    }
    return shash;
}

template <typename FImpl>
void DistillationNoise<FImpl>::save(const std::string filename, const std::string distilname, const unsigned int traj)
{
#ifdef HAVE_HDF5
    if(grid_->IsBoss())
    {
        std::string filestem = filename + "." + std::to_string(traj) + ".h5";
        Hdf5Writer writer(filestem);
        push(writer, distilname);
        auto &group = writer.getGroup();

        ////metadata write
        //hash hits
        std::vector<std::string> shash = generateHash();
        hsize_t         hcount=1, hstride=1, hblock=1, hash_dim = shash.size(), hash_memdim = 1;
        H5NS::DataSpace hash_dataspace(1, &hash_dim), hash_memspace(1, &hash_memdim);
        //maximum string size assuming each unsigned int is hexed to a maximum of 8 chars (ffffffff); better way?
        H5NS::DataType  hash_datatype( H5T_STRING, 8*SHA256_DIGEST_LENGTH );                  
        H5NS::DataSet   hash_dataset = group.createDataSet("NoiseHitHashes", hash_datatype, hash_dataspace);
        for(hsize_t hoffset=0 ; hoffset<shash.size() ; hoffset++)
        {
            hash_dataspace.selectHyperslab(H5S_SELECT_SET, &hcount, &hoffset, &hstride, &hblock);
            hash_dataset.write( shash[hoffset].c_str(), hash_datatype, hash_memspace, hash_dataspace);
        }

        //dilution sizes t,l,s
        std::vector<int> dil_sizes = {dilutionSize(Index::t),dilutionSize(Index::l),dilutionSize(Index::s)};
        hsize_t dilsize_dim = dil_sizes.size();
        H5NS::DataSpace dil_dataspace(1, &dilsize_dim);
        H5NS::Attribute attr_dil = group.createAttribute("DilutionSizes",  Hdf5Type<int>::type(), dil_dataspace);
        attr_dil.write(Hdf5Type<int>::type(), dil_sizes.data());

        //noise hits
        std::vector<hsize_t>    noffset={0,0}, ncount={1,static_cast<hsize_t>(noiseSize_)}, nstride={1,1}, nblock={1,1};
        std::vector<hsize_t>    noise_dim = {static_cast<hsize_t>(noise_.size()), static_cast<hsize_t>(noiseSize_)},
                                    noise_memdim = {1, static_cast<hsize_t>(noiseSize_)};
        H5NS::DataSpace         noise_dataspace(noise_dim.size(), noise_dim.data()),
                                    noise_memspace(noise_memdim.size(), noise_memdim.data());
        H5NS::DataSet           noise_dataset = group.createDataSet("NoiseHits", Hdf5Type<Type>::type(), noise_dataspace);
        for(unsigned int i=0 ; i<noise_.size() ; i++)
        {
            noffset.front() = static_cast<hsize_t>(i);
            noise_dataspace.selectHyperslab(H5S_SELECT_SET, ncount.data(), noffset.data(), nstride.data(), nblock.data());
            noise_dataset.write(&noise_[i].front() , Hdf5Type<Type>::type(), noise_memspace, noise_dataspace);
        }

    }
#else
    HADRONS_ERROR(Implementation, "distillation I/O needs HDF5 library");
#endif
}

template <typename FImpl>
void DistillationNoise<FImpl>::load(const std::string filestem, const std::string distilname)
{
#ifdef HAVE_HDF5
    if(this->grid_->IsBoss())
    {
        Hdf5Reader reader(filestem);
        push(reader, distilname);
        auto &group = reader.getGroup();

        //noise
        std::vector<hsize_t>    memdim = {1, static_cast<hsize_t>(noiseSize_)},
                                    offset={0,0}, count={1,static_cast<hsize_t>(noiseSize_)}, stride={1,1}, block={1,1};
        H5NS::DataSet dataset =  group.openDataSet("NoiseHits");
        H5NS::DataSpace dataspace = dataset.getSpace(),
                            memspace(memdim.size(), memdim.data());
        std::vector<hsize_t> fdim(2); 
        dataspace.getSimpleExtentDims(fdim.data());
        resize(static_cast<int>(fdim[0]));                            //resize noise_ to in file dimension
        for(unsigned int i=0 ; i<noise_.size() ; i++)
        {
            offset.front() = static_cast<hsize_t>(i);
            dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(), stride.data(), block.data());
            dataset.read(&noise_[i].front(), Hdf5Type<Type>::type(), memspace, dataspace);
        }
    }
#else
    HADRONS_ERROR(Implementation, "distillation I/O needs HDF5 library");
#endif
}

/******************************************************************************
 *                 Container for interlaced distillation noise                *
 ******************************************************************************/
template <typename FImpl>
class InterlacedDistillationNoise: public DistillationNoise<FImpl>
{
public:
    typedef typename DistillationNoise<FImpl>::Index Index;
    typedef typename DistillationNoise<FImpl>::LapPack LapPack;
public:
    InterlacedDistillationNoise(GridCartesian *g, GridCartesian *g3d,
                                const LapPack &pack, const unsigned int ti, 
                                const unsigned int li, const unsigned int si, 
                                const unsigned nNoise);
    unsigned int getInterlacing(const Index ind) const;
    using DistillationNoise<FImpl>::dilutionSize;
    virtual int dilutionSize(const Index ind) const;
    void load(const std::string filename, const std::string distilname, const unsigned int traj);
protected:
    virtual void buildMap(void);
private:
    std::array<unsigned int, 3> interlacing_;
};

template <typename FImpl>
InterlacedDistillationNoise<FImpl>::InterlacedDistillationNoise(GridCartesian *g, 
                                                                GridCartesian *g3d, 
                                                                const LapPack &pack,
                                                                const unsigned int ti, 
                                                                const unsigned int li, 
                                                                const unsigned int si,
                                                                const unsigned nNoise)
: interlacing_({ti, li, si}), DistillationNoise<FImpl>(g, g3d, pack, nNoise)
{}

template <typename FImpl>
unsigned int InterlacedDistillationNoise<FImpl>::getInterlacing(const Index ind) const
{
    return interlacing_[ind];
}

template <typename FImpl>
int InterlacedDistillationNoise<FImpl>::dilutionSize(const Index ind) const
{
    return getInterlacing(ind);
}

template <typename FImpl>
void InterlacedDistillationNoise<FImpl>::buildMap(void)
{
    auto                   &map = this->getMap(false);
    std::set<unsigned int> set;

    map[Index::t].clear();
    for (unsigned int it = 0; it < getInterlacing(Index::t); ++it)
    {
        set.clear();
        for (unsigned int t = it; t < this->getNt(); t += getInterlacing(Index::t))
        {
            set.insert(t);
        }
        map[Index::t].push_back(set);
    }
    map[Index::l].clear();
    for (unsigned int il = 0; il < getInterlacing(Index::l); ++il)
    {
        set.clear();
        for (unsigned int l = il; l < this->getNl(); l += getInterlacing(Index::l))
        {
            set.insert(l);
        }
        map[Index::l].push_back(set);
    }
    map[Index::s].clear();
    for (unsigned int is = 0; is < getInterlacing(Index::s); ++is)
    {
        set.clear();
        for (unsigned int s = is; s < this->getNs(); s += getInterlacing(Index::s))
        {
            set.insert(s);
        }
        map[Index::s].push_back(set);
    }
}

template <typename FImpl>
void InterlacedDistillationNoise<FImpl>::load(const std::string filename, const std::string distilname, const unsigned int traj)
{
    std::string filestem = filename + "." + std::to_string(traj) + ".h5";
    DistillationNoise<FImpl>::load(filestem,distilname);
}

/******************************************************************************
 *                 Container for exact distillation                           *
 ******************************************************************************/
template <typename FImpl>
class ExactDistillationPolicy: public DistillationNoise<FImpl>
{
public:
    typedef typename DistillationNoise<FImpl>::Index Index;
    typedef typename DistillationNoise<FImpl>::LapPack LapPack;
    typedef typename DistillationNoise<FImpl>::Type Type;
public:
    ExactDistillationPolicy(GridCartesian *g, GridCartesian *g3d,
                                const LapPack &pack);
    using DistillationNoise<FImpl>::dilutionSize;
    virtual int dilutionSize(const Index ind) const;
    virtual std::vector<std::string> generateHash(void);
protected:
    virtual void buildMap(void);
};

template <typename FImpl>
ExactDistillationPolicy<FImpl>::ExactDistillationPolicy(GridCartesian *g, 
                                                                GridCartesian *g3d, 
                                                                const LapPack &pack)
: DistillationNoise<FImpl>(g, g3d, pack)
{
    const Type shift(1., 0.);
    for (auto &n: this->getNoise())
    for (unsigned int i = 0; i < n.size(); ++i)
    {
        n[i] = shift;
    }
}

template <typename FImpl>
int ExactDistillationPolicy<FImpl>::dilutionSize(const Index ind) const
{
    switch(ind) {
        case Index::t :
            return this->getNt();
        case Index::l :
            return this->getNl();
        case Index::s :
            return this->getNs();
    }
    HADRONS_ERROR(Argument, "'ind' must be a valid index.");
    return 0;
}

template <typename FImpl>
void ExactDistillationPolicy<FImpl>::buildMap(void)
{
    auto                   &map = this->getMap(false);

    map[Index::t].clear();
    for (unsigned int it = 0; it < dilutionSize(Index::t); ++it)
    {
        map[Index::t].push_back({it});
    }
    map[Index::l].clear();
    for (unsigned int il = 0; il < dilutionSize(Index::l); ++il)
    {
        map[Index::l].push_back({il});
    }
    map[Index::s].clear();
    for (unsigned int is = 0; is < dilutionSize(Index::s); ++is)
    {
        map[Index::s].push_back({is});
    }
}

template <typename FImpl>
std::vector<std::string> ExactDistillationPolicy<FImpl>::generateHash(void)
{
    if (this->mapEmpty())
    {
        this->buildMap();
    }
    
    // exact distillation convention
    std::vector<std::string> shash(1);
    shash.front() = "0";    
    return shash;
}

/******************************************************************************
 *              Abstract container for spin color diagonal noise              *
 ******************************************************************************/
template <typename FImpl>
class SpinColorDiagonalNoise: public DilutedNoise
{
public:
    typedef typename FImpl::FermionField    FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
    typedef typename FImpl::ComplexField    ComplexField;
public:
    // constructor/destructor
    SpinColorDiagonalNoise(GridCartesian *g);
    SpinColorDiagonalNoise(GridCartesian *g, const int nNoise);
    virtual ~SpinColorDiagonalNoise(void) = default;
    // access
    std::vector<ComplexField> &         getNoise(void);
    const std::vector<ComplexField> &   getNoise(void) const;
    FermionField &                      getFerm(const int i);
    PropagatorField &                   getProp(const int i);
    virtual void                        resize(const int nNoise);
    virtual int                         size(void) const;
    int                                 fermSize(void) const;
    virtual int                         dilutionSize(void) const = 0;
    GridCartesian                       *getGrid(void) const;
    // generate noise
    void generateNoise(GridParallelRNG &rng);
private:
    void         setFerm(const int i);
    virtual void setProp(const int i) = 0;
    ComplexField                   eta_;
    FermionField                   ferm_;
    GridCartesian                  *grid_;
    std::vector<ComplexField>      noise_;
    PropagatorField                prop_;
protected:
    ComplexField &    getEta(void);
    FermionField &    getFerm(void);
    int               getNd(void) const;
    int               getNsc(void) const;
    PropagatorField & getProp(void);
    void              setPropagator(ComplexField & eta);
};

template <typename FImpl>
class TimeDilutedNoise: public SpinColorDiagonalNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    TimeDilutedNoise(GridCartesian *g);
    TimeDilutedNoise(GridCartesian *g, const int nNoise);
    virtual ~TimeDilutedNoise(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
    Lattice<iScalar<vInteger>> tLat_;
};

template <typename FImpl>
class FullVolumeNoise: public SpinColorDiagonalNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    FullVolumeNoise(GridCartesian *g, const int nNoise);
    virtual ~FullVolumeNoise(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
};

template <typename FImpl>
class CheckerboardNoise: public SpinColorDiagonalNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    CheckerboardNoise(GridCartesian *g, const int nNoise, const int nSparse);
    virtual ~CheckerboardNoise(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
    int nSparse_, nSrc_ec_;
    LatticeInteger coor_, coorTot_;
};

template <typename FImpl>
class SparseNoise: public SpinColorDiagonalNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    SparseNoise(GridCartesian *g, const int nNoise, const int nSparse);
    virtual ~SparseNoise(void) = default;
    int dilutionSize(void) const;
private:
    void setProp(const int i);
    int nSparse_;
    LatticeInteger coor_;
};
/******************************************************************************
 *               SpinColorDiagonalNoise template implementation               *
 ******************************************************************************/
template <typename FImpl>
SpinColorDiagonalNoise<FImpl>::SpinColorDiagonalNoise(GridCartesian *g)
: DilutedNoise(), grid_(g), ferm_(g), prop_(g), eta_(g)
{}

template <typename FImpl>
SpinColorDiagonalNoise<FImpl>::SpinColorDiagonalNoise(GridCartesian *g,
                                                      const int nNoise)
: SpinColorDiagonalNoise(g)
{
    resize(nNoise);
}

template <typename FImpl>
std::vector<typename SpinColorDiagonalNoise<FImpl>::ComplexField> & 
SpinColorDiagonalNoise<FImpl>::getNoise(void)
{
    return noise_;
}

template <typename FImpl>
const std::vector<typename SpinColorDiagonalNoise<FImpl>::ComplexField> & 
SpinColorDiagonalNoise<FImpl>::getNoise(void) const
{
    return noise_;
}

template <typename FImpl>
void SpinColorDiagonalNoise<FImpl>::setFerm(const int i)
{
    int nc  = FImpl::Dimension;
    std::div_t divs;
    divs = std::div(i, nc);

    PropToFerm<FImpl>(ferm_, prop_, divs.quot, divs.rem);
}

template <typename FImpl>
typename SpinColorDiagonalNoise<FImpl>::FermionField & 
SpinColorDiagonalNoise<FImpl>::getFerm(void)
{
    return ferm_;
}

template <typename FImpl>
typename SpinColorDiagonalNoise<FImpl>::FermionField & 
SpinColorDiagonalNoise<FImpl>::getFerm(const int i)
{
    auto nsc   = this->getNsc();
    std::div_t divs;
    divs = std::div(i, nsc);
    setProp(divs.quot);
    setFerm(divs.rem);
    return getFerm();
}

template <typename FImpl>
void SpinColorDiagonalNoise<FImpl>::setPropagator(ComplexField & eta)
{
    prop_ = 1.;
    prop_ = prop_*eta;
}

template <typename FImpl>
typename SpinColorDiagonalNoise<FImpl>::PropagatorField & 
SpinColorDiagonalNoise<FImpl>::getProp(void)
{
    return prop_;
}

template <typename FImpl>
typename SpinColorDiagonalNoise<FImpl>::PropagatorField & 
SpinColorDiagonalNoise<FImpl>::getProp(const int i)
{
    setProp(i);
    return getProp();
}

template <typename FImpl>
int SpinColorDiagonalNoise<FImpl>::size(void) const
{
    return noise_.size();
}

template <typename FImpl>
int SpinColorDiagonalNoise<FImpl>::fermSize(void) const
{
    return dilutionSize()*getNsc();
}

template <typename FImpl>
typename SpinColorDiagonalNoise<FImpl>::ComplexField & 
SpinColorDiagonalNoise<FImpl>::getEta(void)
{
    return eta_;
}

template <typename FImpl>
int SpinColorDiagonalNoise<FImpl>::getNd(void) const
{
    return grid_->GlobalDimensions().size();
}

template <typename FImpl>
int SpinColorDiagonalNoise<FImpl>::getNsc(void) const
{
    return Ns*FImpl::Dimension;
}

template <typename FImpl>
void SpinColorDiagonalNoise<FImpl>::resize(const int nNoise)
{
    noise_.resize(nNoise, grid_);
}

template <typename FImpl>
GridCartesian * SpinColorDiagonalNoise<FImpl>::getGrid(void) const
{
    return grid_;
}

template <typename FImpl>
void SpinColorDiagonalNoise<FImpl>::generateNoise(GridParallelRNG &rng)
{
    Complex        shift(1., 1.);
    for (int n = 0; n < noise_.size(); ++n)
    {
        bernoulli(rng, eta_);
        eta_ = (2.*eta_ - shift)*(1./::sqrt(2.));
        noise_[n] = eta_;
    }
}

/******************************************************************************
 *                  TimeDilutedNoise template implementation                  *
 ******************************************************************************/
template <typename FImpl>
TimeDilutedNoise<FImpl>::
TimeDilutedNoise(GridCartesian *g, int nNoise)
: SpinColorDiagonalNoise<FImpl>(g, nNoise), tLat_(g)
{}

template <typename FImpl>
int TimeDilutedNoise<FImpl>::dilutionSize() const
{
    auto nt = this->getGrid()->GlobalDimensions()[Tp];
    return nt*this->size();
}

template <typename FImpl>
void TimeDilutedNoise<FImpl>::setProp(const int i)
{
    auto eta   = this->getEta();
    auto noise = this->getNoise();
    auto nd    = this->getNd();
    auto nt    = this->getGrid()->GlobalDimensions()[Tp];

    LatticeCoordinate(tLat_, nd - 1);

    std::div_t divs = std::div(i, nt);
    int t = divs.rem;

    eta = where((tLat_ == t), noise[divs.quot], 0.*noise[divs.quot]);
    this->setPropagator(eta);
}

/******************************************************************************
 *                   FullVolumeNoise template implementation                  *
 ******************************************************************************/
template <typename FImpl>
FullVolumeNoise<FImpl>::
FullVolumeNoise(GridCartesian *g, int nNoise)
: SpinColorDiagonalNoise<FImpl>(g, nNoise)
{}

template <typename FImpl>
int FullVolumeNoise<FImpl>::dilutionSize() const
{
    return this->size();
}

template <typename FImpl>
void FullVolumeNoise<FImpl>::setProp(const int i)
{
    auto noise = this->getNoise();
    this->setPropagator(noise[i]);
}

/******************************************************************************
 *                CheckerboardNoise template implementation                   *
 ******************************************************************************/
template <typename FImpl>
CheckerboardNoise<FImpl>::
CheckerboardNoise(GridCartesian *g, int nNoise, int nSparse)
: SpinColorDiagonalNoise<FImpl>(g, nNoise), nSparse_(nSparse),
  coor_(g), coorTot_(g)
{
    if(nNoise%nSparse_==0)
    {
         nSrc_ec_ = nNoise/nSparse_;
    }
    else
    {
         nSrc_ec_ = (nNoise - nNoise%nSparse_)/nSparse_;
    }
}

template <typename FImpl>
int CheckerboardNoise<FImpl>::dilutionSize() const
{
    return this->size();
}

template <typename FImpl>
void CheckerboardNoise<FImpl>::setProp(const int i)
{
    auto eta   = this->getEta();
    auto nd    = this->getNd();
    auto noise = this->getNoise();
    auto nsc   = this->getNsc();
    unsigned int j;
    eta = noise[i];
    j   = i/nSrc_ec_;

    coorTot_ = 0.;
    for(int d = 0; d < nd; ++d) 
    {
        LatticeCoordinate(coor_, d);
        coorTot_ = coorTot_ + coor_;
    }
    coor_ = j;
    coorTot_ = coorTot_ + coor_;
    eta = where(mod(coorTot_,nSparse_), 0.*eta, eta);
    eta *= sqrt(1./nSrc_ec_);
    this->setPropagator(eta);
}

/******************************************************************************
 *                SparseNoise template implementation                   *
 ******************************************************************************/
template <typename FImpl>
SparseNoise<FImpl>::
SparseNoise(GridCartesian *g, int nNoise, int nSparse)
: SpinColorDiagonalNoise<FImpl>(g, nNoise), nSparse_(nSparse), coor_(g)
{}

template <typename FImpl>
int SparseNoise<FImpl>::dilutionSize() const
{
    auto nd  = this->getNd();
    return this->size()*pow(nSparse_, nd);
}

template <typename FImpl>
void SparseNoise<FImpl>::setProp(const int i)
{
    auto eta   = this->getEta();
    auto nd    = this->getNd();
    auto noise = this->getNoise();

    std::div_t divs = std::div(i, pow(nSparse_, nd));
    eta = noise[divs.quot];
    for(int d = 0; d < nd; ++d) 
    {
        LatticeCoordinate(coor_, d);
        eta = where(mod(coor_,nSparse_), 0.*eta, eta);
    }

    for (int d = 0; d < nd; ++d)
    {
        divs = std::div(divs.rem, pow(nSparse_, nd-(d+1)));
        eta = Cshift(eta, d, divs.quot);
    }
    this->setPropagator(eta);
}

END_HADRONS_NAMESPACE

#endif // Hadrons_DilutedNoise_hpp_
