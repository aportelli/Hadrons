/*
 * DilutedNoise.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
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
#include <Hadrons/Modules/MDistil/Distil.hpp> // TODO: for 3D grids, shodul be handled globally

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
    typedef typename ComplexField::scalar_object Type;
    typedef std::array<std::vector<std::set<unsigned int>>, 3> DilutionMap;
    typedef Eigen::TensorMap<Eigen::Tensor<Type, 3, Eigen::RowMajor>> NoiseType;
    enum Index {t = 0, l = 1, s = 2};
public:
    DistillationNoise(GridCartesian *g, const LapPack &pack, const unsigned int nNoise = 1);
    virtual ~DistillationNoise(void) = default;
    std::vector<Vector<Type>> & getNoise(void);
    const std::vector<Vector<Type>> & getNoise(void) const;
    unsigned int dilutionIndex(const unsigned t, const unsigned l, const unsigned s) const;
    std::array<unsigned int, 3> dilutionCoordinates(const unsigned int d) const;
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
    // generate noise
    void generateNoise(GridSerialRNG &rng);
protected:
    virtual void buildMap(void) = 0;
protected:
    DilutionMap               map_;
private:
    GridCartesian             *grid_;
    const LapPack             &pack_;
    FermionField              src_;
    std::vector<Vector<Type>> noise_;
    size_t                    noiseSize_;
};

template <typename FImpl>
DistillationNoise<FImpl>::DistillationNoise(GridCartesian *g, const LapPack &pack, const unsigned int nNoise)
: DilutedNoise(), grid_(g), pack_(pack), src_(grid_)
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
const typename DistillationNoise<FImpl>::FermionField & 
DistillationNoise<FImpl>::makeSource(const unsigned int d, const unsigned int i)
{
    std::unique_ptr<GridCartesian> g3d;
    
    MDistil::MakeLowerDimGrid(g3d, grid_);

    const int         tDir    = grid_->Nd() - 1;
    const int         nt      = grid_->GlobalDimensions()[tDir];
    const int         ntLocal = grid_->LocalDimensions()[tDir];
    const int         tFirst  = grid_->LocalStarts()[tDir];
    ColourVectorField evec3d(g3d.get());
    FermionField      tmp3d(grid_), tmp4d(grid_);
    NoiseType         noise(noise_[i], nt, pack_.eval.size(), Ns);
    auto              c = dilutionCoordinates(d);
    
    src_ = Zero();
    for (int it: map_[Index::t][c[Index::t]])
    { 
        if(it >= tFirst && it < tFirst + ntLocal)
        {
            for (int ik: map_[Index::l][c[Index::l]])
            {
                for (int is: map_[Index::s][c[Index::s]])
                {
                    ExtractSliceLocal(evec3d, pack_.evec[ik], 0, it - tFirst, tDir);
                    evec3d = evec3d*noise(it, ik, is);
                    tmp3d  = Zero();
                    pokeSpin(tmp3d, evec3d, is);
                    tmp4d = Zero();
                    InsertSliceLocal(tmp3d, tmp4d, 0, it - tFirst, tDir);
                    src_ += tmp4d;
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
void DistillationNoise<FImpl>::generateNoise(GridSerialRNG &rng)
{
    constexpr Type shift(1., 1.);
    constexpr double invSqrt2 = 0.7071067812;
    Type      eta;

    for (auto &n: noise_)
    for (unsigned int i = 0; i < n.size(); ++i)
    {
        random(rng, eta);
        n[i] = (2.*eta - shift)*invSqrt2;
    }
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
    InterlacedDistillationNoise(GridCartesian *g, const LapPack &pack,
                                const unsigned int ti, const unsigned int li, 
                                const unsigned int si, const unsigned nNoise = 1);
    unsigned int getInterlacing(const Index ind) const;
    virtual int  dilutionSize(const Index ind) const;
protected:
    virtual void buildMap(void);
private:
    std::array<unsigned int, 3> interlacing_;
};
template <typename FImpl>
InterlacedDistillationNoise<FImpl>::InterlacedDistillationNoise(GridCartesian *g, 
                                                                const LapPack &pack,
                                                                const unsigned int ti, 
                                                                const unsigned int li, 
                                                                const unsigned int si,
                                                                const unsigned nNoise)
: interlacing_({ti, li, si}), DistillationNoise<FImpl>(g, pack, nNoise)
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
    auto                   &map = this->map_;
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
