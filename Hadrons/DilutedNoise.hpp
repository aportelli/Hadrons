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

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Abstract container for spin color diagonal noise              *
 ******************************************************************************/
template <typename FImpl>
class SpinColorDiagonalNoise
{
public:
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;
public:
    // constructor/destructor
    SpinColorDiagonalNoise(GridCartesian *g);
    SpinColorDiagonalNoise(GridCartesian *g, const int nNoise);
    virtual ~SpinColorDiagonalNoise(void) = default;
    // access
    std::vector<LatticeComplex> &       getNoise(void);
    const std::vector<LatticeComplex> & getNoise(void) const;
    FermionField &                      getFerm(const int i);
    PropagatorField &                   getProp(const int i);
    void                                resize(const int nNoise);
    int                                 size(void) const;
    int                                 fermSize(void) const;
    virtual int                         dilutionSize(void) const = 0;
    GridCartesian                       *getGrid(void) const;
    // generate noise
    void generateNoise(GridParallelRNG &rng);
private:
    void         setFerm(const int i);
    virtual void setProp(const int i) = 0;
    LatticeComplex                 eta_;
    FermionField                   ferm_;
    GridCartesian                  *grid_;
    std::vector<LatticeComplex>    noise_;
    PropagatorField                prop_;
protected:
    LatticeComplex &  getEta(void);
    FermionField &    getFerm(void);
    int               getNd(void) const;
    int               getNsc(void) const;
    PropagatorField & getProp(void);
    void              setPropagator(LatticeComplex & eta);
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
: grid_(g), ferm_(g), prop_(g), eta_(g)
{}

template <typename FImpl>
SpinColorDiagonalNoise<FImpl>::SpinColorDiagonalNoise(GridCartesian *g,
                                                      const int nNoise)
: SpinColorDiagonalNoise(g)
{
    resize(nNoise);
}

template <typename FImpl>
std::vector<LatticeComplex> & SpinColorDiagonalNoise<FImpl>::
getNoise(void)
{
    return noise_;
}

template <typename FImpl>
const std::vector<LatticeComplex> & SpinColorDiagonalNoise<FImpl>::
getNoise(void) const
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
void SpinColorDiagonalNoise<FImpl>::setPropagator(LatticeComplex & eta)
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
LatticeComplex & SpinColorDiagonalNoise<FImpl>::getEta(void)
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
