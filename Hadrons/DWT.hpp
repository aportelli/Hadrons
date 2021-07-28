#ifndef GRID_DWT_H_
#define GRID_DWT_H_

#include <Hadrons/Global.hpp>
#include <Hadrons/DWT_filters.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Container class for multi-level wavelet transforms            *
 ******************************************************************************/

// With LatticeWt<Field> w, a component of the wavelet tranform is given
// by w(l, i), where l is the level and i is the component index for that
// level. 0 <= i <= 2^dim - 1 and the binary representation of i gives the mode
// filtering (0: low, 1: high) in each direction.
template<typename Field>
class LatticeWt
{
public:
    typedef std::vector<std::vector<Field>> Wt;
    typedef std::vector<Field>              WtLevel;
public:
    // constructors
    LatticeWt(void) = default;
    LatticeWt(std::vector<GridBase *> gLevel, const unsigned int level = 1);
    // resize
    void resize(std::vector<GridBase *> gLevel, const unsigned int level = 1);
    // access
    const std::vector<GridBase *> &getGrid(void);
    GridBase *                     getGrid(const unsigned int l);
    Wt &            operator()(void);
    const Wt &      operator()(void) const;
    WtLevel &       operator()(const unsigned int l);
    const WtLevel & operator()(const unsigned int l) const;
    Field &         operator()(const unsigned int l, const unsigned int i);
    const Field &   operator()(const unsigned int l, const unsigned int i) const;
private:
    unsigned int            level_;
    std::vector<GridBase *> grid_;
    Wt                      wt_;
};

/******************************************************************************
 *               Class to operate the discrete wavelet transform              *
 ******************************************************************************/
// Base abstract class /////////////////////////////////////////////////////////
template <typename Field>
class DwtBase
{
public:
    // constructor
    DwtBase(std::vector<GridBase *> gLevel, const DwtFilter &filter, 
            const unsigned int level = 1);
    // access
    unsigned int      getLevel(void) const;
    const DwtFilter & getFilter(void) const;
    GridBase *        getFineGrid(void) const;
    GridBase *        getCoarseGrid(const unsigned int l) const;
    // downsampling/upsampling
    static void downsample(Field &out, const Field &in, 
                           const unsigned int dir, const unsigned int n = 2);
    static void upsample(Field &out, const Field &in, 
                         const unsigned int dir, const unsigned int n = 2);
    // DWT interface
    //// single-level forward
    void forward(typename LatticeWt<Field>::WtLevel &out, const Field &in);
    //// multi-level forward
    void forward(LatticeWt<Field> &out, const Field &in);
    //// single-level backward
    void backward(Field &out, const typename LatticeWt<Field>::WtLevel &in);
    //// multi-level backward
    void backward(Field &out, const LatticeWt<Field> &in);
    // grid utilities
    static GridBase * downsampleGrid(GridBase * g, 
                                     const std::vector<unsigned int> n);
    static GridBase * downsampleGrid(GridBase * g, const unsigned int dir, 
                                     const unsigned int n = 2);
    static GridBase * upsampleGrid(GridBase * g, 
                                   const std::vector<unsigned int> n);
    static GridBase * upsampleGrid(GridBase * g, const unsigned int dir, 
                                   const unsigned int n = 2);
    static GridBase * coarsenGrid(GridBase * g, const unsigned int level = 1);
private:
    // transform elementary step, component i in direction mu
    virtual void forwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                             const unsigned int mu, const unsigned int i, 
                             const unsigned l) = 0;
    virtual void backwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                              const unsigned int mu, const unsigned int i, 
                              const unsigned l) = 0;
    // create all the intermediate grids and buffers
    void makeGrids(void);
    void makeBuffers(void);
protected:
    // get buffer grid
    GridBase *   grid(const unsigned int l, const unsigned int dir);
    // figure out level by comparing to the fine grid
    unsigned int levelFromGrid(GridBase *g) const;
private:
    std::vector<GridBase *>                             grid_;
    DwtFilter                                           filter_;
    unsigned int                                        nd_, level_;
    std::vector<std::vector<std::shared_ptr<GridBase>>> bufferGrid_;
    std::vector<std::vector<std::vector<Field>>>        buffer_;
};

// Standard DWT ////////////////////////////////////////////////////////////////
template <typename Field>
class Dwt: public DwtBase<Field>
{
public:
    // reuse base class constructor
    using DwtBase<Field>::DwtBase;
    // filter convolution
    void filterConvolution(Field &out, const Field &in, 
                           const std::vector<Real> &filter,
                           const unsigned int dir, const int offset);
    // transform
    void forwardOp(Field &out, const Field &in, 
                     const std::vector<Real> &filter, const unsigned int dir);
    void backwardOp(Field &out, const Field &in, 
                    const std::vector<Real> &filter, const unsigned int dir);
private:
    virtual void forwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                             const unsigned int mu, const unsigned int i, const unsigned l);
    virtual void backwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                              const unsigned int mu, const unsigned int i, const unsigned l);
};

// Gauge-covariant DWT /////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
class CovariantDwt: public DwtBase<Field>
{
public:
    typedef typename GImpl::LinkField GaugeLinkField;
    typedef typename GImpl::Field     GaugeField;
public:
    // constructor
    CovariantDwt(std::vector<GridBase *> gLevel, const DwtFilter &filter, 
                 const GaugeField &U, const unsigned int level = 1);
    // filter convolution
    void filterConvolution(Field &out, const Field &in, const GaugeLinkField &U,
                           const std::vector<Real> &filter,
                           const unsigned int dir, const int offset);
    // transform
    void forwardOp(Field &out, const Field &in, const GaugeLinkField &Umu,
                   const std::vector<Real> &filter, const unsigned int dir);
    void backwardOp(Field &out, const Field &in, const GaugeLinkField &Umu,
                    const std::vector<Real> &filter, const unsigned int dir);
    // block gauge in direction mu
    static void gaugeBlock(GaugeField &out, const GaugeField &in, const unsigned int mu);
private:
    virtual void forwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                             const unsigned int mu, const unsigned int i, const unsigned l);
    virtual void backwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                              const unsigned int mu, const unsigned int i, const unsigned l);
    // make gauge buffers
    void makeGaugeBuffers(void);
    void fillGaugeBuffers(void);
private:
    bool                                 bufferFilled_{false};
    const GaugeField                     &fineGauge_;
    std::vector<std::vector<GaugeField>> gaugeBuffer_;
};

/******************************************************************************
 *                      LatticeWt template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
LatticeWt<Field>::LatticeWt(std::vector<GridBase *> gLevel, 
                            const unsigned int level)
{
    resize(gLevel, level);
}

// resize //////////////////////////////////////////////////////////////////////
template <typename Field>
void LatticeWt<Field>::resize(std::vector<GridBase *> gLevel, 
                              const unsigned int level)
{
    unsigned int cFactor = 1;

    assert(gLevel.size() >= level + 1);
    grid_  = gLevel;
    level_ = level;
    for (unsigned int mu = 0; mu < grid_[0]->Nd(); ++mu)
    {
        cFactor *= 2;
    }
    wt_.resize(level_);
    for (unsigned int l = 0; l < level_; ++l)
    {
        wt_[l].clear();
        wt_[l].resize(cFactor, grid_[l + 1]);
    }
}

// access //////////////////////////////////////////////////////////////////////
template <typename Field>
const std::vector<GridBase *> & LatticeWt<Field>::getGrid(void)
{
    return grid_;
}

template <typename Field>
GridBase * LatticeWt<Field>::getGrid(const unsigned int l)
{
    return grid_[l];
}

template <typename Field>
typename LatticeWt<Field>::Wt & LatticeWt<Field>::operator()(void)
{
    return wt_;
}

template <typename Field>
const typename LatticeWt<Field>::Wt & LatticeWt<Field>::operator()(void) const
{
    return wt_;
}

template <typename Field>
typename LatticeWt<Field>::WtLevel & LatticeWt<Field>::operator()(const unsigned int l)
{
    return wt_[l];
}

template <typename Field>
const typename LatticeWt<Field>::WtLevel & LatticeWt<Field>::operator()(const unsigned int l) const
{
    return wt_[l];
}

template <typename Field>
Field & LatticeWt<Field>::operator()(const unsigned int l, const unsigned int i)
{
    return wt_[l][i];
}

template <typename Field>
const Field & LatticeWt<Field>::operator()(const unsigned int l, const unsigned int i) const
{
    return wt_[l][i];
}

/******************************************************************************
 *                        DwtBase template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
DwtBase<Field>::DwtBase(std::vector<GridBase *> gLevel, 
                        const DwtFilter &filter, const unsigned int level)
: grid_(gLevel), filter_(filter), level_(level)
{
    nd_ = grid_[0]->Nd();
    makeGrids();
    makeBuffers();
}

// create all the intermediate grids and buffers ///////////////////////////////
template <typename Field>
void DwtBase<Field>::makeGrids(void)
{
    bufferGrid_.clear();
    bufferGrid_.resize(getLevel());
    for (unsigned int l = 0; l < getLevel(); ++l)
    {
        GridBase *gPt = grid_[l];

        for (unsigned int mu = 0; mu < nd_ - 1; ++mu)
        {
            bufferGrid_[l].emplace_back(downsampleGrid(gPt, mu));
            gPt = bufferGrid_[l].back().get();
        }
    }
}

template <typename Field>
void DwtBase<Field>::makeBuffers(void)
{
    buffer_.clear();
    buffer_.resize(getLevel());
    for (unsigned int l = 0; l < getLevel(); ++l)
    {
        unsigned int outSize = 2;

        buffer_[l].resize(nd_ - 1);
        for (unsigned int mu = 0; mu < nd_ - 1; ++mu)
        {
            buffer_[l][mu].resize(outSize, grid(l, mu));
            outSize *= 2;
        }
    }
}

// access //////////////////////////////////////////////////////////////////////
template <typename Field>
unsigned int DwtBase<Field>::getLevel(void) const
{
    return level_;
}

template <typename Field>
const DwtFilter & DwtBase<Field>::getFilter(void) const
{
    return filter_;
}

template <typename Field>
GridBase * DwtBase<Field>::getFineGrid(void) const
{
    return grid_[0];
}

template <typename Field>
GridBase * DwtBase<Field>::getCoarseGrid(const unsigned int l) const
{
    return grid_[l + 1];
}

// get buffer grid /////////////////////////////////////////////////////////////
template <typename Field>
GridBase * DwtBase<Field>::grid(const unsigned int l, const unsigned int dir)
{
    return (dir == nd_ - 1) ? grid_[l + 1] : bufferGrid_[l][dir].get();
}

// figure out level by comparing to the fine grid //////////////////////////////
template <typename Field>
unsigned int DwtBase<Field>::levelFromGrid(GridBase *g) const
{
    unsigned int ratio = grid_[0]->FullDimensions()[0]/g->FullDimensions()[0];
    unsigned int level = 0;

    for (unsigned int mu = 1; mu < nd_; ++mu)
    {
        assert(grid_[0]->FullDimensions()[mu]/g->FullDimensions()[mu] == ratio);
    }
    while (ratio >>= 1) ++level;

    return level;
}

// downsampling/upsampling /////////////////////////////////////////////////////
template <typename Field>
void DwtBase<Field>::downsample(Field &out, const Field &in, 
                                const unsigned int dir, const unsigned int n)
{
    GridBase                                     *dsGrid = out.Grid();
    Grid::Coordinate                             x;
    typename Field::vector_object::scalar_object buf;

    for(int idx=0; idx < dsGrid->lSites(); idx++) 
    {
        dsGrid->LocalIndexToLocalCoor(idx, x);
        x[dir] *= n;
        peekLocalSite(buf, in, x);
        x[dir] /= n;
        pokeLocalSite(buf, out, x);
    }
}

template <typename Field>
void DwtBase<Field>::upsample(Field &out, const Field &in, 
                              const unsigned int dir, const unsigned int n)
{
    GridBase                                     *dsGrid = in.Grid();
    Grid::Coordinate                             x;
    typename Field::vector_object::scalar_object buf;

    out = Zero();
    for(int idx=0; idx < dsGrid->lSites(); idx++) 
    {
        dsGrid->LocalIndexToLocalCoor(idx, x);
        peekLocalSite(buf, in, x);
        x[dir] *= n;
        pokeLocalSite(buf, out, x);
    }
}

// DWT interface ///////////////////////////////////////////////////////////////
// single level forward
template <typename Field>
void DwtBase<Field>::forward(typename LatticeWt<Field>::WtLevel &out, 
                             const Field &in)
{
    GridBase *         gFine = in.Grid();
    unsigned int       nd = gFine->Nd(), l = levelFromGrid(gFine), outSize;
    std::vector<Field> start, *prev, *current;

    outSize = 1;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        outSize *= 2;
    }
    assert(out.size() == outSize);
    start.push_back(in);
    prev    = &start; 
    current = (nd == 1) ? &out : &buffer_[l][0];
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        std::cout << GridLogIterative << "forward  DWT level " << l 
                  << " direction " << mu << " -> " << current->size() << " " 
                  << (*current)[0].Grid()->FullDimensions()
                  << " lattices" << std::endl;
        for (unsigned int i = 0; i < prev->size(); ++i)
        {
            forwardStep(*current, *prev, mu, i, l);
        }
        prev    = current;
        current = (mu == nd - 2) ? &out : &buffer_[l][mu + 1];
    }
}

// multi-level forward
template <typename Field>
void DwtBase<Field>::forward(LatticeWt<Field> &out, const Field &in)
{
    const Field *pt = &in;

    out.resize(grid_, level_);
    for (unsigned int l = 0; l < level_; ++l)
    {
        forward(out(l), *pt);
        pt = &(out(l, 0));
    }
}

// single-level backward
template <typename Field>
void DwtBase<Field>::backward(Field &out, const typename LatticeWt<Field>::WtLevel &in)
{
    GridBase           *gFine = out.Grid();
    unsigned int       nd = gFine->Nd(), l = levelFromGrid(gFine), inSize;
    std::vector<Field> end(1, out.Grid()), *prev, *current;

    inSize = 1;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        inSize *= 2;
    }
    assert(in.size() == inSize);
    prev    = const_cast<std::vector<Field> *>(&in);
    current = (nd == 1) ? &end : &buffer_[l][nd - 2];
    for (unsigned int mup1 = nd; mup1 > 0; --mup1)
    {
        unsigned int mu = mup1 - 1;

        std::cout << GridLogIterative << "backward DWT level " << l 
                  << " direction " << mu << " -> " << current->size() << " " 
                  << (*current)[0].Grid()->FullDimensions()
                  << " lattices" << std::endl;
        for (unsigned int i = 0; i < current->size(); ++i)
        {
            backwardStep(*current, *prev, mu, i, l);
        }
        prev    = current;
        current = (mu == 1) ? &end : &buffer_[l][mu - 2];
    }
    out = end[0];
}

// multi-level backward
template <typename Field>
void DwtBase<Field>::backward(Field &out, const LatticeWt<Field> &in)
{
    LatticeWt<Field> tmp = in;

    assert(level_ == tmp().size());
    for (unsigned int l = level_ - 1; l > 0; --l)
    {
        backward(tmp(l - 1, 0), tmp(l));
    }
    backward(out, tmp(0));
}

// grid utilities //////////////////////////////////////////////////////////////
template <typename Field>
GridBase * DwtBase<Field>::downsampleGrid(GridBase * g, 
                                          const std::vector<unsigned int> n)
{
    auto dim = g->FullDimensions();

    assert(dim.size() == n.size());
    for (unsigned int mu = 0; mu < n.size(); ++mu)
    {
        assert(dim[mu] % n[mu] == 0);
        dim[mu] /= n[mu];
    }
    return new GridCartesian(dim, g->_simd_layout, g->_processors);
}

template <typename Field>
GridBase * DwtBase<Field>::downsampleGrid(GridBase * g, const unsigned int dir, 
                                          const unsigned int n)
{
    std::vector<unsigned int> nVec(g->Nd(), 1);

    assert(dir < g->Nd());
    nVec[dir] = n;

    return downsampleGrid(g, nVec);
}

template <typename Field>
GridBase * DwtBase<Field>::upsampleGrid(GridBase * g, 
                                        const std::vector<unsigned int> n)
{
    auto dim = g->FullDimensions();

    assert(dim.size() == n.size());
    for (unsigned int mu = 0; mu < n.size(); ++mu)
    {
        dim[mu] *= n[mu];
    }
    return new GridCartesian(dim, g->_simd_layout, g->_processors);
}

template <typename Field>
GridBase * DwtBase<Field>::upsampleGrid(GridBase * g, const unsigned int dir, 
                                        const unsigned int n)
{
    std::vector<unsigned int> nVec(g->Nd(), 1);

    assert(dir < g->Nd());
    nVec[dir] = n;

    return upsampleGrid(g, nVec);
}

template <typename Field>
GridBase * DwtBase<Field>::coarsenGrid(GridBase * g, const unsigned int level)
{
    assert(level >= 1);

    std::vector<unsigned int> nVec(g->Nd(), 2*level);
    
    return downsampleGrid(g, nVec);
}


/******************************************************************************
 *                        Dwt template implementation                         *
 ******************************************************************************/
// filter convolution //////////////////////////////////////////////////////////
template <typename Field>
void Dwt<Field>::filterConvolution(Field &out, const Field &in, 
                                   const std::vector<Real> &filter,
                                   const unsigned int dir,
                                   const int offset)
{
    assert(in.Grid()->FullDimensions()[dir] >= filter.size());

    Field buf(in.Grid());

    out = Zero();
    buf = in;
    
    for (int i = 0; i >= offset; i--)
    {
        out += filter[i - offset]*buf;
        buf  = Cshift(buf, dir, -1);
    }
    buf = Cshift(in, dir, 1);
    for (int i = 1; i < filter.size() + offset; ++i)
    {
        out += filter[i - offset]*buf;
        buf  = Cshift(buf, dir, 1);
    }
}

// transform ///////////////////////////////////////////////////////////////////
template <typename Field>
void Dwt<Field>::forwardOp(Field &out, const Field &in, 
                           const std::vector<Real> &filter, 
                           const unsigned int dir)
{
    Field buf(in.Grid());
    int   offset = -filter.size()/2;

    assert(out.Grid()->FullDimensions()[dir] == in.Grid()->FullDimensions()[dir]/2);
    filterConvolution(buf, in, filter, dir, offset);
    this->downsample(out, buf, dir);
}

template <typename Field>
void Dwt<Field>::backwardOp(Field &out, const Field &in, 
                            const std::vector<Real> &filter, 
                            const unsigned int dir)
{
    Field buf(out.Grid());
    int   offset = -filter.size()/2;

    assert(out.Grid()->FullDimensions()[dir] == 2*in.Grid()->FullDimensions()[dir]);
    this->upsample(buf, in, dir);
    filterConvolution(out, buf, filter, dir, -offset - filter.size() + 1);
}

template <typename Field>
void Dwt<Field>::forwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                             const unsigned int mu, const unsigned int i, const unsigned l)
{
    forwardOp(out[2*i], in[i], this->getFilter().fwdL, mu);
    forwardOp(out[2*i + 1], in[i], this->getFilter().fwdH, mu);
}

template <typename Field>
void Dwt<Field>::backwardStep(std::vector<Field> &out, const std::vector<Field> &in, 
                              const unsigned int mu, const unsigned int i, const unsigned l)
{
    Field tmp(out[0].Grid());

    out[i] = Zero();
    backwardOp(tmp, in[2*i], this->getFilter().bwdL, mu);
    out[i] += tmp;
    backwardOp(tmp, in[2*i + 1], this->getFilter().bwdH, mu);
    out[i] += tmp;
}

/******************************************************************************
 *                  CovariantDwt template implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
CovariantDwt<Field, GImpl>::CovariantDwt(std::vector<GridBase *> gLevel, 
                                         const DwtFilter &filter, 
                                         const GaugeField &U, 
                                         const unsigned int level)
: DwtBase<Field>(gLevel, filter, level), fineGauge_(U)
{
    makeGaugeBuffers();
}

// make gauge buffers //////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::makeGaugeBuffers(void)
{
    unsigned int nd = this->getFineGrid()->Nd();

    gaugeBuffer_.resize(this->getLevel());
    for (unsigned int l = 0; l < this->getLevel(); ++l)
    {
        unsigned int outSize = 1;

        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            gaugeBuffer_[l].emplace_back(this->grid(l, mu));
        }
    }
}

template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::fillGaugeBuffers(void)
{
    if (!bufferFilled_)
    {
        const GaugeField *pt  = &fineGauge_;
        unsigned int nd = this->getFineGrid()->Nd();

        for (unsigned int l = 0; l < this->getLevel(); ++l)
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            gaugeBlock(gaugeBuffer_[l][mu], *pt, mu);
            pt = &(gaugeBuffer_[l][mu]);
        }
        bufferFilled_ = true;
    }
}

// block gauge in direction mu /////////////////////////////////////////////////
template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::gaugeBlock(GaugeField &out, const GaugeField &in,
                                                 const unsigned int mu)
{
    GaugeLinkField tmpLink(in.Grid());
    GaugeField     tmp(in.Grid());

    tmpLink = peekLorentz(in, mu);
    tmpLink = GImpl::CovShiftForward(tmpLink, mu, tmpLink);
    tmp     = in;
    pokeLorentz(tmp, tmpLink, mu);
    DwtBase<GaugeField>::downsample(out, tmp, mu);
}

// filter convolution //////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::filterConvolution(Field &out, const Field &in,
                                                   const GaugeLinkField &Umu,
                                                   const std::vector<Real> &filter,
                                                   const unsigned int mu, const int offset)
{
    assert(in.Grid()->FullDimensions()[mu] >= filter.size());

    Field buf(in.Grid());

    out = Zero();
    buf = in;
    for (int i = 0; i >= offset; i--)
    {
        out += filter[i - offset]*buf;
        buf  = GImpl::CovShiftBackward(Umu, mu, buf);
    }
    buf = GImpl::CovShiftForward(Umu, mu, in);
    for (int i = 1; i < filter.size() + offset; ++i)
    {
        out += filter[i - offset]*buf;
        buf  = GImpl::CovShiftForward(Umu, mu, buf);
    }
}

// transform ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::forwardOp(Field &out, const Field &in, const GaugeLinkField &Umu,
                                           const std::vector<Real> &filter, const unsigned int dir)
{
    Field buf(in.Grid());
    int   offset = -filter.size()/2;

    assert(out.Grid()->FullDimensions()[dir] == in.Grid()->FullDimensions()[dir]/2);
    filterConvolution(buf, in, Umu, filter, dir, offset);
    this->downsample(out, buf, dir);
}

template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::backwardOp(Field &out, const Field &in, const GaugeLinkField &Umu,
                                            const std::vector<Real> &filter, const unsigned int dir)
{
    Field buf(out.Grid());
    int   offset = -filter.size()/2;

    assert(out.Grid()->FullDimensions()[dir] == 2*in.Grid()->FullDimensions()[dir]);
    this->upsample(buf, in, dir);
    filterConvolution(out, buf, Umu, filter, dir, -offset - filter.size() + 1);
}

template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::forwardStep(std::vector<Field> &out, 
                                             const std::vector<Field> &in, 
                                             const unsigned int mu, 
                                             const unsigned int i, 
                                             const unsigned l)
{
    const GaugeField *Upt;

    fillGaugeBuffers();
    if (mu == 0)
    {
        if (l == 0)
        {
            Upt = &fineGauge_;
        }
        else
        {
            Upt = &(gaugeBuffer_[l - 1].back());
        }
    }
    else
    {
        Upt = &(gaugeBuffer_[l][mu - 1]);
    }

    GaugeLinkField Umu(Upt->Grid());

    Umu = peekLorentz(*Upt, mu);
    forwardOp(out[2*i], in[i], Umu, this->getFilter().fwdL, mu);
    forwardOp(out[2*i + 1], in[i], Umu, this->getFilter().fwdH, mu);
}

template <typename Field, typename GImpl>
void CovariantDwt<Field, GImpl>::backwardStep(std::vector<Field> &out, 
                                              const std::vector<Field> &in, 
                                              const unsigned int mu, 
                                              const unsigned int i, 
                                              const unsigned l)
{
    const GaugeField *Upt;

    fillGaugeBuffers();
    if (mu == 0)
    {
        if (l == 0)
        {
            Upt = &fineGauge_;
        }
        else
        {
            Upt = &(gaugeBuffer_[l - 1].back());
        }
    }
    else
    {
        Upt = &(gaugeBuffer_[l][mu - 1]);
    }

    GaugeLinkField Umu(Upt->Grid());
    Field          tmp(out[0].Grid());

    Umu    = peekLorentz(*Upt, mu);
    out[i] = Zero();
    backwardOp(tmp, in[2*i], Umu, this->getFilter().bwdL, mu);
    out[i] += tmp;
    backwardOp(tmp, in[2*i + 1], Umu, this->getFilter().bwdH, mu);
    out[i] += tmp;
}

END_HADRONS_NAMESPACE

#endif
