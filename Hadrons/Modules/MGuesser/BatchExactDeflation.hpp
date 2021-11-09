#ifndef Hadrons_MGuesser_BatchExactDeflation_hpp_
#define Hadrons_MGuesser_BatchExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MIO/LoadEigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationPar,
                                    MIO::LoadEigenPackPar, eigenPack,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename Pack, typename GImpl>
class TBatchExactDeflation: public Module<BatchExactDeflationPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    GAUGE_TYPE_ALIASES(GImpl, );
public:
    // constructor
    TBatchExactDeflation(const std::string name);
    // destructor
    virtual ~TBatchExactDeflation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflation, ARG(TBatchExactDeflation<FermionEigenPack<FIMPL>, GIMPL>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationF, ARG(TBatchExactDeflation<FermionEigenPack<FIMPLF>, GIMPLF>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename Pack, typename GImpl>
class BatchExactDeflationGuesser: public LinearFunction<typename Pack::Field>
{
public:
    typedef typename Pack::Field Field;
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    BatchExactDeflationGuesser(const MIO::LoadEigenPackPar &epPar, 
                               const unsigned int evBatchSize,
                               const unsigned int sourceBatchSize,
                               GridBase *grid, GridBase *gridIo,
                               const int traj, 
                               const GaugeLinkField *transform = nullptr)
    : epPar_(epPar)
    , evBatchSize_(evBatchSize)
    , sourceBatchSize_(sourceBatchSize)
    , grid_(grid)
    , gridIo_(gridIo)
    , traj_(traj)
    , transform_(transform)
    {
        epack_.init(evBatchSize_, grid_, gridIo_);
    };

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        unsigned int nBatch     = epPar_.size/evBatchSize_ + (((epPar_.size % evBatchSize_) != 0) ? 1 : 0);
        unsigned int sourceSize = out.size();

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        LOG(Message) << "--- zero guesses" << std::endl;
        for (auto &v: out)
        {
            v = Zero();
        }
        for (unsigned int bv = 0; bv < epPar_.size; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epPar_.size - bv, evBatchSize_);

            LOG(Message) << "--- batch " << bv/evBatchSize_ << std::endl;
            LOG(Message) << "I/O" << std::endl;
            epack_.read(epPar_.filestem, epPar_.multiFile, bv, bv + evBlockSize, traj_);
            if (transform_ != nullptr)
            {
                epack_.gaugeTransform(*transform_);
            }
            LOG(Message) << "project" << std::endl;
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                projAccumulate(in, out, evBlockSize, bs, bs + sourceBlockSize);
            }
        }
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    void projAccumulate(const std::vector<Field> &in, std::vector<Field> &out,
                        const unsigned int evBatchSize,
                        const unsigned int si, const unsigned int sf)
    {
        GridBase *g       = in[0].Grid();
        double   lVol     = g->lSites();
        double   siteSize = sizeof(typename Field::scalar_object);
        double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
        double   nIt      = evBatchSize*(sf - si);
        double   t        = 0.;
        
        t -= usecond();
        for (unsigned int i = 0; i < evBatchSize; ++i)
        for (unsigned int j = si; j < sf; ++j)
        {
            axpy(out[j], 
                 TensorRemove(innerProduct(epack_.evec[i], in[j]))/epack_.eval[i], 
                 epack_.evec[i], out[j]);
        }
        t += usecond();
        // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
        LOG(Message) << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;
    };
private:
    MIO::LoadEigenPackPar epPar_;
    unsigned int          evBatchSize_, sourceBatchSize_;
    GridBase              *grid_, *gridIo_;
    int                   traj_;
    Pack                  epack_;
    const GaugeLinkField  *transform_;
};


/******************************************************************************
 *                     TBatchExactDeflation implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TBatchExactDeflation<Pack, GImpl>::TBatchExactDeflation(const std::string name)
: Module<BatchExactDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TBatchExactDeflation<Pack, GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().eigenPack.gaugeXform.empty())
    {
        in.push_back(par().eigenPack.gaugeXform);
    }
    
    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TBatchExactDeflation<Pack, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Pack, typename GImpl>
DependencyMap TBatchExactDeflation<Pack, GImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    if (!par().eigenPack.gaugeXform.empty())
    {
        dep.insert({par().eigenPack.gaugeXform, getName()});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TBatchExactDeflation<Pack, GImpl>::setup(void)
{
    GridBase       *grid, *gridIo = nullptr, *gridRb = nullptr;
    GaugeLinkField *transform = nullptr;

    LOG(Message) << "Setting batch exact deflation guesser with eigenpack "
                 << "located at '" << par().eigenPack.filestem
                 << "' (" << par().eigenPack.size << " modes), EV batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;
    if (!par().eigenPack.gaugeXform.empty())
    {
        transform = &envGet(GaugeLinkField, par().eigenPack.gaugeXform);
    }
    grid   = getGrid<Field>(par().eigenPack.Ls);
    gridRb = getGrid<Field>(par().eigenPack.redBlack, par().eigenPack.Ls);
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = getGrid<FieldIo>(par().eigenPack.redBlack, par().eigenPack.Ls);
    }
    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationGuesser<Pack, GImpl>),
                     getName(), par().eigenPack.Ls, par().eigenPack, par().evBatchSize,
                     par().sourceBatchSize, gridRb, gridIo, vm().getTrajectory(), transform);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TBatchExactDeflation<Pack, GImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflation_hpp_
