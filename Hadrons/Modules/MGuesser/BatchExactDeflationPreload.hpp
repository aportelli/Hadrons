#ifndef Hadrons_MGuesser_BatchExactDeflationPreload_hpp_
#define Hadrons_MGuesser_BatchExactDeflationPreload_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflationPreloadPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationPreloadPar,
                                    std::string, eigenPack,
                                    unsigned int, epSize,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename FImpl, typename EPack>
class TBatchExactDeflationPreload: public Module<BatchExactDeflationPreloadPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
public:
    // constructor
    TBatchExactDeflationPreload(const std::string name);
    // destructor
    virtual ~TBatchExactDeflationPreload(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflationPreload, ARG(TBatchExactDeflationPreload<FIMPL,BaseFermionEigenPack<FIMPL>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationPreloadF, ARG(TBatchExactDeflationPreload<FIMPL,BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl, typename EPack>
class BatchExactDeflationPreloadGuesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
public:
    BatchExactDeflationPreloadGuesser(const std::vector<EPackField> & evec,const std::vector<RealD> & eval, 
                                        const unsigned int epSize,
                                        const unsigned int evBatchSize,
                                        const unsigned int sourceBatchSize)
    : evec_(evec)
    , eval_(eval)
    , epSize_(epSize)
    , evBatchSize_(evBatchSize)
    , sourceBatchSize_(sourceBatchSize)
    {};

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        bool requireCast = !(std::is_same<Field,EPackField>::value);

        unsigned int nBatch = epSize_/evBatchSize_ + (((epSize_ % evBatchSize_) != 0) ? 1 : 0);
        unsigned int sourceSize = out.size();

        std::vector<Field> evecCast;
        std::vector<RealD> evalCast;
        if (requireCast) {
            evecCast.resize(evBatchSize_, Field(in[0].Grid()) );
            evalCast.resize(evBatchSize_, 0.);
        }

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        if (requireCast) {
            LOG(Message) << "Eigenpack requires precision change" << std::endl;
        }
        for (auto &v: out)
        {
            v = Zero();
        }

        double cast_t = 0.;
        double proj_t = 0.;

        for (unsigned int bv = 0; bv < epSize_; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epSize_ - bv, evBatchSize_);

            cast_t -= usecond();
            if (requireCast) {
                for (unsigned int i = 0; i < evBlockSize; ++i) {
                    precisionChange(evecCast[i],evec_[bv+i]);
                    evalCast[i] = eval_[bv+i];
                }
            }
            cast_t += usecond();

            proj_t -= usecond();
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                if (requireCast) {
                    projAccumulate(in, out, evecCast, evalCast, 0, evBlockSize, bs, bs + sourceBlockSize);
                } else {
                    projAccumulate(in, out, evec_, eval_, bv, bv + evBlockSize, bs, bs + sourceBlockSize);
                }
            }
            proj_t += usecond();
        }

        if (requireCast) {
            LOG(Message) << "Total precision change time " << cast_t/1.e6 << " s" << std::endl;
        }
        LOG(Message) << "Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    void projAccumulateImpl(const std::vector<Field> &in, std::vector<Field> &out,
                        const std::vector<Field>& evec,
                        const std::vector<RealD>& eval,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf,
                        std::true_type)
    {
        GridBase *g       = in[0].Grid();
        double   lVol     = g->lSites();
        double   siteSize = sizeof(typename Field::scalar_object);
        double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
        double   nIt      = (ef-ei)*(sf - si);
        double   t        = 0.;

        t -= usecond();
        for (unsigned int i = ei; i < ef; ++i)
        for (unsigned int j = si; j < sf; ++j)
        {
            axpy(out[j], 
                 TensorRemove(innerProduct(evec[i], in[j]))/eval[i], 
                 evec[i], out[j]);
        }
        t += usecond();
        // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
        LOG(Message) << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;
    };

    template<typename F1, typename F2>
    void projAccumulateImpl(const std::vector<F1> &in, std::vector<F1> &out,
                        const std::vector<F2>& evec,
                        const std::vector<RealD>& eval,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf,
                        std::false_type)
    {
        assert(0 && "Type mismatch");
    };


    template<typename F1, typename F2>
    void projAccumulate(const std::vector<F1> &in, std::vector<F1> &out,
                        const std::vector<F2>& evec,
                        const std::vector<RealD>& eval,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf)
    {
        projAccumulateImpl(in, out, evec, eval, ei, ef, si, sf,
                        std::integral_constant<bool, std::is_same<F1,F2>::value>{});
    }
private:
    const std::vector<EPackField> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int epSize_, evBatchSize_, sourceBatchSize_;
};


/******************************************************************************
 *                     TBatchExactDeflationPreload implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TBatchExactDeflationPreload<FImpl, EPack>::TBatchExactDeflationPreload(const std::string name)
: Module<BatchExactDeflationPreloadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationPreload<FImpl, EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationPreload<FImpl, EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TBatchExactDeflationPreload<FImpl, EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationPreload<FImpl, EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with preloaded eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().epSize << " modes) and batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationPreloadGuesser<FImpl,EPack>),
                     getName(), env().getObjectLs(par().eigenPack),
                     epack.evec, epack.eval, par().epSize, par().evBatchSize,
                     par().sourceBatchSize);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationPreload<FImpl, EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflationPreload_hpp_
