#ifndef Hadrons_MGuesser_BatchExactDeflationPreloadSrcCast_hpp_
#define Hadrons_MGuesser_BatchExactDeflationPreloadSrcCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflationSrcCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationSrcCastPar,
                                    std::string, eigenPack,
                                    unsigned int, epSize,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename FImpl, typename EPack>
class TBatchExactDeflationSrcCast: public Module<BatchExactDeflationSrcCastPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
public:
    // constructor
    TBatchExactDeflationSrcCast(const std::string name);
    // destructor
    virtual ~TBatchExactDeflationSrcCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflationSrcCastEPackF, ARG(TBatchExactDeflationSrcCast<FIMPL, BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl, typename EPack>
class BatchExactDeflationSrcCastGuesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;

    BatchExactDeflationSrcCastGuesser(const std::vector<EPackField> & evec,const std::vector<RealD> & eval, 
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
        assert(in.size() == out.size());

        unsigned int sourceSize = out.size();

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;

        for (auto &v: out)
        {
            v = Zero();
        }

        double cast_t = 0.;
        double proj_t = 0.;

        std::vector<EPackField> inCast( sourceSize , EPackField(evec_[0].Grid()) );
        std::vector<EPackField> outCast( sourceSize , EPackField(evec_[0].Grid()) );

        cast_t -= usecond();
        for (unsigned int i = 0; i < sourceSize; ++i) {
            precisionChange(inCast[i],in[i]);
            outCast[i] = Zero();
        }
        cast_t += usecond();

        for (unsigned int bv = 0; bv < epSize_; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epSize_ - bv, evBatchSize_);

            proj_t -= usecond();
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                BatchDeflationUtils::projAccumulate(inCast, outCast, 
                    evec_, eval_,
                    bv, bv + evBlockSize,
                    bs, bs + sourceBlockSize);
            }
            proj_t += usecond();
        }

        cast_t -= usecond();
        for (unsigned int i = 0; i < sourceSize; ++i) {
            precisionChange(out[i],outCast[i]);
        }
        cast_t += usecond();

        LOG(Message) << "Total precision change time " << cast_t/1.e6 << " s" << std::endl;
        LOG(Message) << "Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }

private:
    const std::vector<EPackField> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int epSize_, evBatchSize_, sourceBatchSize_;
};


/******************************************************************************
 *                     TBatchExactDeflationSrcCast implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TBatchExactDeflationSrcCast<FImpl, EPack>::TBatchExactDeflationSrcCast(const std::string name)
: Module<BatchExactDeflationSrcCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationSrcCast<FImpl, EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationSrcCast<FImpl, EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TBatchExactDeflationSrcCast<FImpl, EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationSrcCast<FImpl, EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with preloaded eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().epSize << " modes) and batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationSrcCastGuesser<FImpl,EPack>),
                     getName(), env().getObjectLs(par().eigenPack),
                     epack.evec, epack.eval, par().epSize, par().evBatchSize,
                     par().sourceBatchSize);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationSrcCast<FImpl, EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflationPreloadSrcCast_hpp_
