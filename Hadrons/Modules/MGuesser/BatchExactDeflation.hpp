#ifndef Hadrons_MGuesser_BatchExactDeflationPreload_hpp_
#define Hadrons_MGuesser_BatchExactDeflationPreload_hpp_

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

class BatchExactDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationPar,
                                    std::string, eigenPack,
                                    unsigned int, epSize,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename FImpl, typename EPack>
class TBatchExactDeflation: public Module<BatchExactDeflationPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
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

MODULE_REGISTER_TMP(BatchExactDeflation, ARG(TBatchExactDeflation<FIMPL,BaseFermionEigenPack<FIMPL>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationF, ARG(TBatchExactDeflation<FIMPLF,BaseFermionEigenPack<FIMPLF>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationEPackF, ARG(TBatchExactDeflation<FIMPL,BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl, typename EPack>
class BatchExactDeflationGuesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
private:
    static constexpr bool requireCast_ = !(std::is_same<Field,EPackField>::value);
    static constexpr std::size_t CastBufferSize_( std::size_t n ){ return requireCast_ ? n : 0; };

    template <typename O>
    typename std::enable_if<std::is_same<Field,O>::value, const std::vector<Field> &>::type
    CopyOrOriginal( const std::vector<Field> &Copy, const std::vector<O> &Original ) { return Original; }
    template <typename O>
    typename std::enable_if<!std::is_same<Field,O>::value, const std::vector<Field> &>::type
    CopyOrOriginal( const std::vector<Field> &Copy, const std::vector<O> &Original ) { return Copy; }
public:
    BatchExactDeflationGuesser(const std::vector<EPackField> & evec,const std::vector<RealD> & eval, 
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
        unsigned int sourceSize = out.size();

        std::vector<Field> evecCast( CastBufferSize_(evBatchSize_) , Field(in[0].Grid()) );
        std::vector<RealD> evalCast( CastBufferSize_(evBatchSize_) , 0. );

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        if (requireCast_) {
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

            if (requireCast_) {
                cast_t -= usecond();
                for (unsigned int i = 0; i < evBlockSize; ++i) {
                    precisionChange(evecCast[i],evec_[bv+i]);
                    evalCast[i] = eval_[bv+i];
                }
                cast_t += usecond();
            }

            proj_t -= usecond();
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                BatchDeflationUtils::projAccumulate(in, out, 
                    CopyOrOriginal(evecCast, evec_),
                    requireCast_ ? evalCast : eval_,
                    requireCast_ ? 0 : bv,
                    requireCast_ ? evBlockSize : bv + evBlockSize,
                    bs, bs + sourceBlockSize);
            }
            proj_t += usecond();
        }

        if (requireCast_) {
            LOG(Message) << "Total precision change time " << cast_t/1.e6 << " s" << std::endl;
        }
        LOG(Message) << "Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    const std::vector<EPackField> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int epSize_, evBatchSize_, sourceBatchSize_;
};


/******************************************************************************
 *                     TBatchExactDeflation implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TBatchExactDeflation<FImpl, EPack>::TBatchExactDeflation(const std::string name)
: Module<BatchExactDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflation<FImpl, EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflation<FImpl, EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TBatchExactDeflation<FImpl, EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflation<FImpl, EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with preloaded eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().epSize << " modes) and batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationGuesser<FImpl,EPack>),
                     getName(), env().getObjectLs(par().eigenPack),
                     epack.evec, epack.eval, par().epSize, par().evBatchSize,
                     par().sourceBatchSize);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflation<FImpl, EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflationPreload_hpp_
