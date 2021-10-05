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
                                    unsigned int, batchSize);
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
                               const unsigned int batchSize,
                               GridBase *grid, GridBase *gridIo,
                               const int traj, 
                               const GaugeLinkField *transform = nullptr)
    : epPar_(epPar)
    , batchSize_(batchSize)
    , grid_(grid)
    , gridIo_(gridIo)
    , traj_(traj)
    , transform_(transform)
    {
        epack_.init(batchSize_, grid_, gridIo_);
    };

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        unsigned int nBatch = epPar_.size/batchSize_ + (((epPar_.size % batchSize_) != 0) ? 1 : 0);

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        LOG(Message) << "--- zero guesses" << std::endl;
        for (auto &v: out)
        {
            v = Zero();
        }
        for (unsigned int b = 0; b < epPar_.size; b += batchSize_)
        {
            unsigned int size = std::min(epPar_.size - b, batchSize_);

            LOG(Message) << "--- batch " << b/batchSize_ << std::endl;
            LOG(Message) << "I/O" << std::endl;
            epack_.read(epPar_.filestem, epPar_.multiFile, b, b + size, traj_);
            if (transform_ != nullptr)
            {
                epack_.gaugeTransform(*transform_);
            }
            LOG(Message) << "project" << std::endl;
            projAccumulate(in, out);
        }
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    void projAccumulate(const std::vector<Field> &in, std::vector<Field> &out)
    {
        for (unsigned int i = 0; i < epack_.evec.size(); ++i)
        for (unsigned int j = 0; j < in.size(); ++j)
        {
            axpy(out[j], 
                 TensorRemove(innerProduct(epack_.evec[i], in[j]))/epack_.eval[i], 
                 epack_.evec[i], out[j]);
        }
    };
private:
    MIO::LoadEigenPackPar epPar_;
    unsigned int          batchSize_;
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
                 << "' (" << par().eigenPack.size << " modes) and batch size " 
                 << par().batchSize << std::endl;
    if (!par().eigenPack.multiFile)
    {
        HADRONS_ERROR(Implementation, "batch deflation not supported yet for single-file eigenpacks");
    }
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
                     getName(), par().eigenPack.Ls, par().eigenPack, par().batchSize,
                     gridRb, gridIo, vm().getTrajectory(), transform);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TBatchExactDeflation<Pack, GImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflation_hpp_
