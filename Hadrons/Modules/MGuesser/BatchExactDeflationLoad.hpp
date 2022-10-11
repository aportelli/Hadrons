#ifndef Hadrons_MGuesser_BatchExactDeflation_hpp_
#define Hadrons_MGuesser_BatchExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MIO/LoadEigenPack.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflationLoadPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationLoadPar,
                                    MIO::LoadEigenPackPar, eigenPack,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename Pack, typename GImpl>
class TBatchExactDeflationLoad: public Module<BatchExactDeflationLoadPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    GAUGE_TYPE_ALIASES(GImpl, );
public:
    // constructor
    TBatchExactDeflationLoad(const std::string name);
    // destructor
    virtual ~TBatchExactDeflationLoad(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflationLoad, ARG(TBatchExactDeflationLoad<FermionEigenPack<FIMPL>, GIMPL>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationLoadF, ARG(TBatchExactDeflationLoad<FermionEigenPack<FIMPLF>, GIMPLF>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflationLoadDIoF, ARG(TBatchExactDeflationLoad<FermionEigenPack<FIMPL, FIMPLF>, GIMPL>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename Pack, typename GImpl>
class BatchExactDeflationLoadGuesser: public LinearFunction<typename Pack::Field>
{
public:
    typedef typename Pack::Field Field;
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    BatchExactDeflationLoadGuesser(const MIO::LoadEigenPackPar &epPar, 
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
        unsigned int sourceSize = out.size();

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        LOG(Message) << "Total Sources: " << sourceSize << std::endl;
        LOG(Message) << "Total Eigenvectors: " << epPar_.size << std::endl;
        LOG(Message) << "Source batch size: " << sourceBatchSize_ << std::endl;
        LOG(Message) << "Eigenvector batch size: " << evBatchSize_ << std::endl;

        double IOAccum = 0; 
        double ProjAccum = 0;
        double GFixAccum = 0;
        TimerArray trA;
 
        LOG(Message) << "--- zero guesses" << std::endl;
        trA.startTimer("Zero");
        for (auto &v: out)
        {
            v = Zero();
        }
        trA.stopTimer("Zero");
        LOG(Message) << "Zeroing the 'out' vector took: " << trA.getDTimer("Zero")/1.0e6 << std::endl;

        for (unsigned int bv = 0; bv < epPar_.size; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epPar_.size - bv, evBatchSize_);

            LOG(Message) << "--- Ev batch " << bv/evBatchSize_ << std::endl;
            LOG(Message) << "--- I/O" << std::endl;
            
            trA.startTimer("IO");
            epack_.read(epPar_.filestem, epPar_.multiFile, bv, bv + evBlockSize, traj_);
            trA.stopTimer("IO");
            IOAccum += trA.getDTimer("IO");
            LOG(Message) << "IO of " << evBlockSize << " eigenvectors took: " << trA.getDTimer("IO")/1.0e6 << std::endl;

            trA.startTimer("GF");
            if (transform_ != nullptr)
            {
                epack_.gaugeTransform(*transform_);
            }
            trA.stopTimer("GF");
            GFixAccum += trA.getDTimer("GF");
            LOG(Message) << "Gauge fixing took: " << trA.getDTimer("GF")/1.0e6 << std::endl;
        
            LOG(Message) << "Project" << std::endl;
            trA.startTimer("Proj");
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);
                LOG(Message) << "--- Source batch " << bs/sourceBatchSize_ << std::endl;
                BatchDeflationUtils::projAccumulate(in, out, 
                                                    epack_.evec, epack_.eval, 
                                                    0, evBlockSize, 
                                                    bs, bs + sourceBlockSize);
            }
            trA.stopTimer("Proj");
            ProjAccum += trA.getDTimer("Proj");
            LOG(Message) << "Projection took: " << trA.getDTimer("Proj")/1.0e6 << std::endl;
        }
        
        LOG(Message) << "=== Timings ===" << std::endl;
        LOG(Message) << "Total Zero time" << trA.getDTimer("Zero")/1.0e6 << std::endl;
        LOG(Message) << "Total IO time: " << IOAccum/1.0e6 << std::endl;
        LOG(Message) << "Total Gauge Fix time: " << GFixAccum/1.0e6 << std::endl;
        LOG(Message) << "Total Project time: " << ProjAccum/1.0e6 << std::endl;
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    MIO::LoadEigenPackPar epPar_;
    unsigned int          evBatchSize_, sourceBatchSize_;
    GridBase              *grid_, *gridIo_;
    int                   traj_;
    Pack                  epack_;
    const GaugeLinkField  *transform_;
};


/******************************************************************************
 *                     TBatchExactDeflationLoad implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TBatchExactDeflationLoad<Pack, GImpl>::TBatchExactDeflationLoad(const std::string name)
: Module<BatchExactDeflationLoadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TBatchExactDeflationLoad<Pack, GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().eigenPack.gaugeXform.empty())
    {
        in.push_back(par().eigenPack.gaugeXform);
    }
    
    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TBatchExactDeflationLoad<Pack, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Pack, typename GImpl>
DependencyMap TBatchExactDeflationLoad<Pack, GImpl>::getObjectDependencies(void)
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
void TBatchExactDeflationLoad<Pack, GImpl>::setup(void)
{
    GridBase       *gridIo = nullptr, *gridRb = nullptr;
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
    gridRb = getGrid<Field>(par().eigenPack.redBlack, par().eigenPack.Ls);
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = getGrid<FieldIo>(par().eigenPack.redBlack, par().eigenPack.Ls);
    }
    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationLoadGuesser<Pack, GImpl>),
                     getName(), par().eigenPack.Ls, par().eigenPack, par().evBatchSize,
                     par().sourceBatchSize, gridRb, gridIo, vm().getTrajectory(), transform);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TBatchExactDeflationLoad<Pack, GImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflation_hpp_
