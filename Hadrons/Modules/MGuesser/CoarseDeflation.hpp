#ifndef Hadrons_MGuesser_CoarseDeflation_hpp_
#define Hadrons_MGuesser_CoarseDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Coarse deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class CoarseDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CoarseDeflationPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename EPack>
class TCoarseDeflation: public Module<CoarseDeflationPar>
{
public:
    typedef typename EPack::Field Field;
    typedef typename EPack::CoarseField CoarseField;
public:
    // constructor
    TCoarseDeflation(const std::string name);
    // destructor
    virtual ~TCoarseDeflation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CoarseDeflation    , ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPL ,HADRONS_DEFAULT_LANCZOS_NBASIS>>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflation250 , ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPL ,250>>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflation400 , ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPL ,400>>), MGuesser);


MODULE_REGISTER_TMP(CoarseDeflationF   , ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflation250F, ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,250>>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflation400F, ARG(TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,400>>), MGuesser);


/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template<class FineField, class CoarseField>
class LocalCoherenceDeflatedGuesser: public LinearFunction<FineField> {
private:
    const std::vector<FineField>   &subspace;
    const std::vector<CoarseField> &evec_coarse;
    const std::vector<RealD>       &eval;
    const unsigned int             epackSize;

public:
    using LinearFunction<FineField>::operator();

    LocalCoherenceDeflatedGuesser(const std::vector<FineField> &subspace_, const std::vector<CoarseField> &evec_coarse_, const std::vector<RealD> &eval_)
    : LocalCoherenceDeflatedGuesser(subspace_,evec_coarse_,eval_,evec_coarse_.size())
    {}

    LocalCoherenceDeflatedGuesser(const std::vector<FineField> &subspace_, const std::vector<CoarseField> &evec_coarse_, const std::vector<RealD> &eval_, unsigned int epackSize_)
    : subspace(subspace_), evec_coarse(evec_coarse_), eval(eval_), epackSize(epackSize_)
    {
        assert(evec_coarse.size()==eval.size());
        assert(epackSize <= evec_coarse.size());
    } 

    virtual void operator()(const FineField &src,FineField &guess) {
        std::vector<FineField> srcVec   = {src};
        std::vector<FineField> guessVec = {guess};

        (*this)(srcVec,guessVec);

        guess = guessVec[0];
    }

    virtual void operator() (const std::vector<FineField> &src, std::vector<FineField> &guess)
    {
        assert(src.size() == guess.size());

        unsigned int sourceSize = src.size();

        unsigned int evBatchSize     = epackSize;
        unsigned int sourceBatchSize = sourceSize;
        // These choices of evBatchSize and sourceBatchSize make the loops
        // below trivial. Left machinery in place in case we want to change
        // this in the future.

        double time_axpy = 0.;
        double time_project = 0.;
        double time_promote = 0.;

        GridBase* coarseGrid = evec_coarse[0].Grid();

        std::vector<CoarseField> src_coarse;   src_coarse.reserve(sourceBatchSize);
        std::vector<CoarseField> guess_coarse; guess_coarse.reserve(sourceBatchSize);
        for (int k=0; k<sourceBatchSize; k++) {
            src_coarse.emplace_back(coarseGrid);
            guess_coarse.emplace_back(coarseGrid);
            guess_coarse[k] = Zero();
        }

        time_project -= usecond();
        batchBlockProject(src_coarse,src,subspace);
        time_project += usecond();

        for (unsigned int bv = 0; bv < epackSize;  bv += evBatchSize) {
        for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize)
        {
            unsigned int evBlockSize     = std::min(epackSize - bv, evBatchSize);
            unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize);

            time_axpy -= usecond();
            BatchDeflationUtils::projAccumulate(src_coarse, guess_coarse, evec_coarse, eval, 
                                                bv, bv + evBlockSize, 
                                                bs, bs + sourceBlockSize);
            time_axpy += usecond();
        }}


        time_promote -= usecond();
        batchBlockPromote(guess_coarse,guess,subspace);
        time_promote += usecond();

        for (int k=0; k<sourceSize; k++)
            guess[k].Checkerboard() = src[k].Checkerboard();

        LOG(Message) << "LocalCoherenceDeflatedGuesser: Total projection time " << time_project/1.e6 << " s" <<  std::endl;
        LOG(Message) << "LocalCoherenceDeflatedGuesser: Total axpy       time " << time_axpy/1.e6 << " s" <<  std::endl;
        LOG(Message) << "LocalCoherenceDeflatedGuesser: Total promotion  time " << time_promote/1.e6 << " s" <<  std::endl;
    }
};

/******************************************************************************
 *                 TExactDeflation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TCoarseDeflation<EPack>::TCoarseDeflation(const std::string name)
: Module<CoarseDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TCoarseDeflation<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename EPack>
std::vector<std::string> TCoarseDeflation<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename EPack>
DependencyMap TCoarseDeflation<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TCoarseDeflation<EPack>::setup(void)
{
    LOG(Message) << "Setting exact deflation guesser with eigenpack '"
                 << par().eigenPack << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ARG(LocalCoherenceDeflatedGuesser<Field, CoarseField>), getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.evecCoarse, epack.evalCoarse);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TCoarseDeflation<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_CoarseDeflation_hpp_
