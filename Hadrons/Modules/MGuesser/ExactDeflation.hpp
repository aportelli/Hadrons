#ifndef Hadrons_MGuesser_ExactDeflation_hpp_
#define Hadrons_MGuesser_ExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Exact deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class ExactDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDeflationPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename EPack>
class TExactDeflation: public Module<ExactDeflationPar>
{
public:
    typedef typename EPack::Field Field;
public:
    // constructor
    TExactDeflation(const std::string name);
    // destructor
    virtual ~TExactDeflation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDeflation, TExactDeflation<BaseFermionEigenPack<FIMPL>>, MGuesser);
MODULE_REGISTER_TMP(ExactDeflationF, TExactDeflation<BaseFermionEigenPack<FIMPLF>>, MGuesser);


/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template<class Field>
class ExactDeflatedGuesser: public LinearFunction<Field> {
private:
  const std::vector<Field> &evec;
  const std::vector<RealD> &eval;
  const unsigned int       epackSize;

public:
    using LinearFunction<Field>::operator();

    ExactDeflatedGuesser(const std::vector<Field> & _evec,const std::vector<RealD> & _eval)
    : ExactDeflatedGuesser(_evec, _eval, _evec.size())
    {}

    ExactDeflatedGuesser(const std::vector<Field> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
    : evec(_evec), eval(_eval), epackSize(_N)
    {
        assert(evec.size()==eval.size());
        assert(epackSize <= evec.size());
    } 

    virtual void operator()(const Field &src,Field &out) {
        out = Zero();
        double proj_t = 0.;

        for (int i=0;i<epackSize;i++) {
            const Field& tmp = evec[i];
            proj_t -= usecond();
            axpy(out,TensorRemove(innerProduct(tmp,src)) / eval[i],tmp,out);
            proj_t += usecond();
        }
        out.Checkerboard() = src.Checkerboard();

        LOG(Message) << "ExactDeflatedGuesser: Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
    }

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        assert(in.size() == out.size());

        unsigned int sourceSize = out.size();

        unsigned int evBatchSize     = 1;
        unsigned int sourceBatchSize = sourceSize;
        // evBatchSize and sourceBatchSize left in case we want a different value in future

        for (auto &v: out)
            v = Zero();

        double proj_t = 0.;

        for (unsigned int bv = 0; bv < epackSize; bv += evBatchSize) {
        for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize)
        {
            unsigned int evBlockSize = std::min(epackSize - bv, evBatchSize);
            unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize);

            proj_t -= usecond();
            BatchDeflationUtils::projAccumulate(in, out, evec, eval, 
                                                bv, bv + evBlockSize, 
                                                bs, bs + sourceBlockSize);
            proj_t += usecond();
        }}

        LOG(Message) << "ExactDeflatedGuesser: Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
    }
};


/******************************************************************************
 *                 TExactDeflation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TExactDeflation<EPack>::TExactDeflation(const std::string name)
: Module<ExactDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TExactDeflation<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename EPack>
std::vector<std::string> TExactDeflation<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename EPack>
DependencyMap TExactDeflation<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflation<EPack>::setup(void)
{
    LOG(Message) << "Setting exact deflation guesser with eigenpack '"
                 << par().eigenPack << "' (" 
                 << par().size << " modes)" << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ExactDeflatedGuesser<Field>, getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflation<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_ExactDeflation_hpp_
