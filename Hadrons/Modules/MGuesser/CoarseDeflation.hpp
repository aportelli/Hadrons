#ifndef Hadrons_MGuesser_CoarseDeflation_hpp_
#define Hadrons_MGuesser_CoarseDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

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
