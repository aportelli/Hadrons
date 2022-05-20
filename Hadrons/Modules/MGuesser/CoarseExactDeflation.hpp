#ifndef Hadrons_MGuesser_CoarseExactDeflation_hpp_
#define Hadrons_MGuesser_CoarseExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Exact deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class CoarseExactDeflationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CoarseExactDeflationPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename EPack, int nbasis>
class TCoarseExactDeflation: public Module<CoarseExactDeflationPar>
{
public:
    typedef typename EPack::Field Field;
public:
    // constructor
    TCoarseExactDeflation(const std::string name);
    // destructor
    virtual ~TCoarseExactDeflation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CoarseExactDeflation, TCoarseExactDeflation<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>>, MGuesser);
MODULE_REGISTER_TMP(CoarseExactDeflationF, TCoarseExactDeflation<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>>, MGuesser);

/******************************************************************************
 *                 TExactDeflation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack, int nbasis>
TCoarseExactDeflation<EPack, nbasis>::TCoarseExactDeflation(const std::string name)
: Module<CoarseExactDeflationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack, int nbasis>
std::vector<std::string> TCoarseExactDeflation<EPack, nbasis>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename EPack, int nbasis>
std::vector<std::string> TCoarseExactDeflation<EPack, nbasis>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename EPack, int nbasis>
DependencyMap TCoarseExactDeflation<EPack, nbasis>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack, int nbasis>
void TCoarseExactDeflation<EPack, nbasis>::setup(void)
{
    LOG(Message) << "Setting exact deflation guesser with eigenpack '"
                 << par().eigenPack << "' (" 
                 << par().size << " modes)" << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>,  LocalCoherenceDeflatedGuesser<Field,nbasis>, getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack, int nbasis>
void TCoarseExactDeflation<EPack, nbasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_CoarseExactDeflation_hpp_
