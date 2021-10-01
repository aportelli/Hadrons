#ifndef Hadrons_MGuesser_ExactDeflation_hpp_
#define Hadrons_MGuesser_ExactDeflation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

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
    envCreateDerived(LinearFunction<Field>, DeflatedGuesser<Field>, getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflation<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_ExactDeflation_hpp_
