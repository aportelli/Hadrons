#ifndef Hadrons_MContraction_DistilMesonFields_hpp_
#define Hadrons_MContraction_DistilMesonFields_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilMesonFields                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DistilMesonFieldsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldsPar,
                                    unsigned int, i);
};

template <typename FImpl>
class TDistilMesonFields: public Module<DistilMesonFieldsPar>
{
public:
    // constructor
    TDistilMesonFields(const std::string name);
    // destructor
    virtual ~TDistilMesonFields(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DistilMesonFields, TDistilMesonFields<FIMPL>, MContraction);

/******************************************************************************
 *                 TDistilMesonFields implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMesonFields<FImpl>::TDistilMesonFields(const std::string name)
: Module<DistilMesonFieldsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonFields<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDistilMesonFields<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonFields<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonFields<FImpl>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DistilMesonFields_hpp_
