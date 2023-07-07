#ifndef Hadrons_MUtilities_Momentum_hpp_
#define Hadrons_MUtilities_Momentum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                simple module creating a constant momentum                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class MomentumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumPar,
                                    std::string, mom,
                                    std::string, output);
};

class TMomentum: public Module<MomentumPar>
{
public:
    // constructor
    TMomentum(const std::string name);
    // destructor
    virtual ~TMomentum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(Momentum, TMomentum, MUtilities);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_Momentum_hpp_
