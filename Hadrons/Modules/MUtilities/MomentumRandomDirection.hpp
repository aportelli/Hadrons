#ifndef Hadrons_MUtilities_MomentumRandomDirection_hpp_
#define Hadrons_MUtilities_MomentumRandomDirection_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        generate spatial twist with given norm and random direction         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class MomentumRandomDirectionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumRandomDirectionPar,
                                    double, norm,
                                    std::string, output);
};

class TMomentumRandomDirection: public Module<MomentumRandomDirectionPar>
{
public:
    // constructor
    TMomentumRandomDirection(const std::string name);
    // destructor
    virtual ~TMomentumRandomDirection(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(MomentumRandomDirection, TMomentumRandomDirection, MUtilities);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_MomentumRandomDirection_hpp_
