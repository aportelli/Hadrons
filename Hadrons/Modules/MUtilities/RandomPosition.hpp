#ifndef Hadrons_MUtilities_RandomPosition_hpp_
#define Hadrons_MUtilities_RandomPosition_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         RandomPosition                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class RandomPositionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomPositionPar,
                                    std::string, output);
};

class TRandomPosition: public Module<RandomPositionPar>
{
public:
    // constructor
    TRandomPosition(const std::string name);
    // destructor
    virtual ~TRandomPosition(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(RandomPosition, TRandomPosition, MUtilities);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_RandomPosition_hpp_
