#ifndef Hadrons_MUtilities_Position_hpp_
#define Hadrons_MUtilities_Position_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Position                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PositionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PositionPar,
                                    std::string, pos,
                                    std::string, output);
};

class TPosition: public Module<PositionPar>
{
public:
    // constructor
    TPosition(const std::string name);
    // destructor
    virtual ~TPosition(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(Position, TPosition, MUtilities);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_Position_hpp_
