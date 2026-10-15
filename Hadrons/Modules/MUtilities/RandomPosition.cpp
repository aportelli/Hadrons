#include <Hadrons/Modules/MUtilities/RandomPosition.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

/******************************************************************************
*                  TRandomPosition implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TRandomPosition::TRandomPosition(const std::string name)
: Module<RandomPositionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TRandomPosition::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TRandomPosition::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TRandomPosition::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
void TRandomPosition::execute(void)
{
    auto &rng = rngSerialHadrons();
    
    std::vector<int> dims = env().getDim();
    std::vector<int> pos(dims.size());
    uint32_t tmp;
    for (int i=0; i < dims.size(); ++i)
    {
        uid(rng, tmp);
        pos[i] = static_cast<int>(tmp % static_cast<uint32_t>(dims[i]));
    }

    LOG(Message) << "Created random position vector " << pos << std::endl;
    saveResult(par().output, "pos", pos);
    auto &out = envGet(HadronsSerializable, getName());
    out = pos;
}
