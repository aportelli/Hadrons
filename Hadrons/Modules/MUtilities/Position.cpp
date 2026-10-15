#include <Hadrons/Modules/MUtilities/Position.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

/******************************************************************************
*                  TPosition implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TPosition::TPosition(const std::string name)
: Module<PositionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TPosition::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TPosition::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TPosition::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
void TPosition::execute(void)
{
    std::vector<int> pos;

    pos = strToVec<int>(par().pos);
    if (pos.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of components");
    }
    LOG(Message) << "Created position vector " << pos << std::endl;
    saveResult(par().output, "pos", pos);
    auto &out = envGet(HadronsSerializable, getName());
    out = pos;
}
