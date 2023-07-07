#include <Hadrons/Modules/MUtilities/Momentum.hpp>
#include <Hadrons/Serialization.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

/******************************************************************************
*                           TMomentum implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TMomentum::TMomentum(const std::string name)
: Module<MomentumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TMomentum::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TMomentum::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TMomentum::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
void TMomentum::execute(void)
{
    std::vector<double> mom;

    mom = strToVec<double>(par().mom);
    if (mom.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of components");
    }
    LOG(Message) << "Created momentum vector " << mom << std::endl;
    saveResult(par().output, "mom", mom);
    auto &out = envGet(HadronsSerializable, getName());
    out = mom;
}
