#include <Hadrons/Modules/MUtilities/MomentumRandomDirection.hpp>
#include <Hadrons/Serialization.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

/******************************************************************************
 *                  TMomentumRandomDirection implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TMomentumRandomDirection::TMomentumRandomDirection(const std::string name)
    : Module<MomentumRandomDirectionPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TMomentumRandomDirection::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

std::vector<std::string> TMomentumRandomDirection::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TMomentumRandomDirection::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
void TMomentumRandomDirection::execute(void)
{
    auto &rng = rngSerial();
    unsigned int nd = env().getNd();
    std::vector<double> r(nd, 0.), mom(nd, 0.);
    double rNorm = 0.;

    for (unsigned int j = 0; j < nd - 1; ++j)
    {
        gaussian(rng, r[j]);
        rNorm += r[j] * r[j];
    }
    rNorm = sqrt(rNorm);
    for (unsigned int j = 0; j < nd; ++j)
    {
        mom[j] = r[j] / rNorm * par().norm;
    }
    LOG(Message) << "Created momentum vector " << mom << " with random direction and norm "
                 << par().norm << std::endl;
    saveResult(par().output, "mom", mom);
    auto &out = envGet(HadronsSerializable, getName());
    out = mom;
}
