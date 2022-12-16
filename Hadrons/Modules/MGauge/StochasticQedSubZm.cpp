#include <Hadrons/Modules/MGauge/StochasticQedSubZm.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

template class Grid::Hadrons::MGauge::TStochasticQedSubZm<vComplex, ZmScheme::qedL>;
template class Grid::Hadrons::MGauge::TStochasticQedSubZm<vComplex, ZmScheme::qedTL>;
