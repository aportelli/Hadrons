#include <Hadrons/Modules/MUtilities/VectorUnpack.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class Grid::Hadrons::MUtilities::TVectorUnpack<FIMPL::ComplexField>;
template class Grid::Hadrons::MUtilities::TVectorUnpack<FIMPL::FermionField>;
template class Grid::Hadrons::MUtilities::TVectorUnpack<FIMPL::PropagatorField>;
