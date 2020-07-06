#include <Hadrons/Modules/MUtilities/RandomField.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class Grid::Hadrons::MUtilities::TRandomField<FIMPL::PropagatorField>;
template class Grid::Hadrons::MUtilities::TRandomField<FIMPL::FermionField>;
template class Grid::Hadrons::MUtilities::TRandomField<FIMPL::ComplexField>;
