#include <Hadrons/Modules/MIO/LoadField.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TLoadField<FIMPL::PropagatorField>;
template class Grid::Hadrons::MIO::TLoadField<FIMPL::PropagatorField, FIMPLF::PropagatorField>;
