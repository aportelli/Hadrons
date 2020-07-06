#include <Hadrons/Modules/MIO/SaveField.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TSaveField<FIMPL::PropagatorField>;
template class Grid::Hadrons::MIO::TSaveField<FIMPL::PropagatorField, FIMPLF::PropagatorField>;
