#include <Hadrons/Modules/MOps/Add.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MOps;

template class HADRONS_NAMESPACE::MOps::TAdd<FIMPL::Field>;
template class HADRONS_NAMESPACE::MOps::TAdd<FIMPL::ComplexField>;
template class HADRONS_NAMESPACE::MOps::TAdd<FIMPL::FermionField>;
template class HADRONS_NAMESPACE::MOps::TAdd<FIMPL::PropagatorField>;
