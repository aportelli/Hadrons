#include <Hadrons/Modules/MSolver/ImplicitlyRestartedLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class HADRONS_NAMESPACE::MSolver::TImplicitlyRestartedLanczos<FIMPL::FermionField>;
template class HADRONS_NAMESPACE::MSolver::TImplicitlyRestartedLanczos<FIMPL::FermionField, FIMPLF::FermionField>;
