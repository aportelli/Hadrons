#include <Hadrons/Modules/MSolver/ImplicitlyRestartedLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TImplicitlyRestartedLanczos<FIMPL::FermionField>;
template class Grid::Hadrons::MSolver::TImplicitlyRestartedLanczos<FIMPL::FermionField, FIMPLF::FermionField>;
