#include <Hadrons/Modules/MGuesser/ExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TExactDeflation<BaseFermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MGuesser::TExactDeflation<BaseFermionEigenPack<FIMPLF>>;
