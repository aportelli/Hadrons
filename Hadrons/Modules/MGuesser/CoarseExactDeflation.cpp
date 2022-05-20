#include <Hadrons/Modules/MGuesser/CoarseExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<BaseFermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<BaseFermionEigenPack<FIMPLF>>;