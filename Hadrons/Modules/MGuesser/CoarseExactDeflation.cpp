#include <Hadrons/Modules/MGuesser/CoarseExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<BaseFermionEigenPack<FIMPL>, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<BaseFermionEigenPack<FIMPLF>, HADRONS_DEFAULT_LANCZOS_NBASIS>;