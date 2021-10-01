#include <Hadrons/Modules/MGuesser/BatchExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TBatchExactDeflation<FermionEigenPack<FIMPL>, GIMPL>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflation<FermionEigenPack<FIMPLF>, GIMPLF>;
