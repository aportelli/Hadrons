#include <Hadrons/Modules/MGuesser/BatchExactDeflationLoad.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TBatchExactDeflationLoad<FermionEigenPack<FIMPL>, GIMPL>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflationLoad<FermionEigenPack<FIMPLF>, GIMPLF>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflationLoad<FermionEigenPack<FIMPL, FIMPLF>, GIMPL>;
