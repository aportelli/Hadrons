#include <Hadrons/Modules/MGuesser/BatchExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TBatchExactDeflation<FIMPL,BaseFermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflation<FIMPLF,BaseFermionEigenPack<FIMPLF>>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflation<FIMPL,BaseFermionEigenPack<FIMPLF>>;
