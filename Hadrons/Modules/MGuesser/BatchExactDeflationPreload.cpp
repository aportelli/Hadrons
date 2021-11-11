#include <Hadrons/Modules/MGuesser/BatchExactDeflationPreload.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TBatchExactDeflationPreload<FIMPL,BaseFermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflationPreload<FIMPL,BaseFermionEigenPack<FIMPLF>>;
