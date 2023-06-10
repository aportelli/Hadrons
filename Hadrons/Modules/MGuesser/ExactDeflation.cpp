#include <Hadrons/Modules/MGuesser/ExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class HADRONS_NAMESPACE::MGuesser::TExactDeflation<BaseFermionEigenPack<FIMPL>>;
template class HADRONS_NAMESPACE::MGuesser::TExactDeflation<BaseFermionEigenPack<FIMPLF>>;
