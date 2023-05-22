#include <Hadrons/Modules/MGuesser/CoarseDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>>;
template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,250>>;
template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,400>>;

template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>>;
template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,250>>;
template class HADRONS_NAMESPACE::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,400>>;
