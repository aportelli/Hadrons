#include <Hadrons/Modules/MGuesser/CoarseDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>>;
template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,250>>;
template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPL,400>>;

template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>>;
template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,250>>;
template class Grid::Hadrons::MGuesser::TCoarseDeflation<CoarseFermionEigenPack<FIMPLF,400>>;
