#include <Hadrons/Modules/MWavelet/Spectrum.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MWavelet;

template class Grid::Hadrons::MWavelet::TSpectrum<FIMPL::FermionField, GIMPL>;
template class Grid::Hadrons::MWavelet::TSpectrum<ColourVectorField<FIMPL>, GIMPL>;