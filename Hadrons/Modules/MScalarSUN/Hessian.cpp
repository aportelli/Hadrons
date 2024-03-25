#include <Hadrons/Modules/MScalarSUN/Hessian.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalarSUN;

template class HADRONS_NAMESPACE::MScalarSUN::THessian<SIMPL>;
template class HADRONS_NAMESPACE::MScalarSUN::THessian<ScalarNxNAdjImplR<2>>;
template class HADRONS_NAMESPACE::MScalarSUN::THessian<ScalarNxNAdjImplR<3>>;
template class HADRONS_NAMESPACE::MScalarSUN::THessian<ScalarNxNAdjImplR<4>>;
template class HADRONS_NAMESPACE::MScalarSUN::THessian<ScalarNxNAdjImplR<5>>;
template class HADRONS_NAMESPACE::MScalarSUN::THessian<ScalarNxNAdjImplR<6>>;
