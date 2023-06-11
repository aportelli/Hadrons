#include <Hadrons/Global.hpp>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonFermion5DImplementation.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5DImplementation.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5Dcache.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandImplementation.h>

#ifndef AVX512
#ifndef QPX
#ifndef A64FX
#ifndef A64FXFIXEDSIZE
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsAsmImplementation.h>
#endif
#endif
#endif
#endif

NAMESPACE_BEGIN(Grid);

template class CayleyFermion5D<LeptonWilsonImplD>; 
template class CayleyFermion5D<LeptonWilsonImplF>; 
template class WilsonFermion5D<LeptonWilsonImplD>; 
template class WilsonFermion5D<LeptonWilsonImplF>; 
template class WilsonKernels<LeptonWilsonImplD>;
template class WilsonKernels<LeptonWilsonImplF>;

NAMESPACE_END(Grid);
