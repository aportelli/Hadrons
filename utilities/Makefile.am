AM_LDFLAGS += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)

bin_PROGRAMS = \
  BatchDeflationBenchmark    \
  HadronsContractor          \
  HadronsContractorBenchmark \
  HadronsXmlRun              \
  HadronsXmlValidate         \
  HadronsFermionEP64To32 

BatchDeflationBenchmark_SOURCES = BatchDeflationBenchmark.cpp
BatchDeflationBenchmark_LDADD = -lHadrons -lGrid
BatchDeflationBenchmark_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsXmlRun_SOURCES = HadronsXmlRun.cpp
HadronsXmlRun_LDADD   = -lHadrons -lGrid
HadronsXmlRun_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsXmlValidate_SOURCES = HadronsXmlValidate.cpp
HadronsXmlValidate_LDADD   = -lHadrons -lGrid
HadronsXmlValidate_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cpp
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = -lHadrons -lGrid
HadronsFermionEP64To32_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsContractor_SOURCES = Contractor.cpp
HadronsContractor_LDADD   = -lHadrons -lGrid
HadronsContractor_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsContractorBenchmark_SOURCES = ContractorBenchmark.cpp
HadronsContractorBenchmark_LDADD   = -lHadrons -lGrid
HadronsContractorBenchmark_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a
