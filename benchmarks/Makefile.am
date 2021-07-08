AM_LDFLAGS  += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)
 
EXTRA_PROGRAMS = \
  Benchmark_meson_field            

CLEANFILES = $(EXTRA_PROGRAMS)

benchmarks-local: $(EXTRA_PROGRAMS)

Benchmark_meson_field_SOURCES=Test_QED_LOCAL.cc
Benchmark_meson_field_LDADD=-lHadrons -lGrid

