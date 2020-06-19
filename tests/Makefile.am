AM_LDFLAGS  += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)
 
EXTRA_PROGRAMS = \
  Test_QED               \
  Test_database          \
  Test_diskvector        \
  Test_distil            \
  Test_free_prop         \
  Test_hadrons_meson_3pt \
  Test_hadrons_spectrum  \
  Test_sigma_to_nucleon

CLEANFILES = $(EXTRA_PROGRAMS)

tests-local: $(EXTRA_PROGRAMS)

Test_QED_SOURCES=Test_QED.cpp
Test_QED_LDADD=-lHadrons -lGrid

Test_diskvector_SOURCES=Test_diskvector.cpp
Test_diskvector_LDADD=-lHadrons -lGrid

Test_database_SOURCES=Test_database.cpp
Test_database_LDADD=-lHadrons -lGrid

Test_distil_SOURCES=Test_distil.cpp
Test_distil_LDADD=-lHadrons -lGrid

Test_free_prop_SOURCES=Test_free_prop.cpp
Test_free_prop_LDADD=-lHadrons -lGrid

Test_hadrons_meson_3pt_SOURCES=Test_hadrons_meson_3pt.cpp
Test_hadrons_meson_3pt_LDADD=-lHadrons -lGrid

Test_hadrons_spectrum_SOURCES=Test_hadrons_spectrum.cpp
Test_hadrons_spectrum_LDADD=-lHadrons -lGrid

Test_sigma_to_nucleon_SOURCES=Test_sigma_to_nucleon.cpp
Test_sigma_to_nucleon_LDADD=-lHadrons -lGrid
