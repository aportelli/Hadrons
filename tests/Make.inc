tests-local: Test_QED_Local Test_QED Test_diskvector Test_exact_distil Test_free_prop Test_hadrons_meson_3pt Test_hadrons_spectrum Test_sigma_to_nucleon Test_stoch_distil Test_xi_to_sigma
EXTRA_PROGRAMS = Test_QED_Local Test_QED Test_diskvector Test_exact_distil Test_free_prop Test_hadrons_meson_3pt Test_hadrons_spectrum Test_sigma_to_nucleon Test_stoch_distil Test_xi_to_sigma

Test_QED_Local_SOURCES=Test_QED.cc
Test_QED_LDADD=-lHadrons -lGrid

Test_QED_SOURCES=Test_QED.cc
Test_QED_LDADD=-lHadrons -lGrid

Test_diskvector_SOURCES=Test_diskvector.cc
Test_diskvector_LDADD=-lHadrons -lGrid

Test_exact_distil_SOURCES=Test_exact_distil.cc
Test_exact_distil_LDADD=-lHadrons -lGrid

Test_free_prop_SOURCES=Test_free_prop.cc
Test_free_prop_LDADD=-lHadrons -lGrid

Test_hadrons_meson_3pt_SOURCES=Test_hadrons_meson_3pt.cc
Test_hadrons_meson_3pt_LDADD=-lHadrons -lGrid

Test_hadrons_spectrum_SOURCES=Test_hadrons_spectrum.cc
Test_hadrons_spectrum_LDADD=-lHadrons -lGrid

Test_sigma_to_nucleon_SOURCES=Test_sigma_to_nucleon.cc
Test_sigma_to_nucleon_LDADD=-lHadrons -lGrid

Test_stoch_distil_SOURCES=Test_stoch_distil.cc
Test_stoch_distil_LDADD=-lHadrons -lGrid

Test_xi_to_sigma_SOURCES=Test_xi_to_sigma.cc
Test_xi_to_sigma_LDADD=-lHadrons -lGrid

CLEANFILES = Test_QED_Local Test_QED Test_diskvector Test_exact_distil Test_free_prop Test_hadrons_meson_3pt Test_hadrons_spectrum Test_sigma_to_nucleon Test_stoch_distil Test_xi_to_sigma

