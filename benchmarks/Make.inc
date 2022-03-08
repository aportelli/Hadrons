benchmarks-local: Benchmark_meson_field
EXTRA_PROGRAMS = Benchmark_meson_field

Benchmark_meson_field_SOURCES=Benchmark_meson_field.cc
Benchmark_meson_field_LDADD=-lHadrons -lGrid

CLEANFILES = Benchmark_meson_field

