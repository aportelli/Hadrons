AM_CXXFLAGS += -I$(top_srcdir)

noinst_LIBRARIES = libsqlite3.a

libsqlite3_a_SOURCES = sqlite/sqlite3.c sqlite/sqlite3.h sqlite/sqlite3ext.h

lib_LIBRARIES = libHadrons.a

include modules.inc

libHadrons_a_SOURCES = \
  Application.cpp     \
  Environment.cpp     \
	Exceptions.cpp      \
  Global.cpp          \
  Module.cpp		      \
	TimerArray.cpp      \
	VirtualMachine.cpp  \
	$(modules_cpp)
	
libHadrons_adir = $(includedir)/Hadrons
nobase_libHadrons_a_HEADERS = \
	A2AVectors.hpp            \
	A2AMatrix.hpp             \
	Application.hpp           \
	DilutedNoise.hpp          \
	DiskVector.hpp            \
	EigenPack.hpp             \
	Environment.hpp           \
	Exceptions.hpp            \
	Factory.hpp               \
	GeneticScheduler.hpp      \
	Global.hpp                \
	Graph.hpp                 \
	Module.hpp                \
	Modules.hpp               \
	ModuleFactory.hpp         \
  NamedTensor.hpp           \
	Solver.hpp                \
	TimerArray.hpp            \
	VirtualMachine.hpp        \
	sqlite/sqlite3.h          \
	sqlite/sqlite3ext.h       \
	$(modules_hpp)

libHadrons_a_LIBADD = libsqlite3.a
