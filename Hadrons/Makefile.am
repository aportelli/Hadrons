AM_CXXFLAGS += -I$(top_srcdir)

lib_LIBRARIES = libHadrons.a

include modules.inc

libHadrons_a_SOURCES = \
	sqlite/sqlite3.c    \
 	Application.cpp     \
	Database.cpp        \
 	Environment.cpp     \
	Exceptions.cpp      \
 	Global.cpp          \
	LeptonImpl.cpp      \
	StatLogger.cpp      \
 	Module.cpp		      \
	TimerArray.cpp      \
	VirtualMachine.cpp  \
	$(modules_cpp)
	
libHadrons_adir = $(includedir)/Hadrons
nobase_libHadrons_a_HEADERS = \
	A2AVectors.hpp            \
	A2AMatrix.hpp             \
	Application.hpp           \
	Database.hpp              \
	DilutedNoise.hpp          \
	DiskVector.hpp            \
	DistilMatrix.hpp          \
	DistillationVectors.hpp   \
	EigenPack.hpp             \
	EmField.hpp               \
	Environment.hpp           \
	Exceptions.hpp            \
	Factory.hpp               \
	FieldIo.hpp               \
	GeneticScheduler.hpp      \
	Global.hpp                \
	Graph.hpp                 \
	LatticeUtilities.hpp      \
	Module.hpp                \
	Modules.hpp               \
	ModuleFactory.hpp         \
	NamedTensor.hpp           \
	Serialization.hpp         \
	Solver.hpp                \
	SqlEntry.hpp              \
	StatLogger.hpp            \
	TimerArray.hpp            \
	VirtualMachine.hpp        \
	sqlite/sqlite3.h          \
	sqlite/sqlite3ext.h       \
	$(modules_hpp)
