#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // parse command line //////////////////////////////////////////////////////
    std::string parameterFileName;
    
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];
    
    // initialise Grid /////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    
    // initialise application //////////////////////////////////////////////////
    Application            application;
    Application::GlobalPar globalPar;
    
    // reading parameters
    {
        XmlReader reader(parameterFileName);

        read(reader, "global", globalPar);

        // read other application-specific parameters here
    }

    // global initialisation
    application.setPar(globalPar);

    // create modules //////////////////////////////////////////////////////////

    // add modules here with application.createModule<...>(...)
    // the one below is just an example
    application.createModule<MGauge::Unit>("gauge");
    
    // execution ///////////////////////////////////////////////////////////////
    try
    {
        application.run();
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }
    
    // epilogue ////////////////////////////////////////////////////////////////
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
