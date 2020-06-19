/*
 * HadronsXmlValidate.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#include <Hadrons/Application.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // parse command line
    std::string parameterFileName;
    
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];

    try
    {
        Application application(parameterFileName);
        auto        &vm = VirtualMachine::getInstance();
        
        application.parseParameterFile(parameterFileName);
        vm.getModuleGraph();
        LOG(Message) << "Application valid (check XML warnings though)" 
                     << std::endl;
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }
    
    return EXIT_SUCCESS;
}
