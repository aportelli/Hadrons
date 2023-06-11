/*
 * Test_field_io.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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
#include <Hadrons/Environment.hpp>
#include <Hadrons/FieldIo.hpp>
#include <Hadrons/Global.hpp>

using namespace Grid;
using namespace Hadrons;

struct Metadata: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                    std::string, foo,
                                    int,         bar);
};

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    auto              &env   = Environment::getInstance();
    GridBase          *grid  = env.getGrid<vComplex>();
    GridBase          *gridf = env.getGrid<vComplexF>();
    auto              &rng   = *env.get4dRng();
    LatticePropagator prop(grid), test(grid);
    Metadata          md, testMd;

    random(rng, prop);

    FieldWriter<LatticePropagator> writer("test_field_io.bin", grid);

    md.foo = "yay!";
    md.bar = 42;
    writer.writeField(prop, md);
    writer.close();

    FieldReader<LatticePropagator> reader("test_field_io.bin", grid);

    reader.readField(test, testMd);
    reader.close();
    test -= prop;
    LOG(Message) << "Norm 2 diff: " << norm2(test) << std::endl;
    LOG(Message) << "Metadata equal: " << (md == testMd) << std::endl;

    FieldWriter<LatticePropagator, LatticePropagatorF> writer32("test_field_io_32.bin", grid, gridf);

    md.foo = "yay!";
    md.bar = 42;
    writer32.writeField(prop, md);
    writer32.close();

    FieldReader<LatticePropagator, LatticePropagatorF> reader32("test_field_io_32.bin", grid, gridf);

    reader32.readField(test, testMd);
    reader32.close();
    test -= prop;
    LOG(Message) << "Norm 2 diff: " << norm2(test) << std::endl;
    LOG(Message) << "Metadata equal: " << (md == testMd) << std::endl;

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
