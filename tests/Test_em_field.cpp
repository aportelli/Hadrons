/*
 * Test_em_field.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
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
#include <Hadrons/EmField.hpp>
#include <Hadrons/Global.hpp>

using namespace Grid;
using namespace Hadrons;

constexpr double tolerance = 1e-13;

template <typename T>
void check(const T &a, const T &b, const double result = 0.)
{
    T      diff(a.Grid());
    double n2diff;

    diff   = a - b;
    n2diff = norm2(diff);
    LOG(Message) << "|a|^2 = " << norm2(a) << " / |b|^2 = " << norm2(b) << std::endl;
    LOG(Message) << "|a-b|^2 = " << n2diff << std::endl;
    LOG(Message) << "Expected: " << result << std::endl;
    if (fabs(n2diff - result) < tolerance)
    {
        LOG(Message) << "PASSED" << std::endl;
    }
    else
    {
        LOG(Message) << "FAILED" << std::endl;
        LOG(Error) << "Regression failed" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    auto                          &env   = Environment::getInstance();
    GridBase                      *grid  = env.getGrid<vComplex>();
    auto                          &rng   = *env.get4dRng();
    EmFieldGenerator              emGen(grid);
    PhotonR                       photonQedl(grid, PhotonR::Gauge::feynman, PhotonR::ZmScheme::qedL);
    PhotonR                       photonQedtl(grid, PhotonR::Gauge::feynman, PhotonR::ZmScheme::qedTL);
    PhotonR                       photonCQedl(grid, PhotonR::Gauge::coulomb, PhotonR::ZmScheme::qedL);
    PhotonR                       photonCQedtl(grid, PhotonR::Gauge::coulomb, PhotonR::ZmScheme::qedTL);
    EmFieldGenerator::GaugeField  a(grid);
    EmFieldGenerator::ScalarField w(grid), v(grid);
    PhotonR::GaugeLinkField       photonW(grid);
    PhotonR::GaugeField           photonA(grid);
    auto                          dim = grid->FullDimensions();
    double                        nt = dim[dim.size() - 1], nl = dim[0];

    LOG(Message) << "============ Regressing Feynman QEDL against Grid" << std::endl;
    rng.SeedUniqueString("bla");
    photonQedl.StochasticWeight(photonW);
    photonQedl.StochasticField(photonA, rng, photonW);
    rng.SeedUniqueString("bla");
    emGen.makeWeightsQedL(w);
    emGen(a, rng, w);
    check(photonA, a);

    LOG(Message) << "============ Regressing Coulomb QEDL against Grid" << std::endl;
    rng.SeedUniqueString("bla");
    photonCQedl.StochasticWeight(photonW);
    photonCQedl.StochasticField(photonA, rng, photonW);
    rng.SeedUniqueString("bla");
    emGen.makeWeightsQedL(w);
    emGen(a, rng, w, &EmFieldGenerator::transverseProjectSpatial);
    check(photonA, a);

    LOG(Message) << "============ Regressing Feynman QEDTL against Grid" << std::endl;
    rng.SeedUniqueString("bli");
    photonQedtl.StochasticWeight(photonW);
    photonQedtl.StochasticField(photonA, rng, photonW);
    rng.SeedUniqueString("bli");
    emGen.makeWeightsQedTL(w);
    emGen(a, rng, w);
    check(photonA, a);

    LOG(Message) << "============ Regressing Coulomb QEDTL against Grid" << std::endl;
    rng.SeedUniqueString("bli");
    photonCQedtl.StochasticWeight(photonW);
    photonCQedtl.StochasticField(photonA, rng, photonW);
    rng.SeedUniqueString("bli");
    emGen.makeWeightsQedTL(w);
    emGen(a, rng, w, &EmFieldGenerator::transverseProjectSpatial);
    check(photonA, a);

    LOG(Message) << "============ Regressing Feynman weights QEDZeta against QEDL" << std::endl;
    double zeta = 0.42, zmn2 = 0.;
    emGen.makeWeightsQedL(w);
    emGen.makeWeightsQedZeta(v, zeta);
    for (unsigned int n0 = 0; n0 < nt; ++n0)
    {
        zmn2 += pow(1/(pow(2.*sin(M_PI/nt*n0), 2.) + pow(1/(zeta*nl), 2.)), 2.);
    }
    check(w, v, zmn2);
}
