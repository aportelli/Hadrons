/*
 * FourFermionFullyConnected.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Julia Kettle <J.R.Kettle-2@sms.ed.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Ryan Abbott <rabbott@mit.edu>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
 * Author: felixerben <46817371+felixerben@users.noreply.github.com>
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

#ifndef Hadrons_MNPR_FourFermionFullyConnected_hpp_
#define Hadrons_MNPR_FourFermionFullyConnected_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MNPR)

class FourFermionFullyConnectedPar : Serializable
{
public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourFermionFullyConnectedPar,
                                        std::string, qIn,
                                        std::string, qOut,
                                        std::string, lIn,
                                        std::string, lOut,
                                        std::string, pIn,
                                        std::string, pOut,
                                        std::string, gamma_basis,
                                        std::string, output);
};


template <typename FImpl>
class TFourFermionFullyConnected : public Module<FourFermionFullyConnectedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaA,
                                        Gamma::Algebra, gammaB,
                                        std::string,  pIn,
                                        std::string,  pOut);
    };
    typedef Correlator<Metadata, SpinColourSpinColourMatrix> Result;

    TFourFermionFullyConnected(const std::string name);
    virtual ~TFourFermionFullyConnected(void) {};

    virtual std::vector<std::string> getInput();
    virtual std::vector<std::string> getOutput();
    virtual std::vector<std::string> getOutputFiles(void);

protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FourFermionFullyConnected, ARG(TFourFermionFullyConnected<FIMPL>), MNPR);

template <typename FImpl>
TFourFermionFullyConnected<FImpl>::TFourFermionFullyConnected(const std::string name)
    : Module<FourFermionFullyConnectedPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TFourFermionFullyConnected<FImpl>::getInput()
{
    std::vector<std::string> in = { par().qIn, par().qOut, par().lIn, par().lOut};

    return in;
}

template <typename FImpl>
std::vector<std::string> TFourFermionFullyConnected<FImpl>::getOutput()
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TFourFermionFullyConnected<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

template <typename FImpl>
void TFourFermionFullyConnected<FImpl>::setup()
{
    LOG(Message) << "Running setup for FourFermionFullyConnected"
        << std::endl;

    envTmpLat(PropagatorField, "quark_bilinear");
    envTmpLat(PropagatorField, "lepton_bilinear");
    envTmpLat(PropagatorField, "tmp_bilinear");
    envTmpLat(ComplexField, "bilinear_phase");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename FImpl>
void TFourFermionFullyConnected<FImpl>::execute()
{
    LOG(Message) << "Computing contractions '" << getName()
        << "' using quark source propagators '" << par().qIn << "' and '" << par().qOut << "'"
        << "' using lepton source propagators '" << par().lIn << "' and '" << par().lOut << "'"
        << std::endl;

    PropagatorField &qIn = envGet(PropagatorField, par().qIn);
    PropagatorField &qOut = envGet(PropagatorField, par().qOut);
    PropagatorField &lIn = envGet(PropagatorField, par().lIn);
    PropagatorField &lOut = envGet(PropagatorField, par().lOut);

    envGetTmp(PropagatorField, quark_bilinear);
    envGetTmp(PropagatorField, lepton_bilinear);
    envGetTmp(PropagatorField, tmp_bilinear);


    std::vector<Result>         result;
    Result                      r;
    Coordinate                  latt_size = GridDefaultLatt();
    std::vector<Real> pIn = strToVec<Real>(par().pIn);
    std::vector<Real> pOut = strToVec<Real>(par().pOut);

    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, bilinear_phase);

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Calculating phases" << std::endl;

    NPRUtils<FImpl>::phase(bilinear_phase,pIn,pOut);

    LOG(Message) << "Computing diagrams" << std::endl;

    r.info.pIn  = par().pIn;
    r.info.pOut = par().pOut;

    auto compute_diagrams = [&](Gamma gamma_A, Gamma gamma_B, bool print = true) {

        r.info.gammaA = gamma_A.g;
        r.info.gammaB = gamma_B.g;

        if (print) {
            LOG(Message) << "Computing diagrams with GammaA = "
                << gamma_A.g << ", " << "GammaB = " << gamma_B.g
                << std::endl;
        }

        // Fully connected diagram
        SpinColourSpinColourMatrix lret;
        quark_bilinear = bilinear_phase * (g5 * adj(qOut) * g5 * gamma_A * qIn);
        lepton_bilinear = bilinear_phase * (g5 * adj(lOut) * g5 * gamma_B * lIn);
        lret = NPRUtils<FImpl>::tensorProdSum(tmp_bilinear, quark_bilinear, lepton_bilinear);

        r.corr.push_back( (1.0 / volume) * lret );
        result.push_back(r);
        r.corr.erase(r.corr.begin());
    };

    std::string gamma_basis = par().gamma_basis;
    if (gamma_basis == "all") {
        for (Gamma gammaA: Gamma::gall) {
            for (Gamma gammaB: Gamma::gall) {
                compute_diagrams(gammaA, gammaB);
            }
        }
    }
    else if (gamma_basis == "diagonal") {
        for (Gamma g: Gamma::gall) {
            compute_diagrams(g, g);
        }
    }
    else if (gamma_basis == "va_av") {
        for (int mu = 0; mu < 4; mu++) {
            Gamma gmu = Gamma::gmu[mu];
            Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
            compute_diagrams(gmu, gmug5);
            compute_diagrams(gmug5, gmu);
        }
    }
    else if (gamma_basis == "diagonal_va" || gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt") {
        for (int mu = 0; mu < 4; mu++) {
            Gamma gmu = Gamma::gmu[mu];
            Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
            compute_diagrams(gmu, gmu);
            compute_diagrams(gmu, gmug5);
            compute_diagrams(gmug5, gmu);
            compute_diagrams(gmug5, gmug5);
        }
        if (gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt") {
            Gamma identity = Gamma(Gamma::Algebra::Identity);

            compute_diagrams(identity, identity);
            compute_diagrams(identity, g5);
            compute_diagrams(g5, identity);
            compute_diagrams(g5, g5);
        }
        if (gamma_basis == "diagonal_va_sp_tt") {
            const std::array<const Gamma, 6> gsigma = {{
                  Gamma(Gamma::Algebra::SigmaXT),
                  Gamma(Gamma::Algebra::SigmaXY),
                  Gamma(Gamma::Algebra::SigmaXZ),
                  Gamma(Gamma::Algebra::SigmaYT),
                  Gamma(Gamma::Algebra::SigmaYZ),
                  Gamma(Gamma::Algebra::SigmaZT)}};

            for (Gamma gammaA: gsigma) {
                    compute_diagrams(gammaA, gammaA);
            }

            compute_diagrams(Gamma(Gamma::Algebra::SigmaXT), Gamma(Gamma::Algebra::SigmaYZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXY), Gamma(Gamma::Algebra::SigmaZT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXZ), Gamma(Gamma::Algebra::SigmaYT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYT), Gamma(Gamma::Algebra::SigmaXZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYZ), Gamma(Gamma::Algebra::SigmaXT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaZT), Gamma(Gamma::Algebra::SigmaXY));
        }
    }
    else {
        LOG(Error) << "Error: unkown gamma_basis: '" << gamma_basis << "'"
            << std::endl;
    }

    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;
    saveResult(par().output, "FourFermionFullyConnected", result);
    auto& out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
