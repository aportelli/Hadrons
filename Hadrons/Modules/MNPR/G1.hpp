/*
 * G1.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Abbott <rabbott@mit.edu>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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
#ifndef Hadrons_MNPR_G1_hpp_
#define Hadrons_MNPR_G1_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         G1                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class G1Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(G1Par,
                                    std::string, qIn,
                                    std::string, qOut,
                                    std::string, pIn,
                                    std::string, pOut,
                                    std::string, gauge,
                                    std::string, output);
};

template <typename FImpl>
class TG1: public Module<G1Par>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                    SpinColourMatrix, fourq_scalar,
                    SpinColourMatrix, fourq_gamma5,
                    SpinColourMatrix, twoq_scalar,
                    SpinColourMatrix, twoq_gamma5);
    };
    // constructor
    TG1(const std::string name);
    // destructor
    virtual ~TG1(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(G1, TG1<FIMPL>, MNPR);

/******************************************************************************
 *                 TG1 implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TG1<FImpl>::TG1(const std::string name)
: Module<G1Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TG1<FImpl>::getInput(void)
{
    std::vector<std::string> in = { par().qOut, par().qIn, par().gauge };

    return in;
}

template <typename FImpl>
std::vector<std::string> TG1<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TG1<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for G1" << std::endl;

    envTmpLat(PropagatorField, "bilinear");

    envTmpLat(ComplexField, "bilinear_phase");

    envTmpLat(GaugeField, "dSdU");
    envTmpLat(ColourMatrixField, "div_field_strength");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TG1<FImpl>::execute(void)
{
    LOG(Message) << "Computing G1 operators '"
        << "' using source propagators '" << par().qIn << "' and '" << par().qOut << "'"
        << std::endl;

    auto &Umu = envGet(GaugeField, par().gauge);

    PropagatorField qIn = envGet(PropagatorField, par().qIn);
    PropagatorField qOut = envGet(PropagatorField, par().qOut);

    envGetTmp(PropagatorField, bilinear);

    Coordinate                  latt_size = GridDefaultLatt();
    std::vector<Real> pIn = strToVec<Real>(par().pIn);
    std::vector<Real> pOut = strToVec<Real>(par().pOut);

    envGetTmp(ComplexField, bilinear_phase);

    envGetTmp(GaugeField, dSdU);
    envGetTmp(ColourMatrixField, div_field_strength);

    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);

    //// Compute volume
    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++)
    {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Calculating phases" << std::endl;

    NPRUtils<FImpl>::phase(bilinear_phase,pIn,pOut);

    LOG(Warning) << "Iwasaki gauge action enforced by this module. Needs re-writing to allow for other actions!" << std::endl;

    IwasakiGaugeAction<FImpl> action(1.0); // Include freedom to choose the gauge action?
    action.deriv(Umu, dSdU);

    for (int mu = 0; mu < Nd; mu++)
    {
        Gamma gmu = Gamma::gmu[mu];

        // The computation of div_field_strength here is taken from Greg's
        // code; the exact discretization is motivated by the Lattice equations
        // of motion (see Greg's thesis for more detail)
	// Greg McGlynn's thesis: https://doi.org/10.7916/D8T72HD7
        ColourMatrixField U = peekLorentz(Umu, mu);
        // From the implementation of IwasakiGaugeAction::deriv(Umu, dSdU) we
        // have div_field_strength(x) = [U_\mu(x) L_\mu(x)]_{Traceless Antihermetian}
        div_field_strength = peekLorentz(dSdU, mu);
        // div_field_strength(x) = 1/2 (U_\mu(x) L_\mu(x) + L_\mu(x - \hat{\mu}) U_\mu(x - \hat{\mu}))_TA
        U = adj(U) * div_field_strength * U;
        div_field_strength = 0.5 * (div_field_strength + Cshift(U, mu, -1));

        // This next line is not necessary since div_field_strength should
        // already be traceless and anti-hermetian
	// TODO: Delete this or assert / enforce it?
        //div_field_strength = Ta(div_field_strength);

        // Scalar G1
        bilinear = g5 * adj(qOut) * g5 * div_field_strength * gmu * qIn;
        bilinear = bilinear_phase * bilinear;
        result.twoq_scalar += (1.0 / volume) * sum(bilinear);

        bilinear = bilinear_phase * bilinear;
        result.fourq_scalar += (1.0 / volume) * sum(bilinear);

        // Pseudoscalar G1
        bilinear = g5 * adj(qOut) * g5 * div_field_strength * gmu * g5 * qIn;
        bilinear = bilinear_phase * bilinear;
        result.twoq_gamma5 += (1.0 / volume) * sum(bilinear);

        bilinear = bilinear_phase * bilinear;
        result.fourq_gamma5 += (1.0 / volume) * sum(bilinear);
    }
    saveResult(par().output, "G1", result);
    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MG1_G1_hpp_
