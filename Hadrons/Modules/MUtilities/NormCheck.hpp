/*
 * NormCheck.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MUtilities_NormCheck_hpp_
#define Hadrons_MUtilities_NormCheck_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MUtilities)

class NormCheckPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NormCheckPar,
                                    std::string, field,
                                    std::string, reference,
                                    std::string, output,
                                    double,      absoluteTolerance,
                                    double,      relativeTolerance,
                                    bool,        failIfExceedTolerance);
};

template <typename Field>
class TNormCheck: public Module<NormCheckPar>
{
public:
    // innerProduct is <field, reference>; swapping the inputs conjugates it.
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,  norm2,
                                        double,  referenceNorm2,
                                        Complex, innerProduct,
                                        double,  differenceNorm2,
                                        double,  relativeDifference,
                                        double,  tolerance,
                                        bool,    isClose);
    };
public:
    TNormCheck(const std::string name);
    virtual ~TNormCheck(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(NormCheckFermion, TNormCheck<FIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckPropagator, TNormCheck<FIMPL::PropagatorField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckComplex, TNormCheck<FIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckColourMatrix, TNormCheck<GIMPL::GaugeLinkField>, MUtilities);

template <typename Field>
TNormCheck<Field>::TNormCheck(const std::string name)
: Module<NormCheckPar>(name)
{}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};

    if (!par().reference.empty())
    {
        in.push_back(par().reference);
    }
    return in;
}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getOutput(void)
{
    return {getName()};
}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getOutputFiles(void)
{
    if (par().output.empty())
    {
        return {};
    }
    return {resultFilename(par().output)};
}

template <typename Field>
void TNormCheck<Field>::setup(void)
{
    if (par().field.empty())
    {
        HADRONS_ERROR(Argument, "field must be specified");
    }
    if (par().absoluteTolerance < 0.0 || par().relativeTolerance < 0.0)
    {
        HADRONS_ERROR(Argument, "tolerances must be non-negative");
    }
    if (par().failIfExceedTolerance && par().reference.empty())
    {
        HADRONS_ERROR(Argument, "reference must be specified when failIfExceedTolerance is set");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename Field>
void TNormCheck<Field>::execute(void)
{
    const auto &field = envGet(Field, par().field);
    Result result{};

    result.norm2 = norm2(field);
    if (!par().reference.empty())
    {
        const auto &reference = envGet(Field, par().reference);
        result.referenceNorm2 = norm2(reference);
        result.innerProduct = innerProduct(field, reference);
        result.differenceNorm2 = norm2(field-reference);

        // Symmetric scale: both fields zero implies zero difference, so 0 is correct.
        const double diff  = std::sqrt(std::max(0.0, result.differenceNorm2));
        const double scale = std::max(std::sqrt(result.norm2), std::sqrt(result.referenceNorm2));
        result.relativeDifference = (scale == 0.0) ? 0.0 : diff/scale;

        // symmetric closeness as in Python's math.isclose:
        // ||field - reference|| <= max(relativeTolerance*max(||field||, ||reference||), absoluteTolerance)
        result.tolerance = std::max(par().relativeTolerance*scale, par().absoluteTolerance);
        // Written as <= so that a NaN difference is not close.
        result.isClose = (diff <= result.tolerance);
    }

    LOG(Message) << "NormCheck '" << par().field << "': norm2 = " << result.norm2 << std::endl;
    if (!par().reference.empty())
    {
        LOG(Message) << "  reference '" << par().reference
                     << "': norm2 = " << result.referenceNorm2
                     << ", difference norm2 = " << result.differenceNorm2
                     << ", relative difference = " << result.relativeDifference
                     << ", tolerance = " << result.tolerance
                     << ", close = " << (result.isClose ? "true" : "false")
                     << std::endl;
    }

    if (par().failIfExceedTolerance && !result.isClose)
    {
        std::ostringstream msg;
        msg << std::scientific << std::setprecision(3)
            << "NormCheck failed: difference norm = "
            << std::sqrt(std::max(0.0, result.differenceNorm2))
            << " (relative " << result.relativeDifference << ")"
            << " exceeds tolerance = " << result.tolerance
            << " = max(relativeTolerance = " << par().relativeTolerance
            << " * max norm, absoluteTolerance = " << par().absoluteTolerance << ")";
        HADRONS_ERROR(Runtime, msg.str());
    }

    saveResult(par().output, "norm_check", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_NormCheck_hpp_
