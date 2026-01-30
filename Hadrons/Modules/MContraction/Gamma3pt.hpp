/*
 * Gamma3pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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

#ifndef Hadrons_MContraction_Gamma3pt_hpp_
#define Hadrons_MContraction_Gamma3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>


// TODO: When traceColour and traceSpin are implemented in Grid, remove these lines
NAMESPACE_BEGIN(Grid);

#define GRID_UNOP(name)   name
#define GRID_DEF_UNOP(op, name)						\
  template <typename T1, typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value,T1>::type * = nullptr> \
  inline auto op(const T1 &arg) ->decltype(LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg)) \
  {									\
    return     LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg); \
  }

GridUnopClass(UnaryTraceColour, traceIndex<ColourIndex>(a));
GridUnopClass(UnaryTraceSpin, traceIndex<SpinIndex>(a));

GRID_DEF_UNOP(traceColour, UnaryTraceColour);
GRID_DEF_UNOP(traceSpin, UnaryTraceSpin);

#undef GRID_UNOP
#undef GRID_DEF_UNOP

NAMESPACE_END(Grid);


BEGIN_HADRONS_NAMESPACE

/*
 * 3pt contraction with gamma matrix insertion.
 *
 * Schematic:
 *
 *                   q2           q3
 *              /----<------*------<----¬
 *             /           gV            \
 *            /                           \
 *   i, gSrc *                            * f, gSnk
 *            \                          /
 *             \                        /
 *              \----------->----------/
 *                          q1
 *
 *      trace(gSnk*q1Snk*(adj(gSrc)*g5)*adj(q2)*(g5*gV)*q3)
 *
 *  options:
 *   - q1: sink smeared propagator, source at i
 *   - q2: propagator, source at i
 *   - q3: propagator, source at f
 *   - gamma: gamma matrices to insert as vector of space-separated strings,
 *                  gSnk     gV   gSrc
 *              {"GammaT GammaX GammaY",
 *               "Gamma5    all Gamma5",
 *               "GammaX GammaY GammaT"}
 *   - tSnk: sink position for propagator q1.
 *   - momProjector: Name of module to do momentum projection.
 *                   If empty, project to zero momentum only.
 *   - output: Name of output file.
 *             If empty, do not save.
 *
 */

 /******************************************************************************
 *                               Gamma3pt                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class Gamma3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Gamma3ptPar,
                                    std::string,                q1,
                                    std::string,                q2,
                                    std::string,                q3,
                                    std::vector<std::string>,   gamma,
                                    unsigned int,               tSnk,
                                    std::string,                momProjector,
                                    std::string,                output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TGamma3pt: public Module<Gamma3ptPar>
{
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    using SpinField = SpinMatrixField1;
    using PhaseField = ComplexField1;
    using MPType = MomentumProject<PhaseField, PhaseField>;

public:
    class GammaTriad: Serializable
    {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(GammaTriad,
                                            Gamma::Algebra, vertex,
                                            Gamma::Algebra, source,
                                            Gamma::Algebra, sink);
    };
    class Result: Serializable
    {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                std::vector<GammaTriad>, gamma,
                std::vector<std::vector<int>>, momentum,
                std::vector<std::vector<std::vector<Complex>>>, corr);
    };

    // constructor
    TGamma3pt(const std::string name);
    // destructor
    virtual ~TGamma3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
        std::vector<GammaTriad> allGammaComb_;          // Vector of all gamma combinations
        std::vector<Gamma::Algebra> allGammaSinks_;     // Vector of all gammas at the sink
        std::vector<Gamma::Algebra> allGammaVertices_;  // Vector of all gammas at the vertex
        std::vector<Gamma::Algebra> allGammaSources_;   // Vector of all gammas at the source
        bool isLegacy_;
};

MODULE_REGISTER_TMP(Gamma3pt, ARG(TGamma3pt<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                       TGamma3pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TGamma3pt<FImpl1, FImpl2, FImpl3>::TGamma3pt(const std::string name)
: Module<Gamma3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().q3};

    if (!par().momProjector.empty())
    {
        in.push_back(par().momProjector);
        in.push_back(par().momProjector + "_momList");
    }

    return in;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::setup(void)
/*
 * Fill allGammaComb_ with all combinations implied by par().gamma.
 *
 * Input example 1:
 *      If strGammas = 'GammaX GammaY Gamma5', the function adds one entry
 *      to allGammaComb_, GammaTriad(source=Gamma5, vertex=GammaY, sink=GammaX).
 *
 * Input example 2:
 *      If strGammas = 'all GammaY all', the function adds one entry per
 *      gamma function combination (cartesian product of every gamma matrix at source
 *      and sink). 'all' can be included either at source, vertex, and/or sink.
 */
{
    isLegacy_ = par().momProjector.empty();

    const unsigned int nd = env().getNd() - 1;

    auto parseInput = [](const std::string strGammas)
    {
        // strGammas = 'sink vertex source' -> listGammas = ['sink', 'vertex', 'source']
        std::vector<std::string> listGammas;
        std::istringstream iss(strGammas);
        std::string gamma;
        while (iss >> gamma)
        {
            listGammas.push_back(gamma);
        }
        // Check the input/parsing is correct
        if (listGammas.size() > 3)
        {
            HADRONS_ERROR(Argument, "Too many gamma matrices! Provide sink vertex source gammas");
        }
        else if (listGammas.size() < 3)
        {
            HADRONS_ERROR(Argument, "Too few gamma matrices! Provide sink vertex source gammas");
        }
        return listGammas;
    };

    auto expandAll = [](const std::string gamma)
    {
        std::vector<Gamma::Algebra> v;
        if (gamma == "all")
        {
            // Do all gamma matrices.
            for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
            {
                v.push_back((Gamma::Algebra)i);
            }
        }
        else
        {
            // Parse individual contractions from input string.
            for (auto g : strToVec<Gamma::Algebra>(gamma))
            {
                v.push_back(g);
            }
        }
        return v;
    };

    auto cartesianProduct = [](
        const std::vector<Gamma::Algebra> &sinks,
        const std::vector<Gamma::Algebra> &vertices,
        const std::vector<Gamma::Algebra> &sources)
    {
        // Cartesian product
        std::vector<GammaTriad> list;
        for (auto gSrc : sources)
        {
            for (auto gSnk : sinks)
            {
                for (auto gVtx : vertices)
                {
                    GammaTriad t;
                    t.source = gSrc;
                    t.vertex = gVtx;
                    t.sink = gSnk;
                    list.push_back(t);
                }
            }
        }
        return list;
    };

    /////////////////////////////////////////////////////////////////////////////
    // Expand 'all' and create full list of gamma matrix combinations to compute
    /////////////////////////////////////////////////////////////////////////////

    allGammaComb_.clear();
    allGammaSinks_.clear();
    allGammaSources_.clear();
    allGammaVertices_.clear();
    for (auto strGammas : par().gamma)
    {
        std::vector<Gamma::Algebra> sinks;
        std::vector<Gamma::Algebra> vertices;
        std::vector<Gamma::Algebra> sources;

        std::vector<std::string> gammaTriad = parseInput(strGammas);

        // Add gammas to sets used in sink, vertex, and source of contraction

        for (const auto &gSnk : expandAll(gammaTriad[0]))
        {
            sinks.push_back(gSnk);
            if (std::find(allGammaSinks_.begin(), allGammaSinks_.end(), gSnk) == allGammaSinks_.end())
            {
                allGammaSinks_.push_back(gSnk);
            }
        }

        for (const auto &gV : expandAll(gammaTriad[1]))
        {
            vertices.push_back(gV);
            if (std::find(allGammaVertices_.begin(), allGammaVertices_.end(), gV) == allGammaVertices_.end())
            {
                allGammaVertices_.push_back(gV);
            }
        }

        for (const auto &gSrc : expandAll(gammaTriad[2]))
        {
            sources.push_back(gSrc);
            if (std::find(allGammaSources_.begin(), allGammaSources_.end(), gSrc) == allGammaSources_.end())
            {
                allGammaSources_.push_back(gSrc);
            }
        }

        for (const auto& e : cartesianProduct(sinks, vertices, sources))
        {
            if (std::find(allGammaComb_.begin(), allGammaComb_.end(), e) == allGammaComb_.end())
            {
                allGammaComb_.push_back(e);
            }
        }
    }

    /////////////////////
    // Allocate variables
    /////////////////////
    {
        envTmpLat(PropagatorField1, "q2Adj");
        envTmpLat(PropagatorField1, "prodSrc");
        envTmpLat(SpinField, "spinField");
        envTmpLat(PhaseField, "contraction");

        envCreate(HadronsSerializable, getName(), 1, 0);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    const unsigned int nt = env().getDim(Tp);

    Gamma g5(Gamma::Algebra::Gamma5);
    auto &out = envGet(HadronsSerializable, getName());
    auto &q1 = envGet(SlicedPropagator1, par().q1);
    auto &q2 = envGet(PropagatorField2, par().q2);
    auto &q3 = envGet(PropagatorField3, par().q3);
    SitePropagator1 q1Snk = q1[par().tSnk];

    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                << " quarks '" << par().q1 << "' (spectator), '" << par().q2 << "', and '"
                << par().q3 << "'" << std::endl;
    LOG(Message) << "   Total: " << allGammaComb_.size() << " correlator(s)" << std::endl;

    // Allocate space for results, then index assign them

    Result result;
    result.gamma = allGammaComb_;
    result.corr.resize(allGammaComb_.size());

    if (isLegacy_)
    {
        result.momentum = {{0, 0, 0}};
    }
    else
    {
        result.momentum = envGet(std::vector<std::vector<int>>, par().momProjector + "_momList");
    }

    envGetTmp(PropagatorField1, q2Adj);
    envGetTmp(PropagatorField1, prodSrc);
    envGetTmp(SpinField, spinField);
    envGetTmp(PhaseField, contraction);
    std::vector<typename PhaseField::scalar_object> fourierModes;

    startTimer("Trace");
    q2Adj = adj(q2);
    stopTimer("Trace");

    // Outer loops reuse partial products to reduce multiplications
    for (auto iSrc : allGammaSources_)
    {
        startTimer("Trace");
        Gamma gSrc(iSrc);
        prodSrc = q1Snk * (adj(gSrc) * g5) * q2Adj;
        stopTimer("Trace");

        for (auto iV : allGammaVertices_)
        {
            startTimer("Trace");
            Gamma gV(iV);
            spinField = traceColour(prodSrc * (g5 * gV) * q3);
            stopTimer("Trace");

            for (auto iSnk : allGammaSinks_)
            {
                // Locate the index of the (iSnk, iV, iSrc) combination
                GammaTriad current;
                current.sink = iSnk;
                current.vertex = iV;
                current.source = iSrc;

                auto loc = std::find(allGammaComb_.begin(), allGammaComb_.end(), current);
                if (loc == allGammaComb_.end())
                {
                    continue;
                }
                unsigned int gammaIdx = static_cast<unsigned int>(std::distance(allGammaComb_.begin(), loc));

                LOG(Message) << "Source: " << iSrc << " | "
                            << "Vertex: " << iV << " | "
                            << "Sink: "   << iSnk   << std::endl;

                startTimer("Trace");
                Gamma gSnk(iSnk);
                contraction = traceSpin(gSnk * spinField);
                stopTimer("Trace");

                startTimer("FT");
                if (isLegacy_)
                {
                    sliceSum(contraction, fourierModes, Tp);
                }
                else
                {
                    auto &momProjector = envGet(MPType, par().momProjector);
                    momProjector.Project(contraction, fourierModes);
                }
                stopTimer("FT");

                startTimer("Result");
                result.corr[gammaIdx].resize(result.momentum.size());
                for (unsigned int m = 0; m < result.momentum.size(); ++m)
                {
                    const auto *first = fourierModes.data() + nt * m;
                    const auto *last  = first + nt;
                    result.corr[gammaIdx][m].assign(first, last);
                }
                stopTimer("Result");
            }
        }
    }
    startTimer("I/O");
    out = result;
    saveResult(par().output, "gamma3pt", result);
    stopTimer("I/O");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Gamma3pt_hpp_
