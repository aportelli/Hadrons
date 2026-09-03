/*
* MomentumProjectorSphere.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
*
* Copyright (C) 2015 - 2023
*
* Author: Teseo San Jose <teseo.sanjose@ed.ac.uk>
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

#ifndef Hadrons_MUtilities_MomentumProjectorSphere_hpp_
#define Hadrons_MUtilities_MomentumProjectorSphere_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              MomentumProjectorSphere                      *
 *****************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class MomentumProjectorSpherePar: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumProjectorSpherePar,
                                        unsigned int, maxFourier);
};

template<typename Field, typename ComplexField>
class TMomentumProjectorSphere: public Module<MomentumProjectorSpherePar>
{
public:
    using MPType = MomentumProject<Field, ComplexField>;
    // constructor
    TMomentumProjectorSphere(const std::string name);
    // destructor
    virtual ~TMomentumProjectorSphere(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    //static methods
    static void fillMomList(
        std::vector<std::vector<int>> &momList, // TODO: Change to Grid::Coordinate?
        const unsigned int maxFourier,
        const unsigned int nd
    );
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MPScalar, ARG(TMomentumProjectorSphere<FIMPL::ComplexField, FIMPL::ComplexField>), MUtilities);

/******************************************************************************
 *                       TMomentumProjectorSphere implementation              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template<typename Field, typename ComplexField>
TMomentumProjectorSphere<Field, ComplexField>::TMomentumProjectorSphere(const std::string name)
: Module<MomentumProjectorSpherePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template<typename Field, typename ComplexField>
std::vector<std::string> TMomentumProjectorSphere<Field, ComplexField>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template<typename Field, typename ComplexField>
std::vector<std::string> TMomentumProjectorSphere<Field, ComplexField>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_momList", getName() + "_phases"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template<typename Field, typename ComplexField>
void TMomentumProjectorSphere<Field, ComplexField>::setup(void)
/*
 * Fill list of momenta to compute. Allocate phases, initialize classes.
 */
{
    envTmpLat(ComplexField, "coor");

    // MomentumProject instance
    env().template createObject<MPType>(getName(), Environment::Storage::standard, 0);

    // List of momenta
    const unsigned int nd = env().getNd() - 1;
    envCreate(std::vector<std::vector<int>>, getName() + "_momList", 1, 0);
    auto &momList = envGet(std::vector<std::vector<int>>, getName() + "_momList");
    fillMomList(momList, par().maxFourier, nd);

    // Phases
    envCreate(std::vector<ComplexField>, getName() + "_phases", 1, momList.size(), envGetGrid(ComplexField));
}

// static /////////////////////////////////////////////////////////////////////
template<typename Field, typename ComplexField>
void TMomentumProjectorSphere<Field, ComplexField>::fillMomList(
    std::vector<std::vector<int>> &momList,
    const unsigned int maxFourier,
    const unsigned int nd)
{
    auto indexToMom = [](std::vector<int> &mom, const unsigned int i, const unsigned int maxFourier, const unsigned int nd)
    {
        // Assign momentum given index following the lexicographic rules. Row major indexing.
        const unsigned int size = 2*maxFourier + 1;

        mom.resize(nd);
        unsigned int buf = i;
        for (int mu = nd - 1; mu >= 0; --mu)
        {
            mom[mu]  = static_cast<int>(buf % size) - static_cast<int>(maxFourier);
            buf     /= size;
        }
    };

    auto euclideanDistance = [](const std::vector<int> &mom)
    {
        // Return Euclidean distance of a momentum 3-vector.
        double mod = 0;
        for (auto p : mom)
        {
            mod += static_cast<double>(p * p);
        }
        return sqrt(mod);
    };

    // Total number of Fourier modes in a cube
    unsigned int nMom = 1;
    const unsigned int momSize = 2*maxFourier + 1;
    for (unsigned int d = 0; d < nd; ++d)
    {
        nMom *= momSize;
    }

    // Gather only modes inside a sphere of radius <= maxFourier
    for (unsigned int m = 0; m < nMom; ++m)
    {
        std::vector<int> mom;
        indexToMom(mom, m, maxFourier, nd);
        if (euclideanDistance(mom) <= static_cast<double>(maxFourier))
        {
            momList.push_back(mom);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template<typename Field, typename ComplexField>
void TMomentumProjectorSphere<Field, ComplexField>::execute(void)
{
    Complex i(0.0,1.0);
    envGetTmp(ComplexField, coor);

    auto &MP = envGet(MPType, getName());
    auto &phases = envGet(std::vector<ComplexField>, getName() + "_phases");
    auto &momList = envGet(std::vector<std::vector<int>>, getName() + "_momList");

    startTimer("Phases");
    for (unsigned int m = 0; m < momList.size(); ++m)
    {
        ComplexField& ph = phases[m];
        ph = Zero();
        const auto &p = momList[m];
        for(unsigned int mu = 0; mu < p.size(); ++mu)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (static_cast<Real>(p[mu]) / static_cast<Real>(env().getDim(mu))) * coor;
        }
        ph = exp(-static_cast<Real>(2*M_PI)*i*ph);
    }
    stopTimer("Phases");

    startTimer("ImportPhases");
    MP.Allocate(momList.size(), env().getGrid());
    MP.ImportMomenta(phases);
    stopTimer("ImportPhases");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif  // Hadrons_MUtilities_MomentumProjectorSphere_hpp_