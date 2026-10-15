#include <Hadrons/Modules/MContraction/QEDBurgerShort.hpp>
#include <cmath>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

template class Grid::Hadrons::MContraction::TQEDBurgerShort<FIMPL, FIMPL::PropagatorField, vComplex>;
template class Grid::Hadrons::MContraction::TQEDBurgerShort<FIMPL, FIMPL::FermionField,    vComplex>;


// culledDiscreteRadialShell
std::vector<Coordinate> QEDBurgerShortFullSiteGenerator::culledDiscreteRadialShell(int squared_radius)
{
    return this->discreteRadialShell(squared_radius, [](const Grid::Coordinate& site){ return true; });
}

std::vector<Coordinate> QEDBurgerShortParitySymmetrySiteGenerator::culledDiscreteRadialShell(int squared_radius)
{
    return this->discreteRadialShell(squared_radius, [this](const Grid::Coordinate& site)
    { 
        // This is a function used to divide the lattice in half, where each half contains sites that are
        // the Nd-reflection of the other half (e.g. (1, 2, -3, 4) and (-1, -2, 3, -4)).
        // We will achieve this by defining a plane-of-separation that imposes this reflection, and then a further
        // Nd-1 separation planes to sort sites that lie in the primary separation plane into either half.

        // The separations are simply done by projecting the site onto the plane normal to see if it is
        // 'above' or 'below' the plane.

        // First do the primary plane-of-separation.
        // We'll use the plane with a normal pointing into the orthant with a fully-positive signature to define
        // the half of the lattice we will use.
        // Any of the planes orthogonal to this one would also work.
        int proj = 0;
        for (int nu=0; nu < this->Nd; ++nu)
            proj += site[nu];
        
        if (proj > 0)
            return true;
        else if (proj < 0)
            return false;

        // The above projection will be 0 if the site lies in the primary separation plane...
        // so in this loop we will decide which half of the lattice these sites will count as being part of.
        // Flipping a single axis of the primary plane-of-separation Nd-1 times guarantees disambiguation.
        for (int mu=1; mu < this->Nd; ++mu)
        {
            int disamb_proj = proj - site[mu]*2;
            
            if (disamb_proj > 0)
                return true;
            else if (disamb_proj < 0)
                return false;
        }

        // If a site gets this far, it is the origin.
        // We will count it, but need to remember that it has a symmetry factor of 1, and not 2.
        return true;
    });
}

std::vector<Coordinate> QEDBurgerShortOrthantSiteGenerator::culledDiscreteRadialShell(int squared_radius)
{
    return this->discreteRadialShell(squared_radius, [this](const Grid::Coordinate& site)
    {
        // Select the site if it is in (or on the boundary of) the orthant with
        // a fully-positive signature.
        // This is an arbitrary orthant to choose, but easy to write the cull condition for.
        for (int mu=0; mu < this->Nd; ++mu)
        {
            if (site[mu] < 0)
                return false;
        }
        return true;
    });
}

std::vector<Coordinate> QEDBurgerShortOctahedral3DSiteGenerator::culledDiscreteRadialShell(int squared_radius)
{
    return this->discreteRadialShell(squared_radius, [this](const Grid::Coordinate& site)
    {
        int dims = this->Nd - 1;
        std::vector<int> desc_site;
        for (int site_i=0; site_i < dims; ++site_i)
            desc_site.push_back(site[site_i]);
        desc_site.push_back(0); // Used to enforce that the site is in a single orthant (i.e. all coords >= 0)

        // Select the Nd! subsection of the orthant.
        // This is defined by (in 3D), any rearrangement of x >= y >= z, which in 3D gives 3! = 6 subsections.
        for (int mu=0; mu < dims; ++mu)
            if (desc_site[mu] < desc_site[mu+1])
                return false;
        return site[dims] >= 0;
    });
}

std::vector<Coordinate> QEDBurgerShortOctahedral4DSiteGenerator::culledDiscreteRadialShell(int squared_radius)
{
    return this->discreteRadialShell(squared_radius, [this](const Grid::Coordinate& site)
    {
        std::vector<int> desc_site;
        for (int c : site)
            desc_site.push_back(c);
        desc_site.push_back(0); // Used to enforce that the site is in a single orthant (i.e. all coords >= 0)

        // Select the Nd! subsection of the orthant.
        // This is defined by (in 4D), any rearrangement of x >= y >= z >= t, which in 4D gives 4! = 24 subsections.
        for (int mu=0; mu < this->Nd; ++mu)
            if (desc_site[mu] < desc_site[mu+1])
                return false;
        return true;
    });
}

// Symmetry Factors
double QEDBurgerShortFullSiteGenerator::symmetryFactor(const Coordinate& site)
{
    return 1;
}

double QEDBurgerShortParitySymmetrySiteGenerator::symmetryFactor(const Coordinate& site)
{
    for (const auto d : site)
    {
        if (d != 0)
            return 2;
    }
    return 1;   
}

double QEDBurgerShortOrthantSiteGenerator::symmetryFactor(const Coordinate& site)
{
    double factor = 1;
    for (int mu=0; mu < this->Nd; ++mu)
    {
        if (site[mu]!=0)
            factor *= 2;
    }
    return factor;
}

double QEDBurgerShortOctahedral3DSiteGenerator::symmetryFactor(const Coordinate& site)
{
    // This is a further reduction over the orthant case.
    // First get the symmetry factor from dividing up a single orthant.
    // The lines of symmetry are along lines of equal coefficients.
    // Therefore, you pick up a duplicate when a site sits on one of these lines because it enters the definition
    // of all q-orthants sharing that line of symmetry. This compounds factorially.
    // e.g. a site at 0,0,1 has 2!=2 duplicates from the 2 permutations of (0,0) and 1 permutation of (1)
    // e.g. a site at 0,1,1 has 2!=2 duplicates from the 2 permutations of (1,1) and 1 permutation of (0)
    // We can convert this into a symmetry factor by considering all 3! q-orthants and dividing out the duplicates.
    double factor = 1.;
    // Init factor as Nd!
    {
        for (int order=1; order <= this->Nd-1; ++order) 
            factor *= order;
    }
    
    // Now find the reduction from duplicates within the q-orthant
    {
        std::vector<double> counts(this->Nd-1, 1.);

        // First count all unique pairs of equal coefficients
        // e.g. this constructs [2,1,1] from [0,0,1]
        for (int i=0;   i < this->Nd-1; ++i)
        for (int j=i+1; j < this->Nd-1; ++j)
        {
            if (abs(site[i]) == abs(site[j])) // Lines of constant coeffs don't care which orthant you live in
                counts[i] += 1.;
        }

        // Multiplying together the 'counts' recovers the factorials
        // i.e. [0,0,1] -> [2,1,1] -> 2! x 1!
        // i.e. [0,1,1] -> [1,2,1] -> 1! x 2!
        for (int i=0; i < this->Nd-1; ++i)
            factor /= counts[i];
    }

    // Now find the reduction from using a single orthant
    for (int mu=0; mu < this->Nd; ++mu)
    {
        if (site[mu]!=0)
            factor *= 2;
    }
    return factor;
}

double QEDBurgerShortOctahedral4DSiteGenerator::symmetryFactor(const Coordinate& site)
{
    // This is a further reduction over the orthant case.
    // First get the symmetry factor from dividing up a single orthant.
    // The lines of symmetry are along lines of equal coefficients.
    // Therefore, you pick up a duplicate when a site sits on one of these lines because it enters the definition
    // of all q-orthants sharing that line of symmetry. This compounds factorially.
    // e.g. a site at 0,0,0,1 has 3!=6 duplicates from the 6 permutations of (0,0,0) and 1 permutation of (1)
    // e.g. a site at 0,0,1,1 has 2!*2! duplicates from the product of the 2 permutations of (0,0) and of (1,1)
    // We can convert this into a symmetry factor by considering all 4! q-orthants and dividing out the duplicates.
    double factor = 1.;
    // Init factor as Nd!
    {
        for (int order=1; order <= this->Nd; ++order) 
            factor *= order;
    }
    
    // Now find the reduction from duplicates within the q-orthant
    {
        std::vector<double> counts(this->Nd, 1.);

        // First count all unique pairs of equal coefficients
        // e.g. this constructs [3,2,1,1] from [0,0,0,1]
        for (int i=0;   i < this->Nd; ++i)
        for (int j=i+1; j < this->Nd; ++j)
        {
            if (abs(site[i]) == abs(site[j])) // Lines of constant coeffs don't care which orthant you live in
                counts[i] += 1.;
        }

        // Multiplying together the 'counts' recovers the factorials
        // i.e. [0,0,0,1] -> [3,2,1,1] -> 3! x 1!
        // i.e. [0,0,1,1] -> [2,1,2,1] -> 2! x 2!
        for (int i=0; i < this->Nd; ++i)
            factor /= counts[i];
    }

    // Now find the reduction from using a single orthant
    for (int mu=0; mu < this->Nd; ++mu)
    {
        if (site[mu]!=0)
            factor *= 2;
    }
    return factor;
}
