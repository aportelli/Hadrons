#ifndef GRID_FWTFILTERS_H_
#define GRID_FWTFILTERS_H_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

struct DwtFilter
{
    const std::vector<Real> fwdL, fwdH, bwdL, bwdH;
};

namespace DwtFilters
{
    extern DwtFilter haar;
    extern DwtFilter db2;
    extern DwtFilter db3;
    extern DwtFilter db4;
    extern DwtFilter db5;
    extern DwtFilter db6;
    extern DwtFilter bior13;
    extern DwtFilter bior15;
    extern DwtFilter bior22;
    extern DwtFilter bior24;
    extern DwtFilter bior31;
    extern DwtFilter bior33;
    extern DwtFilter bior35;

    extern std::map<std::string, const DwtFilter *> fromName;
}

END_HADRONS_NAMESPACE

#endif
