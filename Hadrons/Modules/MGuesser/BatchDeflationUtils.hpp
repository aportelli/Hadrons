#ifndef Hadrons_MGuesser_BatchDeflationUtils_hpp_
#define Hadrons_MGuesser_BatchDeflationUtils_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

namespace BatchDeflationUtils
{
    template <typename Field>
    void projAccumulate(const std::vector<Field> &in, std::vector<Field> &out,
                        const std::vector<Field>& evec,
                        const std::vector<RealD>& eval,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf);
};

template <typename Field>
void BatchDeflationUtils::projAccumulate(const std::vector<Field> &in, std::vector<Field> &out,
                                         const std::vector<Field>& evec,
                                         const std::vector<RealD>& eval,
                                         const unsigned int ei, const unsigned int ef,
                                         const unsigned int si, const unsigned int sf)
{
    GridBase *g       = in[0].Grid();
    double   lVol     = g->lSites();
    double   siteSize = sizeof(typename Field::scalar_object);
    double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
    double   nIt      = (ef - ei)*(sf - si);
    double   t        = 0.;

    t -= usecond();
    for (unsigned int i = ei; i < ef; ++i)
    for (unsigned int j = si; j < sf; ++j)
    {
        axpy(out[j], 
            TensorRemove(innerProduct(evec[i], in[j]))/eval[i], 
            evec[i], out[j]);
    }
    t += usecond();
    // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
    LOG(Message) << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB 
                 << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchDeflationUtils_hpp_
