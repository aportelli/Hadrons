#ifndef Hadrons_MContraction_QEDBurgerShort_hpp_
#define Hadrons_MContraction_QEDBurgerShort_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDBurgerShort                                 *
 ******************************************************************************/
/*
 Universal disconnected "Burger" loop subdiagram
              q
             ___ 
            /   \
           |~~~~~| photon
            \___/
              q

 Tr[q(x,y) * Gamma_{mu} * q(y,x) * Gamma_{nu}] * G^{mu,nu}(x, y)

*/


BEGIN_MODULE_NAMESPACE(MContraction)

class QEDBurgerShortPar: Serializable
{
public:
    GRID_SERIALIZABLE_ENUM(SymmetryMode, undef, none, 0, exact, 1, approximate, 2);

    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDBurgerShortPar,
                                    std::string,  sources,
                                    std::string,  qs,
                                    std::string,  emField,
                                    unsigned int, radius,
                                    SymmetryMode, useLatticeSymmetry,
                                    std::string,  output);
};

template <typename FImpl, typename VType>
class TQEDBurgerShort: public Module<QEDBurgerShortPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    
    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::GaugeField EmField;
    
    // constructor
    TQEDBurgerShort(const std::string name);
    // destructor
    virtual ~TQEDBurgerShort(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void) override;
    virtual std::vector<std::string> getOutput(void) override;
    virtual std::vector<std::string> getOutputFiles(void) override;
    // setup
    virtual void setup(void) override;
    // execution
    virtual void execute(void) override;

    inline std::vector<Coordinate> generateShiftVectors(int rsqmin, int radius, QEDBurgerShortPar::SymmetryMode useLatticeSymmetry) const;
    inline double symmetryFactor(const Coordinate& site, QEDBurgerShortPar::SymmetryMode useLatticeSymmetry) const;
    inline bool isInUsedLatticeHalf(const Coordinate& site) const;
    inline bool isInUsedLatticeOrthant(const Coordinate& site) const;
    inline bool isInUsedLatticeOrthantPair(const Coordinate& site) const;
    inline void fastBurger(const PropagatorField& left, const PropagatorField& right, const typename EmField::scalar_object& pSite, LatticeComplexD& out) const;
    inline void coordCshift(const PropagatorField& field, const Coordinate& coord, PropagatorField& out) const;
};

MODULE_REGISTER_TMP(QEDBurgerShort, ARG(TQEDBurgerShort<FIMPL, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDBurgerShort implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
TQEDBurgerShort<FImpl, VType>::TQEDBurgerShort(const std::string name)
: Module<QEDBurgerShortPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, VType>::getInput(void)
{
    std::vector<std::string> in = {par().sources, par().emField, par().qs};
    
    return in;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, VType>::getOutputFiles(void)
{
    std::vector<std::string> output = {};
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::setup(void)
{
    // Contraction temporaries
    envTmp(FFT,            "fft",            1, env().getGrid());
    envTmp(EmField,        "Gx",             1, env().getGrid());
    envTmp(EmField,        "em_buffer",      1, env().getGrid());
    envTmp(LatticeComplexD, "tmp_cbuffer",   1, env().getGrid());
    envTmp(PropagatorField, "tmp_prop1",     1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_prop2",     1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_q1",        1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_q2",        1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "shifted_quark", 1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "shifted_noise", 1, envGetGrid(PropagatorField));
    // Output
    envCreate(HadronsSerializable, getName(), 1, 0);
    envCreate(HadronsSerializable, getName()+"_full", 1, 0);
    envCreate(HadronsSerializable, getName()+"_bias", 1, 0);
}


template<typename FImpl, typename VType>
std::vector<Coordinate> TQEDBurgerShort<FImpl, VType>::generateShiftVectors(int rsqmin, int radius, QEDBurgerShortPar::SymmetryMode useLatticeSymmetry) const
{
    int rsqmax = radius*radius;

    // Generate list of sites
    std::vector<Coordinate> shifts;
    // Make this depend on Nd rather than hardcode to Nd=4..?
    for (int r0=-radius; r0<=radius; ++r0)
    for (int r1=-radius; r1<=radius; ++r1)
    for (int r2=-radius; r2<=radius; ++r2)
    for (int r3=-radius; r3<=radius; ++r3)
    {
        int rsq = r0*r0+r1*r1+r2*r2+r3*r3;

        if ((rsqmin >= rsq) || (rsq > rsqmax)) continue;

        Coordinate r(Nd);
        r[0]=r0;
        r[1]=r1;
        r[2]=r2;
        r[3]=r3;

        switch (useLatticeSymmetry)
        {
            case QEDBurgerShortPar::SymmetryMode::none:
            {
                shifts.push_back(r);
                break;
            }
            case QEDBurgerShortPar::SymmetryMode::exact:
            {
                if (this->isInUsedLatticeHalf(r))
                    shifts.push_back(r);
                break;
            }
            case QEDBurgerShortPar::SymmetryMode::approximate:
            {
                if (this->isInUsedLatticeOrthant(r))
                    shifts.push_back(r);
                break;
            }
        }
    }
    return shifts;
}

template<typename FImpl, typename VType>
double TQEDBurgerShort<FImpl, VType>::symmetryFactor(const Coordinate& s, QEDBurgerShortPar::SymmetryMode useLatticeSymmetry) const
{
    switch (useLatticeSymmetry)
    {
        case QEDBurgerShortPar::SymmetryMode::none:
        {
            return 1;
        }
        case QEDBurgerShortPar::SymmetryMode::exact:
        {
            // Only returns 1 if all elements of 'site' are 0.
            for (const auto d : s)
            {
                if (d != 0)
                    return 2;
            }
            return 1;
        }
        case QEDBurgerShortPar::SymmetryMode::approximate:
        {
            // A site will have as many duplicates as elements that are not 0, since the vector is reflected
            // in the coordinate axes.
            double factor = 1;
            for (int mu=0; mu < Nd; ++mu)
            {
                if (s[mu]!=0)
                    factor *= 2;
            }
            return factor;
        }
        default:
        {
            HADRONS_ERROR(Argument, "invalid symmetry mode");
        }
    }
}

template<typename FImpl, typename VType>
bool TQEDBurgerShort<FImpl, VType>::isInUsedLatticeHalf(const Coordinate& site) const
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
    for (int nu=0; nu < Nd; ++nu)
        proj += site[nu];
    
    if (proj > 0)
        return true;
    else if (proj < 0)
        return false;

    // The above projection will be 0 if the site lies in the primary separation plane...
    // so in this loop we will decide which half of the lattice these sites will count as being part of.
    // Flipping a single axis of the primary plane-of-separation Nd-1 times guarantees disambiguation.
    for (int mu=1; mu < Nd; ++mu)
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
}

template<typename FImpl, typename VType>
bool TQEDBurgerShort<FImpl, VType>::isInUsedLatticeOrthant(const Coordinate& site) const
{
    // Select the site if it is in (or on the boundary of) the orthant with
    // a fully-positive signature.
    // This is an arbitrary orthant to choose, but easy to write the cull condition for.
    for (int mu=0; mu < Nd; ++mu)
    {
        if (site[mu] < 0)
            return false;
    }
    return true;
}

template<typename FImpl, typename VType>
bool TQEDBurgerShort<FImpl, VType>::isInUsedLatticeOrthantPair(const Coordinate& site) const
{
    // Select the site if it is in (or on the boundary of) the orthant with
    // a fully-positive signature.
    // This is an arbitrary orthant to choose, but easy to write the cull condition for.
    for (int mu=0; mu < 3; ++mu)
    {
        if (site[mu] < 0)
            return false;
    }
    return true;
}

template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::fastBurger(const PropagatorField& left, const PropagatorField& right, const typename EmField::scalar_object& pSite, LatticeComplexD& out) const
{
    Gamma::Algebra Gmu[] = 
    {
        (Gamma::Algebra::GammaX),
        (Gamma::Algebra::GammaY),
        (Gamma::Algebra::GammaZ),
        (Gamma::Algebra::GammaT),
    };

    autoView(vleft,left,AcceleratorRead);
    autoView(vright,right,AcceleratorRead);
    autoView(wvout,out,AcceleratorWrite);

    constexpr int Nc = FImpl::Dimension;
    constexpr int Ns = Grid::Ns;
    const auto& grid = envGetGrid4(PropagatorField);

    // Feynman gauge photon propagator contraction.
    // Tr[G_mu * left * G_mu * right] * Guv[site] ( * i^2, taken care of later )
    accelerator_for(s, grid->oSites(), VType::Nsimd(), 
    {
        const auto& s1 = vleft[s];
        const auto& s2 = vright[s];
        LatticeComplexD::vector_object out_s = Zero();
        for (int mu=0; mu<Nd; ++mu)
        {
            LatticeComplexD::vector_object tmp = Zero();
            const auto& gs1 = Gamma(Gmu[mu])*s1;
            const auto& gs2 = Gamma(Gmu[mu])*s2;
            for(int si=0;si<Ns;si++)
            for(int sj=0;sj<Ns;sj++)
            for(int ci=0;ci<Nc;ci++)
            for(int cj=0;cj<Nc;cj++)
                tmp()()()+=gs1()(si,sj)(ci,cj)*gs2()(sj,si)(cj,ci);
            out_s += tmp*pSite(mu)()();
        }
        // Now also remember the -1 to account for i^2 from the two Aslash insertions.
        wvout[s] = -1.*out_s;
    })
}

template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::coordCshift(const PropagatorField& field, const Coordinate& coord, PropagatorField& out) const
{
    out = Cshift(field,0,coord[0]);
    for (int mu=1;mu<Nd;mu++) out = Cshift(out,mu,coord[mu]);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::execute(void)
{
    // *********** //
    // PREPARATION //
    // *********** //

    Gamma Gmu[] = 
    {
        Gamma(Gamma::Algebra::GammaX),
        Gamma(Gamma::Algebra::GammaY),
        Gamma(Gamma::Algebra::GammaZ),
        Gamma(Gamma::Algebra::GammaT),
    };

    // Get temps
    envGetTmp(LatticeComplexD,           tmp_cbuffer);
    envGetTmp(PropagatorField,           tmp_prop1);
    envGetTmp(PropagatorField,           tmp_prop2);
    envGetTmp(PropagatorField,           tmp_q1);
    envGetTmp(PropagatorField,           tmp_q2);
    envGetTmp(PropagatorField,           shifted_quark);
    envGetTmp(PropagatorField,           shifted_noise);

    // Get env vars
    // Construct photon prop from weights
    const EmField& em_field = envGet(EmField, par().emField);
    envGetTmp(EmField, Gx);
    envGetTmp(EmField, em_buffer);
    envGetTmp(FFT,     fft);
    fft.FFT_all_dim(em_buffer, em_field, FFT::forward);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,0)*peekLorentz(em_buffer,0), 0);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,1)*peekLorentz(em_buffer,1), 1);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,2)*peekLorentz(em_buffer,2), 2);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,3)*peekLorentz(em_buffer,3), 3);
    fft.FFT_all_dim(Gx, em_buffer, FFT::backward);

    // Get qs
    std::vector<PropagatorField*> qs;
    if (envHasType(std::vector<PropagatorField>, par().qs))
    {
        std::vector<PropagatorField>& envQs = envGet(std::vector<PropagatorField>, par().qs);
        for (int i=0; i < envQs.size(); ++i)
            qs.push_back(&(envQs[i]));
    }
    else
        qs = envGet(std::vector<PropagatorField*>, par().qs);
    
    // Get noises
    std::vector<PropagatorField*> noises;
    if (envHasType(std::vector<PropagatorField>, par().sources))
    {
        std::vector<PropagatorField>& envNoises = envGet(std::vector<PropagatorField>, par().sources);
        for (int i=0; i < envNoises.size(); ++i)
            noises.push_back(&(envNoises[i]));
    }
    else
        noises = envGet(std::vector<PropagatorField*>, par().sources);

    // Get other parameters
    std::vector<int> radii;
    for (int i=0; i <= par().radius; ++i)
        radii.push_back(i);
    int rsqmin = -1;
    int Nsrc   = noises.size();

    // ************************* //
    // CONTRACTION ROUTINE START //
    // ************************* //

    // Measure contraction on each site within the radial limit
    LOG(Message) << "Generating and contracting propagators..." << std::endl;
    
    std::vector<std::string> summary_messages;
    std::vector<std::vector<RealD>> burgers;
    std::vector<std::vector<RealD>> burgers_full;
    std::vector<std::vector<RealD>> burgers_bias;
    std::vector<RealD> burger_buffers(Nsrc, 0.0);
    std::vector<RealD> burger_full_buffers(Nsrc, 0.0);
    std::vector<RealD> burger_bias_buffers(Nsrc, 0.0);
    startTimer("Total Contraction Time");
    for (auto radius : radii)
    {   
        int rsqmax = radius*radius;

        const auto shifts = this->generateShiftVectors(rsqmin, radius, par().useLatticeSymmetry);

        LOG(Message) << "Generating propagators for " << shifts.size() << " shifts inside shell from r^2=" << rsqmin << " to (and including) r^2=" << rsqmax << std::endl;
        // Iterate over sites
        for (int site_i=0; site_i < shifts.size(); ++site_i)
        {
            const auto& r          = shifts[site_i];
            double symmetry_factor = this->symmetryFactor(r, par().useLatticeSymmetry);

            // ***************************************************************************************** //
            // To estimate the burger diagram, we average over traces computed from pairs of propagators
            // solved on different noises. The basic way to do this is to loop over two 'noise indices'
            // and compute the trace for each iteration where the indices are not equal.
            // However, this requires O(Nsrc^2) expensive C-shifted products to compute.
            // We can reduce this to O(Nsrc) C-shifted products if we create two individual
            // propagators summed over Nsrc, and compute the diagram using these two propagators. This
            // estimator would however include terms where the same noise is used for both propagators,
            // and so we should compute this part separately in order to subtract it from the total.
            // ***************************************************************************************** //

            // Reset/allocate temps.
            tmp_prop1 = Zero();
            tmp_prop2 = Zero();
            RealD samenoise_contribution = 0.0;

            // Extract gauge field at offset.
            Coordinate gauge_r(Nd);
            Coordinate latt_size = env().getDim();
            gauge_r[0]=(r[0]+latt_size[0])%latt_size[0];
            gauge_r[1]=(r[1]+latt_size[1])%latt_size[1];
            gauge_r[2]=(r[2]+latt_size[2])%latt_size[2];
            gauge_r[3]=(r[3]+latt_size[3])%latt_size[3];
            
            typename EmField::scalar_object pSite;
            peekSite(pSite, Gx, gauge_r);

            auto prev_value = burger_buffers.back();
            auto prev_full  = burger_full_buffers.back();
            auto prev_bias  = burger_bias_buffers.back();
            for (int hit_i=0;hit_i<Nsrc;hit_i++)
            {
                // Calculate S(x+r, x)
                startTimer("Total Contraction Time[C-Shifts]");
                coordCshift(*qs[hit_i], r, shifted_quark);
                stopTimer("Total Contraction Time[C-Shifts]");
                startTimer("Total Contraction Time[Matrix Products]");
                tmp_q1 = shifted_quark*adj(*noises[hit_i]);
                tmp_prop1 += tmp_q1;
                stopTimer("Total Contraction Time[Matrix Products]");

                // Calculate S(x, x+r)
                startTimer("Total Contraction Time[C-Shifts]");
                coordCshift(*noises[hit_i], r, shifted_noise);
                stopTimer("Total Contraction Time[C-Shifts]");
                startTimer("Total Contraction Time[Matrix Products]");
                tmp_q2 = (*qs[hit_i])*adj(shifted_noise);
                tmp_prop2 += tmp_q2;
                stopTimer("Total Contraction Time[Matrix Products]");
                
                // Calculate same-noise contraction
                startTimer("Total Contraction Time[Same-Noise Contraction]");
                fastBurger(tmp_q1, tmp_q2, pSite, tmp_cbuffer);
                samenoise_contribution += toReal(symmetry_factor*sum(tmp_cbuffer));
                stopTimer("Total Contraction Time[Same-Noise Contraction]");
            
                // Calculate full-noise contraction
                startTimer("Total Contraction Time[Full-Noise Contraction]");
                fastBurger(tmp_prop1, tmp_prop2, pSite, tmp_cbuffer);
                auto full_contribution = toReal(symmetry_factor*sum(tmp_cbuffer));

                // Add the result for this site to the total.
                burger_full_buffers[hit_i] += full_contribution;
                burger_bias_buffers[hit_i] += samenoise_contribution;
                burger_buffers     [hit_i] += full_contribution - samenoise_contribution;
                stopTimer("Total Contraction Time[Full-Noise Contraction]");                
            }
            LOG(Message) << "Shift: " << r << " ["
                         << burger_buffers     .back()-prev_value << "/"
                         << burger_full_buffers.back()-prev_full << "/"
                         << burger_bias_buffers.back()-prev_bias << "], "
                         << "Symmetry: " << symmetry_factor << std::endl;
        }

        // There are Nsrc^2 traces that could be computed using N sources for two propagators;
        // we want them all *except* the Nsrc traces using the same noise for both propagators.
        // This leaves us with (Nsrc^2 - Nsrc) traces we have computed and need to normalise for.
        burgers.push_back({});
        burgers_full.push_back({});
        burgers_bias.push_back({});
        for (int hit_i=0; hit_i < Nsrc; ++hit_i)
        {
            int n_hits = hit_i + 1;
            int norm = (n_hits*(n_hits-1));
            auto normed_burger = burger_buffers[hit_i]/norm;
        
            // Handle outputs
            std::stringstream summary_message;
            summary_message << "burger[radius=" << std::to_string(radius) + ", hits=" << n_hits << "] = " << std::setprecision(15) << normed_burger;
            summary_messages.push_back(summary_message.str());
            burgers.back().push_back(normed_burger);
            burgers_full.back().push_back(burger_full_buffers[hit_i]);
            burgers_bias.back().push_back(burger_bias_buffers[hit_i]);
        }
        rsqmin = radius*radius;
    }
    stopTimer("Total Contraction Time");

    for (const auto& summary_message : summary_messages)
        LOG(Message) << summary_message << std::endl;

    saveResult(par().output, "burgershort", burgers);
    auto& out = envGet(HadronsSerializable, getName());
    out = burgers;

    auto& out_full = envGet(HadronsSerializable, getName()+"_full");
    out = burgers_full;

    auto& out_bias = envGet(HadronsSerializable, getName()+"_bias");
    out = burgers_bias;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDBurgerShort_hpp_
