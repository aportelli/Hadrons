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
    GRID_SERIALIZABLE_ENUM(SymmetryMode, undef, 
        none,          0, 
        parity,        1,
        orthant,       2,
        octahedral,    3,
        octahedral3D,  4
    );

    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDBurgerShortPar,
                                    std::string,               sources,
                                    std::string,               qs,
                                    std::vector<std::string>,  photonProps,
                                    unsigned int,              rSq,
                                    std::vector<SymmetryMode>, shellSymmetries,
                                    std::string,               output);
};

struct QEDBurgerShortSiteGenerator
{
public:
    QEDBurgerShortSiteGenerator(int Nd) : Nd{Nd} {};
    virtual ~QEDBurgerShortSiteGenerator() {}
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) = 0;
    virtual double symmetryFactor(const Coordinate& site) = 0;
protected:
    template<typename CullFunction>
    std::vector<Coordinate> discreteRadialShell(int squared_radius, CullFunction cullFn);
protected:
    int Nd;
};

struct QEDBurgerShortFullSiteGenerator : public QEDBurgerShortSiteGenerator
{
    using QEDBurgerShortSiteGenerator::QEDBurgerShortSiteGenerator;
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) override;
    virtual double symmetryFactor(const Coordinate& site) override;
};

struct QEDBurgerShortParitySymmetrySiteGenerator : public QEDBurgerShortSiteGenerator
{
    using QEDBurgerShortSiteGenerator::QEDBurgerShortSiteGenerator;
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) override;
    virtual double symmetryFactor(const Coordinate& site) override;
};

struct QEDBurgerShortOrthantSiteGenerator : public QEDBurgerShortSiteGenerator
{
    using QEDBurgerShortSiteGenerator::QEDBurgerShortSiteGenerator;
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) override;
    virtual double symmetryFactor(const Coordinate& site) override;
};

struct QEDBurgerShortOctahedral3DSiteGenerator : public QEDBurgerShortSiteGenerator
{
    using QEDBurgerShortSiteGenerator::QEDBurgerShortSiteGenerator;
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) override;
    virtual double symmetryFactor(const Coordinate& site) override;
};

struct QEDBurgerShortOctahedral4DSiteGenerator : public QEDBurgerShortSiteGenerator
{
    using QEDBurgerShortSiteGenerator::QEDBurgerShortSiteGenerator;
    virtual std::vector<Coordinate> culledDiscreteRadialShell(int squared_radius) override;
    virtual double symmetryFactor(const Coordinate& site) override;
};

template<typename CullFunction>
std::vector<Coordinate> QEDBurgerShortSiteGenerator::discreteRadialShell(int squared_radius, CullFunction cullFn)
{
    int radius = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(squared_radius))));
    
    // Generate list of sites
    std::vector<Coordinate> shifts;
    // Make this depend on Nd rather than hardcode to Nd=4..?
    for (int r0=-radius; r0<=radius; ++r0)
    for (int r1=-radius; r1<=radius; ++r1)
    for (int r2=-radius; r2<=radius; ++r2)
    for (int r3=-radius; r3<=radius; ++r3)
    {
        int rsq = r0*r0+r1*r1+r2*r2+r3*r3;

        if (rsq != squared_radius) continue;

        Coordinate r(this->Nd);
        r[0]=r0;
        r[1]=r1;
        r[2]=r2;
        r[3]=r3;

        if (cullFn(r))
            shifts.push_back(r);
    }
    return shifts;
}


template <typename FImpl, typename Field, typename VType>
class TQEDBurgerShort: public Module<QEDBurgerShortPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        int,         rSq,
                                        std::string, symmetry,
                                        std::string, photon,
                                        std::vector<Real>, burger,
                                        std::vector<Real>, full,
                                        std::vector<Real>, bias);
    };

    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::ScalarField PhotonProp;
    size_t numPhotonProps;
    
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

    inline std::unique_ptr<QEDBurgerShortSiteGenerator> createSiteGenerator(QEDBurgerShortPar::SymmetryMode symmetry);
    inline void fastBurger(const PropagatorField& left, const PropagatorField& right, const typename PhotonProp::scalar_object& pSite, LatticeComplexD& out) const;
    inline void coordCshift(const Field& field, const Coordinate& coord, Field& out) const;
    inline PropagatorField pointToPointProp(const PropagatorField& left, const PropagatorField& right);
    inline PropagatorField pointToPointProp(const FermionField&    left, const FermionField&    right);
};

MODULE_REGISTER_TMP(QEDBurgerShort,     ARG(TQEDBurgerShort<FIMPL, FIMPL::PropagatorField, vComplex>), MContraction);
MODULE_REGISTER_TMP(QEDBurgerShortFerm, ARG(TQEDBurgerShort<FIMPL, FIMPL::FermionField,    vComplex>), MContraction);

/******************************************************************************
 *                 TQEDBurgerShort implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
TQEDBurgerShort<FImpl, Field, VType>::TQEDBurgerShort(const std::string name)
: Module<QEDBurgerShortPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, Field, VType>::getInput(void)
{
    std::vector<std::string> in = {par().sources, par().qs};
    for (const std::string& photon_prop : par().photonProps)
        in.push_back(photon_prop);
    
    return in;
}

template <typename FImpl, typename Field, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, Field, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_full", getName() + "_bias"};
    
    return out;
}

template <typename FImpl, typename Field, typename VType>
std::vector<std::string> TQEDBurgerShort<FImpl, Field, VType>::getOutputFiles(void)
{
    std::vector<std::string> output = {};
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShort<FImpl, Field, VType>::setup(void)
{
    numPhotonProps = par().photonProps.size();
    if (numPhotonProps==0)
    {
        HADRONS_ERROR(Definition, "No photon props provided to module");
    }

    if (par().shellSymmetries.empty())
    {
        HADRONS_ERROR(Definition, "No shell symmetries were provided");
    }

    for (auto symm : par().shellSymmetries)
    {
        if (symm == QEDBurgerShortPar::SymmetryMode::undef) // Do not remove if-statement braces: Macro expands to multiple lines
        {
            HADRONS_ERROR(Definition, "Encountered an undefined shell symmetry");
        }
    }

    // Contraction temporaries
    envTmp(FFT,                     "fft",           1,                 env().getGrid());
    envTmp(std::vector<PhotonProp>, "Gxs",           1, numPhotonProps, env().getGrid());
    envTmp(LatticeComplexD,         "tmp_cbuffer",   1,                 env().getGrid());
    envTmp(PropagatorField,         "tmp_prop1",     1,                 envGetGrid(PropagatorField));
    envTmp(PropagatorField,         "tmp_prop2",     1,                 envGetGrid(PropagatorField));
    envTmp(PropagatorField,         "tmp_q1",        1,                 envGetGrid(PropagatorField));
    envTmp(PropagatorField,         "tmp_q2",        1,                 envGetGrid(PropagatorField));
    envTmp(Field,                   "shifted_quark", 1,                 envGetGrid(Field));
    envTmp(Field,                   "shifted_noise", 1,                 envGetGrid(Field));
    // Output
    envCreate(HadronsSerializable, getName(), 1, 0);
}


template <typename FImpl, typename Field, typename VType>
std::unique_ptr<QEDBurgerShortSiteGenerator> TQEDBurgerShort<FImpl, Field, VType>::createSiteGenerator(QEDBurgerShortPar::SymmetryMode symmetry)
{
    int Nd = env().getNd();
    switch (symmetry)
    {
        case (QEDBurgerShortPar::SymmetryMode::none):         return std::make_unique<QEDBurgerShortFullSiteGenerator>          (Nd);
        case (QEDBurgerShortPar::SymmetryMode::parity):       return std::make_unique<QEDBurgerShortParitySymmetrySiteGenerator>(Nd);
        case (QEDBurgerShortPar::SymmetryMode::orthant):      return std::make_unique<QEDBurgerShortOrthantSiteGenerator>       (Nd);
        case (QEDBurgerShortPar::SymmetryMode::octahedral3D): return std::make_unique<QEDBurgerShortOctahedral3DSiteGenerator>  (Nd);
        case (QEDBurgerShortPar::SymmetryMode::octahedral):   return std::make_unique<QEDBurgerShortOctahedral4DSiteGenerator>  (Nd);
        default:                                              { HADRONS_ERROR(Argument, "Invalid symmetry mode"); }
    }
}


template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShort<FImpl, Field, VType>::fastBurger(const PropagatorField& left, const PropagatorField& right, const typename PhotonProp::scalar_object& pSite, LatticeComplexD& out) const
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
            out_s += tmp*pSite;
        }
        // Now also remember the -1 to account for i^2 from the two Aslash insertions.
        wvout[s] = -1.*out_s;
    })
}

template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShort<FImpl, Field, VType>::coordCshift(const Field& field, const Coordinate& coord, Field& out) const
{
    out = Cshift(field,0,coord[0]);
    for (int mu=1;mu<Nd;mu++) out = Cshift(out,mu,coord[mu]);
}


template <typename FImpl, typename Field, typename VType>
typename FImpl::PropagatorField TQEDBurgerShort<FImpl, Field, VType>::pointToPointProp(const PropagatorField& left, const PropagatorField& right)
{
    return left * adj(right);
}

template <typename FImpl, typename Field, typename VType>
typename FImpl::PropagatorField TQEDBurgerShort<FImpl, Field, VType>::pointToPointProp(const FermionField&    left, const FermionField&    right)
{
    return outerProduct(left, right);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShort<FImpl, Field, VType>::execute(void)
{
    // *********** //
    // PREPARATION //
    // *********** //
    Gamma::Algebra Gmu[] = 
    {
        (Gamma::Algebra::GammaX),
        (Gamma::Algebra::GammaY),
        (Gamma::Algebra::GammaZ),
        (Gamma::Algebra::GammaT),
    };
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    // Get temps
    envGetTmp(LatticeComplexD, tmp_cbuffer);
    envGetTmp(PropagatorField, tmp_prop1);
    envGetTmp(PropagatorField, tmp_prop2);
    envGetTmp(PropagatorField, tmp_q1);
    envGetTmp(PropagatorField, tmp_q2);
    envGetTmp(Field,           shifted_quark);
    envGetTmp(Field,           shifted_noise);

    // Get env vars
    // Get position-space photon field
    envGetTmp(std::vector<PhotonProp>, Gxs);
    envGetTmp(FFT,     fft);
    for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
    {
        const PhotonProp& photon_prop = envGet(PhotonProp, par().photonProps[photon_prop_idx]);
        fft.FFT_all_dim(Gxs[photon_prop_idx], photon_prop, FFT::backward);
    }

    // Get qs
    std::vector<Field*> qs;
    if (envHasType(std::vector<Field>, par().qs))
    {
        std::vector<Field>& envQs = envGet(std::vector<Field>, par().qs);
        for (int i=0; i < envQs.size(); ++i)
            qs.push_back(&(envQs[i]));
    }
    else
        qs = envGet(std::vector<Field*>, par().qs);
    
    // Get noises
    std::vector<Field*> noises;
    if (envHasType(std::vector<Field>, par().sources))
    {
        std::vector<Field>& envNoises = envGet(std::vector<Field>, par().sources);
        for (int i=0; i < envNoises.size(); ++i)
            noises.push_back(&(envNoises[i]));
    }
    else
        noises = envGet(std::vector<Field*>, par().sources);

    // Get other parameters
    std::vector<size_t> square_radii;
    for (size_t i=0; i <= par().rSq; ++i)
        square_radii.push_back(i);
    int Nsrc   = noises.size();
    int Nd     = env().getNd();

    // ************************* //
    // CONTRACTION ROUTINE START //
    // ************************* //

    // Measure contraction on each site within the radial limit
    LOG(Message) << "Generating and contracting propagators..." << std::endl;
    
    std::vector<std::string> summary_messages;
    std::vector<Result> burgers(numPhotonProps * square_radii.size());
    std::vector<std::vector<RealD>> burger_buffers(numPhotonProps, std::vector<RealD>(Nsrc, 0.0));
    std::vector<std::vector<RealD>> burger_full_buffers(numPhotonProps, std::vector<RealD>(Nsrc, 0.0));
    std::vector<std::vector<RealD>> burger_bias_buffers(numPhotonProps, std::vector<RealD>(Nsrc, 0.0));
    startTimer("Total Contraction Time");
    for (auto square_radius : square_radii)
    {   
        // Use the symmetry for this rSq (or the final defined symmetry if fewer symmetries are given)
        QEDBurgerShortPar::SymmetryMode symmetry_name = par().shellSymmetries[std::min(square_radius, par().shellSymmetries.size()-1)];
        auto site_generator = createSiteGenerator(symmetry_name);

        const auto shifts = site_generator->culledDiscreteRadialShell(square_radius);
        LOG(Message) << "Generating propagators for " << shifts.size() << " shifts on shell r^2=" << square_radius << std::endl;
        
        // Iterate over sites
        for (int site_i=0; site_i < shifts.size(); ++site_i)
        {
            const auto& r          = shifts[site_i];
            double symmetry_factor = site_generator->symmetryFactor(r);

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
            std::vector<RealD> samenoise_contribution(numPhotonProps, 0.0);

            // Extract gauge field at offset.
            Coordinate gauge_r(Nd);
            Coordinate latt_size = env().getDim();
            gauge_r[0]=(r[0]+latt_size[0])%latt_size[0];
            gauge_r[1]=(r[1]+latt_size[1])%latt_size[1];
            gauge_r[2]=(r[2]+latt_size[2])%latt_size[2];
            gauge_r[3]=(r[3]+latt_size[3])%latt_size[3];
            
            std::vector<typename PhotonProp::scalar_object> pSites(numPhotonProps);
            for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
                peekSite(pSites[photon_prop_idx], Gxs[photon_prop_idx], gauge_r);

            // Vars for logging
            std::vector<RealD> prev_values;
            std::vector<RealD> prev_fulls;
            std::vector<RealD> prev_biases;
            for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
            {
                prev_values.push_back(burger_buffers     [photon_prop_idx].back());
                prev_fulls .push_back(burger_full_buffers[photon_prop_idx].back());
                prev_biases.push_back(burger_bias_buffers[photon_prop_idx].back());
            }

            // Now perform the contractions
            for (int hit_i=0;hit_i<Nsrc;hit_i++)
            {
                // Calculate S(x+r, x)
                startTimer("Total Contraction Time[C-Shifts]");
                coordCshift(*qs[hit_i], r, shifted_quark);
                stopTimer("Total Contraction Time[C-Shifts]");
                startTimer("Total Contraction Time[Matrix Products]");
                tmp_q1 = pointToPointProp(shifted_quark, *noises[hit_i]);
                tmp_prop1 += tmp_q1;
                stopTimer("Total Contraction Time[Matrix Products]");

                // Calculate S(x, x+r)
                startTimer("Total Contraction Time[C-Shifts]");
                coordCshift(*noises[hit_i], r, shifted_noise);
                stopTimer("Total Contraction Time[C-Shifts]");
                startTimer("Total Contraction Time[Matrix Products]");
                tmp_q2 = pointToPointProp(*qs[hit_i], shifted_noise);
                tmp_prop2 += tmp_q2;
                stopTimer("Total Contraction Time[Matrix Products]");
                
                // Do the contraction for each photon prop.
                // Only one photon prop enters the contraction, but putting more 
                // than one into the module prevents the expensive c-shifts needing 
                // to be recomputed.
                for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
                {
                    const auto pSite = pSites[photon_prop_idx];
                    
                    // Calculate same-noise contraction
                    startTimer("Total Contraction Time[Same-Noise Contraction]");
                    if constexpr(std::is_same_v<Field, PropagatorField>)
                    {
                        fastBurger(tmp_q1, tmp_q2, pSite, tmp_cbuffer);
                    }
                    else if constexpr(std::is_same_v<Field, FermionField>)
                    {
                        tmp_cbuffer = Zero();
                        for (int mu=0; mu<4; ++mu)
                        {
                            tmp_cbuffer += localInnerProduct(shifted_noise,  closure(Gamma(Gmu[mu])*shifted_quark))
                                         * localInnerProduct(*noises[hit_i], closure(Gamma(Gmu[mu])*(*qs[hit_i])));
                        }
                        tmp_cbuffer *= -1.0 * pSite; // -1.0 accounts for i^2
                    }
                    else
                    {
                        HADRONS_ERROR(Definition, "Inputs are not FermionFields or PropagatorFields");
                    }
                    samenoise_contribution[photon_prop_idx] += toReal(symmetry_factor*sum(tmp_cbuffer));
                    stopTimer("Total Contraction Time[Same-Noise Contraction]");
                
                    // Calculate full-noise contraction
                    startTimer("Total Contraction Time[Full-Noise Contraction]");
                    fastBurger(tmp_prop1, tmp_prop2, pSite, tmp_cbuffer);
                    auto full_contribution = toReal(symmetry_factor*sum(tmp_cbuffer));

                    // Add the result for this site to the total.
                    burger_full_buffers[photon_prop_idx][hit_i] += full_contribution;
                    burger_bias_buffers[photon_prop_idx][hit_i] += samenoise_contribution[photon_prop_idx];
                    burger_buffers     [photon_prop_idx][hit_i] += full_contribution - samenoise_contribution[photon_prop_idx];
                    stopTimer("Total Contraction Time[Full-Noise Contraction]");       
                }         
            }

            for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
            {
                LOG(Message) << par().photonProps[photon_prop_idx] << ": " << "Shift: " << r << " ["
                            << burger_buffers[photon_prop_idx]     .back()-prev_values[photon_prop_idx] << "/"
                            << burger_full_buffers[photon_prop_idx].back()-prev_fulls [photon_prop_idx] << "/"
                            << burger_bias_buffers[photon_prop_idx].back()-prev_biases[photon_prop_idx] << "], "
                            << "Symmetry: " << symmetry_factor << std::endl;
            }
        }

        // There are Nsrc^2 traces that could be computed using N sources for two propagators;
        // we want them all *except* the Nsrc traces using the same noise for both propagators.
        // This leaves us with (Nsrc^2 - Nsrc) traces we have computed and need to normalise for.
        for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
        {
            Result& result      = burgers     [square_radius*numPhotonProps + photon_prop_idx];
            
            result.rSq      = square_radius;
            result.symmetry = symmetry_name;
            result.photon   = par().photonProps[photon_prop_idx];

            result.burger.resize(Nsrc);
            result.full  .resize(Nsrc);
            result.bias  .resize(Nsrc);

            for (int hit_i=0; hit_i < Nsrc; ++hit_i)
            {
                int n_hits = hit_i + 1;
                int norm = n_hits == 1? 1 : (n_hits*(n_hits-1));
                auto normed_burger = burger_buffers[photon_prop_idx][hit_i]/norm;
            
                // Handle outputs
                std::stringstream summary_message;
                summary_message << "burger[radius^2=" << std::to_string(square_radius) + ", hits=" << n_hits << "] = " << std::setprecision(15) << normed_burger;
                summary_messages.push_back(summary_message.str());
                
                result.burger[hit_i] = normed_burger;
                result.full  [hit_i] = burger_full_buffers[photon_prop_idx][hit_i];
                result.bias  [hit_i] = burger_bias_buffers[photon_prop_idx][hit_i];
            }
        }

    }
    stopTimer("Total Contraction Time");

    for (const auto& summary_message : summary_messages)
        LOG(Message) << summary_message << std::endl;

    saveResult(par().output, "burgershort", burgers);
    auto& out = envGet(HadronsSerializable, getName());
    out = burgers;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDBurgerShort_hpp_
