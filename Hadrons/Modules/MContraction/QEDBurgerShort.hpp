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
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDBurgerShortPar,
                                    std::string, sources,
                                    std::string, qs,
                                    std::string, photon,
                                    int,         radius,
                                    std::string, output);
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

    inline void fastBurger(const PropagatorField& left, const PropagatorField& right, const typename EmField::scalar_object& pSite, LatticeComplexD& out);
    inline void coordCshift(const PropagatorField& field, const Coordinate& coord, PropagatorField& out);
private:
    int Ls_;
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
    std::vector<std::string> in = {par().sources, par().photon, par().qs};
    
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
    const auto& sources = envGet(std::vector<PropagatorField>, par().sources);
    // Photon temporaries
    envTmp(FFT,             "fft",           1, env().getGrid());
    envTmp(EmField,         "Gx",            1, env().getGrid());
    envTmp(EmField,         "em_buffer",     1, env().getGrid());
    // Contraction temporaries
    envTmp(LatticeComplexD, "tmp_cbuffer",   1, env().getGrid());
    envTmp(PropagatorField, "tmp_prop1",     1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_prop2",     1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_q1",        1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "tmp_q2",        1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "shifted_quark", 1, envGetGrid(PropagatorField));
    envTmp(PropagatorField, "shifted_noise", 1, envGetGrid(PropagatorField));
    // Output
    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::fastBurger(const PropagatorField& left, const PropagatorField& right, const typename EmField::scalar_object& pSite, LatticeComplexD& out)
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

    constexpr int Nc = 3;
    constexpr int Ns = 4;
    const auto& grid = envGetGrid4(PropagatorField);
    // Tr[G_mu * left * G_mu * right] * Guv[site]
    accelerator_for(s, grid->oSites(), VType::Nsimd(), 
    {
        const auto& s1 = vleft[s];  //coalescedRead(vleft(s));
        const auto& s2 = vright[s]; //coalescedRead(vright(s));
        LatticeComplexD::vector_object out_s = Zero();
        for (int mu=0; mu<4; ++mu)
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
        //coalescedWrite(wvout[s], out_s);
        wvout[s] = out_s;
    })
}

template <typename FImpl, typename VType>
void TQEDBurgerShort<FImpl, VType>::coordCshift(const PropagatorField& field, const Coordinate& coord, PropagatorField& out)
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
    envGetTmp(EmField,                   Gx);
    envGetTmp(EmField,                   em_buffer);
    envGetTmp(LatticeComplexD,           tmp_cbuffer);
    envGetTmp(PropagatorField,           tmp_prop1);
    envGetTmp(PropagatorField,           tmp_prop2);
    envGetTmp(PropagatorField,           tmp_q1);
    envGetTmp(PropagatorField,           tmp_q2);
    envGetTmp(PropagatorField,           shifted_quark);
    envGetTmp(PropagatorField,           shifted_noise);
    envGetTmp(FFT,                       fft);

    // Get env vars
    const EmField& Ax                          = envGet(EmField, par().photon);
    const std::vector<PropagatorField>& qs     = envGet(std::vector<PropagatorField>, par().qs);
    const std::vector<PropagatorField>& noises = envGet(std::vector<PropagatorField>, par().sources);

    // Get other parameters
    int         Nsrc   = noises.size();
    int         radius = par().radius;

    int         rsqmax = radius*radius;

    // Hacky way to create feynman photon gauge field from photon field
    fft.FFT_all_dim(em_buffer,Ax,FFT::forward);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,0)*peekLorentz(em_buffer,0), 0);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,1)*peekLorentz(em_buffer,1), 1);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,2)*peekLorentz(em_buffer,2), 2);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,3)*peekLorentz(em_buffer,3), 3);
    fft.FFT_all_dim(Gx, em_buffer, FFT::backward);

    // ************************* //
    // CONTRACTION ROUTINE START //
    // ************************* //

    // Measure contraction on each site within the radial limit
    LOG(Message) << "Generating and contracting all-to-all propagators..." << std::endl;
    RealD burger = 0.0;

    startTimer("Total Contraction Time");

    // Generate list of sites
    std::vector<Coordinate> rsites;
    for (int r0=-radius; r0<=radius; ++r0)
    for (int r1=-radius; r1<=radius; ++r1)
    for (int r2=-radius; r2<=radius; ++r2)
    for (int r3=-radius; r3<=radius; ++r3)
    {
        int rsq = r0*r0+r1*r1+r2*r2+r3*r3;
        if (rsq > rsqmax) continue;

        Coordinate r(Nd);
        r[0]=r0;
        r[1]=r1;
        r[2]=r2;
        r[3]=r3;
        rsites.push_back(r);
    }

    // Iterate over sites
    for (int i=0; i < rsites.size(); ++i)
    {
        const auto& r = rsites[i];
        LOG(Message) << "[" << i+1 << "/" << rsites.size() << "] Generating all-to-all propagator for site " << r << std::endl;

        
        // ***************************************************************************************** //
        // To estimate the burger diagram, we average over traces computed from pairs of propagators
        // solved on different noises. The basic way to do this is to loop over two 'noise indices'
        // and compute the trace for each iteration where the indices are not equal.
        // However, this requires O(Nsrc^2) expensive C-shifted outer products to compute.
        // We can reduce this to O(Nsrc) C-shifted outer products if we create two individual
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
        for (int i=0;i<Nsrc;i++)
        {
            startTimer("Total Contraction Time[C-Shifts]");
            coordCshift(qs[i], r, shifted_quark);
            stopTimer("Total Contraction Time[C-Shifts]");
            startTimer("Total Contraction Time[Outer Products]");
            tmp_q1 = shifted_quark*adj(noises[i]);
            tmp_prop1 += tmp_q1;
            stopTimer("Total Contraction Time[Outer Products]");

            startTimer("Total Contraction Time[C-Shifts]");
            coordCshift(noises[i], r, shifted_noise);
            stopTimer("Total Contraction Time[C-Shifts]");
            startTimer("Total Contraction Time[Outer Products]");
            tmp_q2 = qs[i]*adj(shifted_noise);
            tmp_prop2 += tmp_q2;
            stopTimer("Total Contraction Time[Outer Products]");

            startTimer("Total Contraction Time[Same-Noise Contraction]");
            fastBurger(tmp_q1, tmp_q2, pSite, tmp_cbuffer);
            samenoise_contribution += toReal(sum(tmp_cbuffer));
            stopTimer("Total Contraction Time[Same-Noise Contraction]");
        }

        // Calculate full correlator.
        startTimer("Total Contraction Time[Full-Noise Contraction]");
        fastBurger(tmp_prop1, tmp_prop2, pSite, tmp_cbuffer);
        auto full_contribution = toReal(sum(tmp_cbuffer));

        // Add the result for this site to the total.
        burger += full_contribution - samenoise_contribution;
        stopTimer("Total Contraction Time[Full-Noise Contraction]");
    }
    stopTimer("Total Contraction Time");

    // There are Nsrc^2 traces that could be computed using N sources for two propagators;
    // we compute them all *except* the Nsrc traces using the same noise for both propagators.
    // This leaves us with (Nsrc^2 - Nsrc) traces we have computed and need to normalise for.
    int norm = (Nsrc*(Nsrc-1));
    burger /= norm*env().getVolume();
    std::cout << "burger = " << std::setprecision(15) << burger << std::endl;
    
    saveResult(par().output, "burgershort", burger);
    auto& out = envGet(HadronsSerializable, getName());
    out = burger;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDBurgerShort_hpp_
