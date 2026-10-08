#ifndef Hadrons_MContraction_QEDBurgerShortPoint_hpp_
#define Hadrons_MContraction_QEDBurgerShortPoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDBurgerShortPoint                                 *
 ******************************************************************************/
/*
 Point-to-point estimator for the universal disconnected "Burger" loop subdiagram
              q
             ___ 
            /   \
           |~~~~~| photon
            \___/
              q

 Tr[q(x,x+r) * Gamma_{mu} * q(x+r,x) * Gamma_{nu}] * G^{mu,nu}(r)

*/

BEGIN_MODULE_NAMESPACE(MContraction)

class QEDBurgerShortPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDBurgerShortPointPar,
                                    std::string,              eta,
                                    std::string,              q,
                                    std::vector<std::string>, photonProps,
                                    std::string,              point,
                                    std::string,              summedPropXtoY,
                                    std::string,              summedPropYtoX);
};


template <typename FImpl, typename Field, typename VType>
class TQEDBurgerShortPoint: public Module<QEDBurgerShortPointPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Complex,     burgerSingle,
                                        Complex,     burgerSummed,
                                        std::string, point,
                                        std::string, photon);
    };

    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::ScalarField PhotonProp;
    size_t numPhotonProps;
public:
    // constructor
    TQEDBurgerShortPoint(const std::string name);
    // destructor
    virtual ~TQEDBurgerShortPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

    // TEMPORARY?
    inline void coordCshift(const Field& field, const Coordinate& coord, Field& out) const;
    inline void fastBurger(const PropagatorField& left, const PropagatorField& right, const typename PhotonProp::scalar_object& pSite, LatticeComplexD& out) const;
    inline void pointToPointProp(PropagatorField& out, const PropagatorField& left, const PropagatorField& right);
    inline void pointToPointProp(PropagatorField& out, const FermionField&    left, const FermionField&    right);
};

MODULE_REGISTER_TMP(QEDBurgerShortPoint, ARG(TQEDBurgerShortPoint<FIMPL, FIMPL::PropagatorField, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDBurgerShortPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
TQEDBurgerShortPoint<FImpl, Field, VType>::TQEDBurgerShortPoint(const std::string name)
: Module<QEDBurgerShortPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
std::vector<std::string> TQEDBurgerShortPoint<FImpl, Field, VType>::getInput(void)
{
    std::vector<std::string> in = {par().eta, par().q};

    for (const std::string& photon_prop : par().photonProps)
        in.push_back(photon_prop);
    
    if (!par().summedPropXtoY.empty() && !par().summedPropYtoX.empty())
    {
        in.push_back(par().summedPropXtoY);
        in.push_back(par().summedPropYtoX);
    }

    return in;
}

template <typename FImpl, typename Field, typename VType>
std::vector<std::string> TQEDBurgerShortPoint<FImpl, Field, VType>::getOutput(void)
{
    std::vector<std::string> out = {
        getName(),
        getName() + "_summedPropXtoY", 
        getName() + "_summedPropYtoX"
    };
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::setup(void)
{
    numPhotonProps = par().photonProps.size();
    if (numPhotonProps==0)
    {
        HADRONS_ERROR(Definition, "No photon props provided to module");
    }

    // Contraction temporaries
    envTmp(LatticeComplexD, "tmp_cbuffer",   1, env().getGrid());
    envTmp(Field,           "shifted_quark", 1, envGetGrid(Field));
    envTmp(Field,           "shifted_eta",   1, envGetGrid(Field));
    // Output
    envCreateLat(PropagatorField, getName() + "_summedPropXtoY");
    envCreateLat(PropagatorField, getName() + "_summedPropYtoX");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////

template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::execute(void)
{
    // Define Gammas
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    // Get temps
    envGetTmp(LatticeComplexD, tmp_cbuffer);
    typename PhotonProp::scalar_object pSite;

    // Output variables
    HadronsSerializable& out  = envGet(HadronsSerializable, getName());
    PropagatorField& propXtoY = envGet(PropagatorField, getName() + "_summedPropXtoY");
    PropagatorField& propYtoX = envGet(PropagatorField, getName() + "_summedPropYtoX");

    // Create the stochastically-estimated point-to-point propagators
    //   S(x+r, eta) eta(x  )
    //   S(x,   eta) eta(x+r)
    const Field& q   = envGet(Field, par().q  );
    const Field& eta = envGet(Field, par().eta);
    Coordinate coord(strToVec<int>(par().point));

    Coordinate photonCoord(coord);
    Coordinate latt_size = env().getDim();
    photonCoord[0]=(coord[0]+latt_size[0])%latt_size[0];
    photonCoord[1]=(coord[1]+latt_size[1])%latt_size[1];
    photonCoord[2]=(coord[2]+latt_size[2])%latt_size[2];
    photonCoord[3]=(coord[3]+latt_size[3])%latt_size[3];
    {
        envGetTmp(Field, shifted_quark);
        coordCshift(q, coord, shifted_quark);
        pointToPointProp(propXtoY, shifted_quark, eta);
    }
    {
        envGetTmp(Field, shifted_eta);
        coordCshift(eta, coord, shifted_eta);
        pointToPointProp(propYtoX, q, shifted_eta);
    }

    // Do the contractions for each photon prop.
    // Only one photon prop enters the contraction, but putting more 
    // than one into the module prevents the expensive c-shifts needing 
    // to be recomputed / a second module to multiply in the photon field site.

    std::vector<Result> results(numPhotonProps);

    // Calculate burger contraction using the same noise hit on both props.
    // These are the "bias" terms that need to be subtracted from the final estimator.
    for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
    {
        const PhotonProp& Gx = envGet(PhotonProp, par().photonProps[photon_prop_idx]);
        peekSite(pSite, Gx, photonCoord);
        
        startTimer("Total Contraction Time [Single Contraction]");
        fastBurger(propXtoY, propYtoX, pSite, tmp_cbuffer);
        // tmp_cbuffer = -traceProduct(propXtoY, propYtoX)*pSite;
        results[photon_prop_idx].burgerSingle = sum(tmp_cbuffer);
        stopTimer("Total Contraction Time [Single Contraction]");
    }

    // Add externally-provided propagators to the point-to-point propagators.
    // These external propagators are intended to be obtained from the "_summedProp..."
    // outputs of previous invocations of this module.
    // Calling this module for successive noise hits therefore accumulates a sum of
    // estimated propagators in the "_summedProp..." variables.
    if (!par().summedPropXtoY.empty() && !par().summedPropYtoX.empty())
    {
        propXtoY += envGet(PropagatorField, par().summedPropXtoY);
        propYtoX += envGet(PropagatorField, par().summedPropYtoX);
    }

    // Calculate burger contraction using the sum of noise-estimated props.
    // These contain all possible pairwise combinations of whatever noises have been used
    // in the prop sum, and need to have all "bias" terms subtracted in the analysis.
    // After bias-subtraction, they will also need to be normalised by (n_hits-1)*n_hits.
    for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
    {
        const PhotonProp& Gx = envGet(PhotonProp, par().photonProps[photon_prop_idx]);
        peekSite(pSite, Gx, photonCoord);
        
        startTimer("Total Contraction Time [Summed Contraction]");
        fastBurger(propXtoY, propYtoX, pSite, tmp_cbuffer);
        results[photon_prop_idx].burgerSummed = sum(tmp_cbuffer);
        stopTimer("Total Contraction Time [Summed Contraction]");       
    }

    // Fill in the remaining 'results' data members.
    for (int photon_prop_idx=0; photon_prop_idx < numPhotonProps; ++photon_prop_idx)
    {
        results[photon_prop_idx].point  = par().point;
        results[photon_prop_idx].photon = par().photonProps[photon_prop_idx];

        LOG(Message) << par().photonProps[photon_prop_idx] << ", " << par().point << ", Burger Single: " << results[photon_prop_idx].burgerSingle << std::endl;
        LOG(Message) << par().photonProps[photon_prop_idx] << ", " << par().point << ", Burger Summed: " << results[photon_prop_idx].burgerSummed << std::endl;
    }

    out = results;
}


template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::fastBurger(const PropagatorField& left, const PropagatorField& right, const typename PhotonProp::scalar_object& pSite, LatticeComplexD& out) const
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
            out_s += tmp;
        }
        // Now also remember the -1 to account for i^2 from the two Aslash insertions.
        wvout[s] = -1.*out_s*pSite;
    })
}


template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::pointToPointProp(PropagatorField& out, const PropagatorField& left, const PropagatorField& right)
{
    out = left * adj(right);
}

template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::pointToPointProp(PropagatorField& out, const FermionField&    left, const FermionField&    right)
{
    out = outerProduct(left, right);
}


template <typename FImpl, typename Field, typename VType>
void TQEDBurgerShortPoint<FImpl, Field, VType>::coordCshift(const Field& field, const Coordinate& coord, Field& out) const
{
    out = Cshift(field,0,coord[0]);
    for (int mu=1;mu<Nd;mu++) out = Cshift(out,mu,coord[mu]);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDBurgerShortPoint_hpp_
