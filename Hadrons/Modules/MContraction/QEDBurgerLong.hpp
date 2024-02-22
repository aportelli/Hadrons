#ifndef Hadrons_MContraction_QEDBurgerLong_hpp_
#define Hadrons_MContraction_QEDBurgerLong_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDBurgerLong                                 *
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

class QEDBurgerLongPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDBurgerLongPar,
                                    int,         radius,
                                    std::string, q,
                                    std::string, photon,
                                    std::string, origin);
};

template <typename FImpl, typename VType>
class TQEDBurgerLong: public Module<QEDBurgerLongPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    
    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::GaugeField EmField;

    // constructor
    TQEDBurgerLong(const std::string name);
    // destructor
    virtual ~TQEDBurgerLong(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(QEDBurgerLong, ARG(TQEDBurgerLong<FIMPL, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDBurgerLong implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
TQEDBurgerLong<FImpl, VType>::TQEDBurgerLong(const std::string name)
: Module<QEDBurgerLongPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename VType>
std::vector<std::string> TQEDBurgerLong<FImpl, VType>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().photon};
    
    return in;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDBurgerLong<FImpl, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDBurgerLong<FImpl, VType>::setup(void)
{
    envTmp(FFT,            "fft",            1, env().getGrid());
    envTmp(EmField,        "Gx",             1, env().getGrid());
    envTmp(EmField,        "em_buffer",      1, env().getGrid());
    envTmp(LatticeComplex, "burger_lattice", 1, env().getGrid());
    envCreate(HadronsSerializable, getName(), 1, 0);
    if (par().radius >= 0)
    {
        envTmp(LatticeInteger, "tmp_intbuffer",    1, env().getGrid());
        envTmp(LatticeInteger, "lsize_buffer",     1, env().getGrid());
        envTmp(LatticeInteger, "lhalfsize_buffer", 1, env().getGrid());
        envTmp(LatticeInteger, "coord_buffer",     1, env().getGrid());
        envTmp(LatticeInteger, "radial_dist_sq",   1, env().getGrid());
        envTmp(LatticeComplex, "tmp_cbuffer",      1, env().getGrid());
    }
}

// execution ///////////////////////////////////////////////////////////////////
#define Burger(q, gamma1, gamma2, g5) trace(gamma1*q*gamma2*g5*adj(q)*g5)
template <typename FImpl, typename VType>
void TQEDBurgerLong<FImpl, VType>::execute(void)
{
    Gamma Gmu[] = 
    {
        Gamma(Gamma::Algebra::GammaX),
        Gamma(Gamma::Algebra::GammaY),
        Gamma(Gamma::Algebra::GammaZ),
        Gamma(Gamma::Algebra::GammaT),
    };

    // Get position-space photon
    std::cout << "Using photon '" << par().photon << "'" << std::endl;
    const EmField& Ax = envGet(EmField, par().photon);

    // Hacky way to create feynman photon gauge field from photon field
    envGetTmp(EmField, Gx);
    envGetTmp(EmField, em_buffer);
    envGetTmp(FFT,     fft);
    fft.FFT_all_dim(em_buffer, Ax, FFT::forward);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,0)*peekLorentz(em_buffer,0), 0);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,1)*peekLorentz(em_buffer,1), 1);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,2)*peekLorentz(em_buffer,2), 2);
    pokeLorentz(em_buffer, peekLorentz(em_buffer,3)*peekLorentz(em_buffer,3), 3);
    fft.FFT_all_dim(Gx, em_buffer, FFT::backward);

    auto& out = envGet(HadronsSerializable, getName());
    if (par().radius >= 0)
    {
        // Compute the distance from the origin for each lattice site.
        // To do this, replace the coordinates of sites larger than half the
        // size of the lattice with their negative equivalents.
        // As long as the radial cut is closer to the origin than the opposing
        // corner of the hypercube, this will ensure we correctly measure the 
        // radial distace for the purposes of splitting the field.
        envGetTmp(LatticeInteger, coord_buffer);
        envGetTmp(LatticeInteger, tmp_intbuffer);
        envGetTmp(LatticeInteger, lsize_buffer);
        envGetTmp(LatticeInteger, lhalfsize_buffer);
        envGetTmp(LatticeInteger, radial_dist_sq);
        radial_dist_sq = Zero();
        Coordinate latt_size   = env().getDim();
        for (int mu=0;mu<Nd;mu++)
        {
            LatticeCoordinate(coord_buffer, mu);
            lhalfsize_buffer = latt_size[mu]/2;
            lsize_buffer     = latt_size[mu];
            tmp_intbuffer    = where(coord_buffer<=lhalfsize_buffer,coord_buffer,lsize_buffer-coord_buffer);
            tmp_intbuffer    *= tmp_intbuffer;
            radial_dist_sq   += tmp_intbuffer;
        }

        // Prepare to shift the propagator back to the source origin
        // in order to measure the contraction distance we will split
        // the diagram at.
        Coordinate origin = strToVec<int>(par().origin);

        // Make a copy since we'll be C-shifting it.
        // This *only* gets close to the short-distance module number
        // when the propagator field is C-shifted, and in that situation
        // the result is the same up to double precision.
        PropagatorField q = envGet(PropagatorField, par().q);
        for (int mu=0;mu<Nd;mu++)
            q = Cshift(q, mu, origin[mu]);

        // // Alternative to shifting the propagator field...
        // for (int mu=0;mu<Nd;mu++)
        // {
        //     radial_dist_sq  = Cshift(radial_dist_sq, mu, -origin[mu]);
        //     Gx              = Cshift(Gx,             mu, -origin[mu]);
        // }

        // Calculate the Burger in Feynman gauge
        // The Feynman gauge photon propagator is delta^{mu,nu}/k^2, so we only need
        // to compute terms where mu == nu.
        Gamma Gamma5(Gamma::Algebra::Gamma5);
        envGetTmp(LatticeComplex, burger_lattice);
        burger_lattice = Zero();
        for (int mu=0; mu<Nd; ++mu)
            burger_lattice += peekLorentz(Gx,mu)*Burger(q, Gmu[mu], Gmu[mu], Gamma5);

        // Get the split-field results.
        envGetTmp(LatticeComplex, tmp_cbuffer);
        tmp_cbuffer     = ComplexD(0.0, 0.0);
        RealD burger    = toReal(sum(where(radial_dist_sq  > static_cast<Integer>(par().radius*par().radius),burger_lattice,tmp_cbuffer)));
         // TODO: Remove this, exists for comparison against SD module
        RealD burger_sd = toReal(sum(where(radial_dist_sq <= static_cast<Integer>(par().radius*par().radius),burger_lattice,tmp_cbuffer)));
        LOG(Message) << "burger ld: " << std::setprecision(15) << burger << std::endl;
        LOG(Message) << "burger sd: " << std::setprecision(15) << burger_sd << std::endl; // For comparison against SD module
        LOG(Message) << "burger:    " << std::setprecision(15) << burger+burger_sd << std::endl;
        
        out = burger;
    }
    else
    {
        PropagatorField q = envGet(PropagatorField, par().q);
        Coordinate origin = strToVec<int>(par().origin);
        // This Cshift shouldn't be necessary? Sum should be invariant...
        // The result is *only* consistent if the PropagatorField or Photon field
        // are C-shifted. In this situation the result is the same up to double precision.
        // for (int mu=0;mu<Nd;mu++)
        //     q = Cshift(q, mu, origin[mu]);
        for (int mu=0;mu<Nd;mu++)
            Gx = Cshift(Gx, mu, -origin[mu]);

        Gamma Gamma5(Gamma::Algebra::Gamma5);
        envGetTmp(LatticeComplex, burger_lattice);
        burger_lattice = Zero();
        for (int mu=0; mu<Nd; ++mu)
            burger_lattice += peekLorentz(Gx,mu)*Burger(q, Gmu[mu], Gmu[mu], Gamma5);
        RealD burger = toReal(sum(burger_lattice));
        LOG(Message) << "burger: " << std::setprecision(15) << burger << std::endl;
        
        out = burger;
    }
}
#undef Burger

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDBurgerLong_hpp_
