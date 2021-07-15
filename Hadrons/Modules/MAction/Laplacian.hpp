#ifndef Hadrons_MAction_Laplacian_hpp_
#define Hadrons_MAction_Laplacian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Covariant Laplacian operator                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class LaplacianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianPar,
                                    std::string, gauge,
                                    std::string, boundary,
                                    std::string, twist,
                                    double,      m2);
};

template <typename Field, typename GImpl>
class Laplacian: public LinearOperatorBase<Field>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef typename GImpl::scalar_type Scalar;
public:
    Laplacian(const GaugeField &U, const double m2, 
              const std::vector<Complex> &boundary, const std::vector<Real> &twist)
    : grid_(U.Grid()), nd_(U.Grid()->Nd())
    , Umu_(U.Grid()), tmp_(U.Grid()), m2_(m2)
    , boundary_(boundary), twist_(twist)
    {
        Lattice<iScalar<vInteger>> coor(grid_);

        assert(boundary.size() == nd_);
        assert(twist.size() == nd_);
        U_.resize(2*nd_, grid_);
        for (unsigned int mu = 0; mu < grid_->Nd(); ++mu)
        {
            unsigned int lmu = grid_->GlobalDimensions()[mu];
            Scalar       ph(real(boundary_[mu]), imag(boundary_[mu]));
            Scalar       tw(cos(twist_[mu]*2.*M_PI/lmu), sin(twist_[mu]*2.*M_PI/lmu));

            LatticeCoordinate(coor, mu);
            Umu_         = PeekIndex<LorentzIndex>(U, mu);
            Umu_         = tw*Umu_;
            tmp_         = where(coor == lmu - 1, ph*Umu_, Umu_);
            U_[mu]       = tmp_;
            Umu_         = adj(Cshift(Umu_, mu, -1));
            Umu_         = where(coor == 0, conjugate(ph)*Umu_, Umu_);
            U_[mu + nd_] = Umu_; 
        }
    }

    virtual void OpDiag(const Field &in, Field &out)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void OpDir(const Field &in, Field &out, int dir, int disp)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void OpDirAll(const Field &in, std::vector<Field> &out)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void Op(const Field &in, Field &out)
    {
        unsigned int nd = grid_->Nd();

        out = (m2_ + 2.*nd)*in;
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            out -= U_[mu]*Cshift(in, mu, 1);
            out -= U_[mu + nd_]*Cshift(in, mu, -1);
        }
    }

    virtual void AdjOp(const Field &in, Field &out)
    {
        Op(in, out);
    }

    virtual void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void HermOp(const Field &in, Field &out)
    {
        Op(in, out);
    }
private:
    std::vector<GaugeLinkField> U_;
    GaugeLinkField              Umu_, tmp_;
    GridBase                    *grid_;
    const unsigned int          nd_;
    double                      m2_;
    std::vector<Complex>        boundary_;
    std::vector<Real>           twist_;
};

template <typename Field, typename GImpl>
class TLaplacian: public Module<LaplacianPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TLaplacian(const std::string name);
    // destructor
    virtual ~TLaplacian(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CovariantLaplacian, 
                    ARG(TLaplacian<ColourVectorField<FIMPL>, GIMPL>), MAction);

/******************************************************************************
 *                      TLaplacian implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
TLaplacian<Field, GImpl>::TLaplacian(const std::string name)
: Module<LaplacianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename GImpl>
std::vector<std::string> TLaplacian<Field, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename Field, typename GImpl>
std::vector<std::string> TLaplacian<Field, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TLaplacian<Field, GImpl>::setup(void)
{
    LOG(Message) << "Setting up covariant Laplacian operator with gauge field '"
                 << par().gauge << "' and m^2= " << par().m2 << std::endl;
    auto                 &U = envGet(GaugeField, par().gauge);
    unsigned int         nd = U.Grid()->Nd();
    std::vector<Complex> boundary(nd, 1.);
    std::vector<Real>    twist(nd, 0.);
    
    if (!par().boundary.empty())
    {
        boundary = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        twist = strToVec<Real>(par().twist);
    }
    envCreateDerived(LinearOperatorBase<Field>, ARG(Laplacian<Field, GImpl>), 
                     getName(), 1, U, par().m2, boundary, twist);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TLaplacian<Field, GImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_Laplacian_hpp_
