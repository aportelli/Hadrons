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
                                    double,      m2);
};

template <typename Field, typename GImpl>
class Laplacian: public LinearOperatorBase<Field>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    Laplacian(const GaugeField &U, const double m2)
    : grid_(U.Grid()), tmp_(U.Grid()), m2_(m2)
    {
        for (unsigned int mu = 0; mu < grid_->Nd(); ++mu)
        {
            U_.push_back(PeekIndex<LorentzIndex>(U, mu)); 
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
            out  -= U_[mu]*Cshift(in, mu, 1);
            tmp_  = adj(U_[mu])*in;
            out  -= Cshift(tmp_, mu, -1);
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
    Field                       tmp_;
    GridBase                    *grid_;
    double                      m2_;
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

MODULE_REGISTER_TMP(FermionLaplacian, 
                    ARG(TLaplacian<FIMPL::FermionField, GIMPL>), MAction);

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
    auto &U = envGet(GaugeField, par().gauge);
    envCreateDerived(LinearOperatorBase<Field>, ARG(Laplacian<Field, GImpl>), 
                     getName(), 1, U, par().m2);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TLaplacian<Field, GImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_Laplacian_hpp_
