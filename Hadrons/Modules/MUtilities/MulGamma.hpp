#ifndef Hadrons_MUtilities_MulGamma_hpp_
#define Hadrons_MUtilities_MulGamma_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MulGamma                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class MulGammaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MulGammaPar,
                                    Gamma::Algebra, gammaleft,
                                    std::string,    q,
                                    Gamma::Algebra, gammaright);
};

template <typename FImpl>
class TMulGamma: public Module<MulGammaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TMulGamma(const std::string name);
    // destructor
    virtual ~TMulGamma(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    enum class Mode
    {
        SINGLE,
        STDVECTOR
    };

    Mode mode;
};

MODULE_REGISTER_TMP(MulGamma, TMulGamma<FIMPL>, MUtilities);

/******************************************************************************
 *                 TMulGamma implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMulGamma<FImpl>::TMulGamma(const std::string name)
: Module<MulGammaPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMulGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMulGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMulGamma<FImpl>::setup(void)
{
    if (envHasType(PropagatorField, par().q))
    {
        LOG(Message) << "SINGLE PROP FIELD" << std::endl;
        envCreateLat(PropagatorField, getName());
        this->mode = Mode::SINGLE;
    }
    else if (envHasType(std::vector<PropagatorField>, par().q))
    {
        LOG(Message) << "VECTOR PROP FIELD" << std::endl;
        std::vector<PropagatorField>& q = envGet(std::vector<PropagatorField>, par().q);
        envCreate(std::vector<PropagatorField>, getName(), 1, q.size(), envGetGrid(PropagatorField));
        this->mode = Mode::STDVECTOR;
    }
    else
    {
        HADRONS_ERROR(Definition, "q is neither a 'PropagatorField' nor a 'std::vector<PropagatorField>'");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMulGamma<FImpl>::execute(void)
{
    Gamma gammaleft (par().gammaleft);
    Gamma gammaright(par().gammaright);
    
    if (this->mode == Mode::SINGLE)
    {
        PropagatorField& q   = envGet(PropagatorField, par().q);
        PropagatorField& out  = envGet(PropagatorField, getName());

        // Optimise out identity multiplications
        if (par().gammaleft == Gamma::Algebra::Identity)
            out = q*gammaright;
        else if (par().gammaright == Gamma::Algebra::Identity)
            out = gammaleft*q;
        else
            out = gammaleft*q*gammaright;
    }
    else if (this->mode == Mode::STDVECTOR)
    {
        std::vector<PropagatorField>& q   = envGet(std::vector<PropagatorField>, par().q);
        std::vector<PropagatorField>& out = envGet(std::vector<PropagatorField>, getName());

        // Optimise out identity multiplications
        if (par().gammaleft == Gamma::Algebra::Identity)
            for (int i=0; i < q.size(); ++i) out[i] = q[i]*gammaright;
        else if (par().gammaright == Gamma::Algebra::Identity)
            for (int i=0; i < q.size(); ++i) out[i] = gammaleft*q[i];
        else
            for (int i=0; i < q.size(); ++i) out[i] = gammaleft*q[i]*gammaright;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_MulGamma_hpp_
