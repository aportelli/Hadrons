#ifndef Hadrons_MFermion_SplitEvenProp_hpp_
#define Hadrons_MFermion_SplitEvenProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SplitEvenProp                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class SplitEvenPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SplitEvenPropPar,
                                    std::string,    q1,
                                    // Must be a propagator solved on g5*src
                                    std::string,    q2g5,
                                    double,         m1,
                                    double,         m2,
                                    // TODO: Remove 'rotated'
                                    bool,           rotated,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl>
class TSplitEvenProp: public Module<SplitEvenPropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    // constructor
    TSplitEvenProp(const std::string name);
    // destructor
    virtual ~TSplitEvenProp(void) {};
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

MODULE_REGISTER_TMP(SplitEvenProp,  TSplitEvenProp<FIMPL>,  MFermion);
MODULE_REGISTER_TMP(ZSplitEvenProp, TSplitEvenProp<ZFIMPL>, MFermion);

/******************************************************************************
 *                 TSplitEvenProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSplitEvenProp<FImpl>::TSplitEvenProp(const std::string name)
: Module<SplitEvenPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSplitEvenProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2g5};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSplitEvenProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSplitEvenProp<FImpl>::setup(void)
{
    if (envHasType(PropagatorField, par().q1))
    {
        if (!envHasType(PropagatorField, par().q2g5))
            {HADRONS_ERROR(Definition, "q1 is a 'PropagatorField' and q2g5 is not");}
        envCreateLat(PropagatorField, getName());
        this->mode = Mode::SINGLE;
    }
    else if (envHasType(std::vector<PropagatorField>, par().q1))
    {
        if (!envHasType(std::vector<PropagatorField>, par().q2g5))
            {HADRONS_ERROR(Definition, "q1 is a 'std::vector<PropagatorField>' and q2g5 is not");}
        std::vector<PropagatorField>& q1 = envGet(std::vector<PropagatorField>, par().q1);
        envCreate(std::vector<PropagatorField>, getName(), 1, q1.size(), envGetGrid(PropagatorField));
        this->mode = Mode::STDVECTOR;
    }
    else
    {
        HADRONS_ERROR(Definition, "q1 and q2g5 are neither 'PropagatorField's nor 'std::vector<PropagatorField>'s");
    }
}

// execution ///////////////////////////////////////////////////////////////////
#define SplitEven(q1, q2g5, g5, gamma, massfactor) massfactor*adj(g5*q2g5)*gamma*q1
#define SplitEvenRotated(q1, q2g5, g5, gamma, massfactor) massfactor*q1*adj(g5*q2g5)*gamma
template <typename FImpl>
void TSplitEvenProp<FImpl>::execute(void)
{
    Gamma g5(Gamma::Algebra::Gamma5);
    Gamma gamma(par().gamma);
    double massfactor = par().m2 - par().m1;
    if (this->mode == Mode::SINGLE)
    {
        PropagatorField& q1   = envGet(PropagatorField, par().q1);
        PropagatorField& q2g5 = envGet(PropagatorField, par().q2g5);
        PropagatorField& out  = envGet(PropagatorField, getName());

        // TODO: Remove 'rotated' - hack for regression
        if (par().rotated) out = SplitEvenRotated(q1, q2g5, g5, gamma, massfactor);
        else               out = SplitEven       (q1, q2g5, g5, gamma, massfactor);
        
    }
    else if (this->mode == Mode::STDVECTOR)
    {
        std::vector<PropagatorField>& q1   = envGet(std::vector<PropagatorField>, par().q1);
        std::vector<PropagatorField>& q2g5 = envGet(std::vector<PropagatorField>, par().q2g5);
        std::vector<PropagatorField>& out  = envGet(std::vector<PropagatorField>, getName());
        for (int i=0; i < q1.size(); ++i)
        {
            // TODO: Remove 'rotated' - hack for regression
            if (par().rotated) out[i] = SplitEvenRotated(q1[i], q2g5[i], g5, gamma, massfactor);
            else               out[i] = SplitEven       (q1[i], q2g5[i], g5, gamma, massfactor);
        }
    }
}
#undef SplitEven
#undef SplitEvenRotated

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SplitEvenProp_hpp_
