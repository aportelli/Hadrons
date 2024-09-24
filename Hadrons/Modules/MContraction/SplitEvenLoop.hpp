#ifndef Hadrons_MContraction_SplitEvenLoop_hpp_
#define Hadrons_MContraction_SplitEvenLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SplitEvenProp                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class SplitEvenLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SplitEvenLoopPar,
                                    std::string,    q1,
                                    std::string,    q2,
                                    double,         m1,
                                    double,         m2);
};

template <typename FImpl>
class TSplitEvenLoop: public Module<SplitEvenLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    // constructor
    TSplitEvenLoop(const std::string name);
    // destructor
    virtual ~TSplitEvenLoop(void) {};
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
        STDVECTOR,
        STDVECTOROFPTR
    };

    Mode mode;
};

MODULE_REGISTER_TMP(SplitEvenLoop,  TSplitEvenLoop<FIMPL>,  MContraction);
MODULE_REGISTER_TMP(ZSplitEvenLoop, TSplitEvenLoop<ZFIMPL>, MContraction);

/******************************************************************************
 *                 TSplitEvenLoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSplitEvenLoop<FImpl>::TSplitEvenLoop(const std::string name)
: Module<SplitEvenLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSplitEvenLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSplitEvenLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSplitEvenLoop<FImpl>::setup(void)
{
    if (envHasType(PropagatorField, par().q1))
    {
        if (!envHasType(PropagatorField, par().q2))
            {HADRONS_ERROR(Definition, "q1 is a 'PropagatorField' and q2 is not");}
        envCreateLat(PropagatorField, getName());
        this->mode = Mode::SINGLE;
    }
    else if (envHasType(std::vector<PropagatorField>, par().q1))
    {
        if (!envHasType(std::vector<PropagatorField>, par().q2))
            {HADRONS_ERROR(Definition, "q1 is a 'std::vector<PropagatorField>' and q2 is not");}
        std::vector<PropagatorField>& q1 = envGet(std::vector<PropagatorField>, par().q1);
        envCreate(std::vector<PropagatorField>, getName(), 1, q1.size(), envGetGrid(PropagatorField));
        this->mode = Mode::STDVECTOR;
    }
    else if (envHasType(std::vector<PropagatorField*>, par().q1))
    {
        if (!envHasType(std::vector<PropagatorField*>, par().q2))
            {HADRONS_ERROR(Definition, "q1 is a 'std::vector<PropagatorField*>' and q2 is not");}
        std::vector<PropagatorField*>& q1 = envGet(std::vector<PropagatorField*>, par().q1);
        envCreate(std::vector<PropagatorField>, getName(), 1, q1.size(), envGetGrid(PropagatorField));
        this->mode = Mode::STDVECTOROFPTR;
    }
    else
    {
        HADRONS_ERROR(Definition, "neither q1 and q2 are 'PropagatorField's, 'std::vector<PropagatorField>' or 'std::vector<PropagatorField*>'s");
    }
}

// execution ///////////////////////////////////////////////////////////////////
#define SplitEvenLoop(q1, q2, g5, massfactor) massfactor*q1*(g5*adj(q2)*g5)
template <typename FImpl>
void TSplitEvenLoop<FImpl>::execute(void)
{
    Gamma g5(Gamma::Algebra::Gamma5);
    double massfactor = par().m2 - par().m1;
    if (this->mode == Mode::SINGLE)
    {
        auto &res  = envGet(PropagatorField, getName());
        auto &q1   = envGet(PropagatorField, par().q1);
        auto &q2   = envGet(PropagatorField, par().q2);
        res = SplitEvenLoop(q1, q2, g5, massfactor);
    }
    else if (this->mode == Mode::STDVECTOR)
    {
        std::vector<PropagatorField>& q1  = envGet(std::vector<PropagatorField>, par().q1);
        std::vector<PropagatorField>& q2  = envGet(std::vector<PropagatorField>, par().q2);
        std::vector<PropagatorField>& out = envGet(std::vector<PropagatorField>, getName());
        for (int i=0; i < q1.size(); ++i)
            out[i] = SplitEvenLoop(q1[i], q2[i], g5, massfactor);
    }
    else if (this->mode == Mode::STDVECTOROFPTR)
    {
        std::vector<PropagatorField*>& q1  = envGet(std::vector<PropagatorField*>, par().q1);
        std::vector<PropagatorField*>& q2  = envGet(std::vector<PropagatorField*>, par().q2);
        std::vector<PropagatorField>&  out = envGet(std::vector<PropagatorField>, getName());
        for (int i=0; i < q1.size(); ++i)
            out[i] = SplitEvenLoop(*q1[i], *q2[i], g5, massfactor);
    }
}
#undef SplitEvenLoop

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_SplitEvenLoop_hpp_
