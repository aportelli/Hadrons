#ifndef Hadrons_MAction_FermionActionModule_hpp_
#define Hadrons_MAction_FermionActionModule_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MAction)

template <typename Par>
class FermionActionModule : public Module<Par>
{
public:
    // constructor
    FermionActionModule(const std::string name) : Module<Par>(name){};

protected:
    template <typename ImplPar>
    void parseBoundary(ImplPar &implParams);
};

template <typename Par>
template <typename ImplPar>
void FermionActionModule<Par>::parseBoundary(ImplPar &implParams)
{
    if (!this->par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(this->par().boundary);
    }
    if (!this->par().twist.empty())
    {
        implParams.twist_n_2pi_L = strToVec<Real>(this->par().twist);
        if (implParams.twist_n_2pi_L.size() <= 1)
        {
            if (envHasType(std::vector<double>, this->par().twist))
            {
                auto &t = envGet(std::vector<double>, this->par().twist);
                implParams.twist_n_2pi_L = t;
            }
            else if (envHasType(HadronsSerializable, this->par().twist))
            {
                auto &t = envGet(HadronsSerializable, this->par().twist);
                implParams.twist_n_2pi_L = t.template get<std::vector<double>>();
            }
            else
            {
                HADRONS_ERROR(Definition, "cannot interpret twist input '" + this->par().twist + "'");
            }
        }
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases
                 << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L
                 << std::endl;
    if (implParams.boundary_phases.size() != this->env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != this->env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_FermionActionModule_hpp_