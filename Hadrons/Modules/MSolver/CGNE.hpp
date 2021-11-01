#ifndef Hadrons_MSolver_CGNE_hpp_
#define Hadrons_MSolver_CGNE_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         CGNE                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class CGNEPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CGNEPar,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , eigenPack);
};

template <typename FImpl, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TCGNE: public Module<CGNEPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCGNE(const std::string name);
    // destructor
    virtual ~TCGNE(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CGNE, TCGNE<FIMPL>, MSolver);

/******************************************************************************
 *                 TCGNE implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TCGNE<FImpl, nBasis>::TCGNE(const std::string name)
: Module<CGNEPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TCGNE<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TCGNE<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TCGNE<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TCGNE<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up normal equation CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;
    
    auto Ls        = env().getObjectLs(par().action);
    auto &mat      = envGet(FMat, par().action);
    auto guesserPt = makeGuesser<FImpl, nBasis>(par().eigenPack);

    auto makeSolver = [&mat, guesserPt, this](bool subGuess) 
    {
        return [&mat, guesserPt, subGuess, this](FermionField &sol,
                                                 const FermionField &source) 
        {
            GridBase                                *g = sol.Grid();
            FermionField                            guess(g), tmp(g);
            MdagMLinearOperator<FMat, FermionField> hermOp(mat);
            ConjugateGradient<FermionField>         cg(par().residual,
                                                       par().maxIteration);

            guess = sol;
            mat.Mdag(source, tmp);
            (*guesserPt)(tmp, sol);
            cg(hermOp, tmp, sol);
            if (subGuess)
            {
                sol -= guess;
            }
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TCGNE<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_CGNE_hpp_
