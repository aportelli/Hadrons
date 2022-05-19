
#include <Grid/Grid.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

#ifndef NBASIS
#define NBASIS 50
#endif

using namespace std;
using namespace Grid;
using namespace Hadrons;


template<typename F>
void ProjAccumRunner(std::vector<typename F::FermionField> &in, std::vector<typename F::FermionField> &out, unsigned int eb, unsigned int sb, unsigned int totSizeE)
{
    GridBase *g       = in[0].Grid();
    EigenPack<typename F::FermionField> Epack(eb, g, g);
    
    std::vector<int> seeds({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    GridSerialRNG            RNG;      RNG.SeedFixedIntegers(seeds);
    GridParallelRNG          RNG5(g);  RNG5.SeedFixedIntegers(seeds5);

    GridStopWatch w1;
    GridTime ProjAccum = GridTime::zero();

    LOG(Message) << "ProjAccumRunner start" << std::endl;
    
    for (int i = 0; i < totSizeE; i += eb)
    {
        unsigned int evBlockSize = std::min(totSizeE - i, eb);
        
        LOG(Message) << "New eigenvector picks" << std::endl;

        for (auto &e: Epack.evec)
        {
            random(RNG5,e);
        }

        for (auto &e: Epack.eval)
        {
            random(RNG,e);

            LOG(Debug) << "eigenval random pick: " << e << std::endl;
        }

        LOG(Message) << "evBlockSize: " << evBlockSize << std::endl;

        for(unsigned int j = 0; j < in.size(); j += sb)
        {
            unsigned int srcBlockSize = std::min((int)in.size() - j, sb);

            LOG(Message) << "srcBlockSize: " << srcBlockSize << std::endl;

            w1.Start();
            BatchDeflationUtils::projAccumulate<typename F::FermionField>(in, out, Epack.evec, Epack.eval, 0, evBlockSize, j, j+1);
            w1.Stop();
            ProjAccum += w1.Elapsed();
            w1.Reset();
        }
    }

    LOG(Message) << "ProjAccumRunner end" << std::endl;
    LOG(Message) << "ProjAccum total: " << ProjAccum << std::endl;
    LOG(Message) << "ProjAccum out norm: " << norm2(out[0]) << std::endl;

}

GridBase * makeGrid(const unsigned int Ls, const bool rb, const bool single=false, const bool coarse=false, std::vector<int> blockSize={1,1,1,1,1})
{
  auto &env = Environment::getInstance();

  if(coarse)
  {
    if(single)
    {
        env.createCoarseGrid<vComplexF>(blockSize, Ls);

        if (rb)
        {
            LOG(Error) << "No Rb Coarse grid in hadrons" << std::endl;
            exit(1);
        }
        else
        {
            if (Ls > 1)
            {   
            return env.getCoarseGrid<vComplexF>(blockSize, Ls);
            }
            else
            {
            return env.getCoarseGrid<vComplexF>(blockSize);
            }
        }

    }
    else
    {
        env.createCoarseGrid(blockSize, Ls);

        if (rb)
        {
            LOG(Error) << "No Rb Coarse grid in hadrons" << std::endl;
            exit(2); 
        }
        else
        {
            if (Ls > 1)
            {   
            return env.getCoarseGrid(blockSize, Ls);
            }
            else
            {
            return env.getCoarseGrid(blockSize);
            }
        }
    }      
  }
  else
  {
    if(single)
    {
        env.createGrid<vComplexF>(Ls);

        if (rb)
        {
            if (Ls > 1)
            {
            return env.getRbGrid<vComplexF>(Ls);
            }
            else
            {
            return env.getRbGrid<vComplexF>();
            }
        }
        else
        {
            if (Ls > 1)
            {   
            return env.getGrid<vComplexF>(Ls);
            }
            else
            {
            return env.getGrid<vComplexF>();
            }
        }

    }
    else
    {
        env.createGrid(Ls);

        if (rb)
        {
            if (Ls > 1)
            {
            return env.getRbGrid(Ls);
            }
            else
            {
            return env.getRbGrid();
            }
        }
        else
        {
            if (Ls > 1)
            {   
            return env.getGrid(Ls);
            }
            else
            {
            return env.getGrid();
            }
        }
    }
  } 
}

template <typename F> 
void scanner(GridBase *g, bool single,
             unsigned int minBatchSizeE, unsigned int maxBatchSizeE, 
             unsigned int minBatchSizeS, unsigned int maxBatchSizeS,
             unsigned int totSizeE, unsigned int totSizeS, unsigned int stepSize,
             unsigned int version)
{
    LOG(Debug) << "Check Grid type" << std::endl;
    LOG(Debug) << " - cb  : " << g->_isCheckerBoarded << std::endl;
    LOG(Debug) << " - fdim: " << g->_fdimensions << std::endl;

    std::vector<typename F::FermionField> srcVec(1,g);
    std::vector<typename F::FermionField> outVec(1,g);

    srcVec.resize(totSizeS,g);
    outVec.resize(totSizeS,g);

    std::vector<int> seeds5({1,2,3,4});
    GridParallelRNG          RNG5(g);
    RNG5.SeedFixedIntegers(seeds5);

    for (int eb = minBatchSizeE; eb <= maxBatchSizeE; eb += stepSize)
    {
        for (int sb = minBatchSizeS; sb <= maxBatchSizeS; sb += stepSize)
        {
            LOG(Message) << "Scan Source batch size: " << sb << std::endl;
            LOG(Message) << "Eigenvector batch size: " << eb << std::endl;
                
            for (auto &s: srcVec)
            {
                random(RNG5,s);
            }

            for (auto &v: outVec)
            {
                v = Zero();
            }

            ProjAccumRunner<F>(srcVec, outVec, eb, sb, totSizeE);

        }
    }
}

template <typename F> 
void scannerCoarse(GridBase *g, GridBase *gc,
             unsigned int minBatchSizeS, unsigned int maxBatchSizeS,
             unsigned int totSizeE, unsigned int totSizeS, 
             unsigned int stepSize)
{
    LOG(Debug) << "Check Grid type" << std::endl;
    LOG(Debug) << " - cb  : " << g->_isCheckerBoarded << std::endl;
    LOG(Debug) << " - fdim: " << g->_fdimensions << std::endl;

    const int nbasis = NBASIS;

    assert(nbasis<totSizeE);

    typedef iVector<vTComplex, nbasis>         CoarseSiteVector;
    typedef Lattice<CoarseSiteVector>          CoarseField;

    unsigned int coarseEvecSize;

    std::vector<typename F::FermionField> srcVec(totSizeS,g);
    std::vector<typename F::FermionField> outVec(totSizeS,g);
    CoarseEigenPack<typename F::FermionField, CoarseField> Epack(1,1,g,gc); 

    std::vector<int> seeds({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    std::vector<int> seeds5C({9,10,11,12});
    GridSerialRNG            RNG;      RNG.SeedFixedIntegers(seeds);
    GridParallelRNG          RNG5(g);  RNG5.SeedFixedIntegers(seeds5);
    GridParallelRNG          RNG5C(gc);  RNG5.SeedFixedIntegers(seeds5C);

    GridStopWatch w1;

    coarseEvecSize = totSizeE - nbasis;

    Epack.resize(nbasis, coarseEvecSize, g, gc);
    
    //for (auto &s: srcVec)
    //{
    //    random(RNG5,s);
    //}

    for (auto &e: Epack.evec)
    {
        random(RNG5,e);
    }

    for (auto &e: Epack.eval)
    {
        random(RNG,e);
    }
        
    for (auto &e: Epack.evecCoarse)
    {
        random(RNG5C,e);
    }

    for (auto &e: Epack.evalCoarse)
    {
        random(RNG,e);
    }

    //for (auto &v: outVec)
    //{
    //    v = Zero();
    //}
    LOG(Message) << "Start LC Deflation with subspace of size " << nbasis << std::endl;
    
    LocalCoherenceDeflatedGuesser<typename F::FermionField, CoarseField> Guesser(Epack.evec, Epack.evecCoarse, Epack.evalCoarse);

    for (int sb = minBatchSizeS; sb <= maxBatchSizeS; sb += stepSize)
    {
        LOG(Message) << "Scan Source size: " << sb << std::endl;

        srcVec.resize(sb,g);
        outVec.resize(sb,g);
            
        for (auto &s: srcVec)
        {
            random(RNG5,s);
        }

        for (auto &v: outVec)
        {
            v = Zero();
        }

        w1.Start();
        Guesser(srcVec, outVec);
        w1.Stop();
        LOG(Message) << "LC Deflation with subspace of size " << nbasis << " took: " << w1.Elapsed() << std::endl;
        w1.Reset();
    }
}

int main(int argc, char *argv[])
{
    unsigned int minBatchSizeE = 0;
    unsigned int minBatchSizeS = 0;
    unsigned int maxBatchSizeE;
    unsigned int maxBatchSizeS;
    unsigned int totSizeE;
    unsigned int totSizeS;
    unsigned int stepSize;
    unsigned int Ls;
    bool single;
    bool coarse;
    bool rb;

    if (argc < 11)
    {
        std::cerr << "usage: " << argv[0] << " <Ls> <RB {0|1}> <coarse {0|1}> <single {0|1}> <minBatchSizeEV> <maxBatchSizeEV> <minBatchSizeSrc> <maxBatchSizeSrc> <Total EV> <Total Src> <Step Size> [Grid options]";
        std::cerr << std::endl;
    }

    Ls = std::stoi(argv[1]);
    rb = (std::string(argv[2]) == "1");
    coarse = (std::string(argv[3]) == "1");
    single = (std::string(argv[4]) == "1");
    minBatchSizeE = std::stoi(argv[5]);
    maxBatchSizeE = std::stoi(argv[6]);
    minBatchSizeS = std::stoi(argv[7]);
    maxBatchSizeS = std::stoi(argv[8]);
    totSizeE      = std::stoi(argv[9]);
    totSizeS      = std::stoi(argv[10]);
    stepSize      = std::stoi(argv[11]);

    if (minBatchSizeE < 1 || minBatchSizeS < 1)
    {
        std::cerr << "--- minBatchSize's have to be greater than 1 ---";
        std::cerr << std::endl;
    }

    Grid_init(&argc, &argv);

    LOG(Message) << "=== Inputs ===" << std::endl;
    LOG(Message) << "Ls: " << Ls << " rb: " << rb << std::endl;
    LOG(Message) << "Total sources: " << totSizeS << " Total Eigenvectors: " << totSizeE << std::endl;
    LOG(Message) << "minBatchSizeE: " << minBatchSizeE << " maxBatchSizeE: " << maxBatchSizeE << std::endl;
    LOG(Message) << "minBatchSizeS: " << minBatchSizeS << " maxBatchSizeS: " << maxBatchSizeS << std::endl;
    LOG(Message) << "Scan Step Size: " << stepSize << std::endl;

    int64_t threads = GridThread::GetThreads();
    auto    mpi     = GridDefaultMpi();

    auto *g = makeGrid(Ls, rb, single);

    LOG(Message) << "Grid is setup to use " << threads << " threads" << std::endl;
    LOG(Message) << "MPI partition " << mpi << std::endl;

    const uint64_t nsimd = g->Nsimd();
    const uint64_t sites = g->oSites();

    LOG(Debug) << "Check Grid type" << std::endl;
    LOG(Debug) << " - cb  : " << g->_isCheckerBoarded << std::endl;
    LOG(Debug) << " - fdim: " << g->_fdimensions << std::endl;

    if (coarse)
    {
        auto *g = makeGrid(Ls, rb, 0, 0);
        auto *gc = makeGrid(Ls, rb, 0, 1,{2,2,2,2});

        scannerCoarse<FIMPL>(g, gc, minBatchSizeS, maxBatchSizeS, totSizeE, totSizeS, stepSize);
    }
    else
    {
        auto *g = makeGrid(Ls, rb, single);

        if(single)
        {
            scanner<FIMPLF>(g, single, minBatchSizeE, maxBatchSizeE, minBatchSizeS, maxBatchSizeS, totSizeE, totSizeS, stepSize);
        }
        else
        {
            scanner<FIMPL>(g, single, minBatchSizeE, maxBatchSizeE, minBatchSizeS, maxBatchSizeS, totSizeE, totSizeS, stepSize);
        }
    }

    LOG(Message) << "---End of Scan---" << std::endl;

    Grid_finalize();
}
