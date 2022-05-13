#include <Grid/Grid.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Modules/MGuesser/BatchDeflationUtils.hpp>

using namespace Grid;
using namespace Hadrons;

template<typename Field>
void ProjAccumRunner(std::vector<Field> &in, std::vector<Field> &out, unsigned int eb, unsigned int sb, unsigned int totSizeE)
{
    GridBase *g       = in[0].Grid();
    EigenPack<Field> Epack(eb, g, g);
    
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
            BatchDeflationUtils::projAccumulate(in, out, Epack.evec, Epack.eval,
                                                0, evBlockSize, j, j + srcBlockSize);
            w1.Stop();
            ProjAccum += w1.Elapsed();
            w1.Reset();
        }
    }

    LOG(Message) << "ProjAccumRunner end" << std::endl;
    LOG(Message) << "ProjAccum total: " << ProjAccum << std::endl;

}

GridBase * makeGrid(const unsigned int Ls = 1, const bool rb = false)
{
  auto &env  = Environment::getInstance();
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

void scanner(unsigned int Ls, bool rb,
             unsigned int minBatchSizeE, unsigned int maxBatchSizeE, 
             unsigned int minBatchSizeS, unsigned int maxBatchSizeS,
             unsigned int totSizeE, unsigned int totSizeS, unsigned int stepSize)
{
    auto *g = makeGrid(Ls, rb);

    LOG(Debug) << "Check Grid type" << std::endl;
    LOG(Debug) << " - cb  : " << g->_isCheckerBoarded << std::endl;
    LOG(Debug) << " - fdim: " << g->_fdimensions << std::endl;

    std::vector<LatticeFermion> srcVec(1,g);
    std::vector<LatticeFermion> outVec(1,g);

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

            ProjAccumRunner<LatticeFermion>(srcVec, outVec, eb, sb, totSizeE);        
        }
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
    bool rb;

    if (argc < 9)
    {
        std::cerr << "usage: " << argv[0] << " <Ls> <RB {0|1}> <minBatchSizeEV> <maxBatchSizeEV> <minBatchSizeSrc> <maxBatchSizeSrc> <Total EV> <Total Src> <Step Size> [Grid options]";
        std::cerr << std::endl;
    }

    Ls = std::stoi(argv[1]);
    rb = (std::string(argv[2]) == "1");
    minBatchSizeE = std::stoi(argv[3]);
    maxBatchSizeE = std::stoi(argv[4]);
    minBatchSizeS = std::stoi(argv[5]);
    maxBatchSizeS = std::stoi(argv[6]);
    totSizeE      = std::stoi(argv[7]);
    totSizeS      = std::stoi(argv[8]);
    stepSize      = std::stoi(argv[9]);

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

    LOG(Message) << "Grid is setup to use " << threads << " threads" << std::endl;
    LOG(Message) << "MPI partition " << mpi << std::endl;

    scanner(Ls, rb, minBatchSizeE, maxBatchSizeE, minBatchSizeS, maxBatchSizeS, totSizeE, totSizeS, stepSize);

    LOG(Message) << "---End of Scan---" << std::endl;

    Grid_finalize();
}
