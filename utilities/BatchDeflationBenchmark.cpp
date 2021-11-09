#include <Grid/Grid.h>

using namespace Grid;

template<typename Field> 
void projAccumulate(const std::vector<Field> &in, std::vector<Field> &out,
                    std::vector<Field> &evec,
                    std::vector<Complex> &eval,
                    const unsigned int evBatchSize,
                    const unsigned int si, const unsigned int sf)
{
    GridBase *g       = in[0].Grid();
    double   lVol     = g->lSites();
    double   siteSize = sizeof(typename Field::scalar_object);
    double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
    double   nIt      = evBatchSize*(sf - si);
    double   t        = 0.;
    
    t -= usecond();
    for (unsigned int i = 0; i < evBatchSize; ++i)
    for (unsigned int j = si; j < sf; ++j)
    {

        axpy(out[j], TensorRemove(innerProduct(evec[i], in[j]))/eval[i], evec[i], out[j]);

    }
    t += usecond();

    // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
    std::cout << GridLogMessage << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;

};

template<typename Field>
void ProjAccumRunner(std::vector<Field> &in, std::vector<Field> &out, unsigned int eb, unsigned int sb)
{
    GridBase *g       = in[0].Grid();

    std::vector<Field> eigenVec(1,g);
    std::vector<Complex> eigenval;
    
    std::vector<int> seeds({1,2,3,4});
    GridSerialRNG            RNG;   RNG.SeedFixedIntegers(seeds);
    std::vector<int> seeds5({5,6,7,8}); // change line up
    GridParallelRNG          RNG5(g);  RNG5.SeedFixedIntegers(seeds5);

    GridStopWatch w1;
    GridTime ProjAccum = w1.Elapsed() - w1.Elapsed();
    
    unsigned int evSize = 4*eb;

    eigenVec.resize(eb,g);
    eigenval.resize(eb);

    std::cout << GridLogMessage << "ProjAccumRunner start" << std::endl;
    
    for (int i = 0; i < evSize; i += eb)
    {

        unsigned int evBlockSize = std::min(evSize - i, eb);

        // make new ev and eval of size evBlockSize to simulate ev IO.
        
        std::cout << GridLogMessage << "New eigenvector picks" << std::endl;

        // touch eigen vectors to bring back to host

        for (auto &e: eigenVec)
        {
            random(RNG5,e);
        }

        for (auto &e: eigenval)
        {
            random(RNG,e);

            std::cout << GridLogDebug << "eigenval random pick" << e << std::endl;
        }

        std::cout << GridLogMessage << "evBlockSize: " << evBlockSize << std::endl;

        for(unsigned int j = 0; j < in.size(); j += sb)
        {
            unsigned int srcBlockSize = std::min((int)in.size() - j, sb);

            std::cout << GridLogMessage << "srcBlockSize: " << srcBlockSize << std::endl;

            w1.Start();
            projAccumulate<Field>(in, out, eigenVec, eigenval, evBlockSize, j, j + srcBlockSize);
            w1.Stop();
            ProjAccum += w1.Elapsed();
            w1.Reset();
        }
    }

    std::cout << GridLogMessage << "ProjAccumRunner end" << std::endl;
    std::cout << GridLogMessage << "ProjAccum total: " << ProjAccum << std::endl;

}

inline void makeGrid(std::shared_ptr<GridBase> &gPt, 
                     const std::shared_ptr<GridCartesian> &gBasePt,
                     const unsigned int Ls = 1, const bool rb = false)
{
  if (rb)
  {
    if (Ls > 1)
    {
      gPt.reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, gBasePt.get()));
    }
    else
    {
      gPt.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(gBasePt.get()));
    }
  }
  else
  {
    if (Ls > 1)
    {
        gPt.reset(SpaceTimeGrid::makeFiveDimGrid(Ls, gBasePt.get()));
    }
    else
    {
        gPt = gBasePt;
    }
  }
}

void scanner(const Coordinate &latt, unsigned int Ls, bool rb,
             unsigned int minBatchSizeE, unsigned int maxBatchSizeE, 
             unsigned int minBatchSizeS, unsigned int maxBatchSizeS)
{

    // set up grid objects

    auto mpi = GridDefaultMpi();
    auto simd = GridDefaultSimd(latt.size(), LatticeFermion::vector_type::Nsimd());

    std::shared_ptr<GridCartesian> gBasePt(SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi));
    std::shared_ptr<GridBase>      gPt;

    makeGrid(gPt, gBasePt, Ls, rb);

    GridBase *g = gPt.get();

    std::vector<LatticeFermion> srcVec(1,g);
    std::vector<LatticeFermion> outVec(1,g);
    
    // resize

    unsigned int SrcFactor = 4;

    srcVec.resize(SrcFactor*maxBatchSizeS,g);
    outVec.resize(SrcFactor*maxBatchSizeS,g);

    // fill vectors with random data

    std::vector<int> seeds5({5,6,7,8});

    GridParallelRNG          RNG5(g);  RNG5.SeedFixedIntegers(seeds5);

    // scan over sourse batch size and eigenvector batch size up to max input

    for (int eb = minBatchSizeE; eb <= maxBatchSizeE; eb++)
    {
        for (int sb = minBatchSizeS; sb <= maxBatchSizeS; sb++)
        {
            std::cout << GridLogMessage << "Scan Source batch size: " << sb << " Eigenvector batch size: " << eb << std::endl;

            // New source vector
                
            for (auto &s: srcVec)
            {
                random(RNG5,s);
            }
         
                // zero the out vector

            for (auto &v: outVec)
            {
                v = Zero();
            }

            ProjAccumRunner<LatticeFermion>(srcVec, outVec, eb, sb);        
            
        }
    }
}

int main(int argc, char *argv[])
{
    // input parameters

    unsigned int minBatchSizeE = 0;
    unsigned int minBatchSizeS = 0;
    unsigned int maxBatchSizeE;
    unsigned int maxBatchSizeS;
    unsigned int Ls;
    bool rb;

    if (argc < 6)
    {
        std::cerr << "usage: " << argv[0] << " <Ls> <RB {0|1}> <minBatchSizeE> <maxBatchSizeE> <minBatchSizeS> <maxBatchSizeS>  [Grid options]";
        std::cerr << std::endl;
    }

    Ls = std::stoi(argv[1]);
    rb = (std::string(argv[2]) == "1");

    minBatchSizeE = std::stoi(argv[3]);
    maxBatchSizeE = std::stoi(argv[4]);
    minBatchSizeS = std::stoi(argv[5]);
    maxBatchSizeS = std::stoi(argv[6]);

    if (minBatchSizeE < 1 || minBatchSizeS < 1)
    {
        std::cerr << "--- minBatchSize's have to be greater than 1 ---";
        std::cerr << std::endl;
    }

    Grid_init(&argc, &argv);

    std::cout << GridLogMessage << "=== Inputs ===" << std::endl;
    std::cout << GridLogMessage << "Ls: " << Ls << " rb: " << rb << std::endl;
    std::cout << GridLogMessage << " minBatchSizeE: " << minBatchSizeE << " maxBatchSizeE: " << maxBatchSizeE << std::endl;
    std::cout << GridLogMessage << " minBatchSizeS: " << minBatchSizeE << " maxBatchSizeS: " << maxBatchSizeE << std::endl;

    int64_t threads = GridThread::GetThreads();
    auto    mpi     = GridDefaultMpi();

    std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;
    std::cout << GridLogMessage << "MPI partition " << mpi << std::endl;

    scanner(GridDefaultLatt(), Ls, rb, minBatchSizeE, maxBatchSizeE, minBatchSizeE, maxBatchSizeS);

    std::cout << GridLogMessage << "---End of Scan---" << std::endl;

    // Print a summary

    std::cout << GridLogMessage << "Grid finalize" << std::endl;

    Grid_finalize();
}