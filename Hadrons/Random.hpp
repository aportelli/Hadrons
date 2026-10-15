#ifndef Hadrons_Random_hpp_
#define Hadrons_Random_hpp_

#include <Grid/lattice/Lattice_rng.h>
#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

class HadronsSerialRNG
{
public:
    typedef uint32_t RngType;

    HadronsSerialRNG()
    {
        _generators.resize(1);
        _uid.resize(1,std::uniform_int_distribution<RngType>() );
    }

    void SeedFixedIntegers(const std::vector<int> &seeds)
    {
        CartesianCommunicator::BroadcastWorld(0,(void *)&seeds[0],sizeof(int)*seeds.size());
        std::seed_seq src(seeds.begin(),seeds.end());
        Seed(src,0);
    }

    void SeedUniqueString(const std::string &s)
    {
      std::vector<int> seeds;
      std::stringstream sha;
      seeds = GridChecksum::sha256_seeds(s);
      for(int i=0;i<seeds.size();i++) { 
        sha << std::hex << seeds[i];
      }
      std::cout << GridLogMessage << "Intialising Hadrons serial RNG with unique string '" 
                << s << "'" << std::endl;
      std::cout << GridLogMessage << "Seed SHA256: " << sha.str() << std::endl;
      SeedFixedIntegers(seeds);
    }

    template <class distribution> 
    inline void fill(RngType &l,std::vector<distribution> &dist)
    {
        dist[0].reset();
        fillScalar(l,dist[0],_generators[0]);
        CartesianCommunicator::BroadcastWorld(0,(void *)&l,sizeof(l));
    }
protected:
    template<class source> void Seed(source &src, int gen)
    {
        _generators[gen] = RngEngine(src);
    }
public:
    std::vector<std::uniform_int_distribution<RngType>> _uid;
protected:
    typedef std::mt19937 RngEngine;
    typedef uint32_t     RngStateType;
    static const int     RngStateCount = std::mt19937::state_size;
    std::vector<RngEngine>                 _generators;
};

inline void uid(HadronsSerialRNG &rng, uint32_t &out) { rng.fill(out, rng._uid); }

END_HADRONS_NAMESPACE


#endif // Hadrons_Random_hpp_
