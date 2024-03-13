#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <thread>
namespace lux
{


class Random {
  public:
	Random(int num) {uIndex = 0; gIndex = 0;numb =num; uniformR.resize(num);}
	~Random() {}




    float createRandUniform(float min, float max) 
    {
        unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();

        static std::default_random_engine generator(seed);
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(generator);
    }

    float createRandGauss(float mean, float stdv)
    {
        unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();

        static std::default_random_engine generator(seed);
        std::normal_distribution<float> distribution(mean, stdv);
        return distribution(generator);
    }
    void populateRand()
    {
        //#pragma omp parallel for
        for (int i = 0; i < numb; i++)
        {
            uniformR[i]=(createRandUniform(0, 1));
            //std::cout << uniformR[i] << '\n';
            //gaussR.push_back(createRandGauss(0, 1));
        }
        std::cout << "rand populated\n";
    }
    float randUniform(float min, float max)
    {
        //std::cout << uniformR.size() << '\n';
        if (uIndex >= uniformR.size())uIndex = 0;
        float rand = (max - min) * uniformR[uIndex] + min;
        uIndex+=1;
        //float rand = (max - min) * createRandUniform(0,1) + min;

        return rand;
    }

    float randGauss(float mean, float stdv)
    {
        if (gIndex == gaussR.size())gIndex = 0;
        return (stdv) * gaussR[gIndex++] + mean;
    }

  private:
  	std::vector<float> uniformR;
	std::vector<float> gaussR;
	int uIndex;
	int gIndex;
    int numb;
};


}

#endif