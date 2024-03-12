#ifndef NOISE_H
#define NOISE_H
#include "Vector.h"
#include <iostream>

namespace lux
{

class NoiseSource
{
  public:
    NoiseSource(){}
    virtual ~NoiseSource(){}

    virtual const float eval(const Vector& pos) const { return 0;}

};
typedef std::shared_ptr<NoiseSource> NoiseSrc;

class PerlinNoise : public NoiseSource
{
  public:
    PerlinNoise() 
    {
      for(int i = 0; i < 256; i++)
      {
        p[i] = permutation[i];
        p[256+i] = permutation[i];
      }
    }
    virtual ~PerlinNoise() {}
    virtual const float eval(const Vector& position) const 
    {
      Vector pos(position.X(), position.Y(), position.Z());
      int X = (int)std::floor(pos.X()) & 255; // FIND UNIT CUBE THAT
      int Y = (int)std::floor(pos.Y()) & 255; // CONTAINS POINT.
      int Z = (int)std::floor(pos.Z()) & 255;
      double x = pos[0];
      double y = pos[1];
      double z = pos[2];
      x -= std::floor(x); // FIND RELATIVE X,Y,Z
      y -= std::floor(y); // OF POINT IN CUBE.
      z -= std::floor(z);
      double u = fade(x), // COMPUTE FADE CURVES
             v = fade(y), // FOR EACH OF X,Y,Z.
             w = fade(z);
      int A = p[X ]+Y, AA = p[A]+Z, AB = p[A+1]+Z, // HASH COORDINATES OF
          B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z; // THE 8 CUBE CORNERS,
      return lerp(w, lerp(v, lerp(u, grad(p[AA ], x , y , z ), // AND ADD
      grad(p[BA ], x-1, y , z )), // BLENDED
      lerp(u, grad(p[AB ], x , y-1, z ), // RESULTS
      grad(p[BB ], x-1, y-1, z ))),// FROM 8
      lerp(v, lerp(u, grad(p[AA+1], x , y , z-1 ), // CORNERS
      grad(p[BA+1], x-1, y , z-1 )), // OF CUBE
      lerp(u, grad(p[AB+1], x , y-1, z-1 ),
      grad(p[BB+1], x-1, y-1, z-1 ))));

    }
    double fade(double t) const { return t * t * t * (t * (t * 6 - 15) + 10); }
    double lerp(double t, double a, double b) const { return a + t * (b - a); }
    double grad(int hash, double x, double y, double z) const
    {
      int h = hash & 15; // CONVERT LO 4 BITS OF HASH CODE
      double u = h<8 ? x : y, // INTO 12 GRADIENT DIRECTIONS.
             v = h<4 ? y : h==12||h==14 ? x : z;
      return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
    }

  private:
    int p[512];
    int permutation[256] = { 151,160,137,91,90,15,
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180};

};


struct NoiseData
{
    NoiseData() :
      frequency (1.0f),
      octaves (1.0f),
      roughness (0.5f),
      fjump(2.0f),
      fade(1),
      fade_radius(1),
      clump(1),
      gamma(1),
    
      fade_x0(Vector(0,0,0)),
      translate (Vector(0,0,0)),
      dP(Vector(0,0,0))
    {} 

    float frequency;
    float octaves;
    float roughness;
    float fjump;
    float fade;
    float fade_radius;
    float clump;
    float gamma;
    Vector fade_x0;
    Vector translate; 
    Vector dP;

};


class Noise
{
  public:
    Noise() {}
    virtual ~Noise(){}

    virtual const float eval( const float x) const {return 0;}
    virtual const float eval(const Vector& x) const { return 0;}

    virtual void setParameters(const NoiseData& param){}
    virtual void getParameters(NoiseData& param) const {}
};
typedef std::shared_ptr<Noise> _Noise;



class FractalSum : public Noise
{
  public:
    FractalSum(NoiseSrc& n) :
      _noise(n),
      frequency (0.666666f),
      octaves (3.f),
      roughness (0.5f),
      fjump(2.f),
      translate (Vector(0,0,0))
    {} 

    virtual ~FractalSum(){}

    virtual const float eval(const Vector& x) const 
    {
        float roughFrac = (1 - roughness) / (1 - std::pow(roughness, octaves));
        float sum = 0;
        for(int i = 0; i < (int)octaves; i++)
        {
            sum += std::pow(roughness, i) * _noise.get()->eval((x-translate)* frequency * fjump);
        }
        return sum * roughFrac;
    }

    virtual void setParameters( const NoiseData& parameters )
    {
        octaves = parameters.octaves;
        fjump = parameters.fjump;
        roughness = parameters.roughness;
        frequency = parameters.frequency; 
        translate = parameters.translate; 
  
    }
  
    
    virtual void getParameters( NoiseData& parameters ) const
    {
        parameters.octaves = octaves;
        parameters.fjump = fjump;
        parameters.roughness = roughness;
        parameters.frequency = frequency; 
        parameters.translate = translate; 

    }

  private:
    NoiseSrc _noise;
    float frequency;
    float octaves;
    float roughness;
    float fjump;
    Vector translate; 
};


class FractalSumFade : public Noise
{
  public:
    FractalSumFade(NoiseSrc& n) :
      _noise(n),
      frequency (0.666666f),
      octaves (3.f),
      roughness (0.5f),
      fjump(2.f),
      fade(1),
      fade_radius(3),
      fade_x0(Vector(0,0,0)),
      translate (Vector(0,0,0))
    {} 

    virtual ~FractalSumFade(){}

    virtual const float eval(const Vector& x) const 
    {
        float roughFrac = (1 - roughness) / (1 - std::pow(roughness, octaves));
        float sum = 0;
        for(int i = 0; i < (int)octaves; i++)
        {
            sum += std::pow(roughness, i) * _noise.get()->eval((x-translate)* frequency * fjump);
        }

        return (sum * roughFrac) * fadeFactor(x);
    }

    virtual void setParameters( const NoiseData& parameters )
    {
        octaves = parameters.octaves;
        fjump = parameters.fjump;
        roughness = parameters.roughness;
        frequency = parameters.frequency; 
        fade = parameters.fade;
        fade_radius = parameters.fade_radius;
        fade_x0 = parameters.fade_x0;
        translate = parameters.translate; 
  
    }
  
    
    virtual void getParameters( NoiseData& parameters ) const
    {
        parameters.octaves = octaves;
        parameters.fjump = fjump;
        parameters.roughness = roughness;
        parameters.frequency = frequency;
        parameters.fade = fade;
        parameters.fade_radius = fade_radius;
        parameters.fade_x0 = fade_x0; 
        parameters.translate = translate; 

    }

    const float fadeFactor(const Vector& x) const
    {
        float q = (1 - ((x - fade_x0).magnitude()/fade_radius) ) / fade;
        if( fade < 0) return 1;
        else if(q > 1) return 1;
        else if(q < 0) return 0;
        else return q;
    }

  private:
    NoiseSrc _noise;
    float frequency;
    float octaves;
    float roughness;
    float fjump;
    float fade;
    float fade_radius;
    Vector fade_x0;
    Vector translate; 
};








}


#endif
