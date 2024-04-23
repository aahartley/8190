
#ifndef __IMPLICITVOLUMESHAPES_H__
#define __IMPLICITVOLUMESHAPES_H__

#include "Volume.h"
#include "FullGrids.h"
#include "SparseGrids.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "Fields.h"
#include "Noise.h"


namespace lux
{


class ConstantVolume : public Volume<float> 
{
  public:

    ConstantVolume( const float v ) ;

   ~ConstantVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const; 

   virtual std::string typelabel() { return "Constant"; }


  private:

    float value;
    Vector gradvalue;
};


class ExpVolume: public Volume<float> 
{
  public:

    ExpVolume( Volume<float>* v );

    ExpVolume( const ScalarField& v );

   ~ExpVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Exp";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;

};



class NegateVolume : public Volume<float> 
{
  public:

    NegateVolume( const ScalarField&  v1);

    ~NegateVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;


  private:

    const ScalarField e1;
};



class AddVolume : public Volume<float> 
{
  public:

    AddVolume( Volume<float> * v1, Volume<float> * v2 );

    AddVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~AddVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Add";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};

class MultiplyVolume : public Volume<float> 
{
  public:

    //MultiplyVolume( Volume<float> * v1, Volume<float> * v2 );

    MultiplyVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~MultiplyVolume(){}


    const float eval( const Vector& P ) const;

    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Mult";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};
class MultiplyVectorVolume : public Volume<float> 
{
  public:


    MultiplyVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~MultiplyVectorVolume(){}


    const float eval( const Vector& P ) const;

    //const Vector grad( const Vector& P ) const; //fix? matrix to vector


    virtual std::string typelabel() 
    { 
       std::string lbl = "MultiplyVectorVolume";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:
  
    VectorField elem1;
    VectorField elem2;
};
  
 class TranslateVolume : public Volume<float> 
 {
   public:
  
     //TranslateVolume( Volume<float> * v, const Vector& s );
  
     TranslateVolume( const ScalarField& v, const Vector& s );
  
     ~TranslateVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
  
     const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Translate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
  
  
   private:
  
     const ScalarField elem;
     Vector translate;
 };

  class ScaleVolume : public Volume<float> 
 {
   public:
  
     //ScaleVolume( Volume<float> * v, const Vector& s );
  
     ScaleVolume( const ScalarField& v, const Vector& s );
  
     ~ScaleVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
  
     const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Scale";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
  
  
   private:
  
     const ScalarField elem;
     Vector scale;
 };

  class RotateVolume : public Volume<float> 
 {
   public:
  
     //RotateVolume( Volume<float> * v, const Vector& s );
  
     RotateVolume( const ScalarField& v, const Vector& axis, float angle );
  
     ~RotateVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
  
     const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Rotate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
  
  
   private:
  
     const ScalarField elem;
     const Vector axis;
     float angle;
 };

class UnionVolume : public Volume<float>
{
  public:
    UnionVolume(const ScalarField& v1, const ScalarField& v2);
    ~UnionVolume() {}

    const float eval( const Vector& P ) const;
     const Vector grad(  const Vector& P ) const;

  private:
    const ScalarField e1, e2;
};

class IntersectionVolume : public Volume<float>
{
  public:
    IntersectionVolume(const ScalarField& v1, const ScalarField& v2);
    ~IntersectionVolume() {}

    const float eval( const Vector& P ) const;
     const Vector grad(  const Vector& P ) const;

  private:
    const ScalarField e1, e2;
};

class CutoutVolume : public Volume<float>
{
  public:
    CutoutVolume(const ScalarField& v1, const ScalarField& v2);
    ~CutoutVolume() {}

    const float eval( const Vector& P ) const;
     const Vector grad(  const Vector& P ) const;

  private:
    const ScalarField e1, e2;
};

class MaskVolume : public Volume<float> 
 {
   public:
  
     //MaskVolume( Volume<float> * v );
  
     MaskVolume( const ScalarField& v );
  
     ~MaskVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
     const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Mask";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   private:
  
     const ScalarField elem;
  
 };
  
 class ClampVolume : public Volume<float> 
 {
   public:
  
    // ClampVolume( Volume<float> * v, float minv, float maxv );
  
     ClampVolume( const ScalarField& v, float minv, float maxv );
  
  
     ~ClampVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
     const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Clamp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   private:
  
    const ScalarField elem;
     float vmin, vmax;
 };

class SphereVolume : public Volume<float> 
{
  public:

    SphereVolume( const Vector& cen, const float rad );

   ~SphereVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Sphere"); }


  private:

    Vector center;
    float radius;
};

class EllipseVolume : public Volume<float> 
{
  public:

    EllipseVolume( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) ;

   ~EllipseVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Ellipse"); }


  private:

    Vector center, normal;
    float radMajor, radMinor;
};

class TorusVolume : public Volume<float> 
{
  public:

    TorusVolume( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) ;

   ~TorusVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Torus"); }


  private:

    Vector center, normal;
    float radMajor, radMinor;
};


class BoxVolume : public Volume<float> 
{
  public:

    BoxVolume( const Vector& cen, const float rad, const float q ) ;

   ~BoxVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Torus"); }


  private:

    Vector center;
    float radius, q;
};


class PlaneVolume : public Volume<float> 
{
  public:

    PlaneVolume( const Vector& cen, const Vector& norm, const Vector& pos) ;

   ~PlaneVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center, normal, position;
};

class ConeVolume : public Volume<float> 
{
  public:

    ConeVolume( const Vector& cen, const Vector& norm, const float h, const float theta) ;

   ~ConeVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center, normal;
    float h, thetaMax;
};

class CylinderVolume : public Volume<float> 
{
  public:

    CylinderVolume( const Vector& cen, const Vector& norm, const float radius) ;

   ~CylinderVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center, normal;
    float radius;
};

class IcosahedronVolume : public Volume<float> 
{
  public:

    IcosahedronVolume() ;

   ~IcosahedronVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Ico"); }


 

};

class SteinerPatchVolume : public Volume<float> 
{
  public:

    SteinerPatchVolume( const Vector& cen) ;

   ~SteinerPatchVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center;
};



 class GriddedGridVolume : public Volume<float> 
 {
   public:
  
     GriddedGridVolume( const ScalarGrid& g );
  
    ~GriddedGridVolume(){}
  
     const float eval( const Vector& P ) const;
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     ScalarGrid scgrid;
     float dx, dy, dz;
 };
 class GriddedSGridVolume : public Volume<float> 
 {
   public:
  
     GriddedSGridVolume( const SScalarGrid& g );
  
    ~GriddedSGridVolume(){}
  
     const float eval( const Vector& P ) const;
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     SScalarGrid scgrid;
     float dx, dy, dz;
 };
  
 class NoiseVolume : public Volume<float> 
 {
   public:
  
     //NoiseVolume( Noise* n, const float d = 0.01 ); 
     NoiseVolume(const _Noise& n); 
  
    ~NoiseVolume(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Noise";
       return lbl;
    }
  
   private:
  
     _Noise noise;
 };

  class SDFVolume : public Volume<float> 
 {
   public:
  
     //NoiseVolume( Noise* n, const float d = 0.01 ); 
     SDFVolume(const ScalarField& e, int n); 
  
    ~SDFVolume(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "SDFVolume";
       return lbl;
    }
  
   private:
     const ScalarField elem;
     int N;
 };


 class PyroclasticSphere : public Volume<float> 
 {
   public:
  
     //NoiseVolume( Noise* n, const float d = 0.01 ); 
     PyroclasticSphere(const Vector& cen, const float r, const float amp, const float gam, const _Noise& n ); 
  
    ~PyroclasticSphere(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "PyroclasticSphere";
       return lbl;
    }
  
   private:
      Vector center;
      float rad;
      float amplitude;
      float gamma;
     _Noise noise;
 };

 class PyroclasticVolume : public Volume<float> 
 {
   public:
  
     PyroclasticVolume(const ScalarField& e, const float amp, const float gam, const _Noise& n, int iter ); 
  
    ~PyroclasticVolume(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "PyroclasticVolume";
       return lbl;
    }
  
   private:
      const ScalarField elem;
      float amplitude;
      float gamma;
      int N;
     _Noise noise;
     ScalarField fspn, warpV;
     VectorField Xnpt;
 };

 class PyroclasticTerrain : public Volume<float> 
 {
   public:
  
     //NoiseVolume( Noise* n, const float d = 0.01 ); 
     PyroclasticTerrain(const Vector& xp, const Vector& norm, const float ampP, const float ampN, const float gamP, const float gamN, const float trans, const _Noise& n ); 
  
    ~PyroclasticTerrain(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "PyroclasticTerrain";
       return lbl;
    }
  
   private:
      Vector Xp;
      Vector normal;
      float posAmp;
      float negAmp;
      float posGam;
      float negGam;
      float transition;
     _Noise noise;
 };

class WarpVolume : public Volume<float>
{
  public:
    WarpVolume( const ScalarField& v, VectorField& map );
    const float eval( const Vector& P ) const;
    const Vector grad(  const Vector& P ) const;  

  private:
    const ScalarField elem;
    VectorField mapX;
};



class DivergenceVVolume : public Volume<float> 
{
  public:

    DivergenceVVolume( const SVectorGrid& u ) ;

   ~DivergenceVVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad( const Vector& P ) const; 

   virtual std::string typelabel() { return "DivergenceVVolume"; }


  private:
    SVectorGrid U;
};

class AdvectVolume : public Volume<float>
{
  public:
    AdvectVolume( const ScalarField& v, const VectorField& u, const float delt ): elem(v), velocity(u), dt(delt) {}
    ~AdvectVolume(){}
    const float eval( const Vector& P ) const
    {
      Vector X = P -  (velocity->eval(P)) * dt;
      return elem->eval(X);
    }
   // const Vector grad( const Vector& P ) const; 

  private:
    const ScalarField elem;
    VectorField velocity;
    float dt;
};

}





#endif



