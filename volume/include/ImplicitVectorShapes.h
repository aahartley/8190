
#ifndef __IMPLICITVECTORSHAPES_H__
#define __IMPLICITVECTORSHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include "Fields.h"
#include "FullGrids.h"
#include "SparseGrids.h"
#include <cmath>
#include <vector>
#include <iostream>

namespace lux
{



class ConstantVectorVolume : public Volume<Vector> 
{
  public:

    ConstantVectorVolume( const Vector& v ) ;

   ~ConstantVectorVolume(){}


    const Vector eval( const Vector& P ) const; 

    //const Matrix grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Vector value;
    Matrix gradvalue;
};

class AddVectorVolume : public Volume<Vector> 
{
  public:

    AddVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) ;

    AddVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~AddVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Add";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:
  
    VectorField elem1, elem2;
};

class SubtractVectorVolume : public Volume<Vector> 
{
  public:


    SubtractVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~SubtractVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "SubtractVectorVolume";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:
  
    VectorField elem1, elem2;
};

class MultiplyVectorScalarVolume : public Volume<Vector> 
{
  public:


    MultiplyVectorScalarVolume( const VectorField& v1, const ScalarField& v2 ) ;

    ~MultiplyVectorScalarVolume(){}


    const Vector eval( const Vector& P ) const;

    //const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "MultiplyVectorScalarVolume";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:
  
    VectorField elem1;
    ScalarField elem2;
};


class DivideVectorVolume : public Volume<Vector> 
{
  public:


    DivideVectorVolume( const VectorField& v1, const ScalarField& v2 ) ;

    ~DivideVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    //const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "DivideVectorVolume";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:
  
    VectorField elem1;
    ScalarField elem2;
};



class GradientVectorVolume : public Volume<Vector> 
{
  public:

    GradientVectorVolume( Volume<float>* v, const float dx = 0.001 ) ;

    GradientVectorVolume( const ScalarField& v, const float dx = 0.001 ) ;

    ~GradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Grad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    float step;
};

 class IdentityVectorVolume : public Volume<Vector> 
 {
   public:
  
     IdentityVectorVolume();
  
     ~IdentityVectorVolume(){}
  
  
     const Vector eval( const Vector& P ) const ;
  
     const Matrix grad( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "Identity";
        return lbl;
     }
  
  
   private:
  
     Matrix gradvalue;
 };
  
class WarpVectorVolume : public Volume<Vector> 
{
  public:

    WarpVectorVolume( const VectorField& v, const VectorField& u ) ;

   ~WarpVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Warp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + warp->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem, warp;
};


 class GriddedSGridVector : public Volume<Vector> 
 {
   public:
  
     GriddedSGridVector( const SVectorGrid& g );
  
    ~GriddedSGridVector(){}
  
     const Vector eval( const Vector& P ) const;
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     SVectorGrid svgrid;
     float dx, dy, dz;
 };
  
class AdvectVectorVolume : public Volume<Vector>
{
  public:
    AdvectVectorVolume( const VectorField& v, const VectorField& u, const float delt ): elem(v), velocity(u), dt(delt) {}
    ~AdvectVectorVolume(){}
    const Vector eval( const Vector& P ) const
    {
      Vector X = P - ( velocity->eval(P) )*dt;
      return elem->eval(X);
    }
  private:
    const VectorField elem;
    VectorField velocity;
    float dt;
};

 class NoiseVectorVolume : public Volume<Vector> 
 {
   public:
  
     //NoiseVolume( Noise* n, const float d = 0.01 ); 
     NoiseVectorVolume(const _Noise& n){noise =n;dx = 0.1;}; 
  
    ~NoiseVectorVolume(){}
  
    const Vector eval( const Vector& P ) const
    {
      float nx = noise->eval(P+Vector(dx,0,0));
      float ny = noise->eval(P+Vector(0,dx,0));
      float nz = noise->eval(P+Vector(0,0,dx));
      Vector g( ny-nz, nz-nx, nx-ny );
      g /= dx;
      return g;
    }
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Noise";
       return lbl;
    }
  
   private:
  
     _Noise noise;
     float dx;
 };

}
#endif
