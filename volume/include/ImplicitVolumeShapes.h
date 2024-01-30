
#ifndef __IMPLICITVOLUMESHAPES_H__
#define __IMPLICITVOLUMESHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <vector>
#include <iostream>



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
  
     RotateVolume( const ScalarField& v, const Vector& s );
  
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
     Vector rotate;
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

    IcosahedronVolume( const Vector& cen) ;

   ~IcosahedronVolume(){}


    const float eval( const Vector& P ) const; 

    //const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center;
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

}





#endif



