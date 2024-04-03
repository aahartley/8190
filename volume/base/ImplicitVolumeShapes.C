
#include "ImplicitVolumeShapes.h"
#include <iostream>
using namespace std;

using namespace lux;


ConstantVolume::ConstantVolume( const float v ) :
       value (v),
       gradvalue (Vector(0,0,0))
    {}


const float ConstantVolume::eval( const Vector& P ) const 
{
   return value; 
}

const Vector ConstantVolume::grad( const Vector& P ) const { return gradvalue; }


ExpVolume::ExpVolume( Volume<float>* v ) :
      elem (v)
    {}

ExpVolume::ExpVolume( const ScalarField& v ) :
      elem (v)
    {}

const float ExpVolume::eval( const Vector& P ) const 
{
   return std::exp( elem->eval(P) ); 
}

const Vector ExpVolume::grad( const Vector& P ) const { return eval(P) * elem->grad(P); }



MultiplyVolume::MultiplyVolume( const ScalarField& v , const ScalarField& v2) :
      e1 (v), e2(v2)
    {}

const float MultiplyVolume::eval( const Vector& P ) const 
{
   return e1->eval(P) * e2->eval(P) ; 
}

// const Vector MultiplyVolume::grad( const Vector& P ) const { return e1->grad(P) * e2->grad(P); }




SphereVolume::SphereVolume( const Vector& cen, const float rad ) :
       center (cen),
       radius (rad)
    {}


const float SphereVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float rad = x.magnitude();
   return radius-rad;
}

const Vector SphereVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   if( x.magnitude() != 0 ){ return -x.unitvector(); }
   return Vector(0,1,0);
}


EllipseVolume::EllipseVolume( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) :
       center (cen),
       radMajor (radMaj),
       radMinor(radMin),
       normal(norm)

    {}


const float EllipseVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float z = x * normal;
   Vector xPerp = x - (z*normal);
   return 1 - ((z*z) / (radMajor*radMajor)) - ((xPerp.magnitude()*xPerp.magnitude()) / (radMinor*radMinor));
}

TorusVolume::TorusVolume( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) :
       center (cen),
       radMajor (radMaj),
       radMinor(radMin),
       normal(norm)

    {}


const float TorusVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   Vector xPerp = x - ((x * normal) * normal);
   return 4*(radMajor*radMajor) * (xPerp.magnitude() * xPerp.magnitude()) -
    std::pow(((x.magnitude() * x.magnitude() ) + (radMajor*radMajor) - (radMinor*radMinor)),2);
}


BoxVolume::BoxVolume( const Vector& cen, const float rad, const float q ) :
       center (cen),
       radius (rad),
       q(q)

    {}


const float BoxVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   return std::pow(radius, 2*q) - (std::pow(fabs(x.X()), 2*q) + std::pow(fabs(x.Y()), 2*q) + std::pow(fabs(x.Z()), 2*q));
}


PlaneVolume::PlaneVolume( const Vector& cen, const Vector& norm, const Vector& pos ) :
       center (cen),
       normal (norm),
       position(pos)

    {}


const float PlaneVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   return -1 * ((x - position) * normal);
}

ConeVolume::ConeVolume( const Vector& cen, const Vector& norm, const float h, const float theta ) :
       center (cen),
       normal (norm),
       h(h),
       thetaMax(theta)

    {}


const float ConeVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   if(x*normal < 0) return x * normal;
   else if(x*normal > h) return h -(x * normal);
   else if ( x*normal < h && x*normal > 0) return (x*normal) - (x.magnitude() * std::cos((thetaMax * M_PI)/180));
   else return -1;
}

CylinderVolume::CylinderVolume( const Vector& cen, const Vector& norm, const float rad) :
       center (cen),
       normal (norm),
       radius(rad)
      

    {}


const float CylinderVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   return radius - (x - ((x*normal)*normal)).magnitude();
}

IcosahedronVolume::IcosahedronVolume( ) 
      
    {}


const float IcosahedronVolume::eval( const Vector& P ) const 
{
   Vector x = P;
   float T = 1.61803399; // golden ratio
   if(x.magnitude() <= 1.8 * M_PI)
   {
      return std::cos(x.X()+T*x.Y()) + std::cos(x.X()-T*x.Y()) +
            std::cos(x.Y()+T*x.Z()) + std::cos(x.Y()-T*x.Z()) +
            std::cos(x.Z()-T*x.X()) + std::cos(x.Z()+T*x.X()) - 2;
   }
   else return -(1.8*M_PI);
}


SteinerPatchVolume::SteinerPatchVolume( const Vector& cen) :
       center (cen)
      
    {}


const float SteinerPatchVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   return -( (x.X()*x.X()) * (x.Y()*x.Y()) + (x.X()*x.X()) * (x.Z()*x.Z()) +
               (x.Y()*x.Y()) * (x.Z()*x.Z()) - x.X()*x.Y()*x.Z() );
}


AddVolume::AddVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float AddVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) + e2->eval(P);
}

const Vector AddVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) + e2->grad(P);
}


NegateVolume::NegateVolume( const ScalarField&  v1 ) :
      e1 (v1)
    {}

const float NegateVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) * -1;
}

const Vector NegateVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) * -1;
}


TranslateVolume::TranslateVolume( const ScalarField& v, const Vector& s ) :
      elem(v),
      translate(s)
   {}

const float TranslateVolume::eval( const Vector& P ) const
{
   return  elem->eval(P-translate); 
}

const Vector TranslateVolume::grad( const Vector& P ) const
{
   return  elem->grad(P-translate);
}

ScaleVolume::ScaleVolume( const ScalarField& v, const Vector& s ) :
      elem(v),
      scale(s)
   {}

const float ScaleVolume::eval( const Vector& P ) const
{
   Vector X(P.X()/scale.X(), P.Y()/scale.Y(), P.Z()/scale.Z());
   return  elem->eval(X); 
}

const Vector ScaleVolume::grad( const Vector& P ) const
{
   Vector X(P.X()/scale.X(), P.Y()/scale.Y(), P.Z()/scale.Z());
   return  elem->grad(X); 
}

RotateVolume::RotateVolume( const ScalarField& v, const Vector& s, float a ) :
      elem(v),
      axis(s),
      angle(a)
   {}

const float RotateVolume::eval( const Vector& P ) const
{
   Vector axisNorm = axis.unitvector();
	float sina = (std::sin(angle*M_PI/180.0));
	float cosa = (std::cos(angle*M_PI/180.0));

   Vector X = P*cosa + axisNorm *(axisNorm*P)*(1.0-cosa)+ (axisNorm^P)*sina;

   return  elem->eval(X); 
}

const Vector RotateVolume::grad( const Vector& P ) const
{
   Vector axisNorm = axis.unitvector();
	float sina = (std::sin(angle*M_PI/180.0));
	float cosa = (std::cos(angle*M_PI/180.0));

   Vector X = P*cosa + axisNorm *(axisNorm*P)*(1.0-cosa)+ (axisNorm^P)*sina;
   return  elem->grad(P); 
}

UnionVolume::UnionVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float UnionVolume::eval( const Vector& P ) const
{
   if(e1->eval(P) > e2->eval(P)) return e1->eval(P);
   else return e2->eval(P);
}
const Vector UnionVolume::grad( const Vector& P ) const
{
   if(e1->eval(P) > e2->eval(P)) return e1->grad(P);
   else return e2->grad(P);
}

IntersectionVolume::IntersectionVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float IntersectionVolume::eval( const Vector& P ) const
{
   if(e1->eval(P) < e2->eval(P)) return e1->eval(P);
   else return e2->eval(P);
}
const Vector IntersectionVolume::grad( const Vector& P ) const
{
   if(e1->eval(P) < e2->eval(P)) return e1->grad(P);
   else return e2->grad(P);

}

CutoutVolume::CutoutVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float CutoutVolume::eval( const Vector& P ) const
{
   if(e1->eval(P) < -e2->eval(P)) return e1->eval(P);
   else return -e2->eval(P);
}

const Vector CutoutVolume::grad( const Vector& P ) const
{
   if(e1->eval(P) < -e2->eval(P)) return e1->grad(P);
   else return -e2->grad(P);
}


MaskVolume::MaskVolume(const ScalarField& v) :
   elem(v)
{}

const float MaskVolume::eval(const Vector& P) const
{
   if (elem->eval(P) <= 0) return 0;
   else return 1;
}
const Vector MaskVolume::grad(const Vector& P) const 
{
   return elem->grad(P);
}

ClampVolume::ClampVolume(const ScalarField& v, float minv, float maxv) :
   elem(v),
   vmin(minv),
   vmax(maxv)
{}

const float ClampVolume::eval(const Vector& P) const 
{
   float eval = elem->eval(P);
   if (eval >= vmin && eval <= vmax) return eval;
   else if (eval < vmin ) return vmin;
   else if (eval > vmax) return vmax;
}
const Vector ClampVolume::grad(const Vector& P) const
{
   return elem->grad(P);
}

GriddedGridVolume::GriddedGridVolume(const ScalarGrid& g) :scgrid(g)
{ }

const float GriddedGridVolume::eval(const Vector& P) const
{
   return scgrid->eval(P);
}

GriddedSGridVolume::GriddedSGridVolume(const SScalarGrid& g) :scgrid(g)
{ }

const float GriddedSGridVolume::eval(const Vector& P) const
{
   return scgrid->eval(P);
}

NoiseVolume::NoiseVolume(const _Noise& n) :noise(n)
{ }

const float NoiseVolume::eval(const Vector& P) const
{
   return noise->eval(P);
}

SDFVolume::SDFVolume(const ScalarField& e, int n) :
  elem(e),
  N(n)
{ }

const float SDFVolume::eval(const Vector& P) const
{
   float f_x = elem->eval(P);
   return (identity()->eval(P) - iteratedNPT(elem, N)->eval(P)).magnitude() * (f_x/std::abs(f_x));

}

PyroclasticSphere::PyroclasticSphere(const Vector& cen, const float r, const float amp, const float gam, const _Noise& n ) :
  center(cen),
  rad(r),
  amplitude(amp),
  gamma(gam),
  noise(n)
{ }

const float PyroclasticSphere::eval(const Vector& P) const
{
   Vector x_closest = rad * (( P - center) / (P- center).magnitude());
   float f_x = rad - (P-center).magnitude();
   return f_x + amplitude * std::pow(std::abs(noise->eval(x_closest)),gamma);
}



PyroclasticVolume::PyroclasticVolume(const ScalarField& e, const float amp, const float gam, const _Noise& n ) :
  elem(e),
  amplitude(amp),
  gamma(gam),
  noise(n)
{   
   fspn = SFNoise(noise);
   Xnpt = iteratedNPT(elem,2);
   warpV = warp(fspn, Xnpt);
    }

const float PyroclasticVolume::eval(const Vector& P) const
{
   float f_x = elem->eval(P);

   //return f_x + amplitude * std::pow(std::abs(noise->eval(iteratedNPT(elem, 4)->eval(P))),gamma);
   return f_x + amplitude * std::pow(std::abs(warpV->eval(P)),gamma);

}

PyroclasticTerrain::PyroclasticTerrain(const Vector& xp, const Vector& norm, const float ampP, const float ampN, const float gamP, 
                                 const float gamN, const float trans, const _Noise& n ) :

  Xp(xp),
  normal(norm),
  posAmp(ampP),
  negAmp(ampN),
  posGam(gamP),
  negGam(gamN),
  transition(trans),
  noise(n)
{ }

const float PyroclasticTerrain::eval(const Vector& P) const
{
   float f_p = -normal * (P - Xp);
   Vector x_closest = P + f_p*normal;
   float fspn = noise->eval(x_closest);
   float gamma;
   if(fspn > transition) gamma = posGam;
   else if( fspn > 0 && fspn < transition) gamma = posGam * (fspn/transition) + 1 - (fspn/transition);
   else if( -transition < fspn && fspn < 0) gamma  = negGam * (std::abs(fspn)/transition) + 1 - (std::abs(fspn)/transition);
   else gamma = negGam;
   if(fspn >= 0)
   {
      //return f_p + posAmp * std::pow(std::abs(fspn), gamma);
      return f_p + posAmp * std::pow(fspn, gamma );
   }
   else
   {
      //return f_p + -negAmp * std::pow(std::abs(fspn),gamma );
      return f_p + -negAmp * std::pow(std::abs(fspn),gamma );
   }
}

WarpVolume::WarpVolume(const ScalarField& sf, VectorField&  map ) :
  elem(sf),
  mapX(map)
{ }

const float WarpVolume::eval(const Vector& P) const
{
      Vector X = mapX->eval(P);
      return elem->eval(X);
}
const Vector WarpVolume::grad( const Vector& P ) const
{
   Vector X = mapX->eval(P);
   Matrix M = mapX->grad(P);
   return M * (elem->grad(X)); 
}
 MultiplyVectorVolume::MultiplyVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}

    const float MultiplyVectorVolume::eval( const Vector& P ) const
    {
      Vector e1 = elem1->eval(P);
      Vector e2 = elem2->eval(P);
      
       return  e1 * e2;
    }


   //  const Vector MultiplyVectorVolume::grad( const Vector& P ) const
   //  {
   //     return  elem1->grad(P) * elem2->grad(P);
   //  }

 DivergenceVVolume::DivergenceVVolume( const SVectorGrid& u ) :
      U(u)
    {}

    const float DivergenceVVolume::eval( const Vector& P ) const
    {
      int i,j,k;
      U->getGridIndex(P, i, j ,k);
      float ux = ( (U->get(i+1, j, k) - U->get(i-1, j, k) ) * Vector(1,0,0) ) / ( 2 * U->dx());
      float uy = ( (U->get(i, j+1, k) - U->get(i, j-1, k) ) * Vector(0,1,0) ) / ( 2 * U->dx()); //dx, must be square gird
      float uz = ( (U->get(i, j, k+1) - U->get(i, j, k-1) ) * Vector(0,0,1) ) / ( 2 * U->dx());
      float divU = ux+ uy + uz;
       return  ((U->dx()*U->dx())/6) * divU;
    }
