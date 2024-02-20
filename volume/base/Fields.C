
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "ImplicitMatrixShapes.h"
#include "ImplicitColors.h"
#include "Fields.h"


using namespace std;
namespace lux
{


const float  evaluate( const ScalarField& v, const Vector& P ) { return v->eval(P); }
const Vector evaluate( const VectorField& v, const Vector& P ) { return v->eval(P); }
const Color  evaluate( const ColorField&  v, const Vector& P ) { return v->eval(P); }
const Matrix evaluate( const MatrixField& v, const Vector& P ) { return v->eval(P); }

ScalarField ScalarField::operator+( const ScalarField& e2 ) { return add(*this,e2); }
ColorField ColorField::operator+( const ColorField& e2 ) { return add(*this,e2); }
VectorField VectorField::operator+( const VectorField& e2 ) { return add(*this,e2); }
MatrixField MatrixField::operator+( const MatrixField& e2 ) { return MatrixField( new AddMatrixVolume(*this,e2) ); }

ScalarField ScalarField::operator-( ) { return negate(*this); }


ScalarField ScalarField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }
ColorField ColorField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }

ScalarField constant( const float v ) { return ScalarField( new ConstantVolume(v) ); }
//VectorField constant( const Vector& v ) { return VectorField(new ConstantVectorVolume(v)); }
ColorField  constant( const Color& v ) { return ColorField(new ConstantColor(v)); }
MatrixField constant( const Matrix& v ) { return MatrixField( new ConstantMatrixVolume(v) ); }




ScalarField exp( const ScalarField& v ) { return ScalarField(new ExpVolume(v)); }
MatrixField exp( const MatrixField& m ) { return MatrixField( new ExpMatrixVolume( m ) ); }



ScalarField add( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new AddVolume(v1,v2));}
VectorField add( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new AddVectorVolume(v1,v2)); }
ColorField  add( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new AddColor(v1,v2)); }
MatrixField add( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new AddMatrixVolume( v1, v2 )); }

ScalarField negate( const ScalarField&  v1) { return ScalarField(new NegateVolume(v1));}


ScalarField multiply( const ScalarField&  v1, const ScalarField& v2 )   { return ScalarField(new MultiplyVolume(v1,v2)); }
ColorField  multiply( const ColorField&  v1, const ScalarField& v2 )   { return ColorField(new MultiplyColor(v1,v2)); }

ScalarField Sphere( const Vector& cen, const float rad ) { return ScalarField( new SphereVolume(cen, rad) ); }
ScalarField Ellipse( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) { return ScalarField( new EllipseVolume(cen, radMaj, radMin, norm) ); }
ScalarField Box( const Vector& cen, const float rad, const float q ) { return ScalarField( new BoxVolume(cen, rad, q) ); }
ScalarField Cone( const Vector& cen, const Vector& norm, const float h, const float theta) { return ScalarField( new ConeVolume(cen, norm, h, theta) ); }
ScalarField Plane( const Vector& cen, const Vector& norm, const Vector& pos) { return ScalarField( new PlaneVolume(cen, norm, pos) ); }
ScalarField Torus( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) { return ScalarField( new TorusVolume(cen, radMaj, radMin, norm) ); }
ScalarField SteinerPatch( const Vector& cen) { return ScalarField( new SteinerPatchVolume(cen) ); } 
ScalarField Icosahedron( ) { return ScalarField( new IcosahedronVolume() ); }
ScalarField Cylinder( const Vector& cen, const Vector& norm, const float radius) { return ScalarField( new CylinderVolume(cen, norm, radius) ); }


ScalarField Union( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new UnionVolume(v1, v2)) ;}
ScalarField Intersection( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new IntersectionVolume(v1, v2)) ;}
ScalarField Cutout( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new CutoutVolume(v1, v2)) ;}

ScalarField translate( const ScalarField& v, const Vector& s ) { return ScalarField( new TranslateVolume(v, s)); }
ScalarField scale( const ScalarField& v, const Vector& s ) { return ScalarField( new ScaleVolume(v, s)); }
ScalarField rotate( const ScalarField& v, const Vector& s , float a) { return ScalarField( new RotateVolume(v, s, a)); }



ScalarField mask( const ScalarField& v ) { return ScalarField( new MaskVolume(v));}
ScalarField clamp( const ScalarField& v, float minv, float maxv ) { return ScalarField ( new ClampVolume(v, minv, maxv));}

}
