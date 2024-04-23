
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

VectorField grad( const ScalarField& e ) {return  VectorField ( new GradientVectorVolume(e));}


ScalarField ScalarField::operator+( const ScalarField& e2 ) { return add(*this,e2); }
ColorField ColorField::operator+( const ColorField& e2 ) { return add(*this,e2); }
VectorField VectorField::operator+( const VectorField& e2 ) { return add(*this,e2); }
MatrixField MatrixField::operator+( const MatrixField& e2 ) { return MatrixField( new AddMatrixVolume(*this,e2) ); }

ScalarField ScalarField::operator-( ) { return negate(*this); }
VectorField VectorField::operator-(const VectorField& v2 ) { return subtract(*this, v2); }

VectorField VectorField::operator/(const ScalarField& v2 ) { return divide(*this, v2); }


ScalarField ScalarField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }
VectorField VectorField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }
ScalarField VectorField::operator*( const VectorField& e2 ) { return multiply(*this,e2); }


ColorField ColorField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }

ScalarField constant( const float v ) { return ScalarField( new ConstantVolume(v) ); }
VectorField constant( const Vector& v ) { return VectorField(new ConstantVectorVolume(v)); }
ColorField  constant( const Color& v ) { return ColorField(new ConstantColor(v)); }
MatrixField constant( const Matrix& v ) { return MatrixField( new ConstantMatrixVolume(v) ); }




ScalarField exp( const ScalarField& v ) { return ScalarField(new ExpVolume(v)); }
MatrixField exp( const MatrixField& m ) { return MatrixField( new ExpMatrixVolume( m ) ); }



ScalarField add( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new AddVolume(v1,v2));}
VectorField add( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new AddVectorVolume(v1,v2)); }
ColorField  add( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new AddColor(v1,v2)); }
MatrixField add( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new AddMatrixVolume( v1, v2 )); }

VectorField subtract( const VectorField&  v1, const VectorField& v2 ){ return VectorField(new SubtractVectorVolume(v1, v2));}


ScalarField negate( const ScalarField&  v1) { return ScalarField(new NegateVolume(v1));}


ScalarField multiply( const ScalarField&  v1, const ScalarField& v2 )   { return ScalarField(new MultiplyVolume(v1,v2)); }
VectorField multiply( const VectorField& v, const ScalarField& u ) { return VectorField(new MultiplyVectorScalarVolume(v, u));}
ScalarField multiply( const VectorField& v, const VectorField& u ) { return ScalarField(new MultiplyVectorVolume(v, u));}
ColorField  multiply( const ColorField&  v1, const ScalarField& v2 )   { return ColorField(new MultiplyColor(v1,v2)); }

VectorField divide( const VectorField& v, const ScalarField& u ){ return VectorField(new DivideVectorVolume(v, u));}


ScalarField Sphere( const Vector& cen, const float rad ) { return ScalarField( new SphereVolume(cen, rad) ); }
ScalarField Ellipse( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) { return ScalarField( new EllipseVolume(cen, radMaj, radMin, norm) ); }
ScalarField Box( const Vector& cen, const float rad, const float q ) { return ScalarField( new BoxVolume(cen, rad, q) ); }
ScalarField Cone( const Vector& cen, const Vector& norm, const float h, const float theta) { return ScalarField( new ConeVolume(cen, norm, h, theta) ); }
ScalarField Plane( const Vector& cen, const Vector& norm, const Vector& pos) { return ScalarField( new PlaneVolume(cen, norm, pos) ); }
ScalarField Torus( const Vector& cen, const float radMaj, const float radMin, const Vector& norm ) { return ScalarField( new TorusVolume(cen, radMaj, radMin, norm) ); }
ScalarField SteinerPatch( const Vector& cen) { return ScalarField( new SteinerPatchVolume(cen) ); } 
ScalarField Icosahedron( ) { return ScalarField( new IcosahedronVolume() ); }
ScalarField Cylinder( const Vector& cen, const Vector& norm, const float radius) { return ScalarField( new CylinderVolume(cen, norm, radius) ); }
ScalarField SFNoise(const _Noise& noise){ return ScalarField( new NoiseVolume(noise));}
ScalarField PyroSphere(const Vector& cen, const float rad, const float amp, const float gam, const _Noise& noise){ return ScalarField( new PyroclasticSphere(cen, rad, amp, gam,noise));}
ScalarField PyroTerrain(const Vector& xp, const Vector& norm, const float ampP, const float ampN, const float gamP, const float gamN, const float trans, const _Noise& n){ return ScalarField( new PyroclasticTerrain(xp, norm, ampP, ampN, gamP,gamN,trans, n));}
ScalarField PyroVolume(const ScalarField& e, const float amp, const float gam, const _Noise& noise, int n ){return ScalarField ( new PyroclasticVolume(e, amp, gam, noise, n));}

VectorField VFNoise(const _Noise& noise){ return VectorField( new NoiseVectorVolume(noise));}


ScalarField Union( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new UnionVolume(v1, v2)) ;}
ScalarField Intersection( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new IntersectionVolume(v1, v2)) ;}
ScalarField Cutout( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new CutoutVolume(v1, v2)) ;}

ScalarField translate( const ScalarField& v, const Vector& s ) { return ScalarField( new TranslateVolume(v, s)); }
ScalarField scale( const ScalarField& v, const Vector& s ) { return ScalarField( new ScaleVolume(v, s)); }
ScalarField rotate( const ScalarField& v, const Vector& s , float a) { return ScalarField( new RotateVolume(v, s, a)); }



ScalarField mask( const ScalarField& v ) { return ScalarField( new MaskVolume(v));}
ScalarField clamp( const ScalarField& v, float minv, float maxv ) { return ScalarField ( new ClampVolume(v, minv, maxv));}

ScalarField warp( const ScalarField& v, VectorField& map ){return ScalarField( new WarpVolume(v, map));}
VectorField warp( const VectorField& v, VectorField& map ){return VectorField( new WarpVectorVolume(v, map));}

VectorField identity(){ return VectorField( new IdentityVectorVolume());}


VectorField iteratedNPT(const ScalarField& f, int N)
{
    VectorField gradient = grad(f);
    ScalarField func = f;
    VectorField NPT = identity() - ((gradient*func)/(gradient*gradient));
    VectorField result = identity();
    for(int i=0;i<N;i++)
    {
        result = warp( result, NPT );
    }
    return result;
}

ScalarField divergenceV ( const SVectorGrid& u) { return ScalarField(new DivergenceVVolume(u));}

void GaussDivFree ( const GridBox& gridB, ScalarField& p_field, VectorField& U_field, const int iter)
{

    SVectorGrid U_grid = makeSGrid(gridB, Vector(0,0,0));
    stamp(U_grid, U_field, 1);

    SScalarGrid p_grid = makeSGrid(gridB, 0);
    stamp(p_grid, p_field, 1);

    ScalarField d_field = divergenceV(U_grid);
    SScalarGrid d_grid = makeSGrid(gridB, 0);
    stamp(d_grid, d_field, 1);

    ProgressMeter pm(1, "div");
    int count = 0;
    while(count < iter)
    {
        for(int j = 0; j < gridB->ny(); j++)
        {
            for(int i = 0; i < gridB->nx(); i++)
            {
                for(int k = 0; k < gridB->nz(); k++)
                {
                    // float p_avg = ( p_grid->eval(p_grid->evalP(i+1, j ,k)) + p_grid->eval(p_grid->evalP(i-1, j ,k)) +
                    //             p_grid->eval(p_grid->evalP(i, j+1 ,k)) + p_grid->eval(p_grid->evalP(i, j-1 ,k)) +
                    //             p_grid->eval(p_grid->evalP(i, j ,k+1)) + p_grid->eval(p_grid->evalP(i, j ,k-1)) ) / 6;
                    float p_avg = ( p_grid->get(i+1, j ,k) + p_grid->get(i-1, j ,k) +
                                    p_grid->get(i, j+1 ,k) + p_grid->get(i, j-1 ,k) +
                                    p_grid->get(i, j ,k+1) + p_grid->get(i, j ,k-1) ) / 6.0;
                    p_grid->set(i, j ,k, p_avg - d_grid->get(i, j, k));
                }
            }
        }
        count++;
    }

    p_field = gridded(p_grid);
    VectorField U = gridded(U_grid);

    U_field = U - grad(p_field);

}

ScalarField advect( const ScalarField& v, const VectorField& u, const float delt ){ return ScalarField(new AdvectVolume(v, u, delt));}
VectorField advect( const VectorField& v, const VectorField& u, const float delt ){ return VectorField( new AdvectVectorVolume(v, u, delt));}

}
