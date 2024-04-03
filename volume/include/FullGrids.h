#ifndef FULLGRIDS_H
#define FULLGRIDS_H
#include "Volume.h"
#include "RectangularGrid.h"
// #include "Vector.h"
// #include "Matrix.h"
// #include "Color.h"
// #include <memory>

namespace lux
{



template< typename T >
class FullGrid : public RectangularGrid
{
  public:
    FullGrid(){}
    FullGrid(T defValue){ default_value = defValue; }
    ~FullGrid()
    {

      if(data == NULL) return;
      delete[] data;
    }

    void init(const Vector&llc, const Vector& urc, const Vector& dims)
    {
      LLC = llc;
      URC = urc;
      dX = dims.X();
      dY = dims.Y();
      dZ = dims.Z();
      Nx = std::ceil(std::fabs(urc.X() - llc.X()) / dX)+1;
      Ny = std::ceil(std::fabs(urc.Y() - llc.Y()) / dY)+1;
      Nz = std::ceil(std::fabs(urc.Z() - llc.Z()) / dZ)+1;
      //URC = LLC + Vector(Nx*dX, Ny*dY, Nz*dZ);
      //RectangularGrid::init(LLC,URC,dims);

      data = new T[Nx*Ny*Nz];

      for(int j = 0; j < Ny; j++)
      {
        for(int i = 0; i < Nx; i++)
        {
          for(int k = 0; k < Nz; k++)
          {
            data[ i + Nx*(j + Ny*k) ] = default_value;
          }
        }
      }
    }

    void setDefVal( const T& def ) { default_value = def; }
    const T& getDefVal() const { return default_value; }

    void set(int i, int j, int k, T value)
    {
      if( in_grid(i,j,k) )
      {
        data[ i + Nx*(j + Ny*k) ] = value;
      }
    }   

    T get(int i, int j, int k) const
    {
      if( in_grid(i,j,k) )
      {
        return data[ i + Nx*(j + Ny*k) ];
      }
      else
      {
        return default_value;
      }   
    }


    const T eval( const Vector& P ) const
    {
      if(!in_grid(P)) return default_value;
      int i, j, k;
      getGridIndex(P, i, j, k);
      float weightX = (1 -((P.X() - (i*dX+LLC.X()))/dX));
      float weightXX = std::abs(((P.X() - (i*dX+LLC.X()))/dX)); //return 1 for grid cell outside of grid
      float weightY = (1 -((P.Y() - (j*dY+LLC.Y()))/dY));
      float weightYY = std::abs(((P.Y() - (j*dY+LLC.Y()))/dY)); //return 1 for grid cell outside of grid
      float weightZ = (1 -((P.Z() - (k*dZ+LLC.Z()))/dZ));
      float weightZZ = std::abs(((P.Z() - (k*dZ+LLC.Z()))/dZ)); //return 1 for grid cell outside of grid
      T accum = default_value * 0.0;
      accum += get(i, j, k) *  weightX * weightY * weightZ;
      accum += get(i+1, j, k) * weightXX * weightY * weightZ;
      accum += get(i, j+1, k) * weightX * weightYY * weightZ;
      accum += get(i, j, k+1) * weightX * weightY * weightZZ;
      accum += get(i+1, j+1, k) * weightXX * weightYY * weightZ;
      accum += get(i+1, j, k+1) * weightXX * weightY * weightZZ;
      accum += get(i, j+1, k+1) * weightX * weightYY * weightZZ;
      accum += get(i+1, j+1, k+1) * weightXX * weightYY * weightZZ;

      return accum;
    }


 private:
    T* data;
    T default_value;
};


 typedef std::shared_ptr<FullGrid<float> > ScalarGridBase;
 typedef std::shared_ptr<FullGrid<Vector> > VectorGridBase;
 typedef std::shared_ptr<FullGrid<Color> > ColorGridBase;
 typedef std::shared_ptr<FullGrid<Matrix> > MatrixGridBase;


  
 class ScalarGrid : public ScalarGridBase
 {
   public:
  
     ScalarGrid( FullGrid<float>* f );

    ~ScalarGrid(){}
  
  
     ScalarGrid operator+=( const ScalarField& e2 );
     ScalarGrid operator-=( const ScalarField& e2 );
     ScalarGrid operator-();
     ScalarGrid operator*=( const ScalarField& e2 );
     ScalarGrid operator/=( const ScalarField& e2 );
 };
  
  
 class VectorGrid : public VectorGridBase
 {
   public:
  
     VectorGrid( FullGrid<Vector>* f );
    ~VectorGrid(){}
  
  
     VectorGrid operator+=( const VectorField& e2 );
     VectorGrid operator-=( const VectorField& e2 );
     VectorGrid operator-();
     VectorGrid operator*=( const ScalarField& e2 );
     VectorGrid operator/=( const ScalarField& e2 );
  
  
 };
  
  
  
  
  
 class ColorGrid : public ColorGridBase
 {
   public:
  
     ColorGrid( FullGrid<Color>* f );
    ~ColorGrid(){}
  
  
     ColorGrid operator+=( const ColorField& e2 );
     ColorGrid operator-=( const ColorField& e2 );
     ColorGrid operator-();
     ColorGrid operator*=( const ScalarField& e2 );
     ColorGrid operator/=( const ScalarField& e2 );
  
  
 };
  
  
  
 class MatrixGrid : public MatrixGridBase
 {
   public:
  
     MatrixGrid( FullGrid<Matrix>* f );
    ~MatrixGrid(){}
  
  
     MatrixGrid operator+=( const MatrixField& e2 );
     MatrixGrid operator-=( const MatrixField& e2 );
     MatrixGrid operator-();
     MatrixGrid operator*=( const ScalarField& e2 );
     MatrixGrid operator/=( const ScalarField& e2 );
  
  
 };
  




}


#endif