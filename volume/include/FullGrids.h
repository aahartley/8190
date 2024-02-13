#ifndef FULLGRIDS_H
#define FULLGRIDS_H
#include "Volume.h"

// #include "Vector.h"
// #include "Matrix.h"
// #include "Color.h"
// #include <memory>

namespace lux
{

class RectangularGrid
 {
   public:
  
     RectangularGrid() :
        Nx   (0),
        Ny   (0),
        Nz   (0),
        dX   (0),
        dY   (0),
        dZ   (0),
        LLC   (Vector(0,0,0)),
        URC   (Vector(0,0,0))

     {}
  
    ~RectangularGrid(){}

    void init( const Vector&llc, const Vector& urc, const Vector& dims );

    
     const int& nx() const { return Nx; }
     const int& ny() const { return Ny; }
     const int& nz() const { return Nz; }
  
     const Vector& llc() const { return LLC; }
     const Vector& urc() const { return URC; }
  
     const float& dx() const { return dX; }
     const float& dy() const { return dY; }
     const float& dz() const { return dZ; }

     const Vector evalP( int i, int j, int k ) const;
  
     const bool in_grid( int i, int j, int k ) const;
     const bool in_grid( const Vector& P ) const;

  
     const bool getBox( const Vector& lower, const Vector& upper, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const;
     const void getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const;
  
   protected:
  
     int Nx, Ny, Nz;
     float dX, dY, dZ;
     Vector LLC, URC;
  
  
 };

typedef std::shared_ptr<RectangularGrid>  GridBoxBase;
  
  
class GridBox : public GridBoxBase
{
  public:
    GridBox( RectangularGrid* f );
    ~GridBox(){}
  
  
     GridBox operator+=( const GridBox& e2 ); // Union
  
     int index( const Vector& P ) const;
     char * __str__();
};
  


template< typename T >
class FullGrid : public RectangularGrid
{
  public:
    FullGrid(){}
    FullGrid(T defValue){ default_value = defValue; }
    ~FullGrid(){delete data;}

    void init(const Vector&llc, const Vector& urc, const Vector& dims)
    {
      LLC = llc;
      URC = urc;
      dX = dims.X();
      dY = dims.Y();
      dZ = dims.Z();
      Nx = std::ceil(std::fabs(urc.X() - llc.X()) / dX);
      Ny = std::ceil(std::fabs(urc.Y() - llc.Y()) / dY);
      Nz = std::ceil(std::fabs(urc.Z() - llc.Z()) / dZ);

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
      T accum = default_value * 0.0;
      accum += get(i, j, k) * (1 - ((P.X()-i)/dX)) * (1 - ((P.Y()-j)/dY )) * (1 - ((P.Z()-k)/dZ));
      accum += get(i+1, j, k) * (((P.X()-i)/dX)) * (1 - ((P.Y()-j)/dY )) * (1 - ((P.Z()-k)/dZ));
      accum += get(i, j+1, k) * (1 - ((P.X()-i)/dX)) * (((P.Y()-j)/dY )) * (1 - ((P.Z()-k)/dZ));
      accum += get(i, j, k+1) * (1 - ((P.X()-i)/dX)) * (1 - ((P.Y()-j)/dY )) * (((P.Z()-k)/dZ));
      accum += get(i+1, j+1, k) * (((P.X()-i)/dX)) * (((P.Y()-j)/dY )) * (1 - ((P.Z()-k)/dZ));
      accum += get(i+1, j, k+1) * (((P.X()-i)/dX)) * (1 - ((P.Y()-j)/dY )) * (((P.Z()-k)/dZ));
      accum += get(i, j+1, k+1) * (1 - ((P.X()-i)/dX)) * (((P.Y()-j)/dY )) * (((P.Z()-k)/dZ));
      accum += get(i+1, j+1, k+1) * (((P.X()-i)/dX)) * (((P.Y()-j)/dY )) * (((P.Z()-k)/dZ));

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