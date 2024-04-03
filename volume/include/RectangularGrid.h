#ifndef RECTANGULARGRID_H
#define RECTANGULARGRID_H

#include "Vector.h"
// #include "Matrix.h"
// #include "Color.h"
#include <memory>

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
     void fix_indices(int& i, int& j, int& k) const;
   
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
  
}

#endif