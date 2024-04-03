#ifndef SPARSEGRIDS_H
#define SPARSEGRIDS_H
#include "Volume.h"
#include "RectangularGrid.h"
// #include "Vector.h"
// #include "Matrix.h"
// #include "Color.h"
// #include <memory>

namespace lux
{



template< typename T >
class SparseGrid : public RectangularGrid
{
  public:
    SparseGrid(){}
    SparseGrid(T defValue, int part){ default_value = defValue; partition_size=part; totalSize =0; nbOccupied=0;}
    ~SparseGrid()
    {

      deallocate();
    }

    void deallocate()
    {
      if(data == NULL) return;

      for(long i = 0; i < totalSize; i++) { delete [] data[i];}
      delete[] data;
      data = NULL;
      nbOccupied = 0;
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
      Mx = (Nx/partition_size);
      My = (Ny/partition_size);
      Mz = (Nz/partition_size);
      if(Mx * partition_size < Nx) Mx++;
      if(My * partition_size < Ny) My++;
      if(Mz * partition_size < Nz) Mz++;
      Nx = Mx * partition_size;
      Ny = My * partition_size;
      Nz = Mz * partition_size;
      totalSize = Mx * My * Mz;
      //URC = LLC + Vector(Nx*dX, Ny*dY, Nz*dZ);
      //RectangularGrid::init(LLC,URC,dims);

      data = new T*[totalSize];

      for(int i = 0; i < totalSize; i++)
      {
        data[i] = NULL;
      }
    }

    void setDefVal( const T& def ) { default_value = def; }
    const T& getDefVal() const { return default_value; }

    void set( int i,int j,int k, const T& val )
    {
        long myIndex = sindex(i, j, k);  
        if(data[myIndex] == NULL)
        {
          data[myIndex] = new T[partition_size * partition_size * partition_size];
          for(int i = 0; i < partition_size * partition_size * partition_size; i++)
          {
              data[myIndex][i] = default_value;
          }
          nbOccupied++;
        }
        int ii = i % partition_size;
        int jj = j % partition_size;
        int kk = k % partition_size;
        int partitionIndex = ii + partition_size*( jj + partition_size * kk );
        data[myIndex][partitionIndex] = val;
    }   

    const T& get(int i,int j,int k) const
    {
       if(!in_grid(i,j,k)){ return(default_value); }
       int myIndex = sindex(i, j, k);
       if(data[myIndex] == NULL)
       {  
          return(default_value); 
       }
       else
       {
          int ii = i % partition_size;
          int jj = j % partition_size;
          int kk = k % partition_size;
          int partitionIndex = ii + partition_size*( jj + partition_size * kk );
          return data[myIndex][partitionIndex];
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

    const long sindex(int i,int j,int k) const
    {
       int ii = i/partition_size;
       int jj = j/partition_size;
       int kk = k/partition_size;
       return (long)ii + (long)Mx*( (long)jj + (long)My*(long)kk );
    }

 private:
    T** data;
    T default_value;
    int partition_size;
    int Mx, My, Mz;
    long totalSize;
    long nbOccupied;
};


 typedef std::shared_ptr<SparseGrid<float> > SScalarGridBase;
 typedef std::shared_ptr<SparseGrid<Vector> > SVectorGridBase;
 typedef std::shared_ptr<SparseGrid<Color> > SColorGridBase;
 typedef std::shared_ptr<SparseGrid<Matrix> > SMatrixGridBase;


  
 class SScalarGrid : public SScalarGridBase
 {
   public:
  
     SScalarGrid( SparseGrid<float>* f );

    ~SScalarGrid(){}
  
  
     SScalarGrid operator+=( const ScalarField& e2 );
     SScalarGrid operator-=( const ScalarField& e2 );
     SScalarGrid operator-();
     SScalarGrid operator*=( const ScalarField& e2 );
     SScalarGrid operator/=( const ScalarField& e2 );
 };
  
  
 class SVectorGrid : public SVectorGridBase
 {
   public:
  
     SVectorGrid( SparseGrid<Vector>* f );
    ~SVectorGrid(){}
  
  
     SVectorGrid operator+=( const VectorField& e2 );
     SVectorGrid operator-=( const VectorField& e2 );
     SVectorGrid operator-();
     SVectorGrid operator*=( const ScalarField& e2 );
     SVectorGrid operator/=( const ScalarField& e2 );
  
  
 };
  
  
  
  
  
 class SColorGrid : public SColorGridBase
 {
   public:
  
     SColorGrid( SparseGrid<Color>* f );
    ~SColorGrid(){}
  
  
     SColorGrid operator+=( const ColorField& e2 );
     SColorGrid operator-=( const ColorField& e2 );
     SColorGrid operator-();
     SColorGrid operator*=( const ScalarField& e2 );
     SColorGrid operator/=( const ScalarField& e2 );
  
  
 };
  
  
  
 class SMatrixGrid : public SMatrixGridBase
 {
   public:
  
     SMatrixGrid( SparseGrid<Matrix>* f );
    ~SMatrixGrid(){}
  
  
     SMatrixGrid operator+=( const MatrixField& e2 );
     SMatrixGrid operator-=( const MatrixField& e2 );
     SMatrixGrid operator-();
     SMatrixGrid operator*=( const ScalarField& e2 );
     SMatrixGrid operator/=( const ScalarField& e2 );
  
  
 };
  




}


#endif