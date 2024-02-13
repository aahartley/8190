#include "FullGrids.h"

using namespace lux;

void RectangularGrid::init(const Vector&llc, const Vector& urc, const Vector& dims)
{
    LLC = llc;
    URC = urc;
    dX = dims.X();
    dY = dims.Y();
    dZ = dims.Z();
    Nx = std::ceil(std::fabs(urc.X() - llc.X()) / dX);
    Ny = std::ceil(std::fabs(urc.Y() - llc.Y()) / dY);
    Nz = std::ceil(std::fabs(urc.Z() - llc.Z()) / dZ);
    URC = LLC + Vector(Nx*dX, Ny*dY, Nz*dZ);
}
const bool RectangularGrid::in_grid(const Vector& P) const
{
    if(!(P.X() >= LLC.X() && P.X() <= URC.X())) return false;
    if(!(P.Y() >= LLC.Y() && P.Y() <= URC.Y())) return false;
    if(!(P.Z() >= LLC.Z() && P.Z() <= URC.Z())) return false;


    return true;
}
const bool RectangularGrid::in_grid(int i, int j, int k) const
{
    if(i < 0 || i >= Nx) return false;
    if(j < 0 || j >= Ny) return false;
    if(k < 0 || k >= Nz) return false;
    
    return true;
}

const void RectangularGrid::getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
{
    ix = std::floor( (P.X() - LLC.X()) / dX );
    iy = std::floor( (P.Y() - LLC.Y()) / dY );
    iz = std::floor( (P.Z() - LLC.Z()) / dZ );
}

const Vector RectangularGrid::evalP(int i, int j, int k) const 
{
    return Vector(LLC.X()+(i*dX), LLC.Y()+(j*dY), LLC.Z()+(k*dZ));
}

GridBox::GridBox(RectangularGrid* f): std::shared_ptr<RectangularGrid>( f ) {}




ScalarGrid::ScalarGrid(FullGrid<float>* fg) :  std::shared_ptr<FullGrid<float> >( fg ) {}
