#include "Grids.h"

#include <iostream>


namespace lux
{

GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx )
{
    RectangularGrid* rg = new RectangularGrid();
    rg->init(llc, urc, dx);
    return GridBox(rg);

}

ScalarGrid makeGrid( const GridBox& gb, float defValue )
{
    FullGrid<float>* fgp = new FullGrid<float>(defValue);
    fgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return ScalarGrid (fgp);

}

ColorGrid makeGrid( const GridBox& gb, Color defValue )
{
    FullGrid<Color>* fgp = new FullGrid<Color>(defValue);
    fgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return ColorGrid (fgp);

}

SScalarGrid makeSGrid( const GridBox& gb, float defValue )
{
    SparseGrid<float>* sgp = new SparseGrid<float>(defValue, 4);
    sgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return SScalarGrid (sgp);

}

SVectorGrid makeSGrid( const GridBox& gb, const Vector& defValue )
{
    SparseGrid<Vector>* sgp = new SparseGrid<Vector>(defValue, 4);
    sgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return SVectorGrid (sgp);

}

SColorGrid makeSGrid( const GridBox& gb, Color defValue )
{
    SparseGrid<Color>* sgp = new SparseGrid<Color>(defValue, 4);
    sgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return SColorGrid (sgp);

}

ScalarField gridded( const ScalarGrid& g )
{
    return ScalarField(new GriddedGridVolume(g));
}


ColorField gridded( const ColorGrid& g )
{
    return ColorField(new GriddedGridColor(g));
}

ScalarField gridded( const SScalarGrid& g )
{
    return ScalarField(new GriddedSGridVolume(g));
}
VectorField gridded( const SVectorGrid& g )
{
    return VectorField(new GriddedSGridVector(g));
}

ColorField gridded( const SColorGrid& g )
{
    return ColorField(new GriddedSGridColor(g));
}

void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples )
{
    ProgressMeter pm(1, "grid_scalar");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, field->eval(grid->evalP(i,j,k)));
                // if(field->eval(grid->evalP(i,j,k))!=0)
                //     std::cout << field->eval(grid->evalP(i,j,k)) << '\n';
            }
        }
    }
}

void stamp( ColorGrid& grid, const ColorField& field, const int nbsamples )
{
    ProgressMeter pm(1, "grid_color");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, field->eval(grid->evalP(i,j,k)));
                // if(field->eval(grid->evalP(i,j,k))!=0)
                //     std::cout << field->eval(grid->evalP(i,j,k)) << '\n';
            }
        }
    }
}

void stampNoise( ScalarGrid& grid, const _Noise& noise)
{
    //ProgressMeter pm(1, "grid_noise");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, std::max(noise->eval(grid->evalP(i,j,k)),grid->get(i, j , k)));
     
            }
        }
    }
}

void stampWisp( ScalarGrid& grid, const Vector& P_4, const float den)
{
    //ProgressMeter pm(1, "grid_wisp");
    int ix,iy,iz;
    grid->getGridIndex(P_4, ix, iy, iz);

    for(int i = ix; i <= ix+1; i++)
    {
        for(int j = iy; j <= iy+1; j++)
        {
            for(int k = iz; k <= iz+1; k++)
            {
                if(grid->in_grid(i,j,k))
                {
                    Vector gp=grid->evalP(i,j,k);
                    float weight = ( 1 - ((Vector(1,0,0) * (gp - P_4))/grid->dx()) ) * ( 1 - ((Vector(0,1,0) * (gp - P_4))/grid->dy()) ) *
                                ( 1 - ((Vector(0,0,1) * (gp - P_4))/grid->dz()) );
                    grid->set(i, j , k, grid->get(i,j,k) + den * weight);
                }

            }
        }
    }



}

void stamp( SScalarGrid& grid, const ScalarField& field, const int nbsamples )
{
    ProgressMeter pm(1, "grid_scalar");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, field->eval(grid->evalP(i,j,k)));
                // if(field->eval(grid->evalP(i,j,k))!=0)
                //     std::cout << field->eval(grid->evalP(i,j,k)) << '\n';
            }
        }
    }
}

void stamp( SVectorGrid& grid, const VectorField& field, const int nbsamples )
{
    ProgressMeter pm(1, "grid_vector");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, field->eval(grid->evalP(i,j,k)));
                // if(field->eval(grid->evalP(i,j,k))!=0)
                //     std::cout << field->eval(grid->evalP(i,j,k)) << '\n';
            }
        }
    }
}

void stamp( SColorGrid& grid, const ColorField& field, const int nbsamples )
{
    ProgressMeter pm(1, "grid_color");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, field->eval(grid->evalP(i,j,k)));
                // if(field->eval(grid->evalP(i,j,k))!=0)
                //     std::cout << field->eval(grid->evalP(i,j,k)) << '\n';
            }
        }
    }
}

void stampNoise( SScalarGrid& grid, const _Noise& noise)
{
    //ProgressMeter pm(1, "grid_noise");
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < grid->ny(); j++)
    {
        for(int i = 0; i < grid->nx(); i++)
        {
            for(int k = 0; k < grid->nz(); k++)
            {
                grid->set(i, j, k, std::max(noise->eval(grid->evalP(i,j,k)),grid->get(i, j , k)));
     
            }
        }
    }
}

void stampWisp( SScalarGrid& grid, const Vector& P_4, const float den)
{
    //ProgressMeter pm(1, "grid_wisp");
    int ix,iy,iz;
    grid->getGridIndex(P_4, ix, iy, iz);

    for(int i = ix; i <= ix+1; i++)
    {
        for(int j = iy; j <= iy+1; j++)
        {
            for(int k = iz; k <= iz+1; k++)
            {
                if(grid->in_grid(i,j,k))
                {
                    Vector gp=grid->evalP(i,j,k);
                    float weight = ( 1 - ((Vector(1,0,0) * (gp - P_4))/grid->dx()) ) * ( 1 - ((Vector(0,1,0) * (gp - P_4))/grid->dy()) ) *
                                ( 1 - ((Vector(0,0,1) * (gp - P_4))/grid->dz()) );
                    grid->set(i, j , k, grid->get(i,j,k) + den * weight);
                }

            }
        }
    }



}

}