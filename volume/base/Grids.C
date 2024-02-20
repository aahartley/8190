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

ScalarField gridded( const ScalarGrid& g )
{
    return ScalarField(new GriddedSGridVolume(g));
}


ColorField gridded( const ColorGrid& g )
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

}