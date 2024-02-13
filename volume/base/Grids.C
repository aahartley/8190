#include "Grids.h"

#include <iostream>


namespace lux
{

GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx )
{
    GridBox grid(new RectangularGrid());
    grid->init(llc, urc, dx);
    return grid;
}

ScalarGrid makeGrid( const GridBox& gb, float defValue )
{
    FullGrid<float>* fgp = new FullGrid<float>(defValue);
    fgp->init(gb->llc(), gb->urc(), Vector(gb->dx(),gb->dy(),gb->dz()) );
    return ScalarGrid (fgp);

}

ScalarField gridded( const ScalarGrid& g )
{
    return ScalarField(new GriddedSGridVolume(g));
}

void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples )
{
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
    std::cout << "stamping done\n";
}

}