#ifndef GRIDS_H
#define GRIDS_H
#include "ImplicitVolumeShapes.h"
#include "FullGrids.h"

namespace lux
{


GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx );

ScalarGrid makeGrid( const GridBox& rg, float defValue );


void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples );

ScalarField gridded( const ScalarGrid& g );







}

#endif