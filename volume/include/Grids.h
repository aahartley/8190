#ifndef GRIDS_H
#define GRIDS_H
#include "ImplicitVolumeShapes.h"
#include "FullGrids.h"
#include "ImplicitColors.h"
#include "ProgressMeter.h"
#include "Noise.h"

namespace lux
{


GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx );

ScalarGrid makeGrid( const GridBox& rg, float defValue );

ColorGrid makeGrid(const GridBox& rg, Color defValue);


void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples );
void stamp( ColorGrid& grid, const ColorField& field, const int nbsamples );

void stampNoise(ScalarGrid& grid, const _Noise& noise);
void stampWisp(ScalarGrid& grid, const Vector& P_4, const float den);


ScalarField gridded( const ScalarGrid& g );
ColorField gridded( const ColorGrid& g );







}

#endif