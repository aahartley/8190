#ifndef GRIDS_H
#define GRIDS_H
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "RectangularGrid.h"
#include "FullGrids.h"
#include "SparseGrids.h"
#include "ImplicitColors.h"
#include "ProgressMeter.h"
#include "Noise.h"

namespace lux
{


GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx );


ScalarGrid makeGrid( const GridBox& rg, float defValue );
ColorGrid makeGrid(const GridBox& rg, Color defValue);

SScalarGrid makeSGrid( const GridBox& rg, float defValue );
SVectorGrid makeSGrid( const GridBox& rg, const Vector& defValue );
SColorGrid makeSGrid(const GridBox& rg, Color defValue);


void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples );
void stamp( ColorGrid& grid, const ColorField& field, const int nbsamples );

void stampNoise(ScalarGrid& grid, const _Noise& noise);
void stampWisp(ScalarGrid& grid, const Vector& P_4, const float den);

ScalarField gridded( const ScalarGrid& g );
ColorField gridded( const ColorGrid& g );



void stamp( SScalarGrid& grid, const ScalarField& field, const int nbsamples );
void stamp( SVectorGrid& grid, const VectorField& field, const int nbsamples );
void stamp( SColorGrid& grid, const ColorField& field, const int nbsamples );

void stampNoise(SScalarGrid& grid, const _Noise& noise);
void stampWisp(SScalarGrid& grid, const Vector& P_4, const float den);

ScalarField gridded( const SScalarGrid& g );
VectorField gridded( const SVectorGrid& g );
ColorField gridded( const SColorGrid& g );




}

#endif