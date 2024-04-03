
#include "ImplicitColors.h"

using namespace lux;

GriddedGridColor::GriddedGridColor(const ColorGrid& g) :scgrid(g)
{ }

const Color GriddedGridColor::eval(const Vector& P) const
{
   return scgrid->eval(P);
}

GriddedSGridColor::GriddedSGridColor(const SColorGrid& g) :scgrid(g)
{ }

const Color GriddedSGridColor::eval(const Vector& P) const
{
   return scgrid->eval(P);
}