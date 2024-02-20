
#include "ImplicitColors.h"

using namespace lux;


GriddedSGridColor::GriddedSGridColor(const ColorGrid& g) :scgrid(g)
{ }

const Color GriddedSGridColor::eval(const Vector& P) const
{
   return scgrid->eval(P);
}