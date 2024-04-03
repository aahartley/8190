#include "FullGrids.h"

using namespace lux;


ScalarGrid::ScalarGrid(FullGrid<float>* fg) :  std::shared_ptr<FullGrid<float> >( fg ) {}
ColorGrid::ColorGrid(FullGrid<Color>* fg) :  std::shared_ptr<FullGrid<Color> >( fg ) {}
