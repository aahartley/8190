#include "SparseGrids.h"

using namespace lux;


SScalarGrid::SScalarGrid(SparseGrid<float>* fg) :  std::shared_ptr<SparseGrid<float> >( fg ) {}
SVectorGrid::SVectorGrid(SparseGrid<Vector>* fg) :  std::shared_ptr<SparseGrid<Vector> >( fg ) {}

SColorGrid::SColorGrid(SparseGrid<Color>* fg) :  std::shared_ptr<SparseGrid<Color> >( fg ) {}
