#ifndef MODELS_H
#define MODELS_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "FullGrids.h"
#include "Grids.h"

namespace lux
{



class Models
{
  public:
    Models(){scalar_volumes_unioned = constant(-1000);  colorfield = constant(Color(0,0,0,0));}
    ~Models(){}

    void addScalarModel(ScalarField& model, Color& color);
    void addHumanoid();
    ScalarField& getMaskedDensityField();
    ScalarField& getClampedDensityField(float min, float max) ;
    ColorField& getColorField() { return colorfield; }

  private:
    ScalarField scalar_volumes_unioned;
    ColorField colorfield;
    ScalarField density;

};




}



#endif