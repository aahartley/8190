#ifndef MODELS_H
#define MODELS_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "Volume.h"

namespace lux
{



class Models
{
  public:
    Models(){ scalar_volumes_unioned = constant(-1000); colorfield = constant(Color(0,0,0,0));}
    ~Models(){}

    void addScalarModel(ScalarField& model, Color& color);
    void addHumanoid();
    ScalarField& getMaskedDensityField() { return density = mask(scalar_volumes_unioned); }
    ScalarField& getClampedDensityField(float min, float max) { return density = clamp(scalar_volumes_unioned, min, max); }
    ColorField& getColorField() { return colorfield; }

  private:
    ScalarField scalar_volumes_unioned;
    ColorField colorfield;
    ScalarField density;

};




}



#endif