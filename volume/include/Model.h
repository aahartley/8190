#ifndef MODEL_H
#define MODEL_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "Volume.h"

namespace lux
{


class Model
{
  public:
    Model(){}
    virtual ~Model() {}

    virtual void display() {}

    ScalarField& getDensityField()  {return densityfield;}
    ColorField& getColorField()  {return colorfield;}

  protected:
    ScalarField densityfield;
    ColorField colorfield;
};


class Humanoid : public Model
{
  public:
    Humanoid();
    virtual ~Humanoid() {}

 
};


}



#endif