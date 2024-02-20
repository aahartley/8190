#ifndef MODELS_H
#define MODELS_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"

#include "FullGrids.h"
#include "Grids.h"
#include "objload.h" //https://github.com/GerhardR/objload
#include <iostream>
#include <algorithm>
#include "objloader.h"

namespace lux
{

// struct Triangle
// {
//   Triangle(){}
//   Triangle(Vector vv1, Vector vv2, Vector vv3) : v1(vv1), v2(vv2), v3(vv3)
//   {
//     e1 = v2 - v1;
//     e2 = v3 - v1;
//     e3 = v3 - v2;
//   }

//   Vector v1, v2, v3, e1, e2, e3;
// };

class Models
{
  public:
    Models(){scalar_volumes_unioned = constant(-1000);  colorfield = constant(Color(0,0,0,0)); gb = makeGridBox(Vector(0,0,0),Vector(0,0,0),Vector(0,0,0));}
    ~Models(){}

    void addScalarModel(ScalarField& model, Color color);
    void addHumanoid();
    void addOBJModel(std::string filepath);

    void setGridBox(GridBox& gB ){gb = gB;}

    float calcDistance(Vector& x_ijk, Triangle& t);

    ScalarField getMaskedDensityField();
    ScalarField getClampedDensityField(float min, float max) ;

    ColorField getColorField() { return colorfield; }
    ScalarField getVolumesUnioned(){return scalar_volumes_unioned;}

    ScalarField getGriddedMaskedDensityField();
    ScalarField getGriddedClampedDensityField(float min, float max) ;

    ColorField getGriddedColorField();
    ScalarField getGriddedVolumesUnioned();
  private:
    ScalarField scalar_volumes_unioned;
    ColorField colorfield;
    ScalarField density;
    GridBox gb = nullptr;
    obj::Model m;
};




}



#endif