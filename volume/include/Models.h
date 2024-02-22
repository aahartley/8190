#ifndef MODELS_H
#define MODELS_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"

#include "FullGrids.h"
#include "Grids.h"
#include <iostream>
#include <cstdio>
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

    void createColorField(ScalarField& model, Color color);
    void createFinalUnion(ScalarField& model);
    ScalarField addHumanoid();
    ScalarField addOBJModel(std::string filepath, Vector llc, Vector urc, Vector dims);

    void scene2();

    void setGridBox(GridBox& gB ){gb = gB;}

    float calcDistance(Vector& x_ijk, Triangle& t);

    ScalarField getMaskedDensityField();
    ScalarField getClampedDensityField(float min, float max) ;

    ColorField getColorField() { return colorfield; }
    ScalarField getVolumesUnioned(){return scalar_volumes_unioned;}

    ScalarField getGriddedMaskedDensityField();
    ScalarField getGriddedClampedDensityField(float min, float max) ;

    ColorField getGriddedColorField();

  private:
    ScalarField scalar_volumes_unioned;
    ColorField colorfield;
    ScalarField density;
    GridBox gb = nullptr;
};




}



#endif