
#include <iostream>
#include "Vector.h"
#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "Volume.h"

using namespace lux;

int main(int argc, char** argv)
{
    std::cout << "test" << std::endl;
    Vector p(0,0,0);
    ScalarField v (new ConstantVolume(1));
    ScalarField v2 = ScalarField(new ConstantVolume(2));
    ScalarField v3 = v + v2;
    std::cout << v3->eval(Vector(0,0,0)) << '\n';
    return 0;
}