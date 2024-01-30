#include "Model.h"

using namespace lux;

void Humanoid::display()
{
    Vector p(0,0,0);
    ScalarField v (new ConstantVolume(1));
    ScalarField v2 = ScalarField(new ConstantVolume(2));
    ScalarField v3 = v + v2;
    
    std::cout << v3->eval(Vector(0,0,0)) << '\n';
    ScalarField e1(new EllipseVolume(p, 2, 3, Vector(0,1,0)));
    std::cout << e1->grad(Vector(0,0,0)).X() << '\n';
}