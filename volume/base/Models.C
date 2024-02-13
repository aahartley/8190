#include "Models.h"

using namespace lux;

void Models::addScalarModel(ScalarField& model, Color& color)
{
    scalar_volumes_unioned = Union(scalar_volumes_unioned, model);

    colorfield = colorfield *mask(-model) + constant(color)*mask(model); 

}

ScalarField& Models::getMaskedDensityField() 
{
    GridBox gb = makeGridBox(Vector(-10,-10,-10), Vector(10,10,10), Vector(0.1,0.1,0.1));
    ScalarGrid grid = makeGrid(gb, 0.f);
    density = mask(scalar_volumes_unioned);
    stamp(grid, density, 1);
    density = gridded(grid);
    return density;
}

ScalarField& Models::getClampedDensityField(float min, float max) 
{
    GridBox gb = makeGridBox(Vector(-15,-15,-15), Vector(15,15,15), Vector(0.1,0.1,0.1));
    ScalarGrid grid = makeGrid(gb, -1000.f);
    density = clamp(scalar_volumes_unioned, min, max);
    stamp(grid, density, 1);
    density = gridded(grid);
    return density;
}

void Models::addHumanoid()
{

    Color red( 1,0.2,0,1); Color green(0,1,0,1); Color black(0,0,0,1); Color blue(0,0,1,1); Color pink(1, 0.75, 0.8, 1); 
    Color white (1,1,1,1);

    //ScalarField c = constant(-1000);
    ScalarField e1 = Sphere(Vector(0,4,0), 2);
    ScalarField e2 = Sphere(Vector(0,2,0),2);
    ScalarField e3 = Sphere(Vector(0,4.7,0),2);
    e3 = scale(e3, Vector(1.3,1,2.5));
    e1 = Cutout(e1, e3);
    ScalarField e4 = Ellipse(Vector(0,-1.9,0), 2.5, 1.5, Vector(0,1,0));
    ScalarField e5 = Ellipse(Vector(-.5, -4, 0.5), .3, .6, Vector(0,1,0)); 
    ScalarField e6 = Ellipse(Vector(.5, -4, 0.5), .3, .6, Vector(0,1,0)); 
    ScalarField e7 = Ellipse(Vector(-2, -1, 1), 0.5, 1.3, Vector(0,1,0));
    ScalarField e8 = Ellipse(Vector(2, -1, 1), 0.5, 1.3, Vector(0,1,0));
    ScalarField e9 = Torus(Vector(0,1.8,1), 1, 0.5, Vector(0,1,0));
    ScalarField e10 = Torus(Vector(0,0.7,1), 1, 0.9, Vector(0,0,1));
    ScalarField e11 = Plane(Vector(0,0.7,0.6), Vector(0,1,0), Vector(0,0.7,0.6));
    e10 = Cutout(e10, e11);
    e10 = scale(e10, Vector(0.4,0.4,0.7));
    ScalarField e12 = Box(Vector(0,13.5, 1), 1, 4);
    e12 = scale(e12, Vector(0.55,0.3,0.7));
    ScalarField e13 = Icosahedron();
    e13 = translate(e13, Vector(18,4,0));
    e13 = scale(e13, Vector(0.1,0.1,0.1));
    ScalarField e14 = SteinerPatch(Vector(3,-2,1));
    e14 = rotate(e14, Vector(1,0,0), 45);
    ScalarField e15 = Cone(Vector(0,9,0.5), Vector(0,-1,0), 0.5, 45);
    e15 = scale(e15, Vector(1.5,0.5,1.5));
    ScalarField e16 = Intersection(e12, e15);
    e16 = translate(e16, Vector(0,-0.3,-0.75));
    ScalarField e17 = Icosahedron();
    e17 = translate(e17, Vector(-18,4,0));
    e17 = scale(e17, Vector(0.1,0.1,0.1));


    addScalarModel(e1, red);
    addScalarModel(e2, green);
    addScalarModel(e9, red);
    addScalarModel(e4, green);
    addScalarModel(e5, green);
    addScalarModel(e6, green);
    addScalarModel(e7, green);
    addScalarModel(e8, green);
    addScalarModel(e10, red);
    addScalarModel(e13, pink);
    addScalarModel(e14, white);
    addScalarModel(e16, red);
    addScalarModel(e17, pink);



    ////ScalarField u1 = Union(c, Union(e2,e9)); //head and eyes
    // ScalarField u1 = e2;
    // u1 = Union(u1, e9);
    // u1 = Union(u1, e1); //horns
    // u1 = Union(u1, e4); // torso
    // u1 = Union(u1, e5); // left foot
    // u1 = Union(u1, e6); // right foot
    // u1 = Union(u1, e7); // left arm
    // u1 = Union(u1, e8); // right arm
    // u1 = Union(u1, e10); // mouth
    // u1 = Union(u1, e13); //tattoo
    // u1 = Union(u1, e14); //stein
    // u1 = Union(u1, e16); // hat
    // u1 = Union(u1, e17); // tattoo

    // //ColorField e1Color = constant(Color(0,0,0,0)); //start with black
    // e1Color = e1Color *mask(-e1) + constant(red)*mask(e1); //horn color
    // e1Color = e1Color *mask(-e2) + constant(green)*mask(e2); //head color
    // e1Color = e1Color *mask(-e4) + constant(green)*mask(e4); //torso color
    // e1Color = e1Color *mask(-e5) + constant(green)*mask(e5); //left feet color
    // e1Color = e1Color *mask(-e6) + constant(green)*mask(e6); //right feet color
    // e1Color = e1Color *mask(-e7) + constant(green)*mask(e7); //left arm color
    // e1Color = e1Color *mask(-e8) + constant(green)*mask(e8); //right arm color
    // e1Color = e1Color *mask(-e9) + constant(red)*mask(e9); //eyes color
    // e1Color = e1Color *mask(-e10) + constant(red)*mask(e10); //mouth color
    // e1Color = e1Color *mask(-e13) + constant(pink)*mask(e13); //ico color
    // e1Color = e1Color *mask(-e14) + constant(Color(1,1,1,1))*mask(e14); //stei color
    // e1Color = e1Color *mask(-e16) + constant(red)*mask(e16); //hat color
    // e1Color = e1Color *mask(-e17) + constant(pink)*mask(e17); //ico color



    // //ScalarField density = mask(u1);
    // ScalarField density = clamp(u1, 0.0, 1.0);

    // densityfield = density;
    // colorfield = e1Color;
}