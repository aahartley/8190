#include "Models.h"

using namespace lux;

void Models::addScalarModel(ScalarField& model, Color color)
{
    scalar_volumes_unioned = Union(scalar_volumes_unioned, model);

    colorfield = colorfield *mask(-model) + constant(color)*mask(model); 

}

ScalarField Models::getMaskedDensityField() 
{

    return density = mask(scalar_volumes_unioned);
}

ScalarField Models::getClampedDensityField(float min, float max) 
{
    return density = clamp(scalar_volumes_unioned, min, max);

}


ScalarField Models::getGriddedMaskedDensityField() 
{

    ScalarGrid dgrid = makeGrid(gb, 0.f);
    density = mask(scalar_volumes_unioned);
    stamp(dgrid, density, 1);
    return gridded(dgrid);
    
}


ScalarField Models::getGriddedClampedDensityField(float min, float max) 
{
    ScalarGrid dgrid = makeGrid(gb, 0.f);
    density = clamp(getGriddedVolumesUnioned(), min, max);
    stamp(dgrid, density, 1);
    return gridded(dgrid);
}


ColorField Models::getGriddedColorField() 
{
    ColorGrid cgrid = makeGrid(gb, Color(0,0,0,0));
    stamp(cgrid, colorfield, 1);
    return gridded(cgrid);
}

ScalarField Models::getGriddedVolumesUnioned()
{
    ScalarGrid mgrid = makeGrid(gb, -10000.f);
    stamp(mgrid,scalar_volumes_unioned , 1);
    return gridded(mgrid);
}

float Models::calcDistance(Vector& x_ijk, Triangle& T)
{
    Vector D = T.v1 - x_ijk;
    float a = T.e1 * T.e1;
    float b = T.e1 * T.e2;
    float c = T.e2 * T.e2;
    float d = T.e1 * D;
    float e = T.e2 * D;
    float f = D * D;

    float s = b * e - c*d;
    float t = b*d - a*e;
    float det = a*c - b*b;
    if(s + t <= det)
    {
        if(s < 0)
        {
            if(t < 0)
            {
                //region 4
                if( d<0)
                {
                    t = 0;
                    if(-d >=a)
                    {
                        s = 1;
                    }
                    else 
                    {
                        s =-d/a;
                    }
                }
                else
                {
                    s =0;
                    if(e>=0)
                    {
                        t = 0;
                    }
                    else if( -e>=c)
                    {
                        t=1;
                    }
                    else
                    {
                        t = -e/c;
                    }
                }
            }
            else
            {
                //region 3
                s = 0;
                if(e >= 0)
                {
                    t =0;
                }
                else if(-e >= c)
                {
                    t =1;
                }
                else
                {
                    t = -e/c;
                }
            }
        }
        else if(t < 0)
        {
            //region 5
            t = 0;
            if(d >=0)
            {
                s =0;
            }
            else if( -d >= a)
            {
                s = 1;
            }
            else
            {
                s = -d/a;
            }
        }
        else
        {
            //region 0
            s /= det;
            t /= det;
        }
    }
    else
    {
        if(s < 0)
        {
            //region 2
            float tmp0 = b +d;
            float tmp1 = c+e;
            if(tmp1 > tmp0)
            {
                float numer = tmp1 - tmp0;
                float denom = a -2 * b+c;
                if(numer >= denom)
                {
                    s = 1;
                }
                else
                {
                    s = numer/denom;
                }
                t = 1-s;
            }
            else
            {
                s = 0;
                if(tmp1 <= 0)
                {
                    t = 1;
                }
                else if( e>=0)
                {
                    t = 0;
                }
                else
                {
                    t = -e/c;
                }
            }
        }
        else if( t < 0)
        {
            //region 6
            float tmp0 = b+e;
            float tmp1 = a+d;
            if(tmp1 > tmp0)
            {
                float numer = tmp1 - tmp0;
                float denom = a -2 * b+c;
                if( numer >= denom)
                {
                    t = 1;
                }
                else
                {
                    t = numer/ denom;
                }
                s = 1 -t;
            }
            else
            {
                t = 0;
                if(tmp1 <= 0)
                {
                    s = 1;
                }
                else if ( d >=0)
                {
                    s =0;
                }
                else
                {
                    s = -d/a;
                }
            }
        }
        else
        {
            //region 1
            float numer = (c+e) - (b+d);
            if(numer <= 0)
            {
                s = 0;
            }
            else
            {
                float denom = a -2 * b+c;
                if( numer >= denom)
                {
                    s =1;
                }
                else
                {
                    s = numer /denom;
                }
            }
            t = 1 -s;
        }
    }
    float distance_squared = (a * (s*s)) + (2 * b * s*t) + (c * (t*t)) + (2*d*s) + (2*e*t) + f;
    return std::sqrt(distance_squared);
    //return distance_squared;
}

void Models::addOBJModel(const std::string filepath)
{
    GridBox gridbox = makeGridBox(Vector(-2,-2,-2),Vector(2,2,2),Vector(0.01,0.01,0.01));
    ScalarGrid modelgrid = makeGrid(gridbox, -10000000);

    m = obj::loadModelFromFile(filepath);

    std::vector<Triangle> triangles = ObjLoader::loadObj(filepath);
    // std::vector<Triangle> triangles;
    // for(auto it = m.faces.cbegin(); it != m.faces.cend(); ++it) //triangles
    // {
    //     for(int i = 0; i < it->second.size(); i+=3) //verts
    //     {
    //         int ind1 = it->second[i]; //vert index
    //         int ind2 = it->second[i+1];
    //         int ind3 = it->second[i+2];

    //         Vector v1(m.vertex[ind1], m.vertex[ind1+1], m.vertex[ind1+2]); //m.vertex stores vertex in uniform order......
    //         Vector v2(m.vertex[ind2], m.vertex[ind2+1], m.vertex[ind2+2]);
    //         Vector v3(m.vertex[ind3], m.vertex[ind3+1], m.vertex[ind3+2]);
    //         std::cout << ind1 << ' ' << ind2 << ' ' << ind3 << ' ' << v1.X() << ' ' << v1.Y() << ' ' << v1.Z() << '\n';
    //         triangles.push_back(Triangle(v1, v2, v3));
    //     }
    //     //std::cout << triangles.size() << '\n';
    // }    
    std::cout << triangles.size() << '\n';
    ProgressMeter pm(triangles.size(), "obj load");
    #pragma omp parallel for
    for(int a = 0; a < triangles.size(); a++)
    {
        int bandwith = 5;
        int i,j,k;
        modelgrid->getGridIndex(triangles[a].v1, i, j, k);
        Vector v1_gridindex(i,j,k);
        modelgrid->getGridIndex(triangles[a].v2, i, j, k);
        Vector v2_gridindex(i,j,k);
        modelgrid->getGridIndex(triangles[a].v3, i, j, k);
        Vector v3_gridindex(i,j,k);
//AABB
        int lowx = std::min({v1_gridindex.X(), v2_gridindex.X(), v3_gridindex.X()}) -1;
        int lowy = std::min({v1_gridindex.Y(), v2_gridindex.Y(), v3_gridindex.Y()}) -1 ;
        int lowz = std::min({v1_gridindex.Z(), v2_gridindex.Z(), v3_gridindex.Z()}) -1 ;
        modelgrid->fix_indices(lowx, lowy, lowz);

        int maxx = std::max({v1_gridindex.X(), v2_gridindex.X(), v3_gridindex.X()}) +1;
        int maxy = std::max({v1_gridindex.Y(), v2_gridindex.Y(), v3_gridindex.Y()}) +1;
        int maxz = std::max({v1_gridindex.Z(), v2_gridindex.Z(), v3_gridindex.Z()}) +1;
        modelgrid->fix_indices(maxx, maxy, maxz);
//BANDWITH
        lowx -= bandwith; lowy -= bandwith; lowz -= bandwith;
        modelgrid->fix_indices(lowx, lowy, lowz);
        maxx += bandwith; maxy += bandwith; maxz += bandwith;
        modelgrid->fix_indices(maxx, maxy, maxz);
        //std::cout << lowx << ' ' << lowy << ' ' << lowz << ' ' << maxx << ' ' << maxy << ' ' << maxz << '\n';
        for(int j = lowy; j <= maxy; j++)
        {
            for(int i = lowx; i <= maxx; i++)
            {
                for(int k = lowz; k<= maxz; k++)
                {
                    Vector x_ijk = modelgrid->evalP(i, j, k);
                    float distance = calcDistance(x_ijk, triangles[a]);
                    if( distance <= std::abs(modelgrid->get(i,j,k)))
                    {
                        modelgrid->set(i, j, k, distance);
                       // std::cout << i << ' ' << j << ' ' << k << ' ' << distance << '\n';
                    }
                }
            }
        }

    }

    for(int j = 0; j < modelgrid->ny(); j++)
    {
        for (int k = 0; k < modelgrid->nz(); k++)
        {
            Vector g_ijk = modelgrid->evalP(0, j, k); //start position of ray
            Vector g_n = Vector(1,0,0); //direction of ray
            std::vector<Vector> intersection_points;
            for(int a = 0; a < triangles.size(); a++)
            {
                Vector t_n = (triangles[a].e2 ^ triangles[a].e1) / (triangles[a].e2 ^ triangles[a].e1).magnitude();
                float d = ((t_n * (triangles[a].v1 - g_ijk)) / (t_n * g_n));
                float u = ( (triangles[a].e2 ^ (g_ijk - triangles[a].v1 + (d*g_n))) * (triangles[a].e2 ^ triangles[a].e1) ) /
                            ((triangles[a].e2 ^ triangles[a].e1).magnitude() * (triangles[a].e2 ^ triangles[a].e1).magnitude());
                float v = ( -1 * (triangles[a].e1 ^ (g_ijk - triangles[a].v1 + (d*g_n))) * (triangles[a].e2 ^ triangles[a].e1) ) /
                            ((triangles[a].e2 ^ triangles[a].e1).magnitude() * (triangles[a].e2 ^ triangles[a].e1).magnitude());
                //find all intersections first
                if(u >= 0 && u <= 1 && v >=0 && v <=1 && (u+v) <= 1 && (u+v) >= 0)
                {
                    //intersection at g_ijk + d * g_n
                    //intersections++;
                    intersection_points.push_back(g_ijk + d*g_n);
                
                }
            }
            // if(intersection_points.size() != 0)
            //     std::cout << "sects points: " << intersection_points.size() << '\n';

            int intersections = 0;
            for(int i = 0; i < modelgrid->nx(); i++)
            {
                float grid_value = modelgrid->get(i, j, k);
                Vector grid_point = modelgrid->evalP(i,j,k);
   
                for(int a = intersections; a < intersection_points.size(); a++) // see if grid point is past intersection
                {
                    // if(grid_point.X() >= intersection_points[a].X())
                    // {
                    //     intersections++;
                    //     break;
                    // }
                    int i_x, i_y, i_z;
                    modelgrid->getGridIndex(intersection_points[a], i_x, i_y, i_z);
                    if(i >= i_x)
                    {
                        intersections++;
                        break;
                    }
                }
                if(intersections % 2 == 0)
                {
                    if(grid_value > 0) grid_value*=-1;
                        modelgrid->set(i, j, k, grid_value);
                }
                else
                {
                    if(grid_value < 0) grid_value*=-1;
                        modelgrid->set(i, j, k, grid_value);
                }                
            }
            // if(intersections != 0)
            //     std::cout << "sects : " << intersections << '\n';
        
        }
 
    }

    ScalarField model = gridded(modelgrid);
    model = scale(model, Vector(6,6,6));
    addScalarModel(model, Color(0.3,0.3,0.3,1));

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
    //addScalarModel(e2, Color(0.1,0.1,0.1,1));
    // addScalarModel(e9, red);
    // addScalarModel(e4, green);
    // addScalarModel(e5, green);
    // addScalarModel(e6, green);
    // addScalarModel(e7, green);
    // addScalarModel(e8, green);
    // addScalarModel(e10, red);
    // addScalarModel(e13, pink);
    // addScalarModel(e14, white);
    // addScalarModel(e16, red);
    // addScalarModel(e17, pink);



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