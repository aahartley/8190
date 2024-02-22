#include "Models.h"

using namespace lux;

void Models::createColorField(ScalarField& model, Color color)
{
    //scalar_volumes_unioned = Union(scalar_volumes_unioned, model);

    colorfield = colorfield *mask(-model) + constant(color)*mask(model); 

}

void Models::createFinalUnion(ScalarField& model)
{
    scalar_volumes_unioned = Union(scalar_volumes_unioned, model);

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
    density = clamp(scalar_volumes_unioned, min, max);
    stamp(dgrid, density, 1);
    return gridded(dgrid);
}


ColorField Models::getGriddedColorField() 
{
    ColorGrid cgrid = makeGrid(gb, Color(0,0,0,0));
    stamp(cgrid, colorfield, 1);
    return gridded(cgrid);
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
    return std::sqrt(std::abs(distance_squared));
    //return distance_squared;
}

ScalarField Models::addOBJModel(const std::string filepath,Vector llc, Vector urc, Vector dims)
{
    //freopen("output.txt","w",stdout);

    std::vector<Triangle> triangles = ObjLoader::loadObj(filepath);
    int bandwith = 2;
    float alx = triangles[0].v1.X(); float aly = triangles[0].v1.Y(); float alz = triangles[0].v1.Z();
    float aux = triangles[0].v1.X(); float auy = triangles[0].v1.Y(); float auz = triangles[0].v1.Z();
    //#pragma omp parallel for
    for(int a = 0; a < triangles.size(); a++)
    {
        Vector v1 = triangles[a].v1;
        Vector v2 = triangles[a].v2;
        Vector v3 = triangles[a].v3;
        if(v1.X() < alx) alx= v1.X();         
        if(v2.X() < alx) alx= v2.X();
        if(v3.X() < alx) alx= v3.X();
        if(v1.Y() < aly) aly= v1.Y();         
        if(v2.Y() < aly) aly= v2.Y();
        if(v3.Y() < aly) aly= v3.Y();
        if(v1.Z() < alz) alz= v1.Z();         
        if(v2.Z() < alz) alz= v2.Z();
        if(v3.Z() < alz) alz= v3.Z();

        if(v1.X() > aux) aux= v1.X();         
        if(v2.X() > aux) aux= v2.X();
        if(v3.X() > aux) aux= v3.X();
        if(v1.Y() > auy) auy= v1.Y();         
        if(v2.Y() > auy) auy= v2.Y();
        if(v3.Y() > auy) auy= v3.Y();
        if(v1.Z() > auz) auz= v1.Z();         
        if(v2.Z() > auz) auz= v2.Z();
        if(v3.Z() > auz) auz= v3.Z();
        // std::cout << triangles[a].v1.X() << ' ' << triangles[a].v1.Y() << ' ' << triangles[a].v1.Z() << ' ' <<
        // triangles[a].v2.X() << ' ' << triangles[a].v2.Y() << ' ' << triangles[a].v2.Z() << ' ' <<
        // triangles[a].v3.X() << ' ' << triangles[a].v3.Y() << ' ' << triangles[a].v3.Z() << '\n';
    }
    // alx = (alx < 0) ? std::floor(alx)-1-bandwith : std::ceil(alx)-1-bandwith;
    // aly = (aly < 0) ? std::floor(aly)-1-bandwith :std::ceil(aly)-1-bandwith;
    // alz = (alz < 0) ? std::floor(alz)-1-bandwith :std::ceil(alz)-1-bandwith;
    // aux = (aux < 0) ? std::floor(aux)+1+bandwith :std::ceil(aux)+1+bandwith;
    // auy = (auy < 0) ? std::floor(auy)+1+bandwith :std::ceil(auy)+1+bandwith;
    // auz = (auz < 0) ? std::floor(auz)+1+bandwith :std::ceil(auz)+1+bandwith;
    std::cout << alx << ' ' << aly << ' '<< alz << ' ' << aux << ' ' << auy << ' ' << auz << '\n';
    //std::cout << dims.X() << ' ' << dims.Y() << ' ' << dims.Z() << '\n';
    std::cout << ((urc.X() - llc.X())/dims.X()) << ' ' << ((urc.Y() - llc.Y())/dims.Y()) << ' '<< ((urc.Z() - llc.Z())/dims.Z()) << '\n';
    //GridBox gridbox = makeGridBox(Vector(alx,aly,alz),Vector(aux,auy,auz),dims);
    GridBox gridbox = makeGridBox(llc,urc,dims);

    ScalarGrid modelgrid = makeGrid(gridbox, -10000000);
    
    std::cout << triangles.size() << '\n';
    ProgressMeter pm(triangles.size(), "obj load");
    //#pragma omp parallel for
    for(int a = 0; a < triangles.size(); a++)
    {
        int ii,jj,kk;
        modelgrid->getGridIndex(triangles[a].v1, ii, jj, kk);
        Vector v1_gridindex(ii,jj,kk);
        modelgrid->getGridIndex(triangles[a].v2, ii, jj, kk);
        Vector v2_gridindex(ii,jj,kk);
        modelgrid->getGridIndex(triangles[a].v3, ii, jj, kk);
        Vector v3_gridindex(ii,jj,kk);
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
                        //#pragma omp critical
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
            std::vector<float> intersection_points;
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
                    intersection_points.push_back(d);
                
                }
            }
            std::sort(intersection_points.begin(), intersection_points.end());
            // if(intersection_points.size() != 0 && intersection_points.size() % 2 != 0)
            //     std::cout << intersection_points.size() << '\n';
            int intersections = 0;
            for(int i = 0; i < modelgrid->nx(); i++)
            {
                float grid_value = modelgrid->get(i, j, k);   
                Vector grid_point = modelgrid->evalP(i,j,k);
                for(int a = intersections; a < intersection_points.size(); a++) // see if grid point is past intersection
                {
                    int i_x, i_y, i_z;
                    Vector intersec_ijk = g_ijk + (intersection_points[a] * g_n);
                    modelgrid->getGridIndex(intersec_ijk, i_x, i_y, i_z);
                    // if(i > i_x)
                    // {
                    //     intersections++;
                    //     break;
                    // }
                    if(grid_point.X() >= intersec_ijk.X())
                    {
                        intersections++;
                        break;
                    }
                }
                if(intersections % 2 == 0)
                {
                    if(grid_value > 0) grid_value*=-1;
                    //#pragma omp critical
                    modelgrid->set(i, j, k, grid_value);
                }
                else
                {
                    if(grid_value < 0) grid_value*=-1;
                    //#pragma omp critical
                    modelgrid->set(i, j, k, grid_value);
                }                
            }
            // if(intersections != 0 && intersections % 2 != 0)
            //     std::cout << "sects : " << intersections << '\n';
        
        }
 
    }

    ScalarField model = gridded(modelgrid);
    // model = translate(model, tsl);
    // model = scale(model, sc);
    //createColorField(model, Color(0.3,0.3,0.3,1));
    return model;
    //createFinalUnion(model);

}

ScalarField Models::addHumanoid()
{

    Color red( 1,0.2,0,1); Color green(0,1,0,1); Color black(0,0,0,1); Color blue(0,0,1,1); Color pink(1, 0.75, 0.8, 1); 
    Color white (1,1,1,1);

    ScalarField c = constant(-1000);
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


    createColorField(e1, red);
    createColorField(e2, green);
    // createColorField(e9, red);
    // createColorField(e4, green);
    // createColorField(e5, green);
    // createColorField(e6, green);
    createColorField(e7, green);
    // createColorField(e8, green);
    // createColorField(e10, red);
    // createColorField(e13, pink);
    // createColorField(e14, white);
    // createColorField(e16, red);
    // createColorField(e17, pink);

    ScalarField u1 =  e2; //head 
    //u1 = Union(u1, e9); //eyes
    u1 = Union(u1, e1); //horns
    // u1 = Union(u1, e4); // torso
    // u1 = Union(u1, e5); // left foot
    // u1 = Union(u1, e6); // right foot
    u1 = Union(u1, e7); // left arm
    // u1 = Union(u1, e8); // right arm
    // u1 = Union(u1, e10); // mouth
    // u1 = Union(u1, e13); //tattoo
    // u1 = Union(u1, e14); //stein
    // u1 = Union(u1, e16); // hat
    // u1 = Union(u1, e17); // tattoo

    ScalarGrid mgrid = makeGrid(gb, -10000.f);
    stamp(mgrid,u1 , 1);
    ScalarField gridsf = gridded(mgrid);
    return gridsf;
    //createFinalUnion(gridsf);

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

void Models::scene2()
{
    ScalarField aj = addOBJModel("models/ajax/smallajax.obj",Vector(-12.2,-25.51,-12.44), Vector(10.2,16.13,9.71), Vector(0.08, 0.08 ,0.08));
    // ScalarField cutter = Sphere(Vector(8,23,9), 25);
    //cutter = scale(cutter, Vector(12.721, 31.483, 12.537));

    //aj = Cutout(aj, cutter);

    // createFinalUnion(cutter);
    // createColorField(cutter, Color(0.,0,1,1));
    //aj = translate(aj, Vector(-35,15,4));
    aj = scale(aj, Vector(0.3,0.3,0.3));
    createColorField(aj, Color(0.3,0.3,0.3,1));

    // ScalarField hum = addHumanoid();

    // //ScalarField bun = addOBJModel("models/bunny/bunny.obj", Vector(-0.0945,0.032,-0.062), Vector(0.061,0.19,0.059), Vector(0.001,0.001,0.001) );
    // ScalarField bun = addOBJModel("models/bunny/bunny.obj", Vector(-0.0945,0.032,-0.062), Vector(0.061,0.19,0.059), Vector(0.00012,0.00012,0.00012));
    // //bun = translate(bun, Vector(-0.15,0.05,0));
    // bun = scale(bun, Vector(30,30,30));
    // createColorField(bun, Color(0.3,0.3,0.3,1));



    //createFinalUnion(hum);
    //createFinalUnion(bun);
    createFinalUnion(aj);   
}