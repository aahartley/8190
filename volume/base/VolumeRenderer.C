#include "VolumeRenderer.h"

using namespace lux;


VolumeRenderer::VolumeRenderer(int frames) : start(0) , end(frames)
{
}

VolumeRenderer::VolumeRenderer(int s, int e) : start(s), end(e)
{
}

void VolumeRenderer::raymarch(double snear, double sfar, double Tmin, double ds, double kappa, ScalarField& density, ColorField& Cm)
{
    //#pragma omp parallel for collapse(2)
    for(int j = 0; j < imgProc->ny(); j++)
 	{
 		for(int i = 0; i < imgProc->nx(); i++)
 		{
  
            Vector direction = camera->view((double)i/(double)(imgProc->nx()-1),(double)j/(double)(imgProc->ny()-1));
            //std::cout << direction.X() << ' ' << direction.Y() << ' ' << direction.Z() << '\n';
            double T = 1;
            Color L(0,0,0,0);
            double s = snear;
            Vector X = camera->eye() + (direction*s);
            //std::cout << X.X() << ' ' << X.Y() << ' ' << X.Z() << '\n';
            while( s < sfar && T > Tmin )
            {
                //std::cout << X.X() << ' ' << X.Y() << ' '<< X.Z() << '\n';
                float den = density->eval(X);
                if( den > 0.0 )
                {
                    float dT = std::exp( -ds * kappa * den );
                    L += Cm->eval(X) * (1-dT) * T / kappa;
                    T *= dT;
                }
                X += direction * ds;
                s += ds;
            }
            L[3] = 1-T; // set the alpha channel to the opacity
            //if(L.X() == 0 &&  L.Y()==0 && L.Z()==0) L[3]=1;
            //std::cout << L.X() << ' ' << L.Y() << ' ' << ' '<< L.Z() << ' ' << L.alpha() << '\n';
            imgProc->set_value(i, j, std::vector<float>{L.X(),L.Y(),L.Z(),L.alpha()});
        }
    }
                
    
}

void VolumeRenderer::display()
{
    Vector p(0,2,0); Vector p2(0,4,0); Vector p3(0,5,0); Vector p4(0,0,0);Vector p5(0,-1.9,0);
    Vector n(1,0,0); Vector n1(0,1,0); Vector n2(0,0,1);
    Color color( 1,0.2,0,1); Color color2(0,1,0,1); Color black(0,0,0,1);


    ScalarField c = constant(-1000);
    //ScalarField e1 = new BoxVolume(p, 1, 3);
    //ScalarField e1 = new EllipseVolume(p,3,1,Vector(1,0,0));
    ScalarField e1 = Sphere(p2, 2);
    ScalarField e2 = Sphere(p,2);
    ScalarField e3 = Sphere(p3,2);
    ScalarField e4 = Ellipse(p5, 2, 1.5, n1);
    ScalarField e5 = Ellipse(Vector(-.5, -3.8, 1), .3, .6, n1);
    ScalarField e6 = Ellipse(Vector(.5, -3.8, 1), .3, .6, n1);
    ScalarField e7 = Ellipse(Vector(-2, -1, 1), 0.5, 1.3, n1);
    ScalarField e8 = Ellipse(Vector(2, -1, 1), 0.5, 1.3, n1);
    ScalarField e9 = Torus(Vector(0,1.8,0.6), 1, 0.5, n1);


    //ScalarField e1 = Cone(p, n, 1, 45);
    //ScalarField e1 = Plane(p, n1, p);
    //ScalarField e1 = Torus(p, 2, 1, n);
    //ScalarField e1 = Cylinder(p, n, 2);
    //ScalarField e1 = SteinerPatch(p);
    //ScalarField e1 = Icosahedron(p);
    //ScalarField e1 = constant(5);
    e1 = Cutout(e1, e3);
    e1 = scale(e1, Vector(2,1,1));
    ScalarField u1 = Union(c, e1); //crown added
    u1 = Union(u1, e2); //head added
    u1 = Union(u1, e4); // torso
    //u1 = Union(u1, e5); // left foot
    //u1 = Union(u1, e6); // right foot
    //u1 = Union(u1, e7); // left arm
    //u1 = Union(u1, e8); // right arm
    u1 = Union(u1, e9); // eyes


    ColorField e1Color = constant(Color(0,0,0,0)); //start with black
    e1Color = e1Color *mask(-e1) + constant(color)*mask(e1); //crown color
    e1Color = e1Color *mask(-e2) + constant(color2)*mask(e2); //head color
    e1Color = e1Color *mask(-e4) + constant(color2)*mask(e4); //torso color
    // e1Color = e1Color *mask(-e5) + constant(color2)*mask(e5); //left feet color
    // e1Color = e1Color *mask(-e6) + constant(color2)*mask(e6); //left feet color
    // e1Color = e1Color *mask(-e7) + constant(color2)*mask(e7); //left arm color
    // e1Color = e1Color *mask(-e8) + constant(color2)*mask(e8); //right arm color
    e1Color = e1Color *mask(-e9) + constant(black)*mask(e9); //eyes color



    ScalarField density = mask(u1);

    ProgressMeter* pm = new ProgressMeter(end,"demo");
    for(int i = start; i < end; i++)
    {
        int last = end - 1;
        if(last == 0) last = 1;
        double angleDegrees = 360.0 * i / (last);
        double angleRadians = angleDegrees * M_PI / 180.0;
        Vector eye(10 * sin(angleRadians), 0, 10 * cos(angleRadians));
        Vector view = -(eye.unitvector());
        camera->setEyeViewUp(eye, view, Vector(0,1,0));
        raymarch(1, 20, 0, 0.9, 1, density, e1Color );
        imgProc->write_image("image_"+std::to_string(i), 'j');
        pm->update();
    }
    delete pm;
}