#include "VolumeRenderer.h"

using namespace lux;


VolumeRenderer::VolumeRenderer(int frames) : start(0) , end(frames)
{
    rotate_table = true;
}

VolumeRenderer::VolumeRenderer(int s, int e) : start(s), end(e)
{
    rotate_table = true;
}

void VolumeRenderer::raymarch(double snear, double sfar, double Tmin, double ds, double kappa, ScalarField& density, ColorField& Cm)
{
    #pragma omp parallel for collapse(2)
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
    ScalarField e2 = Plane(Vector(0,0,0), Vector(0,1,0), Vector(0,0,0));
    ScalarField e1 = Torus(Vector(0,0,0), 2, 1, Vector(0,0,1));
    ScalarField e3 = SteinerPatch(Vector(0,0,0));
    e1 = Cutout(e1, e2);
    ScalarField density = clamp(e3, 0.0,1.0);
    density = mask(e3);
    ColorField colorfield = constant(Color(1,0,0,1));

    ProgressMeter* pm = new ProgressMeter(end,"demo");
    for(int i = start; i < end; i++)
    {
        Vector eye, view;
        double cam_distance = 15;
        if(rotate_table)
        {
            int last  = (end-1 == 0) ? 1 : end - 1;
            double angleDegrees = 360.0 * i / (last);
            double angleRadians = angleDegrees * M_PI / 180.0;
            eye = Vector(cam_distance * sin(angleRadians), 0, cam_distance * cos(angleRadians));
            view = -(eye.unitvector());
        }
        else
        {
            eye = Vector(0,0,cam_distance); view = Vector(0,0,-1);
        }
        camera->setEyeViewUp(eye, view, Vector(0,1,0));
        raymarch(1, 20, 0, 0.9, 1, models[0]->getDensityField(), models[0]->getColorField() );
        //raymarch(1, 20, 0, 0.9, 1, density, colorfield );
        imgProc->write_image("image_"+std::to_string(i), 'j');
        pm->update();

    }
    delete pm;
}