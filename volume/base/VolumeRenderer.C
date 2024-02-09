#include "VolumeRenderer.h"

using namespace lux;



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
            double T = 1;
            Color L(0,0,0,0);
            double s = snear;
            Vector X = camera->eye() + (direction*s);
            while( s < sfar && T > Tmin )
            {
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
            imgProc->set_value(i, j, std::vector<float>{L.X(),L.Y(),L.Z(),L.alpha()});
        }
    }
                
    
}

void VolumeRenderer::generate_frames()
{
    ProgressMeter pm (end,"demo");
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
        raymarch(1, 20, 0, 4, 1, models->getClampedDensityField(0,1), models->getColorField() );
        //imgProc->write_image("image_"+std::to_string(i), 'o');
        //imgProc->write_image("image_"+std::to_string(i), 'j');
        imgProc->write_image("test"+std::to_string(i), 'j');

        pm.update();

    }
}