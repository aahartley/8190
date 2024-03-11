#include "VolumeRenderer.h"

using namespace lux;



VolumeRenderer::VolumeRenderer(int s, int e) : start(s), end(e)
{
    rotate_table = false;
    wedge = true;
}

void VolumeRenderer::addDSM(ScalarField& sf, const ScalarField& density, double ds, double kappa, int index)
{
    GridBox gb = makeGridBox(Vector(-10,-10,-10),Vector(10,10,10),Vector(0.05,0.05,0.05));
    ScalarGrid lgrid = makeGrid(gb, 0.f);

    ProgressMeter pm(1,"dsm"+std::to_string(index));
    #pragma omp parallel for collapse(3)
    for(int j = 0; j < lgrid->ny(); j++)
 	{
 		for(int i = 0; i < lgrid->nx(); i++)
 		{
            for(int k = 0; k < lgrid->nz(); k++ )
            {
                double arg = 0;
                Vector X(lgrid->llc().X()+(i*lgrid->dx()), lgrid->llc().Y()+(j*lgrid->dy()), lgrid->llc().Z()+(k*lgrid->dz()));
                if(density->eval(X) > 0.0)
                {
                    double smax = (lightPos[index]-X).magnitude();
                    Vector direction = (lightPos[index]-X).unitvector();  
                    double s = 0;
                    while(s < smax)
                    {
                        arg += density->eval(X)*ds;
                        X += direction * ds;
                        s += ds;
                    }   
                    lgrid->set(i,j,k,arg);
                }

            }
        }
    }
    sf = exp((constant((float)-kappa)*gridded(lgrid)));
    pm.update();

}

void VolumeRenderer::raymarch(double snear, double sfar, double Tmin, double ds, double kappa, ScalarField& density, ColorField& Cm)
{
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < imgProc->ny(); j++)
 	{
 		for(int i = 0; i < imgProc->nx(); i++)
 		{
  
            Vector direction = camera->view((double)i/(double)(imgProc->nx()-1), (double)j/(double)(imgProc->ny()-1));
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
                    Color Clights = Color(0,0,0,0);
                    for(int a = 0; a < dsmField.size(); a++)
                    {
                        Clights += lightColor[a] * dsmField[a]->eval(X);
                    }
                    L += Cm->eval(X) * Clights * (1-dT) * T / kappa;
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
    ScalarField density;
    ColorField colorfield;
    ScalarField TL, TL2, TL3;
    lightColor = std::vector<Color>{ Color(1,1,1,1), Color(0.25,0.25,0.25,1), Color(0.4,0.4,0.4,1)}; //key, rim, fill
    lightPos = std::vector<Vector>{ Vector(0,9.9,0), Vector(0,-9.9,0), Vector(0,0,-9.9)};
    if(!wedge)
    {
        density = models->getGriddedClampedDensityField(0.0,1.0);
        colorfield =  models->getGriddedColorField();
        //density = models->getClampedDensityField(0.0,0.1);
        //colorfield =  models->getColorField();

        addDSM(TL, density, 0.03, 1, 0);
        addDSM(TL2, density, 0.03, 1, 1);
        addDSM(TL3, density, 0.03, 1, 2);
        dsmField = std::vector<ScalarField>{TL, TL2, TL3};
    }


    ProgressMeter pm (end-start,"demo");
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
        if(wedge)
        {
            //models->addRandPyroSphere();
            //models->addPyroSphere(i);
            models->addWisp(1);
            density = models->getGriddedClampedDensityField(0.0,1.0);
            colorfield =  models->getGriddedColorField();

            addDSM(TL, density, 0.03, 1, 0);
            addDSM(TL2, density, 0.03, 1, 1);
            addDSM(TL3, density, 0.03, 1, 2);
            dsmField = std::vector<ScalarField>{TL, TL2, TL3};
        }
        camera->setEyeViewUp(eye, view, Vector(0,1,0));
        raymarch(1, 20, 0, 0.005, 1, density, colorfield );//fix
        //imgProc->write_image("image_"+std::to_string(i), 'o');
        //imgProc->write_image("image_"+std::to_string(i), 'j');
        imgProc->write_image("test"+std::to_string(i), 'o');
        imgProc->write_image("test"+std::to_string(i), 'j');

        if(wedge) models->reset();
        pm.update();

    }
}
