#include "VolumeRenderer.h"

using namespace lux;



VolumeRenderer::VolumeRenderer(int s, int e) : start(s), end(e)
{
    rotate_table = false;
    wedge = false;
    sim = true;
}

void VolumeRenderer::addDSM(ScalarField& sf, const ScalarField& density, double ds, double kappa, int index)
{
    GridBox gb = makeGridBox(Vector(-10,-10,-10),Vector(10,10,10),Vector(0.05,0.05,0.05));
    SScalarGrid lgrid = makeSGrid(gb, 0.f);

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
    sf = exp((gridded(lgrid)*constant((float)-kappa)));
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
    if(!wedge && !sim)
    {
        density = models->getGriddedClampedDensityField(0.0,1.0);
        colorfield =  models->getGriddedColorField();
        //density = models->getClampedDensityField(0.0,0.1);
        //colorfield =  models->getColorField();

        addDSM(TL, density, 0.03, 1, 0);
        //addDSM(TL2, density, 0.03, 1, 1);
        //addDSM(TL3, density, 0.03, 1, 2);
        dsmField = std::vector<ScalarField>{TL};
        //dsmField = std::vector<ScalarField>{TL, TL2, TL3};
    }
    NoiseData nd;
    NoiseData nd1;
    nd1.frequency = 1000;
    nd1.fjump = 1000;
    nd1.roughness = 10;
    nd1.octaves = 20;
    nd1.gamma = 10;
    nd1.translate = Vector(0,0,0);
    NoiseData nd2;

    std::shared_ptr<PerlinNoise> pn = std::make_shared<PerlinNoise>();
    NoiseSrc ns = pn;
    std::shared_ptr<FractalSum> fs = std::make_shared<FractalSum>(ns);
    _Noise noise = fs;
    noise->setParameters(nd);
    VectorField U = VFNoise(noise);

    ScalarField bun = models->addOBJModel("models/bunny/bunny.obj", Vector(-0.0945,0.032,-0.062), Vector(0.061,0.19,0.059), Vector(0.0005,0.0005,0.0005));
    bun = scale(bun, Vector(30,30,30));
    bun = translate(bun, Vector(0,-3,0));

    // std::shared_ptr<PerlinNoise> pn2 = std::make_shared<PerlinNoise>();
    // NoiseSrc ns2 = pn2;
    // std::shared_ptr<FractalSum> fs2 = std::make_shared<FractalSum>(ns2);
    // _Noise noise2 = fs2;
    // noise2->setParameters(nd1);
    // bun = PyroVolume(bun, 2, 0.3, noise2);

    float dt = 0.05;
    // GridBox gb = makeGridBox(Vector(-3,-3,-3),Vector(3,3,3),Vector(0.05,0.05,0.05));
    // int N = 2;
    // float dtt = dt / N;
    // for(int i = 0; i < N-1; i++)
    // {
    //     bun = advect(bun, U, dtt);
    //     SScalarGrid gridd = makeSGrid(gb, 0);
    //     stamp(gridd,bun,1);
    //     bun = gridded(gridd);
    // }
    // bun = advect(bun, U, dtt);

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
            eye = Vector(cam_distance * sin(angleRadians), 1, cam_distance * cos(angleRadians));
            view = -(eye.unitvector());
        }
        else
        {
            eye = Vector(0,0,cam_distance); view = Vector(0,0,-1);
        }
        if(wedge)
        {
            //models->accumulateNoiseParam(nd1, i, "pyro");
            //models->accumulateNoiseParam(nd1, i, "ifn");
            //models->accumulateNoiseParam(nd2, i, "wisp2");
            //models->addRandPyroSphere(nd1);
            //models->addPyroSphere(nd1);
            //models->addIFNoise(nd1);
            //models->addRandIFNoise(nd1);
            //models->addWisp(nd1, nd2);
            models->addRandWisp(nd1, nd2);
            density = models->getGriddedClampedDensityField(0.0,1.0);
            colorfield =  models->getGriddedColorField();
            //density = models->getClampedDensityField(0.0,0.1);
            //colorfield =  models->getColorField();

            addDSM(TL, density, 0.03, 1, 0);
            addDSM(TL2, density, 0.03, 1, 1);
            addDSM(TL3, density, 0.03, 1, 2);
            dsmField = std::vector<ScalarField>{TL, TL2, TL3};
        }
        if(sim)
        {
            models->sim1(i,nd,nd1,nd2, U, bun, dt);
        
            density = models->getGriddedClampedDensityField(0.0,1.0);
            bun = density;
            colorfield =  models->getGriddedColorField();
            // density = models->getClampedDensityField(0.0,0.1);
            // colorfield =  models->getColorField();

            addDSM(TL, density, 0.01, 1, 0);
            dsmField = std::vector<ScalarField>{TL};
        }
        camera->setEyeViewUp(eye, view, Vector(0,1,0));
        raymarch(1, 20, 0, 0.01, 1, density, colorfield );//0.005
        //imgProc->write_image("image_"+std::to_string(i), 'o');
        //imgProc->write_image("image_"+std::to_string(i), 'j');
        std::string img_num = std::to_string(i);
        std::stringstream ss;
        ss << std::setw(4) << std::setfill('0') << img_num;
        img_num = ss.str();    
        imgProc->write_image("test."+img_num, 'o');
        imgProc->write_image("test."+img_num, 'j');

        if(wedge || models->sim_reset) models->reset();
        pm.update();

    }
}
