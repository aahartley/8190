#ifndef VOLUME_RENDERER_H
#define VOLUME_RENDERER_H

#include <vector>
#include <cmath>
#include "Color.h"
#include "Camera.h"
#include "Model.h"
#include "Fields.h"
#include "ProgressMeter.h"
#include "ImgProc.h"

namespace lux{

class VolumeRenderer
{
  public:
    VolumeRenderer(){}
    VolumeRenderer(int frames);
    VolumeRenderer(int startFrame, int endFrame);
    ~VolumeRenderer(){}

    void addModel(Model* m) { models.push_back(m);}
    void addImgProc(img::ImgProc* imgP) { imgProc = imgP;}
    void addCam(Camera* cam){camera = cam;}
    void display();
    void raymarch(double snear, double sfar, double Tmin, double ds, double kappa, ScalarField& density, ColorField& Cm);
    int getStart(){return start;}

  private:
    std::vector<Model*> models;
    int start, end;
    img::ImgProc* imgProc;
    Camera* camera;
    bool rotate_table;

};

}


#endif