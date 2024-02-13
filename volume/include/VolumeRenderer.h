#ifndef VOLUME_RENDERER_H
#define VOLUME_RENDERER_H

#include <vector>
#include <cmath>
#include "Color.h"
#include "Camera.h"
#include "Models.h"
//#include "Fields.h"
#include "ProgressMeter.h"
#include "ImgProc.h"

namespace lux{

class VolumeRenderer
{
  public:
    VolumeRenderer(){}
    VolumeRenderer(int startFrame, int endFrame);
    ~VolumeRenderer(){}

    void addModels(std::shared_ptr<Models> m) { models = m;}
    void addImgProc(std::shared_ptr<img::ImgProc> imgP) { imgProc = imgP;}
    void addCam(std::shared_ptr<Camera> cam){camera = cam;}
    void generate_frames();
    void raymarch(double snear, double sfar, double Tmin, double ds, double kappa, ScalarField& density, ColorField& Cm);

  private:
    std::shared_ptr<Models> models;
    int start, end;
    std::shared_ptr<img::ImgProc> imgProc;
    std::shared_ptr<Camera> camera;
    bool rotate_table;

};

}


#endif