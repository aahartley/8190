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
    VolumeRenderer(int startFrame, int endFrame);
    ~VolumeRenderer(){}

    void addScalarModel(ScalarModel m) { scalar_models.push_back(m);}
    void addImgProc(std::shared_ptr<img::ImgProc> imgP) { imgProc = imgP;}
    void addCam(std::shared_ptr<Camera> cam){camera = cam;}
    void generate_frames();
    template <typename U>
    void raymarch(double snear, double sfar, double Tmin, double ds, double kappa, U& density, ColorField& Cm);

  private:
    std::vector<ScalarModel> scalar_models;
    int start, end;
    std::shared_ptr<img::ImgProc> imgProc;
    std::shared_ptr<Camera> camera;
    bool rotate_table;

};

}


#endif