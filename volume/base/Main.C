
#include <iostream>
#include <ctype.h>
#include "Vector.h"
#include "VolumeRenderer.h"
#include "ImgProc.h"
#include "Camera.h"

using namespace lux;

int cmdLine(std::vector<std::string>& args, VolumeRenderer** renderer);


int main(int argc, char** argv)
{

    VolumeRenderer* renderer = nullptr;
    std::shared_ptr<img::ImgProc> imgProc (new img::ImgProc());
    std::shared_ptr<Camera> camera ( new Camera());
    std::shared_ptr<Models> models(new Models());
    imgProc->clear(1920, 1080, 4);
    std::vector<std::string> args;
    for(int i = 0; i < argc; i++)
    {
       std::string s(argv[i]);
       args.push_back(s);
    }    
    if(cmdLine(args, &renderer) == 1)
    {
        return 1;
    }
    else
    {
        models->addHumanoid();
        models->addOBJModel("models/bunny/bunny.obj");
        //models->addOBJModel("models/ajax/smallajax.obj");

        renderer->addModels(models);
        renderer->addImgProc(imgProc);
        renderer->addCam(camera);
        renderer->generate_frames();
        //delete renderer;
        return 0;
    }
}

int cmdLine(std::vector<std::string>& args, VolumeRenderer** renderer)
{
    if(args.size() > 1){
        for(unsigned int i = 1; i < args.size(); i++)
        {
            std::string s = args[i];
            if(s == "-frames")
            {
                int start, end;
                if(args.size() == 4)
                {
                    if(std::isdigit(static_cast<unsigned char>(args[i+2].at(0))) && std::isdigit(static_cast<unsigned char>(args[i+1].at(0))))
                    {
                        start = std::stoi(args[i+1]); end = std::stoi(args[i+2]);
                        *renderer = new VolumeRenderer(start, end);
                        return 0;
                    }
                    else
                    {
                        std::cout << "number please\n";
                        return 1;
                    }
                }
              

            }
        }
    }
    std::cout << "error\n";
    return 1;
}