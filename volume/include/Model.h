#ifndef MODEL_H
#define MODEL_H

#include "Fields.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "Volume.h"

namespace lux
{
class ModelBase
{
  public:
    ModelBase(){}
    virtual ~ModelBase(){}
};

template< typename U >
class Model : ModelBase
{
  public:
    Model(){}
    virtual ~Model() {}

    typedef U model_type;

    model_type& getDensityField() {return densityfield;}
    ColorField& getColorField() {return colorfield;}

  protected:
    model_type densityfield;
    ColorField colorfield;
};
typedef std::shared_ptr<Model<ScalarField>> ScalarModelBase;

class ScalarModel : public ScalarModelBase
{
  public:
    ScalarModel() :  std::shared_ptr<Model<ScalarField> >() {}
    ScalarModel(Model<ScalarField>* m) : std::shared_ptr<Model<ScalarField> >( m) {}
    ~ScalarModel() {}

 
};

class Humanoid : public Model<ScalarField>
{
  public:
    Humanoid();
    ~Humanoid(){}
};


}



#endif