

#ifndef __IMPLICITCOLORS_H__
#define __IMPLICITCOLORS_H__


#include "Volume.h"
#include "LinearAlgebra.h"
#include "FullGrids.h"
#include "SparseGrids.h"
namespace lux
{

class ConstantColor : public Volume<Color> 
{
  public:

    ConstantColor( const Color& v ) :
       value (v)
    {}

   ~ConstantColor(){}


    const Color eval( const Vector& P ) const 
    {
       return value; 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Color value;
};





class AddColor : public Volume<Color> 
{
  public:

    AddColor( Volume<Color> * e, Volume<Color>* v ) :
      elem1(e),
      elem2(v)
    {}

    AddColor( const ColorField& e, const ColorField& v ) :
      elem1(e),
      elem2(v)
    {}

   ~AddColor(){}


    const Color eval( const Vector& P ) const 
    {
       return ( elem1->eval(P) + elem2->eval(P) ); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Add";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1, elem2;
};

class MultiplyColor : public Volume<Color> 
{
  public:

    // MultiplyColor( Volume<Color> * e, Volume<Color>* v ) :
    //   elem1(e),
    //   elem2(v)
    // {}

    MultiplyColor( const ColorField& e, const ScalarField& v ) :
      elem1(e),
      elem2(v)
    {}

   ~MultiplyColor(){}


    const Color eval( const Vector& P ) const 
    {
       return ( elem1->eval(P) * elem2->eval(P) ); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1;
    ScalarField elem2;
};

class GriddedGridColor : public Volume<Color> 
 {
   public:
  
     GriddedGridColor( const ColorGrid& g );
  
     ~GriddedGridColor(){}
  
     const Color eval( const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     ColorGrid scgrid;
     float dx, dy, dz;
 };
  
class GriddedSGridColor : public Volume<Color> 
 {
   public:
  
     GriddedSGridColor( const SColorGrid& g );
  
     ~GriddedSGridColor(){}
  
     const Color eval( const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     SColorGrid scgrid;
     float dx, dy, dz;
 };
  



}
#endif
