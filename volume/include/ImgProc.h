//-------------------------------------------------------
//
//  ImgProc.h
//
//  Stores the image pixels and the methods to manipulate
//  them. Handles the IO of images.
//  
//--------------------------------------------------------
#ifndef IMGPROC_H
#define IMGPROC_H

#include <vector>
#include <iostream>
#include <string>
#include <regex>
#include <algorithm>
#include <ctime>

#include <OpenImageIO/imageio.h>

namespace img{

class ImgProc
{
  public:
	
	ImgProc();
	~ImgProc();

	void clear ();
	void clear(int nX, int nY, int nC);

	int nx() const { return Nx;} 
	int ny() const { return Ny;}
	int depth() const { return Nc;}
	int size() const {return Nsize;}
	float* raw() const { return img_data;}

	void set_raw(float* data);
	void set_raw_OI(float* data);

	void value(int i, int j, std::vector<float>& pixel) const;
	void set_value(int i, int j, const std::vector<float>& pixel);

    int read_image(const std::string& s);
	void write_image(std::string fileName, char f) const;

	long index(int i, int j, int c) const;
	long indexNon(int i, int j, int c) const;

	//Affine transformations
	void Flip();
	void Flip(const float* in);
	void Flop();

	ImgProc(const ImgProc& img); //copy constructor
	ImgProc& operator=(const ImgProc& img); //copy assignment

  private:	
  
	int Nx, Ny, Nc;
	long Nsize;
	float* img_data;
};

int read_image(const std::string& s, ImgProc& imgProc);
void write_image(std::string fileName, char f, const ImgProc& imgProc);

//Affine transformations
void Flip(ImgProc& in);
void Flip(const float* in, int nx, int ny, int nc, float* out, int Nx, int Ny, int Nc);
void Flop(ImgProc& in);




}
#endif
