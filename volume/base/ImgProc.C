//-------------------------------------------------------
//
//  ImgProc.C
//
//  Stores the image pixels and the methods to manipulate
//  them. Handles the IO of images.
//  
//--------------------------------------------------------
#include "ImgProc.h"

using namespace OIIO;
using namespace img;

inline int Index(int i, int j, int c,int Nc, int Nx) { return c+(Nc*(i+Nx*j));}


ImgProc::ImgProc() : Nx(0), Ny(0), Nc(0), Nsize(0), img_data(0) {}

ImgProc::~ImgProc()
{
	clear();
}

void ImgProc::clear()
{
	if(img_data !=0) { delete[] img_data; img_data=0; }
	Nx=0;
	Ny=0;
	Nc=0;
	Nsize=0;
}
//clears data, sets size for img
void ImgProc::clear(int nX, int nY, int nC)
{
	clear();
	Nx = nX;
	Ny = nY;
	Nc = nC;
	Nsize = (long)Nx * (long)Ny * (long)Nc;
	img_data = new float[Nsize];
	
	#pragma omp parallel for
	for(long i = 0; i < Nsize; i++) { img_data[i] = 0.0; }
}
//updates pixel color from img_data
void ImgProc::value( int i, int j, std::vector<float>& pixel) const
{
	pixel.clear();
	if( img_data == nullptr ){ return; }
	if( i < 0 || i >= Nx ){ return; }
	if( j < 0 || j >= Ny ){ return; }
	pixel.resize(Nc);
	for(int c = 0; c < Nc; c++)
	{
		pixel[c] = img_data[index(i,j,c)];
	}
	return;
}
//updates img_data from pixel colors
void ImgProc::set_value( int i, int j, const std::vector<float>& pixel)
{
	if( img_data == nullptr ){ return; }
	if( i < 0 || i >= Nx ){ return; }
	if( j < 0 || j >= Ny ){ return; }
	if( Nc > (int)pixel.size() ){ return; }
	//#pragma omp parallel for  //outer loops already in parallel
	for( int c = 0; c < Nc; c++)
	{
		img_data[index(i,j,c)] = pixel[c];
	}
	return;
}

void ImgProc::set_raw(float* in)
{
	if( img_data == nullptr ){ return; }

	#pragma omp parallel for  
	for(int i = 0; i < Nx*Ny*Nc; i++)
	{
		img_data[i] = in[i];

	}
	return;
}
void ImgProc::set_raw_OI(float* data)
{
	Flip(data);
	return;
}

int ImgProc::read_image(const std::string& s)
{
	std::string filename = s;
	std::cout << "\nAttempting to find: " << filename << std::endl;
	auto inp = ImageInput::open(filename);
	if (!inp)
	{
		std::cout << "Couldn't find " << filename << std::endl;
    	return 1;
	}
	const ImageSpec &spec = inp->spec();
	int xres = spec.width;
	int yres = spec.height;
	int nchannels = spec.nchannels;
	auto pixels = std::unique_ptr<float[]>(new float[xres * yres * nchannels]);
	inp->read_image(0, 0, 0, nchannels, TypeDesc::FLOAT, &pixels[0]);
	inp->close();

	clear(xres,yres,nchannels);
	set_raw_OI(&pixels[0]);

	
	std::cout <<"Image loaded\n\n";
	return 0;
}

//function will not be called unless image has been read
void ImgProc::write_image(std::string fileName, char f) const
{
	int xres = Nx;
	int yres = Ny;
	int channels = 0;
	float* pixels = nullptr;

	std::string nfileName;
	std::size_t pos = fileName.find(".");
    std::string fn = fileName.substr(0, pos);
	time_t now = time(0);
  	std::string date = ctime(&now);
	std::replace(date.begin(), date.end(), ':', '_');
	std::replace(date.begin(), date.end(), ' ', '_');
	date.pop_back(); //remove weird character at the end

	if(f =='j') ///jpg
	{
		channels = 3;
		//save new file to the original filepath
		nfileName = fn + ".jpeg";

	}
	else if(f == 'o') //exr
	{
		channels = Nc;
		//nfileName = fn + "exrdemo_" + date + ".exr";
		nfileName = fn + ".exr";

	}
	else return;

	pixels= new float[xres * yres * channels];

	//std::cout << "Writing: " << nfileName <<'\n';

	img::Flip(img_data, Nx, Ny, Nc, pixels, xres, yres, channels);

	std::unique_ptr<ImageOutput> out = ImageOutput::create (nfileName);
	if (! out)
	{
		std::cout<< "Write error \n";
    	return;
	}
	ImageSpec spec(xres, yres, channels, TypeDesc::FLOAT); 
	out->open (nfileName, spec);
	out->write_image (TypeDesc::FLOAT, pixels);
	out->close ();
	//std::cout << "Write successful\n";
	delete [] pixels;
}

long ImgProc::index(int i , int j, int c) const
{
	 return c+(Nc*(i+Nx*j));

}
long ImgProc::indexNon(int i, int j, int c) const
{ 
	return i+(Nx*(j+(Ny*c)));
}

//flip in place
void ImgProc::Flip()
{
	if( img_data == nullptr ){ return; }

	#pragma omp parallel for collapse(2)
	//row to height
	for(int j = 0; j < Ny/2; j++)
	{
		//col to width
		for(int i = 0; i < Nx; i++)
		{				
			for(int k = 0; k < Nc; k++)
			{
				float temp = img_data[index(i,j,k)];

				img_data[index(i,j,k)] = img_data[index(i,Ny-1-j,k)];
				img_data[index(i,Ny-1-j,k)]= temp;

			}
		}
	}

}

//flip as copying
void ImgProc::Flip(const float* data)
{
	if( img_data == nullptr ){ return; }

	#pragma omp parallel for collapse(2)
	for(int j = 0; j < Ny; j++)
	{
		for(int i = 0; i < Nx; i++)
		{
			for(int k = 0; k < Nc; k++)
			{
				img_data[index(i,j,k)] = data[index(i,Ny-1-j,k)];
			}
		}

	}
	return;

}

//flop in place
void ImgProc::Flop()
{
	if( img_data == nullptr ){return;}

	#pragma omp parallel for collapse(2)
	//row to height
	for(int j = 0; j < Ny; j++)
	{
		//col to width
		for(int i = 0; i < Nx/2; i++)
		{				
			for(int k = 0; k < Nc; k++)
			{
				float temp = img_data[index(i,j,k)];

				img_data[index(i,j,k)] = img_data[index(Nx-1-i,j,k)];
				img_data[index(Nx-1-i,j,k)]= temp;

			}
		}
	}
}


ImgProc::ImgProc(const ImgProc& v) :
  Nx (v.Nx), 
  Ny (v.Ny),
  Nc (v.Nc), 
  Nsize (v.Nsize)
{
	img_data = new float[Nsize];
	#pragma omp parallel for
 	for(long i = 0; i < Nsize; i++){ img_data[i] = v.img_data[i]; }
}

ImgProc& ImgProc::operator=(const ImgProc& v)
{
	if(this == &v){ return *this; }
	if(Nx != v.Nx || Ny != v.Ny || Nc != v.Nc)
	{
	    clear();
		Nx = v.Nx;
		Ny = v.Ny;
		Nc = v.Nc;
		Nsize = v.Nsize;
	}
	img_data = new float[Nsize];
	#pragma omp parallel for
	for(long i = 0; i < Nsize; i++){ img_data[i] = v.img_data[i]; }
	return *this;
}





//***************************************************************************************//
//img namespace function definitions  

int img::read_image(const std::string& s, ImgProc& imgProc)
{
	std::string filename = s;
	std::cout << "\nAttempting to find: " << filename << std::endl;
	auto inp = ImageInput::open(filename);
	if (! inp)
	{
		std::cout << "Couldn't find " << filename << std::endl;
    	return 1;
	}
	const ImageSpec &spec = inp->spec();
	int xres = spec.width;
	int yres = spec.height;
	int nchannels = spec.nchannels;
	auto pixels = std::unique_ptr<float[]>(new float[xres * yres * nchannels]);
	inp->read_image(0, 0, 0, nchannels, TypeDesc::FLOAT, &pixels[0]);
	inp->close();

	imgProc.clear(xres,yres,nchannels);
	imgProc.set_raw_OI(&pixels[0]);

	std::cout <<"Image loaded\n\n";
	return 0;
}

//function will not be called unless image has been read
void img::write_image(std::string fileName, char f, const ImgProc& imgProc)
{
	int xres = imgProc.nx();
	int yres = imgProc.ny();
	int channels = 0;
	float* pixels = nullptr;

	std::string nfileName;
	std::size_t pos = fileName.find(".");
    std::string fn = fileName.substr(0,pos);
	time_t now = time(0);
  	std::string date = ctime(&now);
	std::replace(date.begin(), date.end(), ' ', '_');
	std::replace(date.begin(), date.end(), ':', '_');
	date.pop_back(); // remove weird character at the end
	if(f == 'j') ///jpg
	{
		channels=3;
		//save new file to the original filepath
		nfileName = fn + "jpgdemo_" + date + ".jpeg";
	
	}
	else if(f == 'o')
	{
		channels = imgProc.depth();
		nfileName = fn + "exrdemo_" + date + ".exr";
	}
	else return;

	pixels = new float[xres * yres * channels];

	std::cout << "Writing: " << nfileName <<'\n';
	Flip(imgProc.raw(), imgProc.nx(), imgProc.ny(), imgProc.depth(), pixels, xres, yres, channels);

	std::unique_ptr<ImageOutput> out = ImageOutput::create (nfileName);
	if (! out)
	{
		std::cout<< "Write error\n";
    	return;
	}
	ImageSpec spec(xres, yres, channels, TypeDesc::FLOAT); 
	out->open (nfileName, spec);
	out->write_image (TypeDesc::FLOAT, pixels);
	out->close ();
	std::cout << "Write successful\n";
	delete [] pixels;
}

//in place
void img::Flip(ImgProc& in)
{
	#pragma omp parallel for collapse(2)
	//row to height
	for(int j = 0; j < in.ny()/2; j++)
	{
		//col to width
		for(int i = 0; i < in.nx(); i++)
		{				
			for(int k = 0; k < in.depth(); k++)
			{
				float temp = in.raw()[in.index(i,j,k)];

				in.raw()[in.index(i,j,k)] = in.raw()[in.index(i,in.ny()-1-j,k)];
				in.raw()[in.index(i,in.ny()-1-j,k)]= temp;

			}
		}
	}

}

//flip while copying
void img::Flip(const float* in, int nx, int ny, int nc, float* out, int Nx, int Ny, int Nc)
{
	#pragma omp parallel for collapse(2)
	for(int j = 0; j < Ny; j++)
	{
		for(int i = 0; i < Nx; i++)
		{
			for(int k = 0; k < Nc; k++)
			{
				//flip in place
				out[Index(i,j,k,Nc,Nx)] = in[Index(i,ny-1-j,k,nc,nx)];
			}
		}
	}

}

//in place
void img::Flop(ImgProc& in)
{
	#pragma omp parallel for collapse(2)
	//row to height
	for(int j = 0; j < in.ny(); j++)
	{
		//col to width
		for(int i = 0; i < in.nx()/2; i++)
		{				
			for(int k = 0; k <in.depth(); k++)
			{
				float temp = in.raw()[in.index(i,j,k)];

				in.raw()[in.index(i,j,k)] = in.raw()[in.index(in.nx()-1-i,j,k)];
				in.raw()[in.index(in.nx()-1-i,j,k)]= temp;

			}
		}
	}

}




