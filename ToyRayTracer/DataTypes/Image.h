/******************************************************************/
/* Image.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that stores a color buffer.    */
/*     This is particularly useful for offline rendering, and     */
/*     the class includes a very simple routine to save the image */
/*     as a PPM file.                                             */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IMAGE_H
#define IMAGE_H 1

#include <GL/glew.h>
#include <GL/glut.h>

#include "DataTypes/Color.h"
#include "DataTypes/SampleVector.h"

/* An image is a collection of xres*yres Colors.
**
** Beware!!  This allocates memory -> not thread-safe!  However,
**    presumably, you allocate this ONCE at the beginning of the
**    execution.  However, if this is not the case, you should
**    fix this thread-safety issue!
*/
class Image 
{
    Color** image;
    int xres, yres;

    int             channelNum;
    GLuint          texture;
    unsigned char   *canvas_UB;
    float           *buffer;

public:

    SampleVector    *vecPool;

public:

    Image(int xres, int yres);
    Image(int xres, int yres, int channelNum);
    ~Image();

    inline int GetXRes() const                  { return xres; }

    inline int GetYRes() const                  { return yres; }

    inline int GetChanNum() const               { return channelNum; }

    inline GLuint GetTexture()                  { return texture; }

    inline float* & GetBuffer()                 { return buffer; }

    inline unsigned char* GetCanvas()           { return canvas_UB; }

    inline unsigned char& operator()(int x, int y, int offset)  { return canvas_UB[x*channelNum + offset + y*xres*channelNum]; }

    inline float& Buffer(int x, int y, int offset)              { return buffer[x*channelNum + offset + y*xres*channelNum]; }
    //inline Color& operator()(int x, int y)      { return image[y][x]; }

    // update the texture on GPU memory
    void updateTex(); 

    // choose nearest or linear interpolation property for the texture
    void updateTex(const int interpolateFlag);

    void Save(char* file);
};

#endif

