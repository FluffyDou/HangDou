/******************************************************************/
/* Image.cpp                                                      */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that stores a color buffer.    */
/*     This is particularly useful for offline rendering, and     */
/*     the class includes a very simple routine to save the image */
/*     as a PPM file.                                             */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "Image.h"
#include "Color.h"

#include <stdio.h>
#include <stdlib.h>


Image::Image(int xres, int yres): xres(xres), yres(yres)
{
    image =new Color*[yres];
	for (int i=0;i<yres;i++)
		image[i] = new Color[xres];
}

Image::Image(int xres, int yres, int channelNum): xres(xres), yres(yres), channelNum(channelNum)
{
    // create the buffer for a RGBA texture on CPU
    canvas_UB = new unsigned char[xres * yres * channelNum];

    // bind canvas_UB to a GL 2D texture
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, xres, yres, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0); 
    glBindTexture(GL_TEXTURE_2D, 0);
}

Image::~Image()
{
    if(image)
	{
		for (int i=0;i<yres;i++)
			delete[] image[i];
		delete[] image;
    }

    if(canvas_UB)
    {
        delete[] canvas_UB;
    }

    if(buffer)
    {
        delete[] buffer;
    }
}


void Image::updateTex()
{
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, xres, yres, 0, GL_RGBA, GL_UNSIGNED_BYTE, canvas_UB); 
    glBindTexture(GL_TEXTURE_2D, 0);
}

void Image::updateTex(const int interpolateFlag)
{
    glBindTexture(GL_TEXTURE_2D, texture);
    if( interpolateFlag == 0 )
    {
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    }
    else if( interpolateFlag == 1 )
    {
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    }

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, xres, yres, 0, GL_RGBA, GL_UNSIGNED_BYTE, canvas_UB); 
    glBindTexture(GL_TEXTURE_2D, 0);
}

void Image::Save(char* filename)
{
	FILE *file;
	int count=0;
	int r, g, b;

	file = fopen( filename, "w" );
	if (!file) 
		{ 
			fprintf( stderr, "Error: Unable to write to %s!\n", filename ); 
			return; 
		}

	fprintf( file, "P3\n# Image output from Image::save()!\n%d %d\n%d\n", xres, yres, 255 );

	for (int j=yres-1;j>=0;j--)
		for (int i=0;i<xres;i++)
		{
			r = (int)(image[j][i].Red()*255);
			g = (int)(image[j][i].Green()*255);
			b = (int)(image[j][i].Blue()*255);
			r = ( r>255 ? 255 : (r<0 ? 0 : r) );
			g = ( g>255 ? 255 : (g<0 ? 0 : g) );
			b = ( b>255 ? 255 : (b<0 ? 0 : b) );
			fprintf( file, "%d %d %d ", r, g, b );
			if ((++count % 5) == 0)
				fprintf( file, "\n"); 
		}

	fclose( file );
}


