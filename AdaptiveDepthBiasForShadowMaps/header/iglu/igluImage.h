/******************************************************************/
/* igluImage.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* The file defines an image class that reads images of various   */
/*     kinds into a uniform internal structure.                   */
/*                                                                */
/* Chris Wyman (10/03/2011)                                       */
/******************************************************************/

#ifndef IGLU_IMAGE_H
#define IGLU_IMAGE_H

#pragma warning( disable: 4996 )

namespace iglu {

class IGLUImage
{
public:
	// Loads an image from a file
    IGLUImage( char *filename );

	// Loads an image from a memory buffer.  The memory buffer is *not* copied.
	//    IGLUImage can take possession and free() memory (note: not delete!) on exit.
	IGLUImage( unsigned char* image, int width, int height, bool useAlpha=false, bool freeMemory=false );
    ~IGLUImage();

	// Return a pointer to the image data as an unsigned char array.
	inline const unsigned char *ImageData( void ) const   { return m_imgData; }

	// Check if the image was loaded correctly.
	inline bool IsValid() const	                   { return (m_width > 0) && (m_height > 0); }

	// Get size of the texture
	inline int GetWidth()  const                   { return m_width;  }               
	inline int GetHeight() const                   { return m_height; }     

	// Get the OpenGL internal format for this image
	inline unsigned int GetGLFormat() const        { return m_glFormat; }
	inline unsigned int GetGLDatatype() const      { return m_glDatatype; }

	// This is really a test, and might not stay here long
	void SaveAsPPM( char *ppmFile );

	// A pointer to a IGLUImage could have type IGLUImage::Ptr
	typedef IGLUImage *Ptr;

protected:
	int m_width, m_height, m_numComponents;
	unsigned int m_glFormat, m_glDatatype;
	unsigned char *m_imgData;
	bool m_freeMemory;
	char *m_qualifiedFileName;

	bool LoadPPM ( char *filename );
	bool LoadRGB ( char *filename );
	bool LoadBMP ( char *filename );
	bool LoadJPEG( char *filename );
	bool LoadPFM ( char *filename );
	bool LoadPNG ( char *filename );
	bool LoadTGA ( char *filename );
};



// End namespace iglu
}


#endif
