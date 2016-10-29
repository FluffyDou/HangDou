/******************************************************************/
/* igluVideo.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* The file defines an video class that reads various types of    */
/*     video files.  This class is a (currently) crude wrapper    */
/*     around the GPL FFMpeg libraries that allows you to parse   */
/*     a video file essentially one frame at a time.              */
/*                                                                */
/* To use this class in your program, you will need to link to    */
/*     the FFMpeg libraries.  These include the following (though */
/*     not all may be needed):                                    */
/*        Linux:  -lavcodec -lavdevice -lavfilter -lavformat      */
/*                -lavutil -lpostproc -lswresample -lswscale      */
/*        Windows:  avcodec.lib, avdevice.lib, avfilter.lib,      */
/*                  avformat.lib, avutil.lib, postproc.lib,       */
/*                  swresample.lib, swscale.lib                   */
/*     We tested this class against FFMpeg version 0.8.5.  It has */
/*     been known to change APIs quite dramatically, so if you    */
/*     use significantly later (or earlier) libraries, this may   */
/*     not be usable!                                             */
/*                                                                */
/* Chris Wyman (10/06/2011)                                       */
/******************************************************************/


#ifndef IGLU__VIDEO_H
#define IGLU__VIDEO_H

#pragma warning( disable: 4996 )

namespace iglu {

class IGLUVideoIOData;

class IGLUVideo
{
public:
    IGLUVideo( char *filename );
    ~IGLUVideo();

	// Allows you to move around the video.  One of these must be called prior
	//     to calling GetFrameData()!
	bool GetNextFrame( void );
	bool GetPrevFrame( void );
	bool SeekToFrame( unsigned int frameNum );
	bool SeekToTime (        float secFromStart );

	// Return a pointer to the image data as an unsigned char array.
	inline const unsigned char *GetFrameData( void ) const { return IsValid() ? m_imgData : 0; }
	//const unsigned char *GetFrameData( void ) const;

	// Check if the video was at least opened correctly.
	inline bool IsValid() const	                           { return (m_width > 0) && (m_height > 0); }

	// Get size of the video
	inline int GetWidth()  const                           { return m_width;  }               
	inline int GetHeight() const                           { return m_height; }     

	// Get the OpenGL internal format to use for this video
	inline unsigned int GetGLFormat() const                { return m_glFormat; }
	inline unsigned int GetGLDatatype() const              { return m_glDatatype; }

	// This is really a test, and might not stay here long
	void SaveFrameAsPPM( char *ppmFile );

	// A pointer to a IGLUVideo could have type IGLUVideo::Ptr
	typedef IGLUVideo *Ptr;

protected:
	unsigned int m_curFrame;
	int m_width, m_height;
	unsigned int m_glFormat, m_glDatatype;
	unsigned char *m_imgData;
	int m_videoError;
	char *m_qualifiedFileName;

	int OpenVideo( char *filename );

	// An internal containter structure that stores all the needed FFMpeg data structures
	IGLUVideoIOData *m_vidData;
};



// End namespace iglu
}


#endif