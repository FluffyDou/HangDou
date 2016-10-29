/******************************************************************/
/* igluVideoCapture.h                                             */
/* -----------------------                                        */
/*                                                                */
/* Chris Wyman (10/06/2011)                                       */
/******************************************************************/


#ifndef IGLU__VIDEO_CAPTURE_H
#define IGLU__VIDEO_CAPTURE_H

#pragma warning( disable: 4996 )

namespace iglu {

class IGLUVideoIOData;

class IGLUVideoCapture
{
public:
    IGLUVideoCapture();
	~IGLUVideoCapture() {}

	bool StartCapture( char* filename, int width, int height, int codec=0, 
		               int bitrate=1000000, int format=0, int frameRate=30 );
	bool FinishCapture( void );
	bool CaptureFrame( unsigned char *frame, int width, int height );

	unsigned char *GetMemoryPtr( void )            { return m_imgDataRGB; }

	void ListSupportedFormatsAndCodecs( void );

	// A pointer to a IGLUVideo could have type IGLUVideo::Ptr
	typedef IGLUVideoCapture *Ptr;

protected:
	unsigned int m_curFrame;
	int m_width, m_height, m_bitrate, m_framerate;
	
	bool m_isOpen;

	int            m_tmpFrameBufSize;
	unsigned char *m_tmpFrameBuf;
	unsigned char *m_imgData, *m_imgDataRGB;

	// An internal containter structure that stores all the needed FFMpeg data structures
	IGLUVideoIOData *m_vidData;

	bool InitializeVideoFrame( void );
};


enum {
	IGLU_FORMAT_AVI            = 0x000000, // '.avi'
	IGLU_FORMAT_FLV            = 0x000001, // '.flv'
	IGLU_FORMAT_H263           = 0x000002, // '.h263'
	IGLU_FORMAT_H264           = 0x000003, // '.h264'
	IGLU_FORMAT_IPOD           = 0x000004, // '.ipod' 
	IGLU_FORMAT_MPEG4          = 0x000005, // '.mp4'
	IGLU_FORMAT_MJPEG          = 0x000006, // '.mjpeg'
	IGLU_FORMAT_MPEG           = 0x000007, // '.mpeg'
	IGLU_FORMAT_PSP            = 0x000008, // '.psp'
	IGLU_FORMATS_MAX           = 0x000009
};

enum {
	IGLU_CODEC_XVID            = 0x000000, // 'libxvid'
	IGLU_CODEC_MPEG4           = 0x000001, // 'mpeg4'
	IGLU_CODEC_MS_MPEG4        = 0x000002, // 'msmpeg4'
	IGLU_CODEC_H264            = 0x000003, // 'libx264'
	IGLU_CODEC_FLV             = 0x000004, // 'flv' 
	IGLU_CODEC_H263            = 0x000005, // 'h263'
	IGLU_CODEC_MJPEG           = 0x000006, // 'mjpeg'
	IGLU_CODEC_MPEG2           = 0x000007, // 'mpeg2video'
	IGLU_CODEC_WMV8            = 0x000008, // 'WMV2'
	IGLU_CODEC_WMV9            = 0x000009, // 'WMV3'
	IGLU_CODECS_MAX            = 0x00000A
};


// End namespace iglu
}


#endif