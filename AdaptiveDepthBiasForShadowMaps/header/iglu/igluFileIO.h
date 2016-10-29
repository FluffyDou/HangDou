/******************************************************************/
/* igluFileIO.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* A class with (a/some/a buncf of) static methods relating to    */
/*     file IO.                                                   */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_FILE_IO_H
#define IGLU_FILE_IO_H

#include <cstdio>
#include "iglu/igluArray1D.h"

namespace iglu {

class IGLUFileIO {
private:
	IGLUFileIO();
	~IGLUFileIO();

public:
	typedef enum {
		TYPE_UNKNOWN = 0,
		TYPE_OBJ,
		TYPE_OBJAR,
		TYPE_PNG,
		TYPE_BMP,
		TYPE_PPM,
		TYPE_RGB,
		TYPE_TGA,
		TYPE_JPG,
		TYPE_PFM
	} FileType;

public:
	static FileType GetFileType( char *fName );
};





} // End namespace iglu

#endif