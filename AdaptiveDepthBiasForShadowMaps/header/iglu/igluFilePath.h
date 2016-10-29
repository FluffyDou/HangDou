/******************************************************************/
/* igluFilePath.h                                                 */
/* -----------------------                                        */
/*                                                                */
/* A header for interfacing with a simplistic filepath manager    */
/*     and an IGLU version of fopen() for loading files from      */
/*     specified path directories.                                */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_SIMPLE_FILE_PATH_H
#define IGLU_SIMPLE_FILE_PATH_H

#include <cstdio>
#include "iglu/igluArray1D.h"

namespace iglu {

// A generic file path class
class IGLUFilePath {
public:
	// Create a file path.  All searches, by default, only look in the current directory.
	IGLUFilePath();
	~IGLUFilePath();

	// Adds a new search path (<dirName>) into the specified search path
	//    Other than the current directory (""), paths need to end with a 
	//    trailing slash (or, optionally on Windows, backslash).
	void AddPath( const char *dirName, int pathType=COMMON );

	// Checks if a file exists in the specified path.  Returns non-zero if the file exists.
	int FileExists( const char *fileName, int pathType=COMMON );

	// fopen()'s a specified file in the specfied path
	FILE *FileOpen( const char *fileName, const char * fileMode, int pathType=COMMON );

	// Pointer to a fully-qualified filename in the specified path.  
	//    -> If no file exists, NULL is returned.
	//    -> Pointer need not be free'd (it's internal object memory)
	//    -> Data in returned pointer will change when class methods are called
	//       (i.e., copy it if you'll need it long term)
	char *FileInPath( const char *fileName, int pathType=COMMON );

	// Types of file paths.  This allows different open commands to search different
	//    sets of paths.  All searches look in the common path, though searches into
	//    a more specific path look through more specific (e.g., model) paths prior
	//    to the common path.
	enum { COMMON=0, MODEL, IMAGE, VIDEO, SHADER, SCENE, CONFIG, 
		   MAX_PATH_TYPES };

	// A pointer to a IGLUFilePath could have type IGLUFilePath::Ptr
	typedef IGLUFilePath *Ptr;

private:
	IGLUArray1D< char * > m_searchPaths[MAX_PATH_TYPES];
	char                  m_currentFile[1024];

	bool Exists( const char *file );
};


// Specific functions to use the IGLU built in path
void AddSearchPath  ( const char *dirName,  int pathType=IGLUFilePath::COMMON );
bool FileExitsInPath( const char *fileName, int pathType=IGLUFilePath::COMMON );
FILE *FileOpen      ( const char *fileName, const char *fileMode, int pathType=IGLUFilePath::COMMON );

// Returns the filename in the specified path.  Pointer need not be free'd,
//    but may not remain valid after calling other functions in this header.
char *FindFileInPath( const char *fileName, int pathType=IGLUFilePath::COMMON );


} // End namespace iglu

#endif