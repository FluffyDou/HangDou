/******************************************************************/
/* igluFileParser.h                                               */
/* -----------------------                                        */
/*                                                                */
/* This class provides a cleanish parsing interface for various   */
/*   files, allowing easy tokenizing and reading of at least some */
/*   basic IGLU datatype classes directly.                        */
/*                                                                */
/* Chris Wyman (06/09/2010)                                       */
/******************************************************************/

#ifndef __IGLU_FILE_PARSER_H
#define __IGLU_FILE_PARSER_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "igluTextParsing.h"
#include "iglu/igluFilePath.h"

namespace iglu {

class vec2;
class vec3;
class vec4;
class int2;
class int3;
class int4;
class uint2;
class uint3;
class uint4;
class dvec2;
class dvec3;
class dvec4;

class IGLUFileParser {
public:
	IGLUFileParser( char *filename, int searchPathType=IGLUFilePath::COMMON );
	virtual ~IGLUFileParser();

	// Read the next line in the file into an internal buffer.  Discarding blanks
	//    throws ignores blank lines or those containing only comments (starting with '#')
	//    and returns the next line with actual stuff on it.
	// Returns a pointer to the start of the first non-blank character.
	char *ReadNextLine( bool discardBlanks=true );

	// Get a pointer to the as-yet-unprocessed part of the current line
	char *GetCurrentLinePointer( void )						{ return internalBufPtr; }

	// Get a pointer to the beginning of the current line.
	char *GetUnprocessedLinePointer( void )					{ return internalBuf; }

	// Get the next token on the line (i.e., white-space delimited)
	char *GetToken( char *token );

	// Treats the line as a CSV and grabs the next CSV entry
	char *GetCSVEntry( char *token );

	// Get the next token on the line, but force it all to be lower case first
	char *GetLowerCaseToken( char *token );

	// Get the next token on the line, but force it all to be upper case first
	char *GetUpperCaseToken( char *token );

	// Reads a number from the line.  Versions returning a (char *) return a pointer
	//    to the remainder of the current line and store the value in the parameter,
	//    the others return the value directly.
	int		 GetInteger( void );
	unsigned GetUnsigned( void );	
	float	 GetFloat( void );
	vec2	 GetVec2( void );
	vec3	 GetVec3( void );
	vec4	 GetVec4( void );
	int2	 GetIVec2( void );
	int3	 GetIVec3( void );
	int4	 GetIVec4( void );
	uint2	 GetUVec2( void );
	uint3	 GetUVec3( void );
	uint4	 GetUVec4( void );
	dvec2	 GetDVec2( void );
	dvec3	 GetDVec3( void );
	dvec4	 GetDVec4( void );
	double	 GetDouble( void );
	char *	 GetInteger( int *result );
	char *	 GetUnsigned( unsigned *result );
	char *	 GetFloat( float *result );
	char *	 GetDouble( double *result );

	// File line number accessor
	int GetLineNumber( void ) const                         { return lineNum; }

	// Accessors to the underlying file handle
	bool IsFileValid( void ) const							{ return f != NULL; }
	FILE *GetFileHandle( void ) const						{ return f; }

	// Get information about the file
	char *GetFileDirectory( void )							{ return fileDirectory; }
	char *GetFileName( void )								{ return unqualifiedFileName; }
	char *GetQualifiedFileName( void )						{ return fileName; }

	// Get the file's size.  Apparently the approach used to get size could be non-portable to
	//    certain systems, in which case the size will be -1.  This rarely happens.  -1 could
	//    also happen if the file size exceeds a 32-bit limit (currently used in IGLU).
	long GetFileSize( void ) const                          { return fileSize; }

	// Print error messages to standard output.
	void WarningMessage( const char *msg );
	void ErrorMessage( const char *msg );

	// Error messages where the first (msg) parameter includes a %s that
	//    is replaced by the second parameter
	void WarningMessage( const char *msg, const char *paramStr );
	void ErrorMessage( const char *msg, const char *paramStr );

	// Sometimes a line from the file is processed but then we find out that we were
	//    not the correct code to process this line.  Before passing the parser off
	//    to someone else, we should "reset" the line so the other code can start
	//    processing from the beginning of the line.
	void ResetProcessingForCurrentLine( void );

	// A pointer to a IGLUFileParser could have type IGLUFileParser::Ptr
	typedef IGLUFileParser *Ptr;

protected:
	FILE *f;
	bool m_closed;
	int lineNum;
	char *fileName, *unqualifiedFileName, *fileDirectory;
	char internalBuf[ 1024 ], *internalBufPtr;
	long fileSize;

	// A simple call to fgets, storing data internally, and increment our line counter
	char *__ReadLine( void );

	// Derived classes may want to go ahead and close the file when they're ready
	void CloseFile( void );
};


// End namespace iglu
}

#endif


