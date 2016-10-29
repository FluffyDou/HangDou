/******************************************************************/
/* igluCallbackFileParser.h                                       */
/* -----------------------                                        */
/*                                                                */
/* This class provides a callback-based approach for parsing      */
/*    simple text files.  You pass in a filename, setup a list of */
/*    callback functions to execute when encountering lines that  */
/*    with certain tokens, then call Parse().                     */
/*                                                                */
/* Chris Wyman (06/09/2010)                                       */
/******************************************************************/

#ifndef __IGLU_CALLBACK_FILE_PARSER_H
#define __IGLU_CALLBACK_FILE_PARSER_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "igluFileParser.h"
#include "iglu/igluArray1D.h"

namespace iglu {

class IGLUCallbackFileParser : public IGLUFileParser 
{
// Public typedefs
public:
	// A pointer to a IGLUCallbackFileParser could have type IGLUCallbackFileParser::Ptr
	typedef IGLUCallbackFileParser *Ptr;

	// The callbacks for handling parsing of specific data from the file are of type
	//    IGLUCallbackFileParser::Callback, which take a IGLUCallbackFileParser::Ptr to
	//    this object (the parser), which allows one to pull further data from the file.
	typedef void (*Callback)( Ptr );

public:
	// Start parsing by creating a new object with the desired file
	IGLUCallbackFileParser( char *filename, int searchPathType=IGLUFilePath::COMMON );
	virtual ~IGLUCallbackFileParser();

	// Add a callback to the list of token / function pairs
	//   -> Matching with 'token' is not case sensitive.
	void SetCallback( char *token, Callback func );

	// Once the callbacks are setup, parse the file.
	//  -> This also closes the file when done parsing
	void Parse( void );

protected:
	IGLUArray1D< char * >   callbackTokens;
    IGLUArray1D< Callback > callbacks;
};


// End namespace iglu
}
 
#endif


