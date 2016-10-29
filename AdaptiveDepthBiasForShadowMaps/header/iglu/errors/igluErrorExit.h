/*****************************************************************************
** igluErrorExit.h                                                          **
** ------------                                                             **
**                                                                          **
** A header that defines common functions for crashing or printing fatal    **
**    error messages inside of IGLU.                                        **
**                                                                          **
** Chris Wyman (2/23/2012)                                                  **
*****************************************************************************/

#ifndef IGLU__ERROREXIT_H
#define IGLU__ERROREXIT_H

namespace iglu {

// A fatal error message [print error and exit]
void ErrorExit( const char *error,     // Your error message
			    const char *fileName,  // Usually the preprocessor macro __FILE__
				const char *funcName,  // Usually the compiler built-in __FUNCTION__
				int lineNum );         // Usually the preprocessor macro __LINE__

// A non-fatal error message [print warning and return]
void Warning  ( const char *error,     // Your error message
			    const char *fileName,  // Usually the preprocessor macro __FILE__
				const char *funcName,  // Usually the compiler built-in __FUNCTION__
				int lineNum );         // Usually the preprocessor macro __LINE__

// Check if we can start executing OpenGL code yet (i.e., do we have a valid OpenGL context?)
//     -> If yes, silently return true and return to the caller
//     -> If no, print a fatal error about the lack of an OpenGL context.
bool CheckValidOpenGL ( const char *fileName,  // Usually the preprocessor macro __FILE__
		  		        const char *funcName,  // Usually the compiler built-in __FUNCTION__
				        int lineNum );         // Usually the preprocessor macro __LINE__

// Called if there's an error loading file 'fileName' in the function specified.
//     -> Prints an error message and quits
void LoadError ( const char *fileName, const char *inFunction );  

void Warn( const char *error );

} // end namespace iglu.


#endif