/*****************************************************************************
** igluWarning.h                                                            **
** ------------                                                             **
**                                                                          **
** A header that defines common functions for printing warning messages in  **
**    IGLU code.                                                            **
**                                                                          **
** Chris Wyman (2/28/2012)                                                  **
*****************************************************************************/

#ifndef IGLU__WARNING_H
#define IGLU__WARNING_H

namespace iglu {

// A full-on warning message similar to an error, except not fatal!
void Warning( const char *error,     // Your error message
			  const char *fileName,  // Usually the preprocessor macro __FILE__
			  const char *funcName,  // Usually the compiler built-in __FUNCTION__
			  int lineNum );         // Usually the preprocessor macro __LINE__

// A simple one-line status message
void StatusMessage( const char *error );

} // end namespace iglu.


#endif