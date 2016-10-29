/*****************************************************************************
** igluErrors.h                                                             **
** ------------                                                             **
**                                                                          **
** A header that defines numerous common error codes used in various parts  **
**    of the IGLU library.  At some point, there may also be functions to   **
**    convert the symbols into nice human-readable error messages.          **
**                                                                          **
** Chris Wyman (9/26/2011)                                                  **
*****************************************************************************/


#ifndef IGLU__ERRORS_H
#define IGLU__ERRORS_H

namespace iglu {
	
	// A list of error codes used in IGLU
	enum {
		// No error -- should be first
		IGLU_NO_ERROR = 0,

		// File IO errors
		IGLU_ERROR_FILE_OPEN_FAILED,

		// Unsupported file type errors
		IGLU_ERROR_UNSUPPORTED_VIDEO_FILE,
		IGLU_ERROR_UNSUPPORTED_IMAGE_FILE,

		// Memory allocation errors
		IGLU_ERROR_ALLOCATING_TEMP_MEMORY,

		// Video specific error messages
		IGLU_ERROR_NO_VIDEO_STREAM,
		IGLU_ERROR_IN_VIDEO_CODEC,

		// For our internal drawing utilities
		IGLU_ERROR_INITIALIZING_DRAW_ROUTINES,
		IGLU_ERROR_UNKNOWN_FONT,

		// The parameters passed to IGLU functions were invalid
		IGLU_ERROR_INVALID_PARAMETERS,

		// The OpenGL state was not set correctly to perform this operation
		IGLU_ERROR_INVALID_GLSTATE,

		// For model drawing
		IGLU_ERROR_NO_GLSL_VERTEX_SEMANTIC,
		IGLU_ERROR_NO_GLSL_NORMAL_SEMANTIC,
		IGLU_ERROR_NO_GLSL_TEXCOORD_SEMANTIC,
		IGLU_ERROR_MODEL_MISSING_DATA,
		IGLU_ERROR_MODEL_CORRUPT,

		// For GLSL compilation
		IGLU_ERROR_GLSL_COMPILE_FAILED,
		IGLU_ERROR_GLSL_LINK_FAILED,

	};


// End iglu namespace
}

#endif