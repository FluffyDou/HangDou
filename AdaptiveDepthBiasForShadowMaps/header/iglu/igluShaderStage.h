/**************************************************************************
** igluShaderStage.h                                                     **
** -----------------                                                     **
**                                                                       **
** This header implements a class that encapsulates a single GLSL shader **
**   stage (e.g., one glsl file) and associated state.                   **
**                                                                       **
** Chris Wyman (06/14/2011)                                              **
**************************************************************************/

#ifndef __IGLU_SHADER_STAGE_H
#define __IGLU_SHADER_STAGE_H

// Other needed headers
#include <assert.h>

namespace iglu {


// Defines input types for the IGLUShaderStage constructor.
//    Either SHADER_FROM_FILE *or* SHADER_FROM_STRING can be
//    OR'd together with *one* of the other types.  The default
//    is SHADER_FROM_FILE.  If multiple stages are used, the
//    result is undefined (probably lowest enum value used)
enum IGLUShaderStageType {
	IGLU_SHADER_FROM_FILE    = 0x0000,
	IGLU_SHADER_FROM_STRING  = 0x0001,
	IGLU_SHADER_VERTEX       = 0x0002,
	IGLU_SHADER_GEOMETRY     = 0x0004,
	IGLU_SHADER_FRAGMENT     = 0x0008,
	IGLU_SHADER_TESS_CONTROL = 0x0010,
	IGLU_SHADER_TESS_EVAL    = 0x0020,
};

enum IGLUSemanticType {
	// The user tried an unknown semantic type!?
	IGLU_ATTRIB_UNKNOWN   = 5,

	// The user specified a type of known vertex attribute
	IGLU_ATTRIB_VERTEX    = 0,
	IGLU_ATTRIB_NORMAL    = 1,
	IGLU_ATTRIB_TEXCOORD  = 2,
	IGLU_ATTRIB_MATL_ID   = 3,
	IGLU_ATTRIB_OBJECT_ID = 4,

	// The user specified a type of output color buffer
	IGLU_MRT_COLOR0      = 6,
	IGLU_MRT_COLOR1      = 7,
	IGLU_MRT_COLOR2      = 8,
	IGLU_MRT_COLOR3      = 9,
	IGLU_MRT_COLOR4      = 10,
	IGLU_MRT_COLOR5      = 11,
	IGLU_MRT_COLOR6      = 12,
	IGLU_MRT_COLOR7      = 13,
	IGLU_MRT_COLOR8      = 14,
	IGLU_MRT_COLOR9      = 15,
	IGLU_MRT_COLOR10     = 16,
	IGLU_MRT_COLOR11     = 17,
	IGLU_MRT_COLOR12     = 18,
	IGLU_MRT_COLOR13     = 19,
	IGLU_MRT_COLOR14     = 20,
	IGLU_MRT_COLOR15     = 21,
};

// If you add more types to IGLUSemanticType than currently defined (12, or [0-11]) 
//     please change the definition of ___IGLU_SEM_TYPE_COUNT.
#define ___IGLU_SEM_TYPE_COUNT   22


class IGLUShaderStage
{
public:
	// Default constructor & destructor
	IGLUShaderStage( uint type, const char *inputShader, bool verbose=true );
	virtual ~IGLUShaderStage();

	// Returns the GL shader handle/identifier for this shader stage.
	inline GLuint GetShaderID( void ) const                   { return m_shaderID; }

	// Check if this shader is setup or if it needs to be loaded.
	inline bool IsCompiled( void ) const                      { return m_isCompiled; }

	// Reload/reinitializes the shader.  If the shader came from an input string
	//    (instead of a file), this returns trivially.
	int ReloadShader( void );

	// Gets the name for a particular semantic type.
	const char *GetSemanticVariableName( IGLUSemanticType type = IGLU_ATTRIB_VERTEX );

	// Get shader-request GL state settings
	inline uint GetShaderRequestedGLEnables( void ) const     { return m_internalEnables; }
	inline uint GetShaderRequestedGLDisables( void ) const    { return m_internalDisables; }

	// A pointer to a IGLUShaderProgram could have type IGLUShaderProgram::Ptr
	typedef IGLUShaderStage *Ptr;

private:
	bool    m_isCompiled;        // Has this shader been compiled?  (And was it error free?)
	bool    m_verbose;           // Print lots of error messages?	
	bool    m_fromFile;          // Did we load from a file or directly from a string?
	char   *m_stageInput;        // Copy of either the filename/input string
	GLuint  m_shaderID;          // The OpenGL handle (e.g., return from glCreateShader())
	GLenum  m_shaderType;        // Is this a vertex, geom, frag, tess shader?
	uint    m_internalEnables;   // IGLUShaderState flags that the GLSL shader itself asked to enable
	uint    m_internalDisables;  // IGLUShaderState flags that the GLSL shader itself asked to disable

	// IGLU extends GLSL with semantic bindings for in/out vars
	char   *m_semanticNames[___IGLU_SEM_TYPE_COUNT]; 

	// Two helper functions to load the specified file, return it as a string.  
	//    Recursion allows "#include's" inside GLSL (not part of GLSL spec... currently)
	char *ReturnFileAsString( char *filename, int recurseDepth=0 ); 
	char *ReturnProcessedShaderString( char* inputFilename, char* buffer, int recurseDepth );

	// If necessary, loads shader source from a file.  Compiles the source
	int SetupShader( void );

	// Some helper functions that print out errors (if verbose) and return 
	//    appropriate error codes.
	int PrintCompilerError( void );
	int PrintFileIOError( void );

	// Takes a character string and sees if it matches one of our semantic types
	uint GetSemanticType( char *type );

	// Takes a character string and sees if it matches the names of one of our IGLUShaderState flags
	uint GetShaderStateFlag( char *type );
};



// End namespace iglu
}

#endif