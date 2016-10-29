/**************************************************************************
** igluShaderProgram.h                                                   **
** -----------------                                                     **
**                                                                       **
** This header includes a class that can be used as an encapsulation for **
**   igluGLSLPrograms, that provides slightly nicer interfaces for       **
**   setting uniform variables, etc.  It builds on a GLSLProgram, and is **
**   essentially a bit of a glorified smart pointer.                     **
**                                                                       **
** We use a smart pointer for the internal storage to point to the       **
**   underlying GLSL class, so you should be able to assign the shader   **
**   program like pointer, and pass them around without worrying about   **
**   accidentally deleting your shader when temporaries go out of scope. **
**   Please note:  I have not tested this extensively!  You may want to  **
**   avoid abusing this ability.                                         **
**                                                                       **
** NOTE: This class is NOT thread safe.  Do not share among threads!     **
**                                                                       **
** Chris Wyman (06/14/2011)                                              **
**************************************************************************/

#ifndef __IGLU_SHADER_PROGRAM_H
#define __IGLU_SHADER_PROGRAM_H

// Other needed headers
#include "helpers/igluCountedPtr.h"
#include "igluShaderVariable.h"
#include "igluShaderStage.h"
#include <assert.h>

namespace iglu {

class IGLUShaderProgram
{
	friend class IGLUShaderVariable;
public:
	// Default constructor & destructor
	IGLUShaderProgram();
	~IGLUShaderProgram();

	// Constructors that also call Load()
	IGLUShaderProgram( const char *vShaderFile, const char *fShaderFile );
	IGLUShaderProgram( const char *vShaderFile, const char *gShaderFile, const char *fShaderFile );

	// Create a shader based on a set of files
	void Load( const char *vShaderFile, const char *fShaderFile );
	void Load( const char *vShaderFile, const char *gShaderFile, const char *fShaderFile );

	// Create a shader based on a literal string of the shaders
	void CreateFromString( const char *vShader, const char *fShader );
	void CreateFromString( const char *vShader, const char *gShader, const char *fShader );

	// Enable and disable the associated shader
	//    -> When enabled, you can render or set variables
	bool Enable( void );
	bool Disable( void );

	// Setting shader variables outside an Enable/Disable pair is slowish, as the class
	//    needs to Enable() and Diable() the shader in the assignment operators.  These
	//    two methods are a quick Enable()/Disable() that can be used around variable
	//    assignments to help reduce state.
	// Note: These are optional methods only for improving performance.  Do not use inside
	//    an Enable()/Disable() pair!
	bool BeginVariableUpdates( void );
	bool EndVariableUpdates( void );

	// Check if this shader is setup or if it needs to be loaded.
	inline bool IsValid( void ) const                       { return m_isLinked; }

	// Is this shader enabled?
	inline bool IsEnabled( void ) const                     { return m_currentlyEnabled; }

	// Get the program's OpenGL ID/handle
	inline GLuint GetProgramID( void ) const                { return m_programID; }

	// Whenever this shader program is enabled or disabled, you can 
	//     automatically call glEnable() or glDisable() on certain state.
	void SetProgramEnables ( uint stateFlags );
	void SetProgramDisables( uint stateFlags );

	// Link or relink the program.  Usually, you do not need to call this, unless
	//     you change GL state that the OpenGL spec notes "requires a relink" to update
	int Link( void );

	// Reload any shader files associated with this GLSL shader.  If the shader
	//     was created from a string, reload works correctly, but reloads from the *same*
	//     string, so it is effectively a very expensive no-op.
	bool Reload( void );

	// For setting the values of uniform variables, you can use [] to select 
	//    the variable, then the assignment operator to assign the value.
	IGLUShaderVariable &operator[] ( const char *varName );

	// For getting the name of attribute/semantic variables use [], for
	//    instance:  const char *glslVertAttrib = shader[ IGLU_ATTRIB_VERTEX ]
	inline const char *operator[] ( uint attribType ) const { return m_semanticNames[attribType]; }

	// If you are using transform feedback, you need to tell OpenGL which vertex/geometry
	//    shader outputs should be piped to the transform feedback buffer.  There are two
	//    methods here.  The first one if you will only be outputting one variable (which takes
	//    just the string name of the GLSL variable).  The second one should be used if you
	//    output multiple variables from your shaders.  These two methods should NOT be mixed 
	//    for any given IGLUShaderProgram!
	// Note: These should occur once, upon shader instantiation, as they relink (think: slow!)
	// The return value is the return from Link().
	int SetTransformFeedbackVaryings( const char *glslVarName );
	int SetTransformFeedbackVaryings( int numOutputs, const char **glslVarNames, 
	                                  GLenum bufferMode=GL_INTERLEAVED_ATTRIBS );

	// A pointer to a IGLUShaderProgram could have type IGLUShaderProgram::Ptr
	typedef IGLUPtr<IGLUShaderProgram,IGLUShaderVariable&,const char *> Ptr;


private:
	IGLUShaderVariable  m_tmpVar;    // This temp is an issue for multi-threading (if multiple threads touch this class)
	IGLUArray1D<IGLUShaderTexture *> m_activeTex;      // A list of the textures the user has assigned this shader.
	IGLUArray1D<uint>                m_activeTexFlags; // Flags describing properties of the active textures/bound variables/units
	IGLUArray1D<IGLUShaderTexture *> m_activeImage; //A list of textures the user has assigned this shader as images. 
	// A list of shaders that are linked into this program.
	IGLUArray1D<IGLUShaderStage::Ptr> m_shaderStages;

	// Some OpenGL state for this GLSL program
	GLuint m_programID;
	bool   m_isLinked, m_verbose;

	// When the shader is enabled or disabled, certain state may be enabled,
	//     certain state may be disabled, all other enables/disables are left alone.
	uint m_invokeStateEnable, m_invokeStateDisable, m_invokeStateMask;
	uint m_pushedState;

	// IGLU extends GLSL to attach semantic meaning to certain inputs and outputs (in/out)
	//    This is a list of character strings for the program variables matching these
	//    semantic meanings.
	const char *m_semanticNames[___IGLU_SEM_TYPE_COUNT]; 

	// If our program is not enabled and we try to assign shader 
	//     variables, OpenGL fails to assign values somewhat silently
	//     (and might actually change the variables in another program).
	//     We don't want that behavior, so we track if we're enabled, 
	//     and if not, we enable the correct shader to assign a 
	//     variable (then switch back to the currently enabled program)
	bool  m_currentlyEnabled, m_assigningVariables;
	GLint m_prevProgram;
	void PushProgram( void );  // Pushes the current enabled program to enable this one
	void PopProgram( void );   // Pops this one (if appropriate) to go back to the last shader enabled

	// Sometimes, we just need to call glUseProgram(), but don't care about the glEnable() state
	//   (e.g., when we're setting shader variables outside an Enable()/Disable() pair)
	//   -> PartialPushProgram() and PartialPopProgram() do not modify the enable state, but simply
	//      switch to a different program
	void PartialPushProgram( void );
	void PartialPopProgram( void );

	// When we Enable()/Disable() this program this function is called 
	//     to set the appropriate OpenGL state based on the m_invokeState.
	void GetCurrentEnabledState( void );
	void SetInvocationState( int startEnd );
	bool IsStateEnabled( uint bitVector, uint stateFlag );
	bool IsStateDisabled( uint bitVector, uint stateFlag );
	
	// What happens if our link fails?
	int PrintLinkerError( void );

	// Our internal shader stages may have IGLU-specific semantic names.  Copy the appropriate
	//    pointers so we know what they are, and can make them available directly.
	void CopySemanticNames( void );
};


enum IGLUShaderState {
	IGLU_GLSL_DEFAULT_STATE      = 0x0000,
	IGLU_GLSL_BLEND              = 0x0001,
	IGLU_GLSL_DEPTH_TEST         = 0x0002,
	IGLU_GLSL_STENCIL_TEST       = 0x0004,
	IGLU_GLSL_SCISSOR_TEST       = 0x0008,
	IGLU_GLSL_CULL_FACE          = 0x0010,
	IGLU_GLSL_VARY_POINT_SIZE    = 0x0020,
	IGLU_GLSL_RASTERIZE_DISCARD  = 0x0040,
	IGLU_GLSL_LINE_SMOOTH        = 0x0080,
	IGLU_GLSL_MULTISAMPLE        = 0x0100,
	IGLU_GLSL_LOGIC_OP           = 0x0200,
	IGLU_GLSL_COLOR_LOGIC_OP     = 0x0400,
	IGLU_GLSL_SAMPLE_SHADING     = 0x0800,
	IGLU_GLSL_SAMPLE_MASK        = 0x1000,
};

enum IGLUShaderActiveTexFlags {
	IGLU_SHADER_TEX_DEFAULT      = 0x0000,
	IGLU_SHADER_TEX_IS_SHADOW    = 0x0001,
};


// End namespace iglu
}

#endif
