/******************************************************************/
/* igluOBJReader.h                                                */
/* -----------------------                                        */
/*                                                                */
/* Defines a reader for OBJ model files and associated materials  */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_OBJ_READER_H
#define IGLU_OBJ_READER_H

#pragma warning( disable: 4996 )

#include "iglu/igluArray1D.h"
#include "igluModel.h"
#include "iglu/glstate/igluVertexArrayObject.h"

namespace iglu {

enum IGLUModelParams {
	IGLU_OBJ_DEFAULT_STATE      = 0x0000,
	IGLU_OBJ_CENTER             = 0x0001,  // Recenter model around the origin
	IGLU_OBJ_UNITIZE            = 0x0002,  // Resize model so it ranges from [-1..1]
	IGLU_OBJ_COMPACT_STORAGE    = 0x0004,  // Attempt to reuse vertices shared between triangles
	IGLU_OBJ_NO_MATERIALS       = 0x0008,  // Do no load an .mtl file.
	IGLU_OBJ_KEEP_CPU_BUFFER    = 0x0010,  // Keep the temporary copy of the vertex array on the CPU side (created upon load).
	IGLU_OBJ_GUESS_NORMALS      = 0x0020,  // If no file normals, create and use facet normals
	IGLU_OBJ_NO_OBJECTS			= 0x0040,
	IGLU_OBJ_SAVE_BINARY        = 0x0080,  // After loading, save this off as a binary blob (e.g., .objar file)
	IGLU_OBJ_VERBOSE_READ       = 0x0100,  // Outputs debug data when reading
};


struct IGLUOBJTri;
class  IGLUBuffer;

class IGLUOBJReader : public IGLUFileParser, public IGLUModel
{
public:
	// Constructor reads from the file.
	IGLUOBJReader( char *filename, int params=IGLU_OBJ_DEFAULT_STATE );
	virtual ~IGLUOBJReader();

	// Get some data about the object file
	uint GetTriangleCount( void ) const       { return m_objTris.Size(); }

	// Get OpenGL buffers for vertex/triangle/normal data
	IGLUVertexArray::Ptr &GetVertexArray( void )         { return m_vertArr; }

	// Get vertex attrib array data
	uint GetArrayBufferStride( void );

	// Implementation of the necessary IGLUModel virtual methods
	virtual int Draw( IGLUShaderProgram::Ptr &shader );
	virtual int Draw( void );

	virtual bool HasVertices ( void ) const			{ return m_hasVertices; }
	virtual bool HasTexCoords( void ) const         { return m_hasTexCoords; }
	virtual bool HasNormals  ( void ) const         { return m_hasNormals; }
	virtual bool HasMatlID   ( void ) const         { return m_hasMatlID; }
	virtual bool HasObjectID ( void ) const			{ return m_hasObjectID; } 

	// Some meta data about the model bounding box (model space v.s. current)
	inline vec3 GetOriginalMin( void ) const        { return m_origMin; }
	inline vec3 GetOriginalMax( void ) const        { return m_origMax; }
	inline vec3 GetBBoxMin( void ) const            { return m_bboxMin; }
	inline vec3 GetBBoxMax( void ) const            { return m_bboxMax; }


	// A pointer to a IGLUOBJReader could have type IGLUOBJReader::Ptr
	typedef IGLUOBJReader *Ptr;

protected:
	bool m_compactFormat, m_loadMtlFile, m_assignObjects, m_keepCPUArray;
	uint  m_curMatlId;
	uint  m_curObjectId;
	// When drawing we need to bind the vertex buffer to the appropriate 
	//    shader inputs, so we need to know which shader to use.  However, there
	//    is overhead in the setup, so if the same shader is used repeatedly on
	//    the same object, we can reuse the settings from last time.
	GLuint m_shaderID;

	// The original model-space coordinate bounding box
	vec3 m_origMin, m_origMax, m_bboxMin, m_bboxMax;

	// Our buffer(s) storing the vertex array for this geometry
	IGLUVertexArray::Ptr m_vertArr;

	// Basic geometric definitions
	IGLUArray1D<vec3> m_objVerts;
	IGLUArray1D<vec3> m_objNorms;
	IGLUArray1D<vec2> m_objTexCoords;

	// Aux data structure if we don't have built-in normals
	IGLUArray1D<vec3> m_facetNorms;

	// Basic geometric triangles definitions
	IGLUArray1D<IGLUOBJTri *> m_objTris;

	// Locations for material file(s)
	IGLUArray1D<char *> m_objMtlFiles;

	// This is the first .mtl file encountered.  This is used only
	//    to output into a binary VBO object.
	char *m_baseMtlFileName;

	// Array of object names inside the file
	IGLUArray1D<char *> m_objObjectNames;

	// Does the user want us to resize and center the object around the origin?
	//    Some/many models use completely arbitrary coordinates, so it's difficult
	//    to know how big and where they'll be relative to each other
	bool m_resize, m_center;

	// Information about the OBJ file we read
	bool m_hasVertices, m_hasNormals, m_hasTexCoords, m_hasMatlID, m_hasObjectID, m_guessNorms;
	bool m_verboseRead;

	// What is the stride of the data?
	GLuint m_vertStride, m_vertOff, m_normOff, m_texOff, m_matlIdOff, m_objectIdOff;

	// After reading this in, should we save it out as a binary blob?
	bool m_saveAsBinaryVBO;

	// The temporary buffer used to create the GPU vertex array.  Size is in # of floats
	float *m_tmpVertArray;
	int   m_tmpVertArraySz;

private:
	// Indicies in OBJ files may be relative (i.e., negative).  OBJ files also
	//    use a Pascal array indexing scheme (i.e., start from 1).  These methods
	//    fixes these issues to give an array index we can actually use.
	unsigned int GetVertexIndex( int relativeIdx );
	unsigned int GetNormalIndex( int relativeIdx );
	unsigned int GetTextureIndex( int relativeIdx );

	// These methods make sure the vertex array and element array buffers are set up correctly
	void GetArrayBuffer( void );
	void GetElementArrayBuffer( void );
	void GetCompactArrayBuffer( void );

	// Save a binary version of the OBJ file for easy reading later
	bool SaveBinaryVBO( char *outFile );

	// Read appropriate facet tokens
	void Read_V_Token( IGLUOBJTri *tri, int idx, char *token=0 );
	void Read_VT_Token( IGLUOBJTri *tri, int idx, char *token=0 );
	void Read_VN_Token( IGLUOBJTri *tri, int idx, char *token=0 );
	void Read_VTN_Token( IGLUOBJTri *tri, int idx, char *token=0 );

	// Helper function for adding data to a float array for our GPU-packed array
	void AddDataToArray( float *arr, int startIdx, int matlID, int objectID, vec3 *vert, vec3 *norm, vec2 *tex );

	// Helper for centering & resizing the geometry in the array before sending it to OpenGL
	void CenterAndResize( float *arr, int numVerts );

	// This declares an ugly type FnParserPtr that points to one of the Read_???_Token methods
	typedef void (iglu::IGLUOBJReader::*FnParserPtr)(IGLUOBJTri *, int, char *); 

	// Copy from prior triangle facets when triangulating
	void CopyForTriangleFan( IGLUOBJTri *newTri );

	// A method that is called after finding an 'f' line that identifies
	//    the appropriate Read_???_Token() method to call when parsing
	void SelectReadMethod( FnParserPtr *pPtr );

	//Checks that an f line has enough vertices to be parsed as a triangle.
	//For example, it returns false if the f line has only two entries on it
	bool IsValidFLine( FnParserPtr *pPtr);
	
	//Determines if objName has been seen before in this obj file
	int GetObjectID( char *objName );
};


// End namespace iglu
}



#endif

