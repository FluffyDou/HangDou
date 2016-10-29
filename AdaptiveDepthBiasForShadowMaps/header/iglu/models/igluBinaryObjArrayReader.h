/******************************************************************/
/* igluBinaryObjArrayReader.h                                     */
/* -----------------------                                        */
/*                                                                */
/* Defines a reader for our custom binary object files and their  */
/*     associated (.obj-style) materials                          */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_BINARY_OBJ_READER_H
#define IGLU_BINARY_OBJ_READER_H

#pragma warning( disable: 4996 )

#include "iglu/igluArray1D.h"
#include "igluModel.h"
#include "igluOBJReader.h"
#include "iglu/glstate/igluVertexArrayObject.h"

namespace iglu {


struct IGLUOBJTri;
class  IGLUBuffer;

class IGLUBinaryObjArrayReader : public IGLUModel
{
public:
	// Constructor reads from the file.
	IGLUBinaryObjArrayReader( char *filename, int params=IGLU_OBJ_DEFAULT_STATE );
	virtual ~IGLUBinaryObjArrayReader();

	// Implementation of the necessary IGLUModel virtual methods
	virtual int Draw( IGLUShaderProgram::Ptr &shader );
	virtual int Draw( void );

	virtual bool HasVertices ( void ) const			{ return m_hasVertices; }
	virtual bool HasTexCoords( void ) const         { return m_hasTexCoords; }
	virtual bool HasNormals  ( void ) const         { return m_hasNormals; }
	virtual bool HasMatlID   ( void ) const         { return m_hasMatlID; }
	virtual bool HasObjectID ( void ) const			{ return m_hasObjectID; } 

	// A pointer to a IGLUOBJReader could have type IGLUOBJReader::Ptr
	typedef IGLUBinaryObjArrayReader *Ptr;

protected:
	// When drawing we need to bind the vertex buffer to the appropriate 
	//    shader inputs, so we need to know which shader to use.  However, there
	//    is overhead in the setup, so if the same shader is used repeatedly on
	//    the same object, we can reuse the settings from last time.
	GLuint m_shaderID;

	// Our buffer(s) storing the vertex array for this geometry
	IGLUVertexArray::Ptr m_vertArr;

	// Locations for material file(s)
	IGLUArray1D<char *> m_objMtlFiles;

	// Size of the array
	int m_numVertices, m_numElems;

	// Type of OpenGL primitives
	GLenum m_primType;

	// Information about the OBJ file we read
	bool m_hasVertices, m_hasNormals, m_hasTexCoords, m_hasMatlID, m_hasObjectID;

	// What is the stride of the data?
	GLuint m_vertStride, m_vertOff, m_normOff, m_texOff, m_matlIdOff, m_objectIdOff;

	// Load the material file?
	bool m_loadMtlFile;

private:
	bool LoadBinaryBlob( char *file );

	enum {
		HAS_VERTEX   = 0x10,
		HAS_NORMAL   = 0x01,
		HAS_TEXCOORD = 0x02,
		HAS_MATLID   = 0x04,
		HAS_OBJID    = 0x08
	};

	typedef struct {
		uint vboVersionID;       // Magic number / file version header
		uint numVerts;           // Number of vertices in the VBO
		uint numElems;           // Number of elements (i.e., indices in the VBO)
		uint elemType;           // E.g. GL_TRIANGLES
		uint vertStride;         // Number of bytes between subsequent vertices
		uint vertBitfield;       // Binary bitfield describing valid vertex components 
		uint matlOffset;         // In vertex, offset to material ID
		uint objOffset;          // In vertex, offset to object ID
		uint vertOffset;         // In vertex, offset to vertex 
		uint normOffset;         // In vertex, offset to normal 
		uint texOffset;          // In vertex, offset to texture coordinate
		char matlFileName[84];   // Filename containing material information
	} IGLUBinaryVBOHeader;
};


// End namespace iglu
}



#endif

