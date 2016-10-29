/******************************************************************/
/* igluSimpleObject.h                                             */
/* -----------------------                                        */
/*                                                                */
/* Defines a simple object type that developers can stuff full of */
/*    data by passing in vertex attribute arrays.                 */
/*                                                                */
/* Chris Wyman (10/24/2011)                                       */
/******************************************************************/

#ifndef IGLU_SIMPLEOBJECT_H
#define IGLU_SIMPLEOBJECT_H

#include "igluModel.h"
#include "iglu/glstate/igluVertexArrayObject.h"

namespace iglu {

class IGLUSimpleObject : public IGLUModel 
{
public:
	// Simple constructor / destructor
	IGLUSimpleObject();                                     
	virtual ~IGLUSimpleObject();                          

	// Pass in data to this simple object
	//   -> Required data: array of triangle indicies (3 uints per triangle) 
	//                     array of vertex positions (3 floats per vertex)
	//   -> Optional data: array of vertex normals (3 floats per vertex)
	//                     array of vertex texture coords  (2 floats per vertex)
	//                     array of material IDs (1 int per vertex)
	//                     array of object IDs (1 int per vertex)
	void SetGeometry( unsigned int triCount, unsigned int vertCount,
		              unsigned int *idxBuf, float *vPos, 
		              float *vNorm=0, float *vTexCoord=0, int *vMatlID=0, int *vObjID=0 );

	// Draw routines returning either IGLU_NO_ERROR or an error code
	//    -> One enables a specific shader
	//    -> One assumes the developer enabled an appropriate shader
	virtual int Draw( IGLUShaderProgram::Ptr &shader );
	virtual int Draw( void );

	// Specifies what data is available in the renderable primitive
    virtual bool HasVertices ( void ) const         { return m_hasVerts;     }
	virtual bool HasTexCoords( void ) const         { return m_hasTexCoords; }
	virtual bool HasNormals  ( void ) const         { return m_hasNorms;     }
	virtual bool HasMatlID   ( void ) const         { return m_hasMatlID;    }
    virtual bool HasObjectID ( void ) const		    { return m_hasObjID;     }

	// A pointer to a IGLUSimpleObject could have type IGLUSimpleObject::Ptr
	typedef IGLUSimpleObject *Ptr;

private:
	// Data about what this particular object is capable of
	bool m_hasVerts, m_hasTexCoords, m_hasNorms, m_hasMatlID, m_hasObjID;

	// Counts of the underlying object data
	unsigned int m_numTris, m_numVerts;

	// Stride between adjacent vertices in our vertex array;
	int m_vertStride;
	int m_matlIdOff, m_objIdOff, m_vertOff, m_normOff, m_texOff;

	// Our buffer(s) storing the vertex array for this geometry
	IGLUVertexArray::Ptr m_vertArr;

	// An internal utility to help in creating the vertex array
	void AddDataToArray( float *arr, int startIdx, int matlID, int objectID, float *vert, float *norm, float *tex );
};




// End namespace iglu
}



#endif
