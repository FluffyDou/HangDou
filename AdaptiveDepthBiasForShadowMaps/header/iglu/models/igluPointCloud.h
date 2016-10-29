/******************************************************************/
/* igluPointCloud.h                                               */
/* -----------------------                                        */
/*                                                                */
/* Defines a point cloud format.                                  */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_POINTCLOUD_H
#define IGLU_POINTCLOUD_H

#pragma warning( disable: 4996 )

#include "iglu/igluArray1D.h"
#include "igluModel.h"
#include "iglu/glstate/igluVertexArrayObject.h"

namespace iglu {

enum IGLUPointCloudParams {
	IGLU_PCL_HASVERTEX          = 0x0000,  // Vertices should be 3 floats and should be the 1st data in every vertex
	IGLU_PCL_HASNORMAL          = 0x0001,  // Normals should be 3 floats and come after its corresponding vertex data
	IGLU_PCL_HASTEXCOORD        = 0x0002,  // Texture coordinates should be 2 floats and come after its corresponding normal (if included)
	IGLU_PCL_HASMATLID          = 0x0004,  // Material IDs should be 1 float and come after its corresponding tex coord (if included)
	IGLU_PCL_HASOBJID           = 0x0008,  // Object IDs should be 1 float and come after its corresponding material ID (if included)
	IGLU_PCL_WASRESIZED         = 0x1000,
	IGLU_PCL_WASCENTERED        = 0x2000,
	IGLU_PCL_CONVERTTOHALF      = 0x4000,  // If passed to constructor with an array of floats, will convert these to halves internally
};


class IGLUPointCloud : public IGLUModel
{
	// The file header that stores information about the model to load
	struct IGLUCloudHeader;

public:
	// Constructor reads from the file.
	IGLUPointCloud( const char *filename );

	// Constructor reads from an array of floats
	IGLUPointCloud( int numPoints, float *pointData, int formatFlags=IGLU_PCL_HASVERTEX );

	// Constructor reads from an array of half-precision floats (as ushorts) -- NOTE: Experiemental.  May not work!
	IGLUPointCloud( int numPoints, unsigned short *pointData, int formatFlags=IGLU_PCL_HASVERTEX );

	// Destructor
	virtual ~IGLUPointCloud();

	// Get some data about the object file
	uint GetPointCount( void ) const                     { return m_pointCount; }

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

	// Save this to a point cloud file
	bool Save( const char *filename );

	// A pointer to a IGLUOBJReader could have type IGLUOBJReader::Ptr
	typedef IGLUPointCloud *Ptr;

private:
	// How large are the components of each component of our point-cloud data?
	inline int GetElemSz( void ) const              { return m_useHalf ? sizeof(short) : sizeof(float); }

	// Actually go ahead and load the point cloud.
	int LoadCloud( const char *filename );

	// Converts an array of floats to an array of halfs (aka unsigned shorts)
	unsigned short *ConvertToHalves( float *inData, int floatCount );

	// How many points are in the cloud?
	int m_pointCount;

	// What is the cloud/model center?
	vec3 m_cloudCenter;

	// Information about the data at each point in our point cloud
	bool m_hasVertices, m_hasNormals, m_hasTexCoords, m_hasMatlID, m_hasObjectID;

	// Our buffer(s) storing the vertex array for this geometry
	IGLUVertexArray::Ptr m_vertArr;

	// Did the user resize and/or recenter the data before storing it to a point cloud?
	bool m_resize, m_center;

	// Is the data stored in halfs?
	bool m_useHalf;

	// What is the stride of the vertex array?
	GLuint m_vertStride, m_vertOff, m_normOff, m_texOff, m_matlIdOff, m_objectIdOff;

};


// End namespace iglu
}



#endif

