/******************************************************************/
/* igluPinholeCamera.h                                            */
/* -----------------------                                        */
/*                                                                */
/* This class is a fairly straightforward pinhole camera class to */
/*    keep camera data and interactions encapsulated.             */
/*                                                                */
/* Chris Wyman (1/14/2013)                                        */
/******************************************************************/

#ifndef IGLU_PINHOLE_CAMERA_H
#define IGLU_PINHOLE_CAMERA_H

#pragma warning( disable: 4996 )

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include "iglu/igluVector.h"
#include "iglu/igluMatrix4x4.h"
#include "iglu/igluOrthoNormalBasis.h"


namespace iglu {


class IGLUPinholeCamera
{
public:
	// Default constructor (think default GL constructor)
	IGLUPinholeCamera();

	// More general constructor (think gluLookAt())
	IGLUPinholeCamera( vec3 eyePos, vec3 viewVec, vec3 upVec, float fovy, float aspectRatio=(4.0f/3.0f), float zNear=0.1, float zFar=20.0 );

	// Destructor
	~IGLUPinholeCamera() {}

	// Get the current camera projection matrix
	const IGLUMatrix4x4 &GetProjMatrix() const { return m_projMatrix; }

	// Get the current camera view matrix
	const IGLUMatrix4x4 &GetViewMatrix() const { return m_viewMatrix; }

	// Get scalar camera parameters
	float                GetAspect()     const { return m_aspectRatio; }
	float                GetNear()       const { return m_zNear; }
	float                GetFar()        const { return m_zFar; }
	float                GetFovy()       const { return m_fovy; }
	float                GetFovx()       const { return m_fovy*m_aspectRatio; }

	// Get vector camera parameters
	const vec3 &         GetEyePos()     const { return m_eyePos; }
	const vec3 &         GetViewVec()    const { return m_viewVec; }
	const vec3 &         GetUpVec()      const { return m_upVec; }

	// Set scalar camera parameters
	void SetAspect( float aspectRatio );
	void SetNear  ( float zNear );
	void SetFar   ( float zFar );
	void SetFovy  ( float fovy );

	// EXPERIMENTAL:  Movement controls.  
	//   -> Speed is automatically set as a function of distance to far plane.
	//   -> The 'scale' parameters scale the automatically computed speed.
	void MoveForward( float scale=1.0f );
	void MoveBack   ( float scale=1.0f );
	void MoveLeft   ( float scale=1.0f );
	void MoveRight  ( float scale=1.0f );
	void MoveUp     ( float scale=1.0f );
	void MoveDown   ( float scale=1.0f );

	// EXPERIMENTAL:  Given a vector (x,y,z), set the new view vector as
	//     x * rightVec + y * upVec - z * viewViec
	void RotateView ( vec3 relativeVec );

	// A pointer to a IGLUPinholeCamera could have type IGLUPinholeCamera::Ptr
	typedef IGLUPinholeCamera *Ptr;

protected:
	// The camera's current view and projection matrices
	IGLUMatrix4x4 m_projMatrix, m_viewMatrix;

	// Camera parameters
	float m_zNear, m_zFar;
	float m_fovy;
	float m_aspectRatio;

	// Camera vectors
	vec3 m_eyePos;
	vec3 m_upVec, m_viewVec;
	vec3 m_rightVec;

	// Automatically computed translation speed.
	float m_moveSpeed;

	// A common method used to update the matrices;
	void UpdateMatrices( void );
};



// End namespace iglu
}



#endif