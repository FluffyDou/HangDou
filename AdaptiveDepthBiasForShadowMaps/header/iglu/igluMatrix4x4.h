/******************************************************************/
/* igluMatrix4x4.h                                                */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic matrix class that implments most of   */
/*     operations you need to perform on matrices, and interacts  */
/*     correctly with point and vector classes.                   */
/*                                                                */
/* Chris Wyman (02/23/2007)                                       */
/******************************************************************/

#ifndef IGLU_MATRIX4X4_H
#define IGLU_MATRIX4X4_H   1

#include <cmath>
#include <cstring>
#include <cstdio>

namespace iglu {

class vec4;
class vec3;

class IGLUMatrix4x4 {
    float mat[16]; 
	friend class vec4;
	friend class vec3;
public:
	// IGLUMatrix4x4 constructors
	inline IGLUMatrix4x4(float data[16]);         
	inline IGLUMatrix4x4(const IGLUMatrix4x4& v);       
	inline IGLUMatrix4x4(); 

	// Assignment constructors
   	inline IGLUMatrix4x4& operator=(const IGLUMatrix4x4& v);

	// Accessor methods to access the matrix data
	inline float& operator[](const int index)              { return mat[index]; }
	inline float& operator()(const int col, const int row) { return mat[col*4+row]; }

	// Bad accessor method.  This breaks encapsulation, but is particularly useful in
	//    the context of OpenGL when passing matrices to the hardware.  Using this in
	//    other contexts is likely to 
	inline float *GetDataPtr()                             { return mat; }
	inline const float *GetConstDataPtr() const            { return mat; }
	inline void Print() const;
	// Mathematical/assignment operators
	IGLUMatrix4x4& operator*=(const IGLUMatrix4x4& v);
	IGLUMatrix4x4& operator*=(const float s);
	IGLUMatrix4x4& operator+=(const IGLUMatrix4x4& v);
	IGLUMatrix4x4& operator-=(const IGLUMatrix4x4& v);

	// Mathematical matrix operators
	IGLUMatrix4x4 operator()(const IGLUMatrix4x4& v) const;
	IGLUMatrix4x4 operator*(const IGLUMatrix4x4& v)  const;
	IGLUMatrix4x4 operator*(const float s)       const;
	IGLUMatrix4x4 operator+(const IGLUMatrix4x4& v)  const;
	IGLUMatrix4x4 operator-(const IGLUMatrix4x4& v)  const;

	// Friend mathematical matrix operators
	friend IGLUMatrix4x4 operator*(const float s, const IGLUMatrix4x4& v);

	// Mathematic point/vector operators
	vec4 operator()(const vec4& v) const;
	vec4 operator*(const vec4& v)  const;

	// Other interesting functions
	IGLUMatrix4x4 Invert( void ) const;
	IGLUMatrix4x4 Transpose( void ) const;
	float         Determinant( void ) const; // Not yet implemented!

	// Static constants
	static float identityArray[16];
	static float zeroArray[16];

	// Static methods.  Useful for defining commonly used matrices
	static IGLUMatrix4x4 Zero     ( void )          { return IGLUMatrix4x4( zeroArray );     }
	static IGLUMatrix4x4 Identity ( void )          { return IGLUMatrix4x4( identityArray ); }

	// Static methods for defining common transformations.  None of these are optimized for speed!
	static IGLUMatrix4x4 Translate( float x, float y, float z );			 
	static IGLUMatrix4x4 Translate( const vec3& direction );		
	static IGLUMatrix4x4 Scale    ( float x, float y, float z );			 
	static IGLUMatrix4x4 Scale    ( const vec3& scale );
	static IGLUMatrix4x4 Scale    ( float constantScale );					 
	static IGLUMatrix4x4 Rotate   ( float angle, const vec3& axis );		 
	static IGLUMatrix4x4 Perspective  ( float fovy, float aspect, float zNear, float zFar );   
	static IGLUMatrix4x4 Ortho        ( float left, float right, float bottom, float top, float near_=-1, float far_=1 );
	static IGLUMatrix4x4 LookAt       ( const vec3 &eye, const vec3 &at, const vec3 &up );   

	// A scale-and-bias matrix for shadow mapping (e.g., scale 0.5, translate (x,y,z) 0.5, and add a user-spec z-bias.
	static IGLUMatrix4x4 ScaleAndBias( float zBias );

	// A matrix to access the shadow map created from a light with a projection matrix 'lightProj' and light look-at
	//    matrix 'lightView' when seen from an eye position with look at matrix 'eyeView'.
	static IGLUMatrix4x4 ShadowMapAccessMatrix( const IGLUMatrix4x4 &lightProj, const IGLUMatrix4x4 &lightView, 
		                                        const IGLUMatrix4x4 &eyeView, float zBias = -0.005 );

	// Static methods for defining matrices given from outer products (of real-valued vectors)
	static IGLUMatrix4x4 OuterProduct ( const vec3& u );					// Returns u * u^T
	static IGLUMatrix4x4 OuterProduct ( const vec4& u );					// Returns u * u^T
	static IGLUMatrix4x4 OuterProduct ( const vec4& u, const vec4& v );		// Returns u * v^T

	// Check if the matrix is the identity matrix (FIX THIS <-- it's poorly written and not robust!)
	bool IsIdentity( void )  { return (mat[0]==1&&mat[4]==0&&mat[8]==0&&mat[12]==0&&
					   mat[1]==0&&mat[5]==1&&mat[9]==0&&mat[13]==0&&
		                           mat[2]==0&&mat[6]==0&&mat[10]==1&&mat[14]==0&&
					   mat[3]==0&&mat[7]==0&&mat[11]==0&&mat[15]==1); }

	// A pointer to a IGLUMatrix4x4 could have type IGLUMatrix4x4::Ptr
    typedef IGLUMatrix4x4 *Ptr;
};

inline void IGLUMatrix4x4::Print( void ) const
{
	printf("Matrix: %f %f %f %f\n",  mat[0],  mat[1],  mat[2], mat[3]  );
	printf("        %f %f %f %f\n",  mat[4],  mat[5],  mat[6], mat[7]  );
	printf("        %f %f %f %f\n",  mat[8],  mat[9], mat[10], mat[11] );
	printf("        %f %f %f %f\n", mat[12], mat[13], mat[14], mat[15] );
}

inline IGLUMatrix4x4::IGLUMatrix4x4(float data[16]) 
{
	memcpy( mat, data, 16*sizeof(float) );
}

inline IGLUMatrix4x4::IGLUMatrix4x4(const IGLUMatrix4x4& v)      
{
	memcpy( mat, v.mat, 16*sizeof(float) );
}

inline IGLUMatrix4x4::IGLUMatrix4x4()             
{
	memcpy( mat, IGLUMatrix4x4::identityArray, 16*sizeof(float) );
}

inline IGLUMatrix4x4& IGLUMatrix4x4::operator=(const IGLUMatrix4x4& v) 
{
	memcpy( mat, v.mat, 16*sizeof(float) );
    	return *this;
}

// Done namespace iglu
}

#endif

