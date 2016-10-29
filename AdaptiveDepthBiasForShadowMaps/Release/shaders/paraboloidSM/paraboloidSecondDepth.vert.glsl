#version 420 

// Tells the OpenGL/IGLU program to enable the depth test when using this shader
#pragma IGLU_ENABLE  DEPTH_TEST

// Grab the vertex position
layout(location = IGLU_VERTEX)   in vec3 vertex;   
layout(location = IGLU_NORMAL)   in vec3 normal;   
layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)

// The model matrix
uniform mat4x4 model;
//uniform mat4x4 eyeView;

// The location of the light, in world-space
uniform vec4 wsLightPos;

// Information about the near & far planes for projecting the vertex's z-distance
uniform float nearDist, farDist;

// Data passed to the geometry shader for positioning the geometry in the shadow map
out vec4 vertTexCoord;

void main( void )
{	
	// Compute the eye-space position of the vertex
	vec4 wsVert  = model * vec4( vertex, 1.0f );
	
	// Compute the location of the vertex relative to the light
	vec3 toVert = wsVert.xyz - wsLightPos.xyz;
	float dist  = length( toVert );
	toVert = normalize( toVert );
	
	// Compute a perspective NDC z-value (in the range [-1..1]) using the z-row from
	//     standard GL perspective matrices (e.g., gluPerspective).  Also do the perspective divide.
	float z_persp = (-(nearDist+farDist)*dist + (2.0*nearDist*farDist))/(dist*(nearDist-farDist));
	
	// Now compute the paraboloid coordinates (x,y) in [-1,1], at least for valid
	//     points on this hemi-paraboloid map.
	vec2 posMap = vec2(toVert.xy)/(toVert.z+1.0); // + d0 (0,0, 1)
	vec2 negMap = vec2(toVert.xy)/(toVert.z-1.0); // + d1 (0,0,-1)
	
	// Send down the reparameterized vertex to the geometry shader to position it in the image
	vertTexCoord = vec4( dot(posMap,posMap), dot(negMap,negMap), z_persp, 1.0 );
	gl_Position  = vec4( posMap, negMap ); 
}
