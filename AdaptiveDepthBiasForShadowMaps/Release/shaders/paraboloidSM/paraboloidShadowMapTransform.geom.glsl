// Creates triangles in a dual paraboloid shadow map.
//
//    Warning!  Possible error: does not check to make sure geometry generated
//              in on paraboloid will not be drawn into the other half of the window
//              accidentally.  This should only be a problem with large triangles
//              which will have vertices on dramatically different sides of the sphere.

#version 420 
#extension GL_EXT_geometry_shader4 : enable

// Tells the OpenGL/IGLU program to enable the depth test when using this shader
#pragma IGLU_ENABLE  DEPTH_TEST

// Information about this geometry shader.  It's called only once on each
//   input triangle, and creates a set of triangle_strip (really, triangles) 
//   with up to 6 vertices output for each input triangle.
layout (triangles, invocations = 1)       in;
layout (triangle_strip, max_vertices = 6) out;

// Data passed down from the fragment shader
in vec4 vertTexCoord[];

void main( void )
{
	// Check:  Does triange have vertices in positive and/or negative paraboloids?
	bool pos = (vertTexCoord[0].x < 1.0f) || (vertTexCoord[1].x < 1.0f) || (vertTexCoord[2].x < 1.0f);
	bool neg = (vertTexCoord[0].y < 1.0f) || (vertTexCoord[1].y < 1.0f) || (vertTexCoord[2].y < 1.0f);
	
	// If we have vertices in positive paraboloid
	if ( pos )
	{ 
		gl_Layer = 0;
		gl_Position = vec4( gl_PositionIn[0].xy, vertTexCoord[0].z, 1.0 ); 
		EmitVertex();
		gl_Position = vec4( gl_PositionIn[2].xy, vertTexCoord[2].z, 1.0 ); 
		EmitVertex();
		gl_Position = vec4( gl_PositionIn[1].xy, vertTexCoord[1].z, 1.0 ); 
		EmitVertex();
		EndPrimitive();	
	}
	
	// If we have vertices in negative paraboloid
	if ( neg )
	{
		gl_Layer = 1;
		gl_Position = vec4( gl_PositionIn[0].zw, vertTexCoord[0].z, 1.0 ); 
		EmitVertex();
		gl_Position = vec4( gl_PositionIn[1].zw, vertTexCoord[1].z, 1.0 ); 
		EmitVertex();
		gl_Position = vec4( gl_PositionIn[2].zw, vertTexCoord[2].z, 1.0 ); 
		EmitVertex();
		EndPrimitive();
	}
	
}