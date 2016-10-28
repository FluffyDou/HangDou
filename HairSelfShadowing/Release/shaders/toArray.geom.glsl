#version 400 
#extension GL_EXT_geometry_shader4 : enable

layout (triangles, invocations = 1)       in;
layout (triangle_strip, max_vertices = 3) out;

// Our one input...  What layer to store this into
uniform int outputLayer;

in vec2 vertTexCoord[];
out vec2 fragTexCoord;

void main( void )
{
	// Save the layer we want to store this triangle into
	gl_Layer = outputLayer;
	
	// Pass the triangle down to be rasterized...
	fragTexCoord = vertTexCoord[0];
	gl_Position  = gl_PositionIn[0];
	EmitVertex();
	fragTexCoord = vertTexCoord[1];
	gl_Position  = gl_PositionIn[1];
	EmitVertex();
	fragTexCoord = vertTexCoord[2];
	gl_Position  = gl_PositionIn[2];
	EmitVertex();
	EndPrimitive();
}