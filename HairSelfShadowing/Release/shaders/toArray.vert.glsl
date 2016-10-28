#version 400 

layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
layout(location = IGLU_TEXCOORD) in vec2 texCoord; // Sends texture coordinate data (if any) from models here

out vec2 vertTexCoord;

// This is the actual code that is executed.
void main( void )
{  
	vertTexCoord = texCoord;
	gl_Position = vec4( vertex.xyz, 1.0 ); 
}