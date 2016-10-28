#version 420 

layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
//layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)

uniform mat4 model;	 // Transforms: model space -> world space
uniform mat4 view;   // Transforms: world space -> eye space
uniform mat4 proj;   // Transforms: eye space   -> clip space


void main( void )
{  
    // Transform vertex into clip space
	gl_Position = proj * view * model * vec4( vertex, 1.0f );	
}