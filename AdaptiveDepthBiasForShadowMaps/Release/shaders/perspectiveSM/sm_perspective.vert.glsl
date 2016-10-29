#version 420 

layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)

uniform mat4 model;	 // Transforms: model space -> world space
uniform mat4 view;   // Transforms: world space -> eye space
uniform mat4 proj;   // Transforms: eye space   -> clip space

out vec4 esFragNormal;
//out vec4 lsFragNormal;
//out vec4  esFragPos;
//out float fragMatlID;


void main( void )
{  
	// Store the eye-space position of this vertex, send to frag shader
	//esFragPos = view * model * vec4( vertex, 1.0f );
	
	// Transform vertex normal to eye-space, send to frag shader
	esFragNormal = inverse( transpose( view * model ) ) * vec4( normal, 0.0 );	

    // Light space fragment normal --- the view is light view
    //lsFragNormal = inverse( transpose( view * model ) ) * vec4( normal, 0.0 );	    
    
    // Transform vertex into clip space
	gl_Position = proj * ( view * model * vec4( vertex, 1.0f ) );	
	
	// Pass down the material ID
	//fragMatlID   = matlID;
}