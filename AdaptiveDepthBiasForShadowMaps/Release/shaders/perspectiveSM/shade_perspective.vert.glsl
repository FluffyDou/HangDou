#version 420 

layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)
layout(location = IGLU_TEXCOORD) in vec2 texcoord; // Sends texture as albedo

uniform mat4 model;	     // Transforms: model space -> world space
uniform mat4 view;       // Transforms: world space -> eye space
uniform mat4 proj;       // Transforms: eye space   -> clip space

uniform mat4 lightView;  // For shadow map looking up

// The buffer containing our material information
uniform samplerBuffer  matlInfoTex;

out vec3  fragTexcoord;
out vec4  esFragNormal; 
out vec4  lsFragNormal; 
out vec4  esFragPos;
out vec4  wsFragPos;
out float fragMatlID;


void main( void )
{  
    wsFragPos = model * vec4( vertex, 1.0f );

	// Store the eye-space position of this vertex, send to frag shader
	esFragPos = view * wsFragPos;
	
	// Transform vertex normal to eye-space, send to frag shader
	esFragNormal = inverse( transpose( view * model ) ) * vec4( normal, 0.0 );	
    // Light space fragment normal
    //lsFragNormal = inverse( ( lightView * model ) ) * vec4( normal, 0.0 );	
    lsFragNormal = ( lightView * model ) * vec4( normal, 0.0 );	
    
    // Transform vertex into clip space
	gl_Position = proj * esFragPos;	
	
	// Pass down the material ID
	fragMatlID   = matlID;
    
    // look up the material and output the diffuse color
	int matlIdx  = 4*int(matlID) + 1;
    // Get texture index, we are not worried that it interpolate
    fragTexcoord =  vec3(texcoord, texelFetch(matlInfoTex, matlIdx + 2).r);     
}