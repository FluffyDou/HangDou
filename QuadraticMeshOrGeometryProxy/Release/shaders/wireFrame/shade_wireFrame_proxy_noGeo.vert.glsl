#version 450 
//#extension GL_NV_conservative_raster : enable
layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)
//layout(location = IGLU_TEXCOORD) in vec2 texcoord; // Sends texture as albedo

uniform mat4 model;	     // Transforms: model space -> world space
uniform mat4 view;       // Transforms: world space -> eye space
uniform mat4 proj;       // Transforms: eye space   -> clip space

//uniform mat4 lightView;  // For shadow map looking up

// The buffer containing our material information
//uniform samplerBuffer  matlInfoTex;

//out vec3  fragTexcoord;
//out float fragMatlID;

//out float w1;
//out float w2;
//out float w3;

out vec2 w;

// This is w1*x1
out vec3 vertFragPos1;
// This is w2*x2
out vec3 vertFragPos2;
// This is actually "w1*x1+w2*x2+w3*x3"
out vec3 vertFragPos3;

//flat out int testID;

//flat          out  vec4 wsPoint0;
//noperspective out  vec4 wsPoint1;

void main( void )
{  
    w = vec2(gl_VertexID%3 == 0 ? 1.0 : 0.0, gl_VertexID%3 == 1 ? 1.0 : 0.0);
    //testID = (gl_VertexID);
    //w1 = gl_VertexID%3 == 0 ? 1.0 : 0.0;
    //w2 = gl_VertexID%3 == 1 ? 1.0 : 0.0;
    //w3 = gl_VertexID%3 == 2 ? 1.0 : 0.0;  
    vec4 vertFragPos = view * model * vec4( vertex, 1.0f );

    vertFragPos1 = gl_VertexID%3 == 0 ? vertFragPos.xyz : vec3(0);
    vertFragPos2 = gl_VertexID%3 == 1 ? vertFragPos.xyz : vec3(0);
    //vertFragPos3 = gl_VertexID%3 == 2 ? vertFragPos.xyz : vec3(0);      
    
    vertFragPos3 = vertFragPos.xyz;
    
    // Transform vertex into clip space
	gl_Position = proj * vertFragPos;    
    
	// Pass down the material ID
	//fragMatlID   = matlID;
    
    // look up the material and output the diffuse color
	//int matlIdx  = 4*int(matlID) + 1;
    // Get texture index, we are not worried that it interpolate
    //fragTexcoord =  vec3(texcoord, texelFetch(matlInfoTex, matlIdx + 2).r);     
}