#version 450
//#extension GL_NV_conservative_raster : enable
layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
//layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)
//layout(location = IGLU_TEXCOORD) in vec2 texcoord; // Sends texture as albedo

layout(binding = 0, offset = 0) uniform atomic_uint ac;

//layout(binding = 0, r32i) coherent restrict volatile uniform iimage1D imageCounter;

uniform mat4 model;	     // Transforms: model space -> world space
uniform mat4 view;       // Transforms: world space -> eye space
uniform mat4 proj;       // Transforms: eye space   -> clip space

uniform int vertexID;


//uniform mat4 lightView;  // For shadow map looking up

// The buffer containing our material information
//uniform samplerBuffer  matlInfoTex;

//out vec3  fragTexcoord;
//out float fragMatlID;

//out float w1;
//out float w2;
//out float w3;

//out vec2 w;

//out vec2 vertFragPos1;
//out vec2 vertFragPos2;
//out vec2 vertFragPos3;

out vec4 tmpColor;
out float dummy;
//flat          out  vec4 wsPoint0;
//noperspective out  vec4 wsPoint1;

void main( void )
{  
    //memoryBarrier();
    //ivec4 vID = imageLoad(imageCounter, 0);
    //imageStore(imageCounter, 0, ivec4(vID.x+1, 0, 0, 0));
    
    uint counter = atomicCounterIncrement(ac);
    //counter %= 3;
    //w = vec2(gl_VertexID%3 == 0 ? 1.0 : 0.0, gl_VertexID%3 == 1 ? 1.0 : 0.0);

    //w1 = gl_VertexID%3 == 0 ? 1.0 : 0.0;
    //w2 = gl_VertexID%3 == 1 ? 1.0 : 0.0;
    //w3 = gl_VertexID%3 == 2 ? 1.0 : 0.0;  
    dummy = float(gl_VertexID);
    
    tmpColor = counter == vertexID ? vec4(1,0,0,1) : vec4(1);
    
    // Transform vertex into clip space
	gl_Position = proj * view * model * vec4( vertex, 1.0f );
	//vec2 fragPos = gl_Position.xy/gl_Position.w;
    //vertFragPos1 = gl_VertexID%3 == 0 ? fragPos : vec2(0);
    //vertFragPos2 = gl_VertexID%3 == 1 ? fragPos : vec2(0);
    //vertFragPos3 = gl_VertexID%3 == 2 ? fragPos : vec2(0);      
    
	// Pass down the material ID
	//fragMatlID   = matlID;
    
    // look up the material and output the diffuse color
	//int matlIdx  = 4*int(matlID) + 1;
    // Get texture index, we are not worried that it interpolate
    //fragTexcoord =  vec3(texcoord, texelFetch(matlInfoTex, matlIdx + 2).r);     
}