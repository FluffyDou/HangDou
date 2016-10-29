#version 450 
//#extension GL_NV_conservative_raster : enable
layout(location = IGLU_VERTEX)   in vec3 vertex;   // Sends vertex data from models here
layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)
//layout(location = IGLU_TEXCOORD) in vec2 texcoord; // Sends texture as albedo

layout(binding = 0, offset = 0) uniform atomic_uint ac;

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

//uniform int vertexID;
//out float[100] st;
//out vec4 tmpColor;
out vec2 w;

out vec2 vertFragPos1;
out vec2 vertFragPos2;
out vec2 vertFragPos3;

//flat          out vec3 p1;
//noperspective out vec3 p2;
//smooth        out vec3 p3;
//smooth        out vec3 interpolatedNormal;

//flat          out  vec4 wsPoint0;
//noperspective out  vec4 wsPoint1;

void main( void )
{  

    uint counter = atomicCounterIncrement(ac);
    counter %= 3;
    
    //st[counter] = 1.0;
    
    //uint counter = gl_VertexID;
    //tmpColor = gl_VertexID == vertexID ? vec4(1,0,0,1) : vec4(1);

    w = vec2(counter == 0 ? 1.0 : 0.0, counter == 1 ? 1.0 : 0.0);

    //w1 = gl_VertexID%3 == 0 ? 1.0 : 0.0;
    //w2 = gl_VertexID%3 == 1 ? 1.0 : 0.0;
    //w3 = gl_VertexID%3 == 2 ? 1.0 : 0.0;  
    
    vec4 viewPos = view * model * vec4( vertex, 1.0f );
    
    //p1 = vertex.xyz;
    //p2 = vertex.xyz;
    //p3 = vertex.xyz;
    
    //interpolatedNormal = normal;
    
    // Transform vertex into clip space
	gl_Position = proj * viewPos;
	vec2 fragPos = gl_Position.xy/gl_Position.w;
    vertFragPos1 = counter == 0 ? fragPos : vec2(0);
    vertFragPos2 = counter == 1 ? fragPos : vec2(0);
    vertFragPos3 = counter == 2 ? fragPos : vec2(0);      
    
	// Pass down the material ID
	//fragMatlID   = matlID;
    
    // look up the material and output the diffuse color
	//int matlIdx  = 4*int(matlID) + 1;
    // Get texture index, we are not worried that it interpolate
    //fragTexcoord =  vec3(texcoord, texelFetch(matlInfoTex, matlIdx + 2).r);     
}