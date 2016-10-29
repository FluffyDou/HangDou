#version 450
//#extension GL_EXT_geometry_shader4 : enable
//#extension GL_NV_geometry_shader_passthrough : enable

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

//uniform mat4 proj;
//uniform mat4 view;
//uniform mat4 model;
//uniform vec3 eye_position;
//uniform int  radiusScale;


in vertexData
{
    //vec3  fragTexcoord;
    //vec4  esFragNormal; 
    //vec4  lsFragNormal;
    //vec4  esFragPos;
    //vec4  wsFragPos;
    //vec4  lsFragPos;
    //float fragMatlID;
    vec2 vertFragPos;
		
} vert[];


// parameters for ray intersection test
//out	vec3  fragTexcoord;
//out	vec4  esFragNormal; 
//out	vec4  lsFragNormal;
//out	vec4  lsFragFlatNormal;
//out	vec4  esFragPos;
//out	vec4  wsFragPos;

out	vec2 vertFragPos1;
out	vec2 vertFragPos2;
out	vec2 vertFragPos3;
// Redeclare gl_PerVertex to pass through "gl_Position".
//layout(passthrough) in gl_PerVertex
//{
//  vec4 gl_Position;
//};

// Declare "Inputs" with "passthrough" to automatically copy members.
//layout(passthrough) in Inputs 
//{
//    vec3 vertFragPos;
//};


void main()
{
    vertFragPos1 = vert[0].vertFragPos;
    vertFragPos2 = vert[1].vertFragPos;
    vertFragPos3 = vert[2].vertFragPos;
    
    gl_Position  = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    gl_Position = gl_in[2].gl_Position;
    EmitVertex();   

    EndPrimitive();
}
