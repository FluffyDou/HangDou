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
    vec3 vertFragPos;
		
} vert[];


// parameters for ray intersection test
//out	vec3  fragTexcoord;
//out	vec4  esFragNormal; 
//out	vec4  lsFragNormal;
//out	vec4  lsFragFlatNormal;
//out	vec4  esFragPos;
//out	vec4  wsFragPos;
out vec3 p0;
out	vec3 p1;
out	vec3 p2;
out	vec3 p3;
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
    p1 = vert[0].vertFragPos;
    p2 = vert[1].vertFragPos;
    p3 = vert[2].vertFragPos;
    
    p0 = vert[0].vertFragPos;
    gl_Position  = gl_in[0].gl_Position;
    EmitVertex();

    p0 = vert[1].vertFragPos;
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    p0 = vert[2].vertFragPos;
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();   

    EndPrimitive();
}