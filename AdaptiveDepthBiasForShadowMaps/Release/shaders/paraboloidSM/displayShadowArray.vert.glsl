#version 420

#pragma IGLU_ENABLE  DEPTH_TEST

layout(location = 0) in vec3 vertex;
layout(location = 2) in vec2 texcoord;
out vec2 fragTexCoord;

void main( void ) 
{
    gl_Position = vec4( vertex.xyz, 1.0f );
    fragTexCoord = texcoord.xy;
}