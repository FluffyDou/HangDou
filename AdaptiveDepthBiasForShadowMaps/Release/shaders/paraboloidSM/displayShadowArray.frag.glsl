#version 420
uniform sampler2DArray inputTex;
uniform float layer;

in vec2 fragTexCoord;
out vec4 result;
void main( void ) 
{
    vec4 texColor = texture( inputTex, vec3( fragTexCoord.xy, layer) );
    result = texColor;
}