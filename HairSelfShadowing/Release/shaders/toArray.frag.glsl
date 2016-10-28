#version 400

uniform sampler2D inputTex;
in vec2 fragTexCoord;
out vec4 result;

// This is the actual code that is executed.
void main( void )
{
	// Output the result
	result = texture( inputTex, vec2( fragTexCoord.x,
	                                  1.0f-fragTexCoord.y ) );
}