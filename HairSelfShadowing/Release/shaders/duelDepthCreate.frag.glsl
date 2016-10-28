/*The frag shader for creating the duel depth map*/

#version 420

layout(location = 0) out vec4 result;

void main( void ) {

    // we set the blend func for color channel to be min and alpha channel to be max in c code
	result = vec4(0, 0, gl_FragCoord.z, gl_FragCoord.z);
  
}