// This is a shader to "create" a shadow map.  It's actually
//    really not important, because we don't care about the
//    color parts of the shadow mask (only the depth).  The
//    depth value is automatically set to the z-buffer when
//    it exists (and z-testing is enabled).

// Of course for more complex techniques, one might use the
//    color buffers that correspond to the z-buffer to store
//    useful stuff (e.g., during g-buffering), but that's for a
//    different example.

#version 420
in vec4 esFragNormal;

out vec4 result;

void main( void ) 
{
    // Store the light space frag normal
	//result = vec4(normalize(esFragNormal.xyz), 1.0);
    result = vec4(0);
}