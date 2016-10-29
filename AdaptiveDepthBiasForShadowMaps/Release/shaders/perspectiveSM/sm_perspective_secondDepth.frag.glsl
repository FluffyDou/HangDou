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

#define  EPSILON 0.0001

in vec4 smCoord;

out vec4 result;

uniform sampler2D firstDepth; // First nearest depth texture

void main( void ) 
{                 
    vec4 smTexCoord = smCoord/smCoord.w;                        // [-1, 1]
    smTexCoord  = 0.5 * smTexCoord + vec4(0.5, 0.5, 0.5, 0.0);  // [ 0, 1]

    // First nearest depth value
    float fDepth = texture( firstDepth, smTexCoord.xy ).x;
    
    if( gl_FragCoord.z <= fDepth)
        discard;
    
    result = vec4(0);
}