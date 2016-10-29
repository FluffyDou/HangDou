#version 420 

// Our shadow map texture.  Layer 0 = positive z-hemisphere, layer 1 = negative
uniform sampler2DArray firstDepth; 

in vec3 smTexCoord; // [-1,1]

out vec4 outColor;

void main( void )
{
    vec3 parabTexCoord = vec3( smTexCoord.z!=0 ? vec3(0.5*smTexCoord.xy+vec2(0.5), 0) : vec3(0.5*smTexCoord.xy+vec2(0.5), 1));
    // First nearest depth value
    float fDepth = texture( firstDepth, parabTexCoord ).x;
    
    if( gl_FragCoord.z <= fDepth )
        discard;
    
    outColor = vec4(0);
}

