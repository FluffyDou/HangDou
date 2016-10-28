#version 420

#define BITLENGTH 95.0

uniform vec2  screenRes; // texture length
uniform float znear;     // z_near
uniform float zfar;      // z_far

// texture storing the depth_near & depth_far in [0,1]
layout(binding = 0) uniform sampler2D depthRange;

layout(location = 0) out uvec4 result;

void main( void ) {

	// get the texture coordinate of the current fragment location
    vec2 texCoord      = vec2(gl_FragCoord.x-0.5, gl_FragCoord.y-0.5) / screenRes;
	// get the duel depth for the current fragment location
    vec4 duelDepth     = texture( depthRange, texCoord );
	// pack the depth_near and depth_far values into a float number
    uint packedDepthUI = packUnorm2x16( vec2(duelDepth.z, duelDepth.w) );           
    
    // get the real depth of the fragment and light near and light far
	float realNear = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*duelDepth.z    - 1.0) );
	float realFar  = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*duelDepth.w    - 1.0) );
	float realZi   = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*gl_FragCoord.z - 1.0) );	
	
    int bitFlag = int( floor( BITLENGTH * ((realZi-realNear)/(realFar-realNear)) ) ); // to tell which bit in [0,95]
    int bitLocation = bitFlag / 32; // to tell in which channel to set the bit
    
    // set the depth bit according to the fragment's depth
    uvec3 bitBar = uvec3(0u, 0u, 0u);
    bitBar.x = bitLocation == 0 ? 1u << uint(31 - bitFlag) : 0u;
    bitBar.y = bitLocation == 1 ? 1u << uint(63 - bitFlag) : 0u;
    bitBar.z = bitLocation == 2 ? 1u << uint(95 - bitFlag) : 0u;   
    
    // the alpha channel gets the packed duel depth and the color channel stores frag occupancy
	result = uvec4(bitBar, packedDepthUI);
    
	//result = gl_FragCoord.x < 0.5*screen.x ? vec4(1.0, 0.0, 0.0, 1.0) : vec4(0.0, 1.0, 0.0, 1.0)
	//result = vec4(duelDepth.z, duelDepth.w, 0.0, packedDepth2x16F);
	//result = vec4(0.5, 0.5, 0.0, 1.0);
}