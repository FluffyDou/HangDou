#version 420

#define BITLENGTH 95.0

// Information from the vertex shader about the eye-space location
//    of this fragment we're working on
in vec4  esFragPos;
in vec4  lsFragPos;
in vec4  esFragNormal;
in float fragMatlID;

// Where to store the color we compute
out vec4 result;

uniform float useLambertian;  // Asked to use a Lambertian shading scheme?
uniform float zBias;          // avoiding self-shadowing
uniform vec4 esLightPos;      // Location of the light, relative to the eye

uniform float znear; // light z_near
uniform float zfar;  // light z_far

uniform float intenseBias; // brightness bias when transparent light depends directly on the number of frag before
uniform float fragAlpha;   // alpha value for each fragment(alpha == 0 means totally opaque fragment)
uniform int   shadowFlag;  // use normal shadow map or occupancy shadow map
uniform int   layer;       // choose between frag number trans or alpha trans

// the duel depth map texture.    
uniform usampler2D shadowTex;

// The buffer containing our material information
uniform samplerBuffer  matlInfoTex;


// This is the actual code that is executed.
void main( void )
{	
	// Look up data about this particular fragment's material
	int  matlIdx     = 4 * int(fragMatlID+0.5);
	vec4 matlInfo   = texelFetch( matlInfoTex, matlIdx+1 );  // Get diffuse color

	// Transform from eye-space to shadow map texture coordinates
	vec4 smTexCoord = lsFragPos/lsFragPos.w;
    // map x y from [-1,1] to [0,1]
    smTexCoord    = (smTexCoord + vec4(1.0))/2.0;
	smTexCoord.z += zBias;
	
    // sample in the occupancey map
	uvec4 alphaBar = texture(shadowTex, smTexCoord.xy);
    
    // unpack the duel depth
    uint packedDepth = alphaBar.w;        
    vec2 duelDepth   = unpackUnorm2x16(packedDepth);
    
	// get the real depth of the fragment and light near and light far
	float realNear = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*duelDepth.x  - 1.0) );
	float realFar  = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*duelDepth.y  - 1.0) );
	float realZi   = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*smTexCoord.z - 1.0) );	
	
    int bitFlag = int( floor( BITLENGTH * ((realZi-realNear)/(realFar-realNear)) ) ); // to tell which bit in [0,95]
    int bitLocation = bitFlag / 32; // to tell in which channel to set the bit
	
    // count the number of set bits before the current bit ---- how many frags lie before the current frag
	int fragNum = 0;
    fragNum = bitLocation == 0 ?  bitCount(alphaBar.x >> uint(32 - bitFlag)) : fragNum;
    fragNum = bitLocation == 1 ?  bitCount(alphaBar.x) + bitCount(alphaBar.y >> uint(64 - bitFlag)) : fragNum;
    fragNum = bitLocation == 2 ?  bitCount(alphaBar.x) + bitCount(alphaBar.y) + bitCount(alphaBar.z >> uint(96 - bitFlag)) : fragNum;	
    fragNum = bitLocation  > 2 ?  bitCount(alphaBar.x) + bitCount(alphaBar.y) + bitCount(alphaBar.z) : fragNum;
	
    // two types of energy attenuation functions
	float alphaTrans = fragAlpha < 1 ? pow(1.0-fragAlpha, fragNum) : 0.0;	
    float fragTrans = 1.0 - float(fragNum)/intenseBias;	
	
    // choose a type of energy attenuation function
    float accumTrans = layer == 0 ? fragTrans : alphaTrans;
    
	// isLit based on an occupancey map
    // duelDepth == 0 means that no occluder drawn in the light view there
    // transpaency is 0 means totally opaque
	float isLitTrans = length(alphaBar) == 0 ? 1.0 : accumTrans;	

	// isLit based on a normal shadow map
	float isLitNorm = duelDepth.x > smTexCoord.z || length(alphaBar) == 0 ? 1.0 : 0.0;
	
    // choose to use normal depth compare or energy attenuation
	float isLit = shadowFlag == 0 ? isLitTrans : isLitNorm;	
    
	// if the frag coord in light space is out of light near and light far range, set the isLit to be 0
	isLit = smTexCoord.z > 1 || smTexCoord.z < 0 ? 0.0 : isLit;	
	
	// Find simple vectors for basic shading
	vec3 toLight = normalize( esLightPos.xyz - esFragPos.xyz);
	vec3 toEye   = normalize( -esFragPos.xyz );
	
	// Figure out the Lambertian color (if shadowed, set to 0)
	vec3 lambert = isLit * useLambertian * max( 0.0f, dot( normalize(esFragNormal.xyz), toLight ) ) * matlInfo.xyz;	
	// Figure out the ambient color (if shadowed, set to 0)
	vec3 ambient = isLit * (1.0-useLambertian) * matlInfo.xyz;
	
	// The ambient & lambertian terms we computed above are mutually
	//    independent, so we simply add them together, giving us a correct
	//    answer, no matter which display approach the user has chosen.
	result = vec4( lambert + ambient, 1.0f );
}