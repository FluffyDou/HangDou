#version 420
// This is not a very robust implementation of omni shadow map. The dual parabola always point to (0,1,0) and (0,-1,0)
// in world space. We do this only to test adaptive depth bias applied to omni shadow map

// Tells the OpenGL/IGLU program to enable the depth test when using this shader
#pragma IGLU_ENABLE  DEPTH_TEST

#define EPSILON  0.0001

// Information from the vertex shader about the eye-space location
//    of this fragment we're working on
in vec3  fragTexcoord;
in vec4  esFragPos;
in vec4  wsFragPos;
in vec4  wsFragNormal;
in vec4  esFragNormal;
in float fragMatlID;

// Where to store the color we compute
out vec4 result;

uniform int   adaptiveFlag;      // Use adaptive bias or not

uniform float constantBias;      // Constant depth bias based on Scene Scale
uniform float useLambertian;     // Asked to use a lambertian shading scheme
uniform float usePhong;          // Asked to use a phong shading scheme
uniform float smBufferRes;       // Shadow map resolution
uniform float nearDist, farDist; // Light near & far planes
uniform float lightIntense;      // Light intensity
uniform float constAmbient;      // Constant ambient value
uniform float phongAlpha;        // Phong alpha value

// Bias for shadow map
uniform float zBias;

// Location of the light, relative to the eye
uniform vec4 wsLightPos;
uniform vec4 esLightPos;

// Our shadow map texture.  Layer 0 = positive z-hemisphere, layer 1 = negative
uniform sampler2DArray shadowTex;

// The buffer containing our material information
uniform samplerBuffer  matlInfoTex;

// This is the actual code that is executed.
void main( void )
{	
    /****************** Obtain the fragment's color *******************/
	// Look up data about this particular fragment's material
	int  matlIdx     = 4 * int(fragMatlID+0.5);
	vec4 matlInfo   = texelFetch( matlInfoTex, matlIdx+1 );  // Get diffuse color
    
    /***************** Visibility Check ******************************/
    
    /** Define tangent plane **/
    vec3 n = normalize(wsFragNormal.xyz); // World space frag normal
    
    /** Locate shadow map texel center **/
    //// Compute the coordinate in paraboloid space --- the length is the same as world space
	vec3  toVert = wsFragPos.xyz - wsLightPos.xyz;
	float dist   = length( toVert );
	toVert       = normalize( toVert );
	
	// Compute the coordinates in the + & - hemi-paraboloids
	vec2 posMap = vec2(toVert.xy)/(toVert.z+1);
	vec2 negMap = vec2(toVert.xy)/(toVert.z-1);
	
	// Select if we're using the + or - hemi-paraboloid for omni shadow map.
	bool usePos = dot(posMap,posMap) < dot(negMap,negMap);
	
	// Now map this into texture space (i.e., [0..1])
	posMap = 0.5*posMap + 0.5;
	negMap = 0.5*negMap + 0.5;

    // Select valid texture coordinate --- on pos or neg hemi-sphere
    vec2 smTexCoord = usePos ? posMap : negMap;
    
    // Locate the texel grid center
    vec2  index = floor( vec2( smTexCoord * smBufferRes) );
    float delta = 1.0 / smBufferRes;
    // Normalized eye space grid center in [0,1]
    vec2  nlsGridCenter = delta*(index + vec2(0.5));    
    // The light space is just a traslation from the world space.
    // In our implementation, the two hemi paraboloids always face positive and negative z axis
    nlsGridCenter = 2.0*nlsGridCenter - vec2(1.0); // From [0,1] to [-1,1]
    
    // Compute the z value in light space on paraboloid
    float nlsGridCenterZ = 0.5*(1.0 - dot(nlsGridCenter, nlsGridCenter));
    // Compute the paraboloid normal in light space, which is the same in world space since there is no rotation between two spaces
    vec3  nlsGridCenterN = normalize( vec3(nlsGridCenter/nlsGridCenterZ, 1.0/nlsGridCenterZ) );
    nlsGridCenterN = usePos ? nlsGridCenterN : -nlsGridCenterN; 
    // Choose the vector based on correct hem-sphere (d0)
    vec3  d = usePos ? vec3(0,0,1) : vec3(0,0,-1);
    
    /** Define light ray **/
    // Compute the light ray direction of texel center
    vec3  v = -normalize( reflect(d, nlsGridCenterN) );
    // Locate the shadow map grid center point in world space: translate by light offset
    vec3  p = wsLightPos.xyz + vec3(nlsGridCenter, nlsGridCenterZ);
    
    /** Plane ray intersection **/
    // Apply world space line plane intersection
    float ws_t_hit = dot(n, wsFragPos.xyz-p) / dot(n, v);
    vec3  ws_hit_p = p + ws_t_hit * v;
    // Locate potential occluder
    float potentialDist = distance(ws_hit_p, wsLightPos.xyz);
    
	// Map depth of potential occluder (potentialDist) to [0,1] to compare with shadow map's value
	potentialDist = 0.5 * ((nearDist+farDist)*(-potentialDist)+(2.0*nearDist*farDist))/(potentialDist*(nearDist-farDist)) + 0.5;
    // Depth of constant depth bias
    dist          = 0.5 * ((nearDist+farDist)*(-dist)+(2.0*nearDist*farDist))/(dist*(nearDist-farDist)) + 0.5;    
    
    /** Compute Adaptive Epsilon **/ 
    // Normalized depth value in shadow map
    float SMDepth = texture( shadowTex, vec3( usePos ? vec3(posMap, 0) : vec3(negMap, 1)) ).x;
    float A = -(farDist + nearDist) / (farDist - nearDist);
    float B = -2.0*farDist*nearDist / (farDist - nearDist);
    float adaptiveDepthBias = 0.5*pow(1.0 - A - 2.0*SMDepth, 2)*constantBias / B;
    
    //float actualDepth = potentialDist;
    float actualDepth = min(potentialDist,dist);
    //float actualDepth = potentialDist;
    float actualBias  = adaptiveDepthBias;    
    
    // Constant depth bias VS adaptive depth bias
    actualDepth = adaptiveFlag != 0 ? actualDepth : dist;
    actualBias  = adaptiveFlag != 0 ? actualBias  : zBias;     
    
    /** Check the visibility **/
    float isLit = SMDepth > actualDepth + actualBias ? 1.0 : 0.0;
    
    // Filter out out invalid value
    isLit = clamp(dist,0.0,1.0) != dist || clamp(smTexCoord,0.0,1.0) != smTexCoord ? 0.0 : isLit;   
    
    /***************** Shade the Fragment ******************************/
	// Find simple vectors for basic shading
	vec3 toLight = normalize( esLightPos.xyz - esFragPos.xyz);
	vec3 toEye   = normalize( -esFragPos.xyz );
    
    // Add a little coarse ambient for shadowed area
    isLit = isLit == 0 && useLambertian != 0 ? 0.1 : isLit;    
    
	// Figure out the Lambertian color (if shadowed, set to 0)
    float weight = max( 0.0f, dot( normalize(esFragNormal.xyz), toLight ) );
	vec3 lambert = isLit * useLambertian * weight * matlInfo.xyz;
    // Add a little ambient when lambertian
	lambert += constAmbient * matlInfo.xyz;
    
    // Phong color
    float phongValue = pow(dot(toEye, normalize(reflect( -toLight, normalize(esFragNormal.xyz)))), phongAlpha);
    vec3  phong = isLit * usePhong * max(0.0, phongValue) * matlInfo.xyz;    
    
	// Figure out the ambient color (if shadowed, set to 0)
	vec3 ambient = isLit * (1.0-useLambertian) * matlInfo.xyz;

    // Final shading
    result = vec4( lightIntense*(lambert + ambient + phong), 1.0f );
}