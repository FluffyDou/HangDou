#version 420

#define EPSILON  0.0001

// Texture coordinate
in vec3 fragTexcoord;

in vec4  esFragPos;
in vec4  wsFragPos;
in vec4  esFragNormal;
in vec4  lsFragNormal;
in float fragMatlID;

// Output color
out vec4 result;

uniform int   adaptiveFlag;   // Use adaptive bias or not
uniform float constantBias;   // Constant depth bias based on Scene Scale
uniform float lightIntense;   // Light intensity
uniform float viewBound;      // View left/right/top/bottom (quad view)
uniform float smBufferRes;    // Shadow map resolution
uniform float zBias;
uniform float useLambertian;  // Asked to use a Lambertian shading scheme?
uniform float usePhong;       // Asked to use a phong shading scheme?
uniform float phongAlpha;     // Phong alpha value
uniform float constAmbient;   // Constant ambient value
uniform float lightNear;      // Light near distance
uniform vec4  esLightPos;      // Location of the light, relative to the eye
uniform mat4  lightView;       // Light view matrix
uniform mat4  lightProj;       // Light projection matrix

uniform sampler2D shadowTex;  // Our shadow map texture.  

// The buffer containing our material information
uniform samplerBuffer  matlInfoTex;

void main( void )
{
    /****************** Obtain the fragment's color *******************/
	// Look up data about this particular fragment's material
	int  matlIdx    = 4 * int(fragMatlID+0.5);
	vec4 matlInfo   = texelFetch( matlInfoTex, matlIdx+1 );  // Get diffuse color
    
    
    /***************** Visibility Check ******************************/
    // Note: Actually the light rays can be precomputed.
    // The light ray can therefore be obtained by looking up shadow map
    // texel center, followed by a texture lookup.
    
    /** Define tangent plane **/
    // Light space frag normal
    vec3 n = normalize(lsFragNormal.xyz);
    
    /** Locate shadow map texel center **/
	// Transform from eye-space to shadow map texture coordinates
    vec4 lsFragPos  = lightView * wsFragPos;    
    vec4 smTexCoord = lightProj * lsFragPos;
    smTexCoord /= smTexCoord.w;
    smTexCoord  = 0.5 * smTexCoord + vec4(0.5, 0.5, 0.5, 0.0);
    // Locate corresponding light space shadow map grid center
    vec2  index = floor( vec2( smTexCoord.xy * smBufferRes) );
    float delta = 1.0 / smBufferRes;
    // Normalized coordinate in [0,1]
    vec2  nlsGridCenter = delta*(index + vec2(0.5)); // Normalized eye space grid center --- [0,1]
    // Unnormalized coordinate in [-lightLeft,lightLeft]
    vec2 lsGridCenter = viewBound*( 2.0*nlsGridCenter - vec2(1.0) );
    
    /** Define light ray **/
    // Light ray direction in light space
    vec3 lsGridLineDir = normalize( vec3(lsGridCenter, -lightNear) ); // Light space grid line direction    
    
    /** Plane ray intersection **/
    // Locate the potential occluder for the shading fragment
    float ls_t_hit = dot(n, lsFragPos.xyz) / dot(n, lsGridLineDir);
    vec3  ls_hit_p = ls_t_hit * lsGridLineDir;
    
    /** Compute Adaptive Epsilon **/
    // Normalized depth value in shadow map
    float SMDepth = texture( shadowTex, smTexCoord.xy ).x;
    // A and B are computed bnased on light near and far planes.
    // They can be retrieved directly from light projection matrix
    float A = lightProj[2][2];
    float B = lightProj[3][2];    
    float adaptiveDepthBias = 0.5*pow(1.0 - A - 2.0*SMDepth, 2)*constantBias / B; 
	
    // Use the intersection point as new look up point
    vec4 lsPotentialoccluder = lightProj * vec4(ls_hit_p, 1.0);
    lsPotentialoccluder      = lsPotentialoccluder/lsPotentialoccluder.w;
    lsPotentialoccluder      = 0.5 * lsPotentialoccluder + vec4(0.5, 0.5, 0.5, 0.0);
   
    float actualDepth = min(lsPotentialoccluder.z, smTexCoord.z);
    float actualBias  = adaptiveDepthBias;
        
    // Constant depth bias VS adaptive depth bias
    actualDepth = adaptiveFlag != 0 ? actualDepth : smTexCoord.z;
    actualBias  = adaptiveFlag != 0 ? actualBias  : zBias;    
    
    /** Check the visibility **/
    float isLit = SMDepth < actualDepth + actualBias ? 0.0 : 1.0;    
    // Shadow area, if the surface faces the light direction. EPSILON is a small number.
    isLit = dot(lsFragNormal.xyz,lsFragPos.xyz) > EPSILON ? 0.0 : isLit;    
    // Set region out of light frustum to be dark
    isLit = clamp(smTexCoord.xyz,0.0,1.0) != smTexCoord.xyz ? 0.0 : isLit;
    
    /***************** Shade the Fragment ******************************/
    
	// Find simple vectors for Lambertian shading
	vec3 toLight = normalize( esLightPos.xyz - esFragPos.xyz);
	vec3 toEye   = normalize( -esFragPos.xyz );
	
    // Add a little naive ambient for shadowed area
    isLit = isLit == 0 && useLambertian != 0 ? 0.1 : isLit;
    
	// Figure out the Lambertian color (if shadowed, set to 0)
    float weight = max( 0.0f, dot( normalize(esFragNormal.xyz), toLight ) );
	vec3 lambert = isLit * useLambertian * weight * matlInfo.xyz;
    // Add constant ambient
    lambert += constAmbient * matlInfo.xyz;
    
    // Phong color
    float phongValue = pow(dot(toEye, normalize(reflect( -toLight, normalize(esFragNormal.xyz)))), phongAlpha);
    vec3  phong = isLit * usePhong * max(0.0, phongValue) * matlInfo.xyz;
    
	// Figure out the ambient color for no lambertian case
	vec3 ambient = 0.5 * isLit * (1.0-useLambertian) * matlInfo.xyz;
    
    // Final shading
	result = vec4( lightIntense*(lambert + ambient + phong), 1.0f );
}
