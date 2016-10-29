#version 450
//#extension GL_NV_conservative_raster : enable

#define EPSILON  1e-5 // 0.000001
uniform float lineWidth;
uniform vec3 wireColor = vec3(0);
uniform vec3 fillColor = vec3(1);

uniform float winSize;

// Texture coordinate
//in vec3 fragTexcoord;
//in float fragMatlID;

//in vec4 tmpColor;
//in float[100] st;
//in float w1;
//in float w2;
//in float w3;

//flat          in vec3 p1;
//noperspective in vec3 p2;
//smooth        in vec3 p3;
//smooth        in vec3 interpolatedNormal;

in vec2 w;

in vec2 vertFragPos1;
in vec2 vertFragPos2;
in vec2 vertFragPos3;

//flat          in  vec4 wsPoint0;
//noperspective in  vec4 wsPoint1;

// Output color
out vec4 result;

//uniform float smBufferRes;    // Shadow map resolution
//uniform float zBias;
//uniform float useLambertian;  // Asked to use a Lambertian shading scheme?
//uniform float usePhong;       // Asked to use a phong shading scheme?
//uniform float phongAlpha;     // Phong alpha value
//uniform float constAmbient;   // Constant ambient value

//uniform sampler2D shadowTex;  // Our shadow map texture.  

// The buffer containing our material information
//uniform samplerBuffer  matlInfoTex;
// The texture buffer holding the texutre for obj model
//uniform sampler2DArray textureArray;


void main( void )
{
    //vec3 v1 = normalize(p2-p1);
    //vec3 v2 = normalize(p3-p1);
    //vec3 flatNormal = cross(v1, v2);
    //flatNormal = length(flatNormal) < EPSILON ? normalize(interpolatedNormal) : normalize(flatNormal);
    //flatNormal = dot(flatNormal, interpolatedNormal) > 0 ? flatNormal : -flatNormal;
    //flatNormal = normalize(flatNormal);
    
    /****************** Obtain the fragment's color *******************/
	// Look up data about this particular fragment's material
	//int  matlIdx    = 4 * int(fragMatlID+0.5);
	//vec4 matlInfo   = texelFetch( matlInfoTex, matlIdx+1 );  // Get diffuse color
    
    //vec4 albedo     = texture(textureArray, fragTexcoord);   // Get texture color
    //albedo = length(albedo.xyz) == 0 ? vec4(1.0, 1.0, 1.0, 1.0) : albedo;
	
    /** Compute Adaptive Epsilon **/
    // Normalized depth value in shadow map
    //float SMDepth = texture( shadowTex, smTexCoord.xy ).x;

    /** Check the visibility **/
    //float isLit = SMDepth < actualDepth + actualBias ? 0.0 : 1.0;    
    // Shadow area, if the surface faces the light direction. EPSILON is a small number.
    //isLit = dot(lsFragNormal.xyz,lsFragPos.xyz) > EPSILON ? 0.0 : isLit;    
    //isLit = dot(lsFragFlatNormal.xyz,lsFragPos.xyz) > EPSILON ? 0.0 : isLit;
    
    // Set region out of light frustum to be dark
    //isLit = clamp(smTexCoord.xyz,0.0,1.0) != smTexCoord.xyz ? 0.0 : isLit;
    
    /***************** Shade the Fragment ******************************/
    
	// Find simple vectors for Lambertian shading
	//vec3 toLight = normalize( esLightPos.xyz - esFragPos.xyz);
	//vec3 toEye   = normalize( -esFragPos.xyz );
	
    // Add a little naive ambient for shadowed area
    //isLit = isLit == 0 && useLambertian != 0 ? 0.1 : isLit;
    
	// Figure out the Lambertian color (if shadowed, set to 0)
    //float weight = max( 0.0f, dot( normalize(esFragNormal.xyz), toLight ) );
	//vec3 lambert = isLit * useLambertian * weight * matlInfo.xyz * albedo.xyz;
    // Add constant ambient
    //lambert += constAmbient * matlInfo.xyz * albedo.xyz;
    
    // Phong color
    //float phongValue = pow(dot(toEye, normalize(reflect( -toLight, normalize(esFragNormal.xyz)))), phongAlpha);
    //vec3  phong = isLit * usePhong * max(0.0, phongValue) * matlInfo.xyz * albedo.xyz;
    
	// Figure out the ambient color for no lambertian case
	//vec3 ambient = 0.5 * isLit * (1.0-useLambertian) * matlInfo.xyz * albedo.xyz;
    
    // Final shading
	//result = vec4( lightIntense*(lambert + ambient + phong), 1.0f );
    //result = length(lsFragFlatNormal) < EPSILON ? vec4(1,0,0,1) : result;
    //result = result.z == 0.1 ? vec4(normalize(esFragNormal.xyz),1.0): vec4(normalize(esFragNormal.xyz),1.0);    
    //result = result.z != 0.177 ? vec4(1,0,0,1): vec4(0,1,0,1);
    
    // Compute the pixel wise distance to three rasterized edges
    //vec2  p0 = gl_FragCoord.xy;
    //float w3 = 1.0 - w.x - w.y;
    //vec2 p1 = (0.5 * (w.x > EPSILON ? vertFragPos1/w.x : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    //vec2 p2 = (0.5 * (w.y > EPSILON ? vertFragPos2/w.y : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    //vec2 p3 = (0.5 * (w3  > EPSILON ? vertFragPos3/w3  : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    vec2 p1 = (0.5 * vertFragPos1/w.x + vec2(0.5)) * winSize + vec2(0.5);
    vec2 p2 = (0.5 * vertFragPos2/w.y + vec2(0.5)) * winSize + vec2(0.5);
    vec2 p3 = (0.5 * vertFragPos3/(1.0-w.x-w.y) + vec2(0.5)) * winSize + vec2(0.5);
    float d1 = abs( ((p2.y-p1.y)*gl_FragCoord.x - (p2.x-p1.x)*gl_FragCoord.y + p2.x*p1.y - p2.y*p1.x)/length(p2-p1) );
    float d2 = abs( ((p3.y-p1.y)*gl_FragCoord.x - (p3.x-p1.x)*gl_FragCoord.y + p3.x*p1.y - p3.y*p1.x)/length(p3-p1) );
    float d3 = abs( ((p2.y-p3.y)*gl_FragCoord.x - (p2.x-p3.x)*gl_FragCoord.y + p2.x*p3.y - p2.y*p3.x)/length(p2-p3) );
    
    // Wireframe rendering line:
    //float thickness = 1.5;
    //result = min(d1, min(d2, d3)) < thickness ? vec4(1.0) : vec4(0.0,0.0,0.0,1.0);
    //result = w.x < EPSILON || w.y < EPSILON || w3 < EPSILON ? vec4(1,0,0,1) : result;
    float d = min(d1, min(d2, d3));
    float I = exp2(-lineWidth*d*d);
    result  = vec4( I*wireColor + (1.0 - I)*fillColor, 1.0);
    
    //result = result.x != -0.98763 ? vec4(2.0*flatNormal,1) : vec4(2.0*flatNormal,1);
}







