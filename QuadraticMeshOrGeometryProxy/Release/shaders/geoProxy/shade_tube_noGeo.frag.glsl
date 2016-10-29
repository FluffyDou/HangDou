#version 450
//#extension GL_NV_conservative_raster : enable

#define EPSILON        1e-20
#define EPSILON_SPHERE 1e-20
#define BIG_FLOAT      1e+20
uniform vec3  wireColor = vec3(0,1,0);
uniform vec3  fillColor = vec3(1);

uniform float cylinderRadiusSquare;
uniform float sphereRadiusSquare;
uniform float winSize;
uniform float viewBound;
uniform float znear;

uniform mat4 proj;

//uniform int shapeFlag;

// Texture coordinate
//in vec3 fragTexcoord;
//in float fragMatlID;

//in float w1;
//in float w2;
//in float w3;
in float w1;

//in vec3 vertFragPos0;
in vec3 vertFragPos1;
in vec3 vertFragPos2;
//in vec3 vertFragPos3;

//flat          in  vec4 wsPoint0;
//noperspective in  vec4 wsPoint1;

// Output color
out vec4 result;
//layout( depth_unchanged )out float gl_FragDepth;

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

// Intersect with three same radius spheres
void SphereIntersection(in vec3 center, in vec3 rayDir, out vec3 sNormal, out vec3 hitP, out float t, out bool discardFlag )
{
    t = BIG_FLOAT;
    float v = dot(center, rayDir);
    float d = sphereRadiusSquare - ( dot(center,center) - v*v);
    // to see if we intersect with the sphere
    discardFlag = (d < EPSILON_SPHERE) ? true : (t = v-sqrt(d), false);
    
    hitP    = t*rayDir;
    sNormal = normalize(hitP - center);    
}

// Intersect the cylinder
void CylinderIntersection(in vec3 p1, in vec3 p2, in vec3 rayDir, out vec3 cNormal, out vec3 hitP, out float t, out bool discardFlag)
{
    //// Intersect with the first cylinder
    vec3 AB = p2 - p1;
    //vec3 AO = -p0;
    vec3 AOxAB = cross(-p1,AB); // cross product
    vec3 VxAB  = cross(rayDir,AB); // cross product
    float ab2  = dot(AB,AB); // dot product
    float a    = dot(VxAB,VxAB); // dot product
    float b    = dot(VxAB,AOxAB); // dot product --- 2.0*dot(VxAB,AOxAB)
    float c    = dot(AOxAB,AOxAB) - (cylinderRadiusSquare * ab2);
    float d    = b*b - a*c; // b*b - 4.0*a*c

    // to see if we get intersection on the cylinder
	d = d < EPSILON ? 0.0 : sqrt(d);
    // Todo: check "a"
    t = min((-b + d)/a, (-b - d)/a);
    //t = t < EPSILON ? BIG_FLOAT : t;
	//discardFlag = d < EPSILON || t < EPSILON || dot(t*rayDir - p1, AB) < 0.0 ? t = BIG_FLOAT, true : false;
    discardFlag = d < EPSILON || t < EPSILON ? t = BIG_FLOAT, true : false;
    
    hitP = t*rayDir;
    vec3 axis = normalize(AB);

    // It turns out we do need to check that the hit point NOT falls on the cylinder extension
    float projOnCylinder = dot(hitP-p1, axis);
    discardFlag = projOnCylinder < 0.0 || projOnCylinder > length(AB) ? t = BIG_FLOAT, true : discardFlag;
    
    cNormal = normalize( hitP-p1 - projOnCylinder*axis );    
}

void main( void )
{
    /****************** Obtain the fragment's color *******************/
	// Look up data about this particular fragment's material
	//int  matlIdx    = 4 * int(fragMatlID+0.5);
	//vec4 matlInfo   = texelFetch( matlInfoTex, matlIdx+1 );  // Get diffuse color
    
    //vec4 albedo     = texture(textureArray, fragTexcoord);   // Get texture color
    //albedo = length(albedo.xyz) == 0 ? vec4(1.0, 1.0, 1.0, 1.0) : albedo;
    
	// Figure out the Lambertian color (if shadowed, set to 0)
    //float weight = max( 0.0f, dot( normalize(esFragNormal.xyz), toLight ) );
	//vec3 lambert = useLambertian * weight * matlInfo.xyz * albedo.xyz;
    // Add constant ambient
    //lambert += constAmbient * matlInfo.xyz * albedo.xyz;
    
    // Phong color
    //float phongValue = pow(dot(toEye, normalize(reflect( -toLight, normalize(esFragNormal.xyz)))), phongAlpha);
    //vec3  phong = usePhong * max(0.0, phongValue) * matlInfo.xyz * albedo.xyz;
    
	// Figure out the ambient color for no lambertian case
	//vec3 ambient = 0.5 * (1.0-useLambertian) * matlInfo.xyz * albedo.xyz;
    
    // Final shading
	//result = vec4( lightIntense*(lambert + ambient + phong), 1.0f );
    //result = length(lsFragFlatNormal) < EPSILON ? vec4(1,0,0,1) : result;
    //result = result.z == 0.1 ? vec4(normalize(esFragNormal.xyz),1.0): vec4(normalize(esFragNormal.xyz),1.0);    
    
    // Compute the pixel wise distance to three rasterized edges
    float w2 = 1.0 - w1;
    //vec3 p3 = wz > EPSILON ? (vertFragPos3-vertFragPos2-vertFragPos1)/wz  : vec3(0);
    vec3 p1 = vertFragPos1/w1;
    vec3 p2 = (vertFragPos2-vertFragPos1)/w2;
    // Find the closes vertex and only intersect with that one
    
    // calculate the direction of the intersection ray based on the screen coordinate
    vec2 vsGridCenter = vec2( gl_FragCoord.x, gl_FragCoord.y) / winSize; // normalized view space, [0,1]
    vsGridCenter = viewBound*( 2.0*vsGridCenter - vec2(1.0) ); // view space
    vec3 rayDir  = normalize( vec3(vsGridCenter, -znear) );
    
    float t;
    // If we need to discard the fragment due to its corresponding location on the cylinder
    bool discardFlag;
        
    // surface normal and the hit point on the cylinder
    vec3 hitP, normal;
    
    /*************** cylinder intersection *****************/    
    //if(1 == shapeFlag)
    if(w1 > EPSILON && w1 < 1.0 - EPSILON)   
        CylinderIntersection(p1, p2, rayDir, normal, hitP, t, discardFlag);                 
    else
        SphereIntersection(vertFragPos2, rayDir, normal, hitP, t, discardFlag);
	//ThreeSphereIntersection(p1, p2, p3, rayDir, normal_sphere, hitP_sphere, t_sphere, discardFlagSphere);
    
    if(discardFlag) discard;    
    //if(discardFlagCylinder && discardFlagSphere) discard;    
    
    // If we reach here, at least one of the intersection is valid.
    // Well, by default, we use cylinder intersection result
    // Sphere always win when ...... (do not let cylinder stick out of the sphere --- once in, it is there)
    // After this, we might not need cylinder body hit check
	//normal_cylinder = (t_sphere < t_cylinder && !discardFlagSphere) || discardFlagCylinder ?
	//(hitP_cylinder = hitP_sphere, normal_sphere) : normal_cylinder;     

    // compute the surface lambertian color --- the extra 0.5 serves as ambient color
	vec3 cColor = wireColor * (max( 0.0, dot(normal, -rayDir) ) + 0.15);    
    
    //// Todo: implement cylinder edge antialiasing like what normal wireframe does
    result = vec4(cColor, 1.0);
    // Simple display of wireframe
    //result = w.x < 0.01 || w.y < 0.01 || wz < 0.01 ? vec4(1.0, 0.0,0.0,1.0) : result;
    
    //  compute and write in the "real cylinder depth"
    vec4 tempFragPos = proj * vec4(hitP, 1.0);
    
    //gl_FragDepth = 0.5*tempFragPos.z/tempFragPos.w + 0.5;
    gl_FragDepth = discardFlag ?
                   result = vec4(1), gl_FragCoord.z : 
                   0.5*tempFragPos.z/tempFragPos.w + 0.5;
}







