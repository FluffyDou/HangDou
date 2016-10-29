#version 450
//#extension GL_NV_conservative_raster : enable

#define EPSILON        1e-15
#define EPSILON_SPHERE 1e-15
#define BIG_FLOAT      1e+10
uniform vec3  wireColor = vec3(0,1,0); // vec3(1);
uniform vec4  fillColor = vec4(1); // vec4(0,0,0,1);

//uniform float znear;
//uniform float viewBound;
uniform float cylinderRadiusSquare;
uniform float sphereRadiusSquare;
//uniform float winSize;

uniform mat4 proj;

// Texture coordinate
//in vec3 fragTexcoord;
//in float fragMatlID;
//flat in int testID;
//in float w1;
//in float w2;
//in float w3;
in vec2 w;

//in vec3 vertFragPos0;
in vec3 vertFragPos1;
in vec3 vertFragPos2;
in vec3 vertFragPos3;

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
void ThreeSphereIntersection(in vec3 center1, in vec3 center2, in vec3 center3, in vec3 rayDir, out vec3 sNormal, out vec3 hitP, out float t, out bool discardFlag )
{
    float rSquare = sphereRadiusSquare;
    t = BIG_FLOAT;
    vec3 center = center1;
    float v = dot(center1, rayDir);
    float d = rSquare - ( dot(center1,center1) - v*v);
    // to see if we intersect with the sphere
    discardFlag = (d < EPSILON_SPHERE) ? true : (t = v-sqrt(d), false);
    
    v = dot(center2, rayDir);
    d = rSquare - ( dot(center2,center2) - v*v);
    float tempT = v-sqrt(d);
    // to see if we intersect with the sphere
    discardFlag = (d < EPSILON_SPHERE) ? discardFlag : (t=t>tempT?center=center2,tempT:t, false);
    
    v = dot(center3, rayDir);
    d = rSquare - ( dot(center3,center3) - v*v);
    tempT = v-sqrt(d);
    // to see if we intersect with the sphere
    discardFlag = (d < EPSILON_SPHERE) ? discardFlag : (t=t>tempT?center=center3,tempT:t, false);
    
    hitP    = t*rayDir;
    sNormal = normalize(hitP - center);    
}

// Intersect the cylinder
void ThreeCylinderIntersection(in vec3 p1, in vec3 p2, in vec3 p3, in vec3 rayDir, out vec3 cNormal, out vec3 hitP, out float t, out bool discardFlag)
{
    float rSquare = cylinderRadiusSquare;
    bool flag2, flag3;
    vec3 p = p1;
    //// Intersect with the first cylinder
    vec3 AB = p2 - p1;
    //vec3 AO = -p0;
    vec3 AOxAB = cross(-p1,AB); // cross product
    vec3 VxAB  = cross(rayDir,AB); // cross product
    float ab2  = dot(AB,AB); // dot product
    float a    = dot(VxAB,VxAB); // dot product
    float b    = dot(VxAB,AOxAB); // dot product --- 2.0*dot(VxAB,AOxAB)
    float c    = dot(AOxAB,AOxAB) - (rSquare * ab2);
    float d    = b*b - a*c; // b*b - 4.0*a*c

    // to see if we get intersection on the cylinder
	d = d < EPSILON ? 0.0 : sqrt(d);
    // Todo: check "a"
    t = min((-b + d)/a, (-b - d)/a);
    //t = t < EPSILON ? BIG_FLOAT : t;
	//discardFlag = d < EPSILON || t < EPSILON || dot(t*rayDir - p1, AB) < 0.0 ? t = BIG_FLOAT, true : false;
    discardFlag = d < EPSILON || t < EPSILON ? t = BIG_FLOAT, true : false;
    
    //// Intersect with the second cylinder
    vec3 AB2 = p3 - p2;
    AOxAB = cross(-p2,AB2); // cross product
    VxAB  = cross(rayDir,AB2); // cross product
    ab2   = dot(AB2,AB2); // dot product
    a     = dot(VxAB,VxAB); // dot product
    b     = dot(VxAB,AOxAB); // dot product --- 2.0*dot(VxAB,AOxAB)
    c     = dot(AOxAB,AOxAB) - (rSquare * ab2);
    d     = b*b - a*c; // b*b - 4.0*a*c
    // to see if we get intersection on the cylinder
	d = d < EPSILON ? 0.0 : sqrt(d);
    float t2 = min((-b + d)/a, (-b - d)/a);
    //t2 = t2 < EPSILON ? BIG_FLOAT : t2;
	//flag2 = d < EPSILON || t2 < EPSILON || dot(t2*rayDir - p2, AB2) < 0.0 ? t2 = BIG_FLOAT, true : false;
    flag2 = d < EPSILON || t2 < EPSILON ? t2 = BIG_FLOAT, true : false;
    
    //// Intersect with the thrid cylinder
    vec3 AB3 = p1 - p3;
    AOxAB = cross(-p3,AB3); // cross product
    VxAB  = cross(rayDir,AB3); // cross product
    ab2   = dot(AB3,AB3); // dot product
    a     = dot(VxAB,VxAB); // dot product
    b     = dot(VxAB,AOxAB); // dot product --- 2.0*dot(VxAB,AOxAB)
    c     = dot(AOxAB,AOxAB) - (rSquare * ab2);
    d     = b*b - a*c; // b*b - 4.0*a*c
    // to see if we get intersection on the cylinder
	d = d < EPSILON ? 0.0 : sqrt(d);
    float t3 = min((-b + d)/a, (-b - d)/a);
    //t3 = t3 < EPSILON ? BIG_FLOAT : t3;
	//flag3 = d < EPSILON || t3 < EPSILON || dot(t3*rayDir - p3, AB3) < 0.0 ? t3 = BIG_FLOAT, true : false;     
    flag3 = d < EPSILON || t3 < EPSILON ? t3 = BIG_FLOAT, true : false;
    
    //// Pick the closest valid hit
    // "AB" is cylinder axis (unormalized)
    // "p" is the cylinder bottom end 
    AB = !flag2 && t2 < t ? discardFlag = flag2, t = t2, p = p2, AB2 : AB;
    AB = !flag3 && t3 < t ? discardFlag = flag3, t = t3, p = p3, AB3 : AB;
    
    hitP = t*rayDir;
    vec3 axis = normalize(AB);

    // It turns out we do need to check that the hit point NOT falls on the cylinder extension
    float projOnCylinder = dot(hitP-p, axis);
    discardFlag = projOnCylinder < 0.0 || projOnCylinder > length(AB) ? t = BIG_FLOAT, true : discardFlag;
    
    cNormal = normalize( hitP-p - projOnCylinder*axis );    
}


// Intersect the cylinder
void CylinderIntersection(in vec3 p0, in vec3 p1, in vec3 rayDir, out vec3 cNormal, out vec3 hitP, out float t, out bool discardFlag)
{
    vec3 AB = p1 - p0;
    //vec3 AO = -p0;
    vec3 AOxAB = cross(-p0,AB); // cross product
    vec3 VxAB  = cross(rayDir,AB); // cross product
    float ab2  = dot(AB,AB); // dot product
    float a    = dot(VxAB,VxAB); // dot product
    float b    = dot(VxAB,AOxAB); // dot product --- 2.0*dot(VxAB,AOxAB)
    float c    = dot(AOxAB,AOxAB) - (cylinderRadiusSquare * ab2);
    float d    = b*b - a*c; // b*b - 4.0*a*c

    // to see if we get intersection on the cylinder
	d = d < EPSILON ? 0.0 : sqrt(d);
    
	//a = 2.0*a;
	//float t1 = (-b + d)/a; // hit point on the body
	//float t2 = (-b - d)/a; // hit point on the body
    // Todo: check "a"
    t = min((-b + d)/a, (-b - d)/a);
    
	discardFlag = d < EPSILON ? true : false;
    //discardFlag = false;
    
    hitP = t*rayDir;
    vec3 axis = normalize(AB);
    cNormal = normalize( hitP-p0 - dot(hitP-p0, axis)*axis );
}
/*
void CylinderIntersection(in vec3 p0, in vec3 p1, in vec3 rayDir, out vec3 cNormal, out vec3 hitP, out float t, out int discardFlag)
{
    // generate cylinder space
	vec3 axis    = p1 - p0; // p0 is the origin for the cylinder space
	vec3 vector  = cross( axis, vec3(1.0, 0.0, 0.0) );
	vector       = length(vector) > 0.01 ? vector : cross( axis, vec3(0.0, 1.0, 0.0) );
	vec3 uVec    = normalize( cross(axis, vector) );
	vec3 vVec    = normalize( cross(axis, uVec) );
    float height = length(axis);
	axis         = normalize(axis);
	//float height = distance(geo.p1, geo.p0);
	// get the eye coordinate and the ray direction value in the cylinder space
	vec3 temp = -p0;
	vec3 orig = vec3( dot(temp, uVec),   dot(temp, vVec),   dot(temp, axis) );
    vec3 dir  = vec3( dot(rayDir, uVec), dot(rayDir, vVec), dot(rayDir, axis) );
	// find the location on the cylinder through ray intersection     
	float temp_a = (geo.p2.y - geo.p2.x)*(geo.p2.y - geo.p2.x)/(height*height);
	float temp_b = 2.0*geo.p2.x*(geo.p2.y - geo.p2.x)/height;
	float temp_c = geo.p2.x*geo.p2.x;
	float a = dir.x*dir.x + dir.y*dir.y - dir.z*dir.z*temp_a; // if the view ray is parallel with the axis, then a = 0
	float b = 2.0*( orig.x*dir.x + orig.y*dir.y - temp_a*orig.z*dir.z) - temp_b*dir.z;
	float c = orig.x*orig.x + orig.y*orig.y - temp_a*orig.z*orig.z - temp_b*orig.z - temp_c;
	float d = b*b - 4.0*a*c;
	        
	discardFlag = d < EPSILON ? 1 : 0; // to see if we get intersection on the cylinder
		
	d = d < EPSILON ? 0.0 : sqrt(d);
	a = 2.0*a;
	        
	float t1 = (-b + d)/a;               // hit point on the body
	float t2 = (-b - d)/a;               // hit point on the body
    float t3 = (-orig.z)/dir.z;          // hit point on the bottom
    float t4 = (height - orig.z)/dir.z;  // hit point on the top
        
    float d3 = length( (orig+t3*dir).xy );
    float d4 = length( (orig+t4*dir).xy );         
	float z1 = orig.z + t1*dir.z;
	float z2 = orig.z + t2*dir.z;
	        
    discardFlag = ( (z1>height || z1<0) && (z2>height || z2<0) && d3 > geo.p2.x && d4 > geo.p2.y ) ? 1 : discardFlag;
    t = 1000000000.0;
    // find the closet hit on the platform body
    t = (t1 > 0 && t1 < t && z1 < height && z1 > 0) ? t1 : t;               
    t = (t2 > 0 && t2 < t && z2 < height && z2 > 0) ? t2 : t;        
    t = (t3 > 0 && t3 < t && d3 < geo.p2.x)         ? t3 : t; // if we comment out this line, the platform will miss the bottom
    t = (t4 > 0 && t4 < t && d4 < geo.p2.y)         ? t4 : t; // if we comment out this line, the platform will miss the top

    vec3 hitP_view = t*rayDir; // get the hit point in the view space
    vec3 hitP_cylin = orig + t*dir; // get the hit point in the cylinder space 
        
    // compute the normal of platform body                             
    vec3 tip = geo.p2.x > geo.p2.y ? (geo.p0 + (height*geo.p2.x / (geo.p2.x - geo.p2.y))*axis) 
                                   : (geo.p1 - (height*geo.p2.y / (geo.p2.y - geo.p2.x))*axis);        
    vec3 tempNormal = geo.p2.x > geo.p2.y ? 
    ( normalize( dot( normalize(tip - hitP_view), (geo.p1 - hitP_view) )*normalize(tip - hitP_view) + hitP_view - geo.p1 ) ) :
    ( normalize( dot( normalize(tip - hitP_view), (geo.p0 - hitP_view) )*normalize(tip - hitP_view) + hitP_view - geo.p0 ) ) ;
    
    // if the trancated cone degenerate to a cylinder, we get the normal in the cylinder way
    tempNormal = abs( geo.p2.x - geo.p2.y ) < 1e-3f ? normalize( hitP_view-geo.p0 - dot(hitP_view-geo.p0, axis)*axis ) : tempNormal;        
    // find the normal of the hit point --- the hit point can be on the top lid, bottom lid or on the body        
    cNormal = axis;        
    cNormal = ( hitP_cylin.z > EPSILON1 && hitP_cylin.z + EPSILON1 < height ) ? tempNormal : cNormal;     
    cNormal = ( hitP_cylin.z < EPSILON1 )                                     ? -axis      : cNormal;     
           
    hitP = hitP_view;   
        
    discardFlag = ( dot(hitP-geo.lastVert, geo.lastNorm) < 0 || dot(hitP-geo.nextVert, geo.nextNorm) < 0 ) ? 1 : discardFlag; 
                       
}
*/

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
    //result = result.z != 0.177 ? vec4(1,0,0,1): vec4(0,1,0,1);
    
    // Compute the pixel wise distance to three rasterized edges
    float wz = 1.0 - (w.x + w.y);
    //vec3 p1 = w.x > EPSILON ? vertFragPos1/w.x : vec3(0);
    //vec3 p2 = w.y > EPSILON ? vertFragPos2/w.y : vec3(0);
    //vec3 p3 = wz  > EPSILON ? (vertFragPos3-vertFragPos2-vertFragPos1)/wz  : vec3(0);
    vec3 p1 = vertFragPos1/w.x;
    vec3 p2 = vertFragPos2/w.y;
    vec3 p3 = (vertFragPos3-vertFragPos2-vertFragPos1)/wz;
    // Find the closes vertex and only intersect with that one
    
    //vec2 p1 = (0.5 * (w.x > EPSILON ? vertFragPos1/w.x : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    //vec2 p2 = (0.5 * (w.y > EPSILON ? vertFragPos2/w.y : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    //vec2 p3 = (0.5 * (w3  > EPSILON ? vertFragPos3/w3  : vec2(0)) + vec2(0.5)) * winSize + vec2(0.5);
    
    // calculate the direction of the intersection ray based on the screen coordinate
    //vec2 vsGridCenter = vec2( gl_FragCoord.x, gl_FragCoord.y) / winSize; // normalized view space, [0,1]
    //vsGridCenter = viewBound*( 2.0*vsGridCenter - vec2(1.0) ); // view space
    //vec3 rayDir  = normalize( vec3(vsGridCenter, -znear) );
    vec3 rayDir = normalize(vertFragPos3);
    
    float t_cylinder;
    float t_sphere;
    // If we need to discard the fragment due to its corresponding location on the cylinder
    bool discardFlagCylinder;
    bool discardFlagSphere;
        
    // surface normal and the hit point on the cylinder
    vec3 normal_cylinder, hitP_cylinder;
    vec3 normal_sphere, hitP_sphere;  
    
    /*************** cylinder intersection *****************/       
    //CylinderIntersection(p1, p2, rayDir, normal_cylinder, hitP_cylinder, t_cylinder, discardFlagCylinder );
    
    ThreeCylinderIntersection( p1, p2, p3, rayDir, normal_cylinder, hitP_cylinder, t_cylinder, discardFlagCylinder ); 
	ThreeSphereIntersection(p1, p2, p3, rayDir, normal_sphere, hitP_sphere, t_sphere, discardFlagSphere);
    
    //if(discardFlagCylinder && discardFlagSphere) discard;    
    
    // If we reach here, at least one of the intersection is valid.
    // Well, by default, we use cylinder intersection result
    // Sphere always win when ...... (do not let cylinder stick out of the sphere --- once in, it is there)
    // After this, we might not need cylinder body hit check
	normal_cylinder = (t_sphere < t_cylinder && !discardFlagSphere) || discardFlagCylinder ?
	(hitP_cylinder = hitP_sphere, normal_sphere) : normal_cylinder;     

    // compute the surface lambertian color --- the extra 0.5 serves as ambient color
	vec3 cColor = wireColor * (max( 0.0, dot(normal_cylinder, -rayDir) ) + 0.15);    
    
    //// Todo: implement cylinder edge antialiasing like what normal wireframe does
    result = vec4( cColor, 1.0);
    // Simple display of wireframe
    //result = w.x < 0.01 || w.y < 0.01 || wz < 0.01 ? vec4(1.0, 0.0,0.0,1.0) : result;
    
    //  compute and write in the "real cylinder depth"
    vec4 tempFragPos = proj * vec4(hitP_cylinder, 1.0);
    
    //gl_FragDepth = 0.5*tempFragPos.z/tempFragPos.w + 0.5;
    gl_FragDepth = discardFlagCylinder && discardFlagSphere ? 
                           result = fillColor, gl_FragCoord.z : 
                           0.5*tempFragPos.z/tempFragPos.w + 0.5;   
}







