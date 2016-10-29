#version 450
//#extension GL_NV_conservative_raster : enable

#define EPSILON  1e-15 // 0.000001
uniform float lineWidth;
uniform vec3 wireColor = vec3(0);
uniform vec3 fillColor = vec3(1);

uniform float winSize;

// Texture coordinate
//in vec3 fragTexcoord;
//in float fragMatlID;

in vec2 vertFragPos1;
in vec2 vertFragPos2;
in vec2 vertFragPos3;

//flat          in  vec4 wsPoint0;
//noperspective in  vec4 wsPoint1;

// Output color
out vec4 result;

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
    // Compute the pixel wise distance to three rasterized edges
    //vec2 p0 = gl_FragCoord.xy;
    vec2 p1 = (0.5 * vertFragPos1 + vec2(0.5)) * winSize + vec2(0.5);
    vec2 p2 = (0.5 * vertFragPos2 + vec2(0.5)) * winSize + vec2(0.5);
    vec2 p3 = (0.5 * vertFragPos3 + vec2(0.5)) * winSize + vec2(0.5);
    
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
}
