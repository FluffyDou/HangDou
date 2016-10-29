#version 450 
//#extension GL_NV_conservative_raster : enable

layout(location = IGLU_VERTEX)     in vec3 vertex;   // Sends vertex data from models here
//layout(location = IGLU_NORMAL)   in vec3 normal;   // Sends normal data (if any) from models here
//layout(location = IGLU_MATL_ID)  in float matlID;  // Sends ID for current material from model (if any)
//layout(location = IGLU_TEXCOORD) in vec2 texcoord; // Sends texture as albedo

uniform mat4 model;	     // Transforms: model space -> world space
uniform mat4 view;       // Transforms: world space -> eye space
//uniform mat4 proj;       // Transforms: eye space   -> clip space

/** real cylinder test **/
out  vertexData
{
	// for rasterize position
	vec3 center;	
    
} outData;

// render a cylinder like quad from a line
void main(void)
{
        // pass in the real position after transform into geometry shader
        outData.center  = (view * model * vec4(vertex, 1.0)).xyz;                
        
}