#version 450
//#extension GL_EXT_geometry_shader4 : enable
//#extension GL_NV_geometry_shader_passthrough : enable

layout (lines) in;
layout (triangle_strip, max_vertices = 4) out;

uniform float radius;
uniform mat4 proj;

in  vertexData
{
	vec3 center;
} vert[];

//out vec3 vertFragPos0;
out vec3 vertFragPos1;
out vec3 vertFragPos2;

// generate three quads to represent a cylinder --- 4 triangles for each quad
struct vertData{

	vec3  center;
	vec3  up;
	vec3  down;
	vec3  upFront;
	vec3  downFront;
	vec3  upBack;      
	vec3  downBack;
	//float radius;

} vert0, vert1;


void main()
{
    vert0.center = vert[0].center;
    vert1.center = vert[1].center;
	
	vec3 axis     = normalize( vert1.center - vert0.center ); // cylinder axis vector
    vec3 extDir   = normalize( cross( vert0.center, -axis ) ); // generated vertices direction	
	vec3 frontDir = normalize( cross( extDir, -axis ) );
	
    // Extend the quad region --- for debug
    //float tempScale  = 1.0f;
    
	// extend the center a little bit for inclined plane
	vert0.center = vert0.center - radius*axis;
	vert1.center = vert1.center + radius*axis;	
    
	// bottom
	vert0.up         = vert0.center + radius * extDir;
	vert0.down       = vert0.center - radius * extDir;
	vert0.upFront    = vert0.up     + radius * frontDir;
	vert0.upBack     = vert0.up     - radius * frontDir;	
	vert0.downFront  = vert0.down   + radius * frontDir;  
	vert0.downBack   = vert0.down   - radius * frontDir;   
	// top
	vert1.up         = vert1.center + radius * extDir;
	vert1.down       = vert1.center - radius * extDir;
	vert1.upFront    = vert1.up     + radius * frontDir;
	vert1.upBack     = vert1.up     - radius * frontDir;	
	vert1.downFront  = vert1.down   + radius * frontDir;  
	vert1.downBack   = vert1.down   - radius * frontDir;
		
	/*********** generate the front quad *************/
	
    // cylinder parameters
    vertFragPos1 = vert[0].center;
    vertFragPos2 = vert[1].center;     
        
    // vert1 is the end we see cap(top), vert0 is the end we see curve(bottom)
    // --- the end close to eye always shows the cap        
    float dis0 = dot( vert[0].center, vert[0].center );
    float dis1 = dot( vert[1].center, vert[1].center);
      	
    if(dis0 > dis1)	
    {
        //vertFragPos0 = vert0.upFront;
        gl_Position = proj * vec4( vert0.upFront, 1.0 );
        EmitVertex();               
        gl_Position = proj * vec4( vert0.downFront, 1.0 );
        EmitVertex();     
        gl_Position = proj * vec4( vert1.upBack, 1.0 );
        EmitVertex();
        gl_Position = proj * vec4( vert1.downBack, 1.0 );
        EmitVertex();
    }
    else
    {
        gl_Position = proj * vec4( vert0.upBack, 1.0 );
        EmitVertex();
        gl_Position = proj * vec4( vert0.downBack, 1.0 );
        EmitVertex();
        gl_Position = proj * vec4( vert1.upFront, 1.0 );
        EmitVertex();
        gl_Position = proj * vec4( vert1.downFront, 1.0 );
        EmitVertex();         
    }
    
    EndPrimitive();

}
