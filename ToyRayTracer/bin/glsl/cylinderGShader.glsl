#version 400
//#extension GL_ARB_geometry_shader4 : enable

layout (lines) in;
layout (triangle_strip, max_vertices = 8) out;

uniform vec3 eye_position;

in  vertexData
{
	// for rasterize position
	vec4 up_coord;
	vec4 down_coord;
	vec4 center_point;
	float line_id;
	float click_id;
		
} vert[];


out float edgeColorR;
out float line_id;
out float click_id;

void main()
{
/*	int i;
	
	for(i = 0; i < gl_VerticesIn; i++)
	{
		gl_Position = gl_PositionIn[i];
		EmitVertex();
	}
	EndPrimitive();*/
	
	// draw up quad.
	gl_Position = vert[0].up_coord;
	line_id     = vert[0].line_id;
	click_id    = vert[0].click_id;
	edgeColorR  = 0.1;		// cylinder edge color
	EmitVertex();
	gl_Position = vert[0].center_point;
	line_id     = vert[0].line_id;
	click_id    = vert[0].click_id;
	edgeColorR  = 1.0;		// cylinder center color
	EmitVertex();
	gl_Position = vert[1].up_coord;
	line_id     = vert[1].line_id;
	click_id    = vert[1].click_id;
	edgeColorR  = 0.1;		// cylinder edge color
	EmitVertex();
	gl_Position = vert[1].center_point;
	line_id     = vert[1].line_id;
	click_id    = vert[1].click_id;
	edgeColorR  = 1.0;		// cylinder center color
	EmitVertex();

	EndPrimitive();
	
	// draw down quad
	gl_Position = vert[0].down_coord;
	line_id     = vert[0].line_id;
	click_id    = vert[0].click_id;
	edgeColorR  = 0.1;
	EmitVertex();
	gl_Position = vert[0].center_point;
	line_id     = vert[0].line_id;
	click_id    = vert[0].click_id;
	edgeColorR  = 1.0;
	EmitVertex();
	gl_Position = vert[1].down_coord;
	line_id     = vert[1].line_id;
	click_id    = vert[1].click_id;
	edgeColorR  = 0.1;
	EmitVertex();
	gl_Position = vert[1].center_point;
	line_id     = vert[1].line_id;
	click_id    = vert[1].click_id;
	edgeColorR  = 1.0;
	EmitVertex();

	EndPrimitive();

}
