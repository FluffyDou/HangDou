#version 330

layout(location = 0) in vec2 vertex_coord;
layout(location = 1) in vec2 tex_coord;

out vec2 fragTex_coord;


void main(void)
{
	gl_Position   = vec4(vertex_coord, 0.0, 1.0);
	fragTex_coord = tex_coord;
}
