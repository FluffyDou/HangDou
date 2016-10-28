#version 400

in      vec3 vertex_coord;

uniform mat4 transform;

void main(void)
{
	gl_Position =  transform * vec4(vertex_coord, 1.0);
}
