#version 330

uniform sampler2D renderTex;

in  vec2 fragTex_coord;

out vec4 fragColor;


void main(void)
{

	fragColor = texture(renderTex, fragTex_coord);

}
