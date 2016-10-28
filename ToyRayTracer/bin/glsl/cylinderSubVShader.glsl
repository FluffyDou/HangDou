#version 400

//uniform mat4 transform;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 eye_position;
uniform vec3 small_box_center;
uniform float bounding_radius;

in vec3 vertex_coord;
in vec3 tangent;
in float line_id;
in float click_id;

out  vertexData 
{
	// for rasterize position
	vec4 up_coord;
	vec4 down_coord;
	vec4 center_point;
	
	int flag; // specify if draw the line or not in the G shader
	float line_id;
	float click_id;
	
} outData;


void main(void)
{
	vec3 vertex_new  = ( model * vec4(vertex_coord, 1.0) ).xyz; // the real vertex position after transform
	vec3 tmpVector   = vertex_coord + normalize(tangent);       // temp vector use to computer the tangent value after the transform
	vec3 tangent_new = ( model * vec4(tmpVector, 1.0) ).xyz - vertex_new;            // new tangent value after the transform
	vec3 ExtendDir   = normalize( cross( vertex_new - eye_position, tangent_new ) ); //generated quad direction
	
	// 1 if the vertex is inside the small box, 0 otherwise
	outData.flag = length( vertex_new - small_box_center) >  bounding_radius ? 0 : 1; 
	
	// get up and down extension points for the generalized quad
	outData.up_coord      = projection * view * vec4( vertex_new + 0.1 * ExtendDir, 1.0 );
	outData.down_coord    = projection * view * vec4( vertex_new - 0.1 * ExtendDir, 1.0 );
	outData.center_point  = projection * view * vec4( vertex_new, 1.0 );	
	
	outData.line_id  = line_id;
	outData.click_id = click_id;

}











