#version 400

//uniform mat4 transform;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;

in  vec3 vertex_coord;
in  vec3 normal;

out vec3 normalDirection;


void main(void)
{
        gl_Position     = projection * view * model * vec4(vertex_coord, 1.0);
        
//        normalDirection = normalize( vec3(model * vec4(normal, 1.0) ) ); 
        normalDirection = normalize(normal);

}











