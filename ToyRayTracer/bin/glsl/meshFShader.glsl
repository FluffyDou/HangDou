#version 400

uniform float znear;
uniform float zfar;

in  vec3 normalDirection;

out vec4 fragData_Color;
out vec4 fragData_Depth;

struct lightSource
{
        vec4 position;
        vec4 diffuse;
};

lightSource light0 = lightSource(
//        vec4(30.0,  42.5,  0.0, 1.0),
        vec4(-1.0,  1.0,  -1.0, 1.0),
        vec4( 1.0,  1.0,   1.0, 1.0)
);

lightSource light1 = lightSource(
//        vec4(30.0,  42.5,  0.0, 1.0),
        vec4(1.0,  -1.0,  1.0, 1.0),
        vec4(1.0,   1.0,  1.0, 1.0)
);

vec4 scene_ambient = vec4(0.4, 0.4, 0.4, 1.0);
 
struct material
{
        vec4 ambient;
        vec4 diffuse;
};

material frontMaterial = material(
        vec4(0.6, 0.6, 0.6, 1.0),
        vec4(0.6, 0.6, 0.6, 1.0)
);
 

void main()
{
        // direction light
        vec3 light0Direction = normalize( vec3(light0.position) );
        vec3 light1Direction = normalize( vec3(light1.position) );
        
        vec3 ambientLighting = vec3(scene_ambient) * vec3(frontMaterial.ambient);

        vec3 diffuseReflectionLight0 = vec3(light0.diffuse) * vec3(frontMaterial.diffuse)
                                 * max( 0.0, dot(normalDirection, light0Direction) );

        vec3 diffuseReflectionLight1 = vec3(light1.diffuse) * vec3(frontMaterial.diffuse)
                                 * max( 0.0, dot(normalDirection, light1Direction) );

        fragData_Color = vec4( (ambientLighting + diffuseReflectionLight0 + diffuseReflectionLight1).x, 0.0, 0.0, 0.0);
        
       	// out put depth to color attachment 1
	float realDepth = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*gl_FragCoord.z - 1.0) );
	fragData_Depth = vec4( 0.0, 0.0, 0.0, realDepth);

}
