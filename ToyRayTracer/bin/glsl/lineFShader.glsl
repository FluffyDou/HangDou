#version 400

uniform float znear;
uniform float zfar;

out vec4 fragData_Color;
out vec4 fragData_Depth;

void main()
{
       	// out put depth to color attachment 1
	float realDepth = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*gl_FragCoord.z - 1.0) );

        // leave the alpha to be zero, sicne the alpha value is used as line id in the cylinder shader
	fragData_Color = vec4(1.0, 0.0, 0.0, 0.0);
	fragData_Depth = vec4( 1.0, 0.0, 0.0, realDepth);
}
