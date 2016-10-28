#version 400

uniform float znear;
uniform float zfar;

in float edgeColorR;
in float line_id;
in float click_id;

out vec4 fragData_Color;
out vec4 fragData_Depth;

void main()
{
	// click_id == 0 means not clicked
	float alpha      = (click_id != 0.0) ? 0.0 : line_id;
	float edgeColorG = click_id != 0.0 ? edgeColorR : 0.0;

	// out put color to color attachment 0
	fragData_Color = vec4( edgeColorR, edgeColorG, 0.0, alpha );
	
	// out put depth to color attachment 1
	float realDepth = 2.0*zfar*znear / ( zfar+znear - (zfar-znear)*(2.0*gl_FragCoord.z - 1.0) );
	fragData_Depth = vec4( edgeColorR, edgeColorG, 0.0, realDepth);

}
