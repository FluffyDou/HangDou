#version 420

uniform usampler2D depthTex;
uniform float      intenseBias;
uniform int        layer;

in vec2 fragTexCoord;

layout(location = 0) out vec4 result;

void main( void ) 
{
    uvec4 texColor = texture( depthTex, fragTexCoord.xy );
    
    // unpack the duel depth
    uint packedDepth = texColor.w;        
    vec2 duelDepth   = unpackUnorm2x16(packedDepth);
    
    //float renderColor = layer==0 ? duelDepth.x : duelDepth.y;
    //float renderColor = layer==0 ? float(texColor.x) : float(texColor.y);
	float renderColor = float(bitCount(texColor.x)+bitCount(texColor.y)+bitCount(texColor.z)) / intenseBias;
    renderColor = layer==1 ? duelDepth.x : renderColor;
    renderColor = layer==2 ? duelDepth.y : renderColor;
    
    result = vec4(renderColor, renderColor, renderColor, 1.0);
    //result = texColor.x != 0 ? vec4(1.0, 0.0, 0.0, 1.0) : vec4(0.0, 1.0, 0.0, 1.0);
    //result = vec4( 2*texColor.xyz - 1, 1.0 ); // map the intense from [0.9, 1] to [0, 1]    
	//result = vec4( 10*texColor.xyz - 9, 1.0 ); // map the intense from [0.9, 1] to [0, 1]
}