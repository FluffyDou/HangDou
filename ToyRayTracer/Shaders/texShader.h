#ifndef MYTEXTURESHADER_H
#define MYTEXTURESHADER_H

#include "iglu.h"

class texShader{

    iglu::IGLUShaderProgram::Ptr prog;
    GLint  attribute_verCoord;
    GLint  attribute_texCoord;
    GLint  uniform_texSampler;
    GLuint quad_verCoord;
    GLuint quad_texCoord;


public:

    texShader();
    texShader(char* vshader, char* fshader);
    ~texShader();

    inline void Enable()       { prog->Enable(); }
    inline void Disable()      { prog->Disable(); }

    void RenderTexture( GLuint texture );
    void testRenderTexture( GLuint texture );
    void InitAttribute( GLuint texture );
    void updateUniformTex( GLuint texture );
};


#endif