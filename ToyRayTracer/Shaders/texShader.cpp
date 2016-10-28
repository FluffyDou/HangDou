#include "Shaders/texShader.h"

inline texShader::texShader()
{
    prog               = 0;
    quad_verCoord      = 0;
    quad_texCoord      = 0;
    attribute_verCoord = 0;
    attribute_texCoord = 0;
}

texShader::texShader(char* vshader, char* fshader)
{
    //////// create VBO of the quad where we bind the render texture

    // quad's coordinate
    GLfloat verCoord[] = { -1.0, -1.0,  1.0,-1.0,   1.0,1.0,  -1.0,1.0 };
    glGenBuffers(1, &quad_verCoord);
    glBindBuffer(GL_ARRAY_BUFFER, quad_verCoord);  	
    glBufferData(GL_ARRAY_BUFFER, sizeof(verCoord), verCoord, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // texture coordinate
    GLfloat texCoord[] = { 0.0,0.0,  1.0,0.0,  1.0,1.0,  0.0,1.0 };
    glGenBuffers(1, &quad_texCoord);
    glBindBuffer(GL_ARRAY_BUFFER, quad_texCoord);  	
    glBufferData(GL_ARRAY_BUFFER, sizeof(texCoord), texCoord, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);


    ////// initialize the shader and locate attribute variables //////

    prog = new iglu::IGLUShaderProgram( vshader, fshader );
    attribute_verCoord = 0;
    attribute_texCoord = 1;
    //texShaderProg->SetProgramEnables( iglu::IGLU_GLSL_DEPTH_TEST ); 

    /*attribute_verCoord = glGetAttribLocation(prog->GetProgramID(), "vertex_coord" );
    //attribute_texCoord = glGetAttribLocation(prog->GetProgramID(),    "tex_coord" );

    //std::cout << "a_ver is "<< attribute_verCoord <<", a_tex is "<<attribute_texCoord<<"\n";
    //
    //// check if the shader is valid
    //if(! prog->IsValid() || attribute_verCoord < 0 || attribute_texCoord < 0)
    //{
    //    printf("The texture shader is not valid..\n");
    //}*/

}


void texShader::InitAttribute(GLuint texture)
{
    // set attribute array pointer
    glBindBuffer( GL_ARRAY_BUFFER, quad_verCoord );
    glVertexAttribPointer(
        attribute_verCoord,
        2,
        GL_FLOAT,
        GL_FALSE,
        0,
        0
        );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );

    glBindBuffer( GL_ARRAY_BUFFER, quad_texCoord );
    glVertexAttribPointer(
        attribute_texCoord,
        2,
        GL_FLOAT,
        GL_FALSE,
        0,
        0
        );

    glBindBuffer(GL_ARRAY_BUFFER, 0);

}


void texShader::updateUniformTex( GLuint texture )
{
    // set uniform location
    uniform_texSampler = glGetUniformLocation( prog->GetProgramID(), "renderTex" );
    glBindTexture( GL_TEXTURE_2D, texture);
    glUniform1i( uniform_texSampler, 0 );
    glBindTexture( GL_TEXTURE_2D, 0);
}

void texShader::RenderTexture( GLuint texture )
{
    
    glActiveTexture( GL_TEXTURE0 );
    glBindTexture( GL_TEXTURE_2D, texture);

    glEnableVertexAttribArray( attribute_verCoord );
    glEnableVertexAttribArray( attribute_texCoord );

    glDrawArrays(GL_QUADS, 0, 4);
    
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisableVertexAttribArray( attribute_verCoord );
    glDisableVertexAttribArray( attribute_texCoord );
}


//void texShader::testRenderTexture( GLuint texture )
//{
//    GLint loc1 = glGetUniformLocation( prog->GetProgramID(), "renderTex" );
//    glActiveTexture( GL_TEXTURE0 );
//    glBindTexture( GL_TEXTURE_2D, texture);
//    glUniform1i( loc1, 0 );
//
//    glBindBuffer( GL_ARRAY_BUFFER, quad_verCoord );
//    glEnableVertexAttribArray( 0 );
//    glVertexAttribPointer(
//        0,
//        2,
//        GL_FLOAT,
//        GL_FALSE,
//        0,
//        0
//    );
//    glBindBuffer( GL_ARRAY_BUFFER, 0 );
//
//    glBindBuffer( GL_ARRAY_BUFFER, quad_texCoord );
//    glEnableVertexAttribArray( 1 );
//    glVertexAttribPointer(
//        1,
//        2,
//        GL_FLOAT,
//        GL_FALSE,
//        0,
//        0
//    );
//
//    glDrawArrays(GL_QUADS, 0, 4);
//
//    glBindTexture(GL_TEXTURE_2D, 0);
//    glBindBuffer(GL_ARRAY_BUFFER, 0);
//    glDisableVertexAttribArray(quad_verCoord);
//    glDisableVertexAttribArray(quad_texCoord);
//}