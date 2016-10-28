#include "GLObjects/frameBuf.h"

inline frameBuf::frameBuf()
{
    rto = 0;
    fbo = 0;
}

frameBuf::frameBuf(int screenWidth, int screenHight)
{
    // generate texture object
    glGenTextures(1, &rto);
    glBindTexture(GL_TEXTURE_2D, rto);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, screenWidth, screenHight, 0, GL_RGBA, GL_FLOAT, 0); 
    glBindTexture(GL_TEXTURE_2D, 0);
    
    // generate frame buffer object
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    // bind the texture to fbo
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rto, 0);
    
    // check if the frame is valid
    GLenum status = glCheckFramebufferStatusEXT( GL_FRAMEBUFFER );
    assert(glCheckFramebufferStatus( GL_FRAMEBUFFER ) == GL_FRAMEBUFFER_COMPLETE );
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}