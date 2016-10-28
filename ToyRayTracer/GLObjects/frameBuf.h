#ifndef MYFRAMEBUFFER
#define MYFRAMEBUFFER

#include "iglu.h"

class frameBuf{

     GLuint rto;
     GLuint fbo;

public:
    frameBuf();
    frameBuf(int screenWidth, int screenHeight);
    ~frameBuf();

    void Enalble()   { glBindFramebuffer(GL_FRAMEBUFFER, fbo); }
    void Disable()   { glBindFramebuffer(GL_FRAMEBUFFER, 0); }

};


#endif