/******************************************************************/
/* igluFLTKGLWindow.h                                             */
/* -----------------------                                        */
/*                                                                */
/* A Fl_Gl_Window has one big drawback...  You can't create a     */
/*    specific OpenGL context type with it.  This class is        */
/*    is *identical* to a FLTK 1.3.2 Fl_Gl_Window in interfaces   */
/*    and code, except in how the OpenGL context is created.      */
/*                                                                */
/*                                                                */
/* Chris Wyman (02/21/2012)                                       */
/******************************************************************/

#ifndef Fl_BetterGl_Window_H
#define Fl_BetterGl_Window_H

#include <FL/Fl_Gl_Window.H>

class Fl_BetterGl_Window : public Fl_Gl_Window {
public:
	// This overrides the Fl_Gl_Window's OpenGL context creation to
	//   give us more flexibility
	void make_current();

	~Fl_BetterGl_Window() {}
	Fl_BetterGl_Window(int W, int H, const char *l=0) : Fl_Gl_Window(W,H,l) {}
	Fl_BetterGl_Window(int X, int Y, int W, int H, const char *l=0) : Fl_Gl_Window(X,Y,W,H,l) {}

protected:
};

#endif

