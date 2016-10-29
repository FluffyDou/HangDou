/******************************************************************/
/* igluWidgetWindow.h                                             */
/* -----------------------                                        */
/*                                                                */
/* A class encapuslating an FLTK window that can contain other    */
/*     windows & widgets.                                         */
/*                                                                */
/*                                                                */
/* Chris Wyman (02/21/2012)                                       */
/******************************************************************/

#ifndef IGLU_WIDGET_WINDOW
#define IGLU_WIDGET_WINDOW

#include "igluWindow.h"
#include <FL/Fl_Double_Window.H>

class Fl_Scroll;
class Fl_Widget;

namespace iglu {

class IGLUInt;
class IGLUBool;
class IGLUFloat;

class IGLUWidgetWindow : public Fl_Double_Window {
public:
	IGLUWidgetWindow( int width, int height, const char *title=0 );
    virtual ~IGLUWidgetWindow();

	// Actually create the window
	void CreateWindow( int argc=0, char **argv=0 );

	// Add widgets
	void AddWidget( IGLUFloat *uiFloat, IGLUCallback *callback=0 );
	void AddWidget( IGLUInt *uiInt,     IGLUCallback *callback=0 );
	void AddWidget( IGLUBool *uiBool,   IGLUCallback *callback=0 );
	void AddWidgetSpacer( void )             { m_currentY += m_labelOffset; }

	// Add a button bound to an IGLUBool
	void AddButton( IGLUBool *uiBool,   IGLUCallback *callback=0 );

	// A pointer to a IGLUWidgetWindow could have type IGLUWidgetWindow::Ptr
	typedef IGLUWidgetWindow *Ptr;

protected:
	virtual int handle( int );   // This is FLTK's event handler callback

private:
	// Variables for basic window operation
	bool m_isRunning;
	bool m_noResize;
	unsigned int m_igluWindowFlags;

	// Offsets between sliders
	int m_labelOffset, m_betweenOffset;

	// Position in the scroller for the next widget
	unsigned int m_currentY;

	// Container for all widgets in the window
	Fl_Scroll *m_scroller;

	// Array of widgets...
	IGLUArray1D<Fl_Widget *> widgets;
};






} // end namespace iglu

#endif