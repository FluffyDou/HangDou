/******************************************************************/
/* igluWindow.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* A class encapuslating the a simple FLTK OpenGL window.         */
/*                                                                */
/*                                                                */
/* Chris Wyman (02/21/2012)                                       */
/******************************************************************/

#ifndef IGLU_WINDOW
#define IGLU_WINDOW

#include "iglu/window/igluFLTKGLWindow.h"

#ifdef CreateWindow
#undef CreateWindow
#endif

namespace iglu {

class IGLUFrameRate;
class IGLUWidgetWindow;

class IGLUWindow : 
#if defined(WIN32)
			public Fl_BetterGl_Window {
#else
			public Fl_Gl_Window {
#endif

public:
	IGLUWindow( int width, int height, const char *title=0 );
    virtual ~IGLUWindow();

	// Actually create the window
	void CreateWindow( int argc=0, char **argv=0 );

	// If the user wants to tell us about a widget window, we can explicitly control showing it
	void SetWidgetWindow( IGLUWidgetWindow *widgets );

	// If any preprocessing needs to occur after glutCreateWindow() but prior
	//    to glutMainLoop(), this can go inside a user-specified preprocess function.
	void SetPreprocessOnGLInit( void (*preprocess)( void ) );

	// Setup any window flags
	void SetWindowProperties( unsigned int flags );

	// Toggle on/off the automatic framerate display
	inline void DisplayFramerate( bool doDisplay )				{ m_displayFPS = doDisplay; }

	// Allow the user to setup GLUT callback functions of their own.
	void SetDisplayCallback( void (*newDisplay)( void ) );
	void SetIdleCallback   ( void (*newIdle)( void ) );
	void SetKeyboardCallback( void (*newKeyboard)( unsigned char, int, int ) );
	void SetReshapeCallback( void (*newReshape)( int, int ) );
	void SetMouseButtonCallback( void (*newButton)( int, int, int, int ) );
	void SetActiveMotionCallback( void (*newActive)( int, int ) );
	void SetPassiveMotionCallback( void (*newPassive)( int, int ) );
	void SetEntryCallback( void (*newEntry)( int ) );      // Doesn't currently work quite as expected.
	void SetVisibilityCallback( void (*newVisible)( int ) );
	void SetSpecialKeyCallback( void (*newSpecial) ( int, int, int ) );

	// After all initializing all windows, call Run().  This is *not* window-specific.
	//    This is equivalent to calling glutMainLoop(), and starts processing inputs
	//    for *all* windows!
	static void Run( void );

	// This is a trivial function that can be used as a no-op idle function.
	//    Usually this is used when IGLU_WINDOW_REDRAW_ON_IDLE is a window property.
	static void NullIdle( void );

	// A pointer to a IGLUWindow could have type IGLUWindow::Ptr
	typedef IGLUWindow *Ptr;
protected:
	virtual void draw( void );   // This is FLTK's redraw callback
	virtual int handle( int );   // This is FLTK's event handler callback

private:
	// Variables for basic window operation
	bool m_isRunning, m_isInitialized;
	bool m_redrawOnIdle, m_fltkIdleAdded;
	bool m_noIdle, m_noResize, m_displayFPS;
	int m_currentWidth, m_currentHeight;

	// Variables used to setup the GLUT window (not useful later)
	int m_initialPosX, m_initialPosY;
	unsigned int m_igluWindowFlags, m_fltkMode;

	bool m_overPlusButton;

	// Access to screen widgets
	int m_currentButtonDown;

	IGLUFrameRate *frameRate;

	// Do we control a widget window ?
	IGLUWidgetWindow *m_widgetWindow;

	// Called the first time we see we're ready to draw...
	void Initialize( void );

	// To call an idle function, we'll have a separate internal class idle
	//    that calls the function specified by SetIdleCallback().  This will
	//    get called by DefaultIdleFunction(), which takes a 'this' pointer
	//    to the particular class instantiation.
	void CallIdleFunction( void );
	static void DefaultIdleFunction   ( void *igluWinPtr );

	// Variables storing user-specified window preprocess
	void (*m_preprocessOnGLInit)( void );

	// Variables storing externally-specified callback functions
	void (*m_callbackDisplay)( void );
	void (*m_callbackIdle)( void );
	void (*m_callbackKeyboard)( unsigned char, int, int );
	void (*m_callbackReshape)( int, int );
	void (*m_callbackButton)( int, int, int, int );
	void (*m_callbackActiveMove)( int, int );
	void (*m_callbackPassiveMove)( int, int );
	void (*m_callbackEntry)( int );
	void (*m_callbackVisible)( int );
	void (*m_callbackSpecial)( int, int, int );

	// We're drawing a few custom buttons.  Do some event handling in here.
	void CheckButtonCoverage( int, int );
};


#define IGLU_EVENT_DOWN		      GLUT_DOWN
#define IGLU_EVENT_UP             GLUT_UP
#define IGLU_EVENT_LEFT_BUTTON    GLUT_LEFT_BUTTON
#define IGLU_EVENT_RIGHT_BUTTON   GLUT_RIGHT_BUTTON
#define IGLU_EVENT_MIDDLE_BUTTON  GLUT_MIDDLE_BUTTON
#define IGLU_EVENT_ENTERED        GLUT_ENTERED
#define IGLU_EVENT_LEFT           GLUT_LEFT
#define IGLU_EVENT_VISIBLE        GLUT_VISIBLE
#define IGLU_EVENT_NOT_VISIBLE    GLUT_NOT_VISIBLE

// NOTE:  These are FLTK keys, NOT GLUT keys -- their values are different!
#define IGLU_KEY_F1			1
#define IGLU_KEY_F2			2
#define IGLU_KEY_F3			3
#define IGLU_KEY_F4			4
#define IGLU_KEY_F5			5
#define IGLU_KEY_F6			6
#define IGLU_KEY_F7			7
#define IGLU_KEY_F8			8
#define IGLU_KEY_F9			9
#define IGLU_KEY_F10		10
#define IGLU_KEY_F11		11
#define IGLU_KEY_F12		12
#define IGLU_KEY_LEFT		FL_Left
#define IGLU_KEY_UP			FL_Up
#define IGLU_KEY_RIGHT		FL_Right
#define IGLU_KEY_DOWN		FL_Down
#define IGLU_KEY_PAGE_UP	FL_Page_Up
#define IGLU_KEY_PAGE_DOWN	FL_Page_Down
#define IGLU_KEY_HOME		FL_Home
#define IGLU_KEY_END		FL_End
#define IGLU_KEY_INSERT		FL_Insert




} // end namespace iglu

#endif