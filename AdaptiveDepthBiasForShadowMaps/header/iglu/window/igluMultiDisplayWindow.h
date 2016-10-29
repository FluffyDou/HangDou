/*******************************************************************/
/* igluMultiDisplayWindow.h                                        */
/* -----------------------                                         */
/*                                                                 */
/* Often, when prototyping OpenGL code, you really want to compare */
/*   multiple rendering techniques, or switch to partially-done or */
/*   debug display modes.  This can be done either by cluttering   */
/*   up your display and computation logic in a single window, or  */
/*   by having seperate display "modes" for each type of display.  */
/* This class is designed to support multiple display modes, to    */
/*   allow for the second programming style.  Ideally, these modes */
/*   will all be (relatively) simple and easy to read, improving   */
/*   code reuse, legibility, and portability (and reducing costs   */
/*   of debugging complex logic).                                  */
/*                                                                 */
/*                                                                 */
/* Chris Wyman (04/05/2012)                                        */
/*******************************************************************/

#ifndef IGLU_MULTI_DISPLAY_WINDOW
#define IGLU_MULTI_DISPLAY_WINDOW

// This include needed to ensure we have consistent window #defines and enums.
#include "igluWindow.h"

namespace iglu {

class IGLUFrameRate;
class IGLUMultiDisplayWidgetWindow;
class IGLUDisplayMode;
class IGLUInt;

class IGLUMultiDisplayWindow :
#if defined(WIN32)
			public Fl_BetterGl_Window {
#else
			public Fl_Gl_Window {
#endif

public:
	IGLUMultiDisplayWindow( int width, int height, const char *title="",
		                    int widgetW=600, int widgetH=300, const char *widgetTitle=0 );
    virtual ~IGLUMultiDisplayWindow();

	// Actually create the window
	void CreateWindow( int argc=0, char **argv=0 );

	// If any preprocessing needs to occur after glutCreateWindow() but prior
	//    to glutMainLoop(), this can go inside a user-specified preprocess function.
	//    The callback can either take no parameters *or* a pointer to this window object.
	void SetPreprocessOnGLInit( void (*preprocess)( void ) );
	void SetPreprocessOnGLInit( void (*preprocess)( void * ) );

	// Setup any window flags
	void SetWindowProperties( unsigned int flags );

	// Adds a display mode
	void AddDisplayMode( IGLUDisplayMode *mode );

	// Add widgets into the UI's widget window for the correct display mode.  
	//    Optionally, add a callback to be called whenever the variable is changed. 
	void AddWidget( int displayMode, IGLUFloat *uiFloat, IGLUCallback *callback=0 );
	void AddWidget( int displayMode, IGLUInt   *uiInt,   IGLUCallback *callback=0 );
	void AddWidget( int displayMode, IGLUBool  *uiBool,  IGLUCallback *callback=0 );
	void AddButton( int displayMode, IGLUBool  *uiBool,  IGLUCallback *callback=0 );

	// Add widgets to the common pane of the UI window (seen for all display modes)
	void AddCommonWidgetSpacer( void );
	void AddCommonWidget( IGLUBool  *uiBool,  IGLUCallback *callback=0 );
	void AddCommonWidget( IGLUInt   *uiInt,   IGLUCallback *callback=0 );
	void AddCommonWidget( IGLUFloat *uiFloat, IGLUCallback *callback=0 );
	void AddCommonButton( IGLUBool  *uiBool,  IGLUCallback *callback=0 );

	// Add widgets into the UI's widget window for the correct display mode.  
	//    Optionally, add a callback to be called whenever the variable is changed. 
	void AddWidget( int displayMode, IGLUFloat *uiFloat, void (*changeCallback)( void * ) );
	void AddWidget( int displayMode, IGLUInt *uiInt,     void (*changeCallback)( void * ) );
	void AddWidget( int displayMode, IGLUBool *uiBool,   void (*changeCallback)( void * ) );
	void AddWidgetSpacer( int displayMode ); 

	// IGLUVariables can have a more complex callback, which takes as input some specific extra value
	//    (i.e., a 'this' pointer for some class) in addition to the IGLUVariable pointer.  These callbacks
	//    can be set up using these methods.
	void AddWidget( int displayMode, IGLUFloat *uiFloat, void (*changeCallback)( void *, void * ), void *callbackData );
	void AddWidget( int displayMode, IGLUInt *uiInt,     void (*changeCallback)( void *, void * ), void *callbackData );
	void AddWidget( int displayMode, IGLUBool *uiBool,   void (*changeCallback)( void *, void * ), void *callbackData );

	// Allow the user to setup GLUT callback functions of their own.
	void SetIdleCallback   ( void (*newIdle)( void ) );
	void SetKeyboardCallback( void (*newKeyboard)( unsigned char, int, int ) );
	void SetReshapeCallback( void (*newReshape)( int, int ) );
	void SetMouseButtonCallback( void (*newButton)( int, int, int, int ) );
	void SetActiveMotionCallback( void (*newActive)( int, int ) );
	void SetPassiveMotionCallback( void (*newPassive)( int, int ) );
	void SetEntryCallback( void (*newEntry)( int ) );      // Doesn't currently work quite as expected.
	void SetVisibilityCallback( void (*newVisible)( int ) );
	void SetSpecialKeyCallback( void (*newSpecial) ( int, int, int ) );

	// Sometimes the on-screen UI is annoying (e.g., screen/video captures) and you may want it disabled
	inline void EnableWindowUI( void )                  { m_displayWindowUI = true; }
	inline void DisableWindowUI( void )                 { m_displayWindowUI = false; }
	inline bool IsUIOnscreen( void ) const              { return m_displayWindowUI; }

	// If you need the widget window pointer, you can grab it here.  Be careful how you
	//    use this pointer.  The IGLUMultiDisplayWindow mostly manages basic operations,
	//    so if you want to do something clever, make doubly sure you won't break IGLU
	//    or add the funcitonality into this window class!
	inline IGLUMultiDisplayWidgetWindow *GetWidgetWindow( void ) const { return m_widgetWindow; }

	// After all initializing all windows, call Run().  This is *not* window-specific.
	//    This is equivalent to calling glutMainLoop(), and starts processing inputs
	//    for *all* windows!
	static void Run( void );

	// This is a trivial function that can be used as a no-op idle function.
	//    Usually this is used when IGLU_WINDOW_REDRAW_ON_IDLE is a window property.
	static void NullIdle( void );

	// A pointer to a IGLUWindow could have type IGLUWindow::Ptr
	typedef IGLUMultiDisplayWindow *Ptr;
protected:
	virtual void draw( void );   // This is FLTK's redraw callback
	virtual int handle( int );   // This is FLTK's event handler callback

private:
	// Variables for basic window operation
	bool m_isRunning, m_isInitialized;
	bool m_redrawOnIdle, m_fltkIdleAdded;
	bool m_noIdle, m_noResize, m_displayFPS;
	int m_currentWidth, m_currentHeight;
	bool m_displayWindowUI;

	// Information about the display modes available in this window
	IGLUInt *m_currentDisplayMode;
	IGLUArray1D< IGLUDisplayMode * > m_displayModes;

	// Variables used to setup the GLUT window (not useful later)
	int m_initialPosX, m_initialPosY;
	unsigned int m_igluWindowFlags, m_fltkMode;

	bool m_overPlusButton;

	// Access to screen widgets
	int m_currentButtonDown;

	IGLUFrameRate *frameRate;

	// Do we control a widget window ?
	IGLUMultiDisplayWidgetWindow *m_widgetWindow;

	// Called the first time we see we're ready to draw...
	void Initialize( void );

	// To call an idle function, we'll have a separate internal class idle
	//    that calls the function specified by SetIdleCallback().  This will
	//    get called by DefaultIdleFunction(), which takes a 'this' pointer
	//    to the particular class instantiation.
	void CallIdleFunction( void );
	static void DefaultIdleFunction   ( void *igluWinPtr );

	// Variables storing user-specified window preprocess
	void (*m_preprocessOnGLInit)       ( void );
	void (*m_preprocessOnGLInitWithPtr)( void * );

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

	// This callback is executed when a window's display mode counter changes
	void Callback_UpdateDisplayMode( IGLUVariable *igluInt );
};


} // end namespace iglu

#endif