/******************************************************************/
/* igluMainWindow.h                                               */
/* -----------------------                                        */
/*                                                                */
/* A class encapuslating the main GLUT window for a simple OpenGL */
/*    program.                                                    */
/*                                                                */
/* Standard GLUT callbacks are available to the the programmer    */
/*    using the Set*Callback() methods.                           */
/*                                                                */
/* Chris Wyman (12/06/2009)                                       */
/******************************************************************/

#ifndef IGLU_MAIN_WINDOW
#define IGLU_MAIN_WINDOW


class IGLUMainWindow {
public:
	
	// Called at the beginning of main(), with command line parameters as inputs
	static void Initialize( int *argc, char **argv );

	// If any preprocessing needs to occur after glutCreateWindow() but prior
	//    to glutMainLoop(), this can go inside a user-specified preprocess function.
	static void SetPreprocessOnGLInit( void (*preprocess)( void ) );

	// Setup any window flags
	static void SetWindowProperties( unsigned int flags );

	// After all initializing approproate callbacks, open a GL window and run!
	//    This function calls the user-preprocess function, if needed.
	static void Run( const char *windowTitle, int windowWidth, int windowHeight );

	// Allow the user to setup GLUT callback functions of their own.
	static void SetDisplayCallback( void (*newDisplay)( void ) );
	static void SetReshapeCallback( void (*newReshape)( int, int ) );
	static void SetKeyboardCallback( void (*newKeyboard)( unsigned char, int, int ) );
	static void SetSpecialKeyCallback( void (*newSpecial) ( int, int, int ) );
	static void SetMouseButtonCallback( void (*newButton)( int, int, int, int ) );
	static void SetActiveMotionCallback( void (*newActive)( int, int ) );
	static void SetPassiveMotionCallback( void (*newPassive)( int, int ) );
	static void SetEntryCallback( void (*newEntry)( int ) );
	static void SetVisibilityCallback( void (*newVisible)( int ) );
	static void SetIdleCallback( void (*newIdle)( void ) );

private:
	// Variables for basic window operation
	static bool m_isRunning, m_isInitialized, m_isWindowInitialized;
	static int m_currentWidth, m_currentHeight;

	// Variables used to setup the GLUT window (not useful later)
	static int m_initialPosX, m_initialPosY;
	static unsigned int m_igluWindowFlags;

	// Access to screen widgets
	static int m_currentButtonDown;

	// This window has a number of callbacks that override the user's callbacks,
	//    but these then call the specified user functions.  This allows for
	static void GlobalDisplayCallback( void );
	static void GlobalReshapeCallback( int width, int height );
	static void GlobalKeyboardCallback( unsigned char key, int x, int y );
	static void GlobalSpecialKeyCallback( int key, int x, int y );
	static void GlobalMouseButtonCallback( int button, int state, int x, int y );
	static void GlobalActiveMotionCallback( int x, int y );
	static void GlobalPassiveMotionCallback( int x, int y );
	static void GlobalEntryCallback( int state );      // Called with window focus changes.
	static void GlobalVisibilityCallback( int state );  // Called upon minimization/maximization
	static void GlobalIdleFunction ( void );

	// Variables storing user-specified window preprocess, which occurs
	//    after glutCreateWindow() but before glutMainLoop()
	static void (*m_preprocessOnGLInit)( void );

	// Variables storing externally-specifier GLUT callback functions
	static void (*m_callbackDisplay)( void );
	static void (*m_callbackReshape)( int, int );
	static void (*m_callbackKeyboard)( unsigned char, int, int );
	static void (*m_callbackSpecial)( int, int, int );
	static void (*m_callbackButton)( int, int, int, int );
	static void (*m_callbackActiveMove)( int, int );
	static void (*m_callbackPassiveMove)( int, int );
	static void (*m_callbackEntry)( int );
	static void (*m_callbackVisible)( int );
	static void (*m_callbackIdle)( void );
};



#define IGLU_WINDOW_NO_FLAGS			0x00000000
#define IGLU_WINDOW_NO_RESIZE			0x00000001
#define IGLU_WINDOW_NO_IDLE             0x00000002
#define IGLU_WINDOW_REDRAW_ON_IDLE      0x00000004
#define IGLU_WINDOW_DOUBLE		        0x00000008
#define IGLU_WINDOW_DEPTH		        0x00000010
#define IGLU_WINDOW_STENCIL     		0x00000020
#define IGLU_WINDOW_BLEND               0x00000040
#define IGLU_WINDOW_MULTISAMPLE         0x00000080

#define IGLU_WINDOW_W_FRAMERATE         0x10000000


#endif