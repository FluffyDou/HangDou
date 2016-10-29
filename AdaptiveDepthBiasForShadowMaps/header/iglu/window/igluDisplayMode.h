/*******************************************************************/
/* igluDisplayMode.h                                               */
/* -----------------------                                         */
/*                                                                 */
/* A base class for deriving various "display modes" for windows   */
/*     created with IGLUMultiDisplayWindow.  The idea of a display */
/*     mode is that it encapsulates the drawing of a fairly simple */
/*     and straightforward screen, plus any special initialization */
/*     needed for display mode-specific variables.                 */
/* Additional bits of data can/will/are eventually be encapsulated */
/*     here, including adding widgets to the display mode specific */
/*     UI, specific mouse/keyboard callbacks for this display,     */
/*     changes to the window title, or other user-displayable help */
/*     text.                                                       */
/*                                                                 */
/*                                                                 */
/* Chris Wyman (04/05/2012)                                        */
/*******************************************************************/

#ifndef IGLU_DISPLAY_MODE
#define IGLU_DISPLAY_MODE

namespace iglu {

class IGLUDisplayMode {
public:
	IGLUDisplayMode() {}
	virtual ~IGLUDisplayMode() {}

	// Initialization specific to this display mode (occurs after GL context created)
	virtual void Initialize( int trackNum ) = 0;

	// The display mode's callback.  This controls redrawing the GL viewport 
	virtual void Display() = 0;

	// A method that is called every time this mode becomes active (e.g., a switch
	//   from one display mode to another)
	virtual void UponSwitch() {};

	// Returns a name to be displayed in the window's title bar
	virtual const char *GetWindowTitle( void ) 
		{ return "*** Unknown Display ***"; }

	// Returns a string to display in a help window / about box.
	virtual const char *GetHelpDescriptor( void )
		{ return "(No display description provided!)"; }

	// A pointer to a IGLUDisplayMode could have type IGLUDisplayMode::Ptr
	typedef IGLUDisplayMode *Ptr;

private:
};




} // end namespace iglu

#endif