/*******************************************************************/
/* igluMultiDisplayWidgetWindow.h                                  */
/* ------------------------------                                  */
/*                                                                 */
/* NOTE: This is not really designed as a stand-alone class, but   */
/*    is rather tightly integrated with the IGLUMultiDisplayWindow */
/*    class.  Feel free to use it, but it may change rather in     */
/*    unexpected and strange ways as IGLU evolves to help improve  */
/*    the ease of using the IGLUMultiDisplayWindow.  Some changes  */
/*    may be unintuitive unless you understand IGLU internals.     */
/*                                                                 */
/* Basically, this class is designed as a two-pane widget window,  */
/*    where one uses built-in IGLUMultiDisplayWindow variables to  */
/*    control which display mode is being used, while the other    */
/*    pane contains display mode-specific IGLUVariable widgets.    */
/*    The second pane is tabbed (or otherwise displayed) so that   */
/*    essentially each display mode can think about having its own */
/*    IGLUWidgetWindow.                                            */
/*                                                                 */
/*                                                                 */
/* Chris Wyman (02/21/2012)                                        */
/*******************************************************************/

#ifndef IGLU_MULTI_DISPLAY_WIDGET_WINDOW
#define IGLU_MULTI_DISPLAY_WIDGET_WINDOW

#include "igluMultiDisplayWindow.h"

class Fl_Scroll;
class Fl_Widget;
class Fl_Hor_Slider;
class Fl_Box;
class Fl_Text_Display;

namespace iglu {

class IGLUInt;
class IGLUBool;
class IGLUFloat;

class IGLUMultiDisplayWidgetWindow : public Fl_Double_Window {
public:
	IGLUMultiDisplayWidgetWindow( int width, int height, const char *title=0 );
    virtual ~IGLUMultiDisplayWidgetWindow();

	// Actually create the window
	void CreateWindow( int argc=0, char **argv=0 );

	// A new right pane should be created for each display mode
	void CreateNewRightPane( void );
	void SetCurrentRightPane( int paneNum );

	// Add widgets (for individual display window panes)
	void AddWidget( uint paneID, IGLUFloat *uiFloat );
	void AddWidget( uint paneID, IGLUInt *uiInt );
	void AddWidget( uint paneID, IGLUBool *uiBool );
	void AddWidgetSpacer( uint paneID )                { m_currentY[paneID] += m_spacerOffset; }

	// Add a button bound to an IGLUBool
	void AddButton( uint paneID, IGLUBool *uiBool );

	// Add common widgets and buttons (displayed for all panes)
	void AddCommonWidget( IGLUFloat *uiFloat );
	void AddCommonWidget( IGLUInt *uiInt );
	void AddCommonWidget( IGLUBool *uiBool );
	void AddCommonWidgetSpacer( void )                 { m_currentLeftY += m_spacerOffset; }
	void AddCommonButton( IGLUBool *uiBool );

	// Setup the left-pane widgets that control the display mode interface
	void SetDisplayModeText( const char *modeName );
	void SetDisplayCounter( IGLUInt *uiInt );
	void SetDisplayDescriptionText( const char *modeDscr );

	// A pointer to a IGLUWidgetWindow could have type IGLUWidgetWindow::Ptr
	typedef IGLUWidgetWindow *Ptr;

protected:
	virtual int handle( int );   // This is FLTK's event handler callback

private:
	// Setup the left pane widgets & display
	void SetupLeftPane( void );

	// Variables for basic window operation
	bool m_isRunning;
	bool m_noResize;
	unsigned int m_igluWindowFlags;

	// Offsets between sliders
	int m_labelOffset, m_betweenOffset, m_spacerOffset;

	// Position in the right panes for the next widget
	IGLUArray1D< uint > m_currentY;

	// Position in the left pane for the next common widget
	uint m_currentLeftY;

	// Which right-hand side pane are we showing?
	int m_curRightPane;

	// Container for all widgets in the window
	IGLUArray1D<Fl_Scroll *> m_rightPane;
	Fl_Group  *m_leftPane;

	// Widgets needed when changing the left-pane display
	Fl_Box          *m_displayModeName;
	Fl_Text_Display *m_displayModeDescriptor;
	Fl_Hor_Slider   *m_displayModeSelector;
	IGLUInt         *m_displayCounter;

	// Array of widgets...
	IGLUArray1D<Fl_Widget *> widgets;

	// These get called when m_displayModeSelector changes.  First the static
	//    member (the actual callback).  This is a wrapper around UpdateSelectedMode()
	static void SelectedModeCallback( Fl_Widget *, void * );
	void UpdateSelectedMode( Fl_Widget * );
};






} // end namespace iglu

#endif