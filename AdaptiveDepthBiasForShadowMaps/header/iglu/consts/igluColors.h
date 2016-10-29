/******************************************************************/
/* igluColors.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* Some common colors predefined to avoid vec4(x,y,z,w) calls for */
/*     commonly used colors.                                      */
/*                                                                */
/* Chris Wyman (02/09/2012)                                       */
/******************************************************************/

#ifndef IGLU_CONST_COLORS_H
#define IGLU_CONST_COLORS_H

#include "../vectors/vec4.h"

namespace iglu {

class color
{
public:
	// Pure colors
	static const vec4 Black, White, Red, Green, Blue, Purple, Yellow, Cyan, Orange;

	// Grayscale colors
	static const vec4 Gray00, Gray10, Gray20, Gray30, Gray40, Gray50;
	static const vec4 Gray60, Gray70, Gray80, Gray90, Gray100;

	// Blue colors
	static const vec4 MidnightBlue, NavyBlue, CornflowerBlue, SkyBlue, PowderBlue;
	static const vec4 Turquoise, RoyalBlue, SlateBlue, DodgerBlue;

	// Yellow colors
	static const vec4 Goldenrod, Gold, PaleGoldenrod, DarkGoldenrod;

	// Brown colors
	static const vec4 Wheat, SaddleBrown, Beige, Khaki, DarkKhaki, Tan, Sienna;

	// Green colors
	static const vec4 Aquamarine, SeaGreen, PaleGreen, LimeGreen, ForestGreen;

	// Pinks & purples
	static const vec4 HotPink, Pink, Maroon, Plum;
};


} // end namespace

#endif