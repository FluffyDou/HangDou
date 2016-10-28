
#ifndef RGBCOLOR_H
#define RGBCOLOR_H

#include "DataTypes/vec3.h"
// A Color is a RGB-triple of floats. 

class Color 
{
    float d[3];
public:
	// Constructors & destructors
	inline Color();                          // Default constructor
    inline Color(float r, float g, float b); // color from three floats
	inline Color(const Color& copy);         // Copy constructor 
    inline Color(const vec3 color);          // initialize color from vector 3

	// Mathematical operations
    inline Color operator*(const Color& c) const  { return Color(d[0]*c.d[0], d[1]*c.d[1], d[2]*c.d[2]); }
	inline Color operator*(float s) const            { return Color(d[0]*s, d[1]*s, d[2]*s); }
	inline Color& operator*=(float s) 	            { d[0]*=s; d[1]*=s; d[2]*=s; return *this; }
    inline Color operator+(const Color& c) const  { return Color(d[0]+c.d[0], d[1]+c.d[1], d[2]+c.d[2]); }
	inline Color& operator+=(const Color& c) 	    { d[0]+=c.d[0]; d[1]+=c.d[1]; d[2]+=c.d[2]; return *this; }
    inline Color operator-(const Color& c) const 	{ return Color(d[0]-c.d[0], d[1]-c.d[1], d[2]-c.d[2]); }
	inline Color& operator-=(const Color& c) 	    { d[0]-=c.d[0]; d[1]-=c.d[1]; d[2]-=c.d[2]; return *this; }

	// Accessor functions
	inline float Red() const                            { return d[0]; }
    inline float Green() const                          { return d[1]; }
    inline float Blue() const                           { return d[2]; }

    inline vec3 ColorVec() const                        { return vec3(d[0],d[1],d[2]); }

	// Returns a grayscale ('luminance') value roughly corresponding to the color
    inline float Luminance() const { return (float)(0.213*d[0] + 0.175f*d[1] +  0.072*d[2]); }

	// Returns the maximum (or minimum component)
    inline float MaxComponent() const { float temp = (d[1] > d[0]? d[1] : d[0]); return (d[2] > temp? d[2] : temp); }
	inline float MinComponent() const { float temp = (d[1] > d[0]? d[0] : d[1]); return (d[2] > temp? temp : d[2]); }

	// Static methods.  Useful for defining commonly used colors
	static Color Black( void )  { return Color(0,0,0); }
	static Color White( void )  { return Color(1,1,1); }
};




inline Color::Color() 
{ 
	d[0] = d[1] = d[2] = 0; 
}

inline Color::Color(vec3 color) 
{
    d[0] = color.X();
    d[1] = color.Y();
    d[2] = color.Z();
}

inline Color::Color(float r, float g, float b) 
{ 
	d[0] = r; 
	d[1] = g; 
	d[2] = b; 
}
	
inline Color::Color(const Color& copy) 
{ 
	d[0] = copy.d[0]; 
	d[1] = copy.d[1]; 
	d[2] = copy.d[2]; 
}
    


#endif

