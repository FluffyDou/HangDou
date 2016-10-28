/******************************************************************/
/* Light.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a light class that stores basic information   */
/*     about the light (position and color).                      */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/


#ifndef LIGHT_H
#define LIGHT_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DataTypes/datatypes.h"

class Scene;

class Light
{
  vec3  pos;
  Color intense;
public:
    Light();
	Light( const vec3 &pos, const Color &lcolor ): pos(pos), intense(lcolor) {}
    Light( FILE *f, Scene *s ) {};

    ~Light();

	inline vec3  &GetPosition( void )   { return pos; }
	inline Color &GetColor( void )      { return intense; }
    inline vec3  GetIntense( void )    { return vec3( intense.Red(), intense.Green(), intense.Blue() ); }

    void SetPosition(vec3 newPos);
};


#endif

