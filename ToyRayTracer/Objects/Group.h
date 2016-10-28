/******************************************************************/
/* Group.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines an object container (i.e., an object that     */
/*     contains numerous others).  The intersection technique     */
/*     simply loops through all objects that have been added to   */
/*     the group, calls its Intersect() routine, and returns the  */
/*     closest object.  Obviously, this is NOT the most efficient */
/*     container class possible.                                  */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef GROUP_H
#define GROUP_H

#include "Objects/Object.h"
#include "Core/Ray.h"
#include "DataTypes/datatypes.h"

class Group : public Object {
private:
	Array1D<Object *> objs;
	int size;
public:

	// Set up a default (empty) group.
	Group();   

	// Currently, one cannot read a group of objects from a file.  If you wish to
	//    define how this would be done, you can change this.
	Group( FILE *f, Scene *s ) { 
		printf("Constructor Group::Group( FILE *f ) called.  This is not implemented!"); 
	}

	// Free all the memory inside this group.
	virtual ~Group();

	// Add an object to the group.
	void Add( Object *obj );

	// Get an object from the group.  (Note: no bounds checking is done)
	Object *Get( int i );
	const Object *Get( int i ) const;

	// Get the number of objects in the group.
	inline int GetSize( void ) const   { return size; }

	// Here's the actual intersection routine.  Given a ray, intersect it with all
	//    the objects in this group (one at a time).  This could/should be done in a
	//    much more intelligent way (in classes inheriting from Group()).
	virtual void Intersect( Ray &ray );

};


#endif

