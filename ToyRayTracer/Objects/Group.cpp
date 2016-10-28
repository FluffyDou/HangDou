/******************************************************************/
/* Group.cpp                                                      */
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

#include "Group.h"
#include <assert.h>

using namespace std;

Group::Group() : size(0)
{
}

Group::~Group()
{
}

void Group::Add(Object* obj)
{
  objs.Add(obj);
  size++;
}

Object *Group::Get( int i ) 
{
  return objs[i];
}

const Object* Group::Get(int i) const
{
    return objs[i];
}

void Group::Intersect( Ray &ray ) 
{
	Object **ptr = objs.GetData();
	for(int i=0; i<size; i++)
		ptr[i]->Intersect( ray );

}

