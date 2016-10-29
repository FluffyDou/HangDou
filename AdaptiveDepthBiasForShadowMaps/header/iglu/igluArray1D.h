/******************************************************************/
/* IGLUArray1D.h                                                  */
/* -----------------------                                        */
/*                                                                */
/* The file defines an very basic (templated) dynamic 1D array    */
/*    class.  You can actually use the standard template library  */
/*    for this, without a performance hit when running optimized  */
/*    versions of your code.  But the STL vector class sucks when */
/*    running in Debug mode.  (And you'll waste time trying to    */
/*    optimize away the 10-30% of the time spent in its methods)  */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IGLU_ARRAY1D_H
#define IGLU_ARRAY1D_H

namespace iglu {


template<class T> class IGLUArray1D {
  T* objs;
  unsigned int size, memAlloc;
public:

  // Allocate a new 1D array
  IGLUArray1D();
  ~IGLUArray1D(); 

  // Get access to the elements that have already been added to the array.
  //   NOTE: no array bounds checks are performed.
  inline const T& operator[](unsigned int n) const { return objs[n]; }
  inline T& operator[](unsigned int n) { return objs[n]; }
    
  // Find out how many elements have been added to the array
  inline unsigned int Size() const { return size; }

  // Set the array size
  void SetSize( unsigned int n );

  // Add a new array element to the end of the array (returns index of added element)
  unsigned int Add(const T&);

  // Get a pointer to the array elements
  T* GetData( void );

  // A pointer to a IGLUArray1D<T> could have type IGLUArray1D<T>::Ptr
  typedef IGLUArray1D<T> *Ptr;

};




// Setup an empty array
template<class T>
IGLUArray1D<T>::IGLUArray1D() : memAlloc(0), size(0), objs(0)
{
}	

// Free an array
template<class T>
IGLUArray1D<T>::~IGLUArray1D()
{
  if (objs) delete [] objs;
}


// Add an element to an array
template<class T>
unsigned int IGLUArray1D<T>::Add( const T& obj )
{
  if (size+1 > memAlloc || !objs)
  {
	  // If there's not enough space in the array, allocate the 
	  //    larger of 10 or 2*currSize new elements for the array
	  unsigned int growSize = memAlloc*2 > 10 ? memAlloc*2 : 10;	
	  unsigned int newSize = memAlloc+growSize;
	  T *newObjs = new T[newSize];      // allocate new memeory
	  for (unsigned int i=0;objs && i<size;i++)  // copy the data from the old array
		  newObjs[i] = objs[i];
	  if (objs) delete [] objs;         // get rid of the old data
	  objs = newObjs;                   // set the array ptr to the new array
	  memAlloc = newSize;               // remember how big our new allocation is
  }
  objs[size++] = obj;                   // add our new element to the array
  return size-1;
}

template<class T>
void IGLUArray1D<T>::SetSize( unsigned int n )
{
	unsigned int newSize = n > memAlloc ? n : memAlloc;
	if (newSize <= memAlloc) return;
	T *newObjs = new T[newSize];      // allocate new memeory
	if (objs)
	{
		for (unsigned int i=0;objs && i<size;i++)  // copy the data from the old array
			newObjs[i] = objs[i];
		delete [] objs;
	}
	memAlloc = newSize;
	size = newSize;
	objs = newObjs;
}

// Get a pointer to the array data
template<class T>
T* IGLUArray1D<T>::GetData( void )
{
  return objs;
}



}

#endif
