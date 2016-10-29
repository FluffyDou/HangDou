/******************************************************************/
/* igluOBJMaterial.h                                              */
/* -----------------------                                        */
/*                                                                */
/* Defines a reader for OBJ material files                        */
/*                                                                */
/* Chris Wyman (10/20/2011)                                       */
/******************************************************************/

#ifndef IGLU_OBJ_MATERIAL_H
#define IGLU_OBJ_MATERIAL_H

namespace iglu {

class IGLUOBJMaterial 
{
public:
	IGLUOBJMaterial();
	~IGLUOBJMaterial();

	/////////////////////////////////////////////////////////////
	// MTL characteristics that this class respects

	char *m_matlName;            // Name specified in MTL file
	vec3 m_amb, m_dif, m_spec;   // Ka, Kd, Ks from MTL file
    float m_shininess;           // Ns from MTL file
	float m_dissolve;            // d from MTL file
	float m_idxRefract;          // Ni from MTL file
	int m_ambTexId, m_difTexId, m_specTexId; // IDs for textures.


	/////////////////////////////////////////////////////////////
	// MTL characteristics that this class does not (yet) respect

	int m_illumModel;            // illum from MTL file
	vec3 m_trans;                // Tf from MTL file
	float m_reflSharpness;       // sharpness from MTL file
};

class IGLUOBJMaterialReader : public IGLUFileParser
{
public:
	// Read & create a bunch of materials from a .mtl file
	IGLUOBJMaterialReader( char *filename );

	// Arrays storing lists of textures & materials used by all OBJ materials loaded thus far.
	//    These values are set during the constructor IGLUOBJMaterialReader().  Each new object
	//    (i.e., file that is parsed) adds to any previous values stored in these static arrays.
	static IGLUArray1D<IGLUTexture2D::Ptr>   s_matlTexture;
	static IGLUArray1D<IGLUOBJMaterial *>    s_matl;

	//////////////////////////////////////////////////////////////////////////
	// The following static method can be used to create a default material
	//   --> This should be called automatically by any instantiation of a
	//       IGLUOBJMaterialReader object.  But this public method can be used
	//       to start creating materials manually.
	static void AddDefaultMaterial( void );

	//////////////////////////////////////////////////////////////////////////
	// The following static methods can interact with the prior two arrays  
	//   --> Lookup a material by name
	//////////////////////////////////////////////////////////////////////////
	static IGLUOBJMaterial *GetNamedMaterial( char *mtlName );
	static int              GetNamedMaterialId( char *mtlName );

	//////////////////////////////////////////////////////////////////////////
	// Once ALL OBJ & MTL files have been processed, the following method 
	// should be called, to create the textures that are used in IGLU shaders
	//////////////////////////////////////////////////////////////////////////
	static bool FinalizeMaterialsForRendering( void );

	//////////////////////////////////////////////////////////////////////////
	// After calling the FinalizeMaterialsForRendering() method, the following
	// data can be accessed.
	//////////////////////////////////////////////////////////////////////////

	// A created texture buffer that stores the material properties.  Essentially,
	//    this texture stores all the s_matl information in a texture buffer that
	//    can be indexed into by a material ID in shaders.
	static IGLUTextureBuffer::Ptr            s_matlCoefBuf;   

	// A created texture array that stores the material textures.  Essentially,
	//    this is an array of all the textures in s_matlTexture, resized so they
	//    all fit into a a tex array with fixed width & height.  The layer number
	//    is the index into the s_matlTexture array.
	static IGLUTexture::Ptr                  s_matlTexArray;    

	// Check the array s_matlTexture to see if a texture from file <filename> 
	//    has already been added to our list.  If not, return -1, else return
	//    array index.
	static int GetTextureID( const char *filename );
};


}

#endif

