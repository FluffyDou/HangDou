/*****************************************************************************
** iglu.h                                                                   **
** ------                                                                   **
**                                                                          **
** Main header file for the (I)owa Open(GL) (U)tility functions.            **
**                                                                          **
** Chris Wyman (9/26/2011)                                                  **
*****************************************************************************/

#ifndef IGLU_H
#define IGLU_H

// Some #define's to set the version of IGLU.  Currently, this is a very early 
//     release with somewhat limited and untested functionality.
// If you update the code and check it back in, please update at least the
//     build number (by incrementing by 1).
#define IGLU_MAJOR_VERSION  0
#define IGLU_MINOR_VERSION  4
#define IGLU_BUILD_NUMBER   0

// We need some OpenGL Extension handler.  We use GLEW.  If you prefer
//    glext, glee, or some other extension library/headers, you should
//    be able to insert it here (instead).
#include <GL/glew.h>

// Common error codes that may be used by various IGLU classes/methods/functions
#include "iglu/errors/igluErrors.h"
#include "iglu/errors/igluErrorExit.h"
#include "iglu/errors/igluWarning.h"

// Simple, basic math operations (not or non-portably defined elsewhere)
#include "iglu/igluMath.h"

// Vector datatypes (igluVector.h #includes specific files)
#include "iglu/igluVector.h"

// Matrix datatypes
#include "iglu/igluMatrix4x4.h"

// Some vector-based utilities
#include "iglu/igluOrthoNormalBasis.h"

// Bounding region datatypes
#include "iglu/igluInterval.h"
#include "iglu/igluRange.h"

// Array datatypes
#include "iglu/igluArray1D.h"

// Simple OpenGL Context state
#include "iglu/glstate/igluGLContext.h"

// IGLU constants 
#include "iglu/consts/igluColors.h"

// IGLU interface callback class
#include "iglu/interface/igluCallback.h"

// IGLU Variables
#include "iglu/variables/igluVariable.h"
#include "iglu/variables/igluInt.h"
#include "iglu/variables/igluFloat.h"
#include "iglu/variables/igluBool.h"

// OpenGL buffer encapsulations
#include "iglu/igluBuffer.h"
#include "iglu/buffers/igluUniformBuffer.h"

// Encapsulation for timing routines
#include "iglu/igluGPUTimer.h"
#include "iglu/igluCPUTimer.h"
#include "iglu/igluFrameRate.h"

// Random number generation
#include "iglu/igluRandom.h"
#include "iglu/igluDRand48.h"

// Sampling patterns and quasi-random numbers
#include "iglu/igluSampling.h"

// Cameras
#include "iglu/camera/igluPinholeCamera.h"

// OpenGL GLSL shader encapsulations
#include "iglu/igluShaderVariable.h"  
#include "iglu/igluShaderStage.h"  
#include "iglu/igluShaderProgram.h" 

// Capture utilities
#include "iglu/igluFrameGrab.h"

// Image input/video IO utilities
#include "iglu/igluImage.h"
#include "iglu/igluVideo.h"
#include "iglu/igluVideoCapture.h"

// Texturing utilities
#include "iglu/igluTexture2D.h"
#include "iglu/igluPixelBufferTexture2D.h"
#include "iglu/igluTextureLightprobeCubemap.h"
#include "iglu/igluTextureBuffer.h"
#include "iglu/igluRandomTexture2D.h"
#include "iglu/igluRandomTexture3D.h"
#include "iglu/igluVideoTexture2D.h"

// Render-to-texture utilities
#include "iglu/igluRenderTexture.h"
#include "iglu/igluRenderTexture2D.h"
#include "iglu/igluRenderTexture2DArray.h"
#include "iglu/igluRenderTexture2DMultisample.h"
#include "iglu/igluFramebuffer.h"

// Drawing utilities
#include "iglu/igluDrawingUtils.h"
#include "iglu/igluGeomUtils.h"

// File parsing utilities
#include "iglu/igluParsing.h"
#include "iglu/igluModels.h"

// User interface utilities
#include "iglu/interactors/igluMouseInteractor.h"
#include "iglu/interactors/igluTrackball.h"
#include "iglu/interactors/igluMouseTranslator.h"
#include "iglu/interactors/iglu2DMouseRotator.h"
#include "iglu/interactors/igluCameraRotator.h"
#include "iglu/interface/igluMainWindow.h"

// OpenGL state utility classes
#include "iglu/glstate/igluTransformFeedback.h"
#include "iglu/glstate/igluVertexArrayObject.h"

// IGLU OpenGL Windowing utilities
#include "iglu/window/igluWindow.h"
#include "iglu/window/igluMultiDisplayWindow.h"
#include "iglu/window/igluDisplayMode.h"
#include "iglu/window/igluWidgetWindow.h"
#include "iglu/window/igluMultiDisplayWidgetWindow.h"

// IGLU's simple file path & file IO utilities
#include "iglu/igluFilePath.h"
#include "iglu/igluFileIO.h"

// Ok.  Maybe there's some macro definitions of min & max...  Let's get rid of those
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#endif