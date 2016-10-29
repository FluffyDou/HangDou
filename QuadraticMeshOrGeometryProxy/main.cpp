#include <GL/glew.h>
#include <stdio.h>

// All headers are automatically included from "iglu.h"
#include <GL/glut.h> // This is only for track ball. Get rid of it from iglu when you have time.
#include "iglu.h"

#include "header/uiViewInteraction.h"

using namespace iglu;

ivec2 winSize = ivec2(1024, 1024);
float constScale = 1000.0;
// These are handles to the two windows, should you need to modify their
//   behavior after window creation.
IGLUWindow::Ptr         myWin = 0;
IGLUWidgetWindow::Ptr   uiWin = 0;

// All the variables in the control panel
IGLUBool                printLightInfo     ( false,                                     "Print light info" );
IGLUBool                displayGPUTime     ( false,                                     "Display GPU timing" );
IGLUBool                displayFPS         ( false,                                     "Display framerate" );
IGLUBool				useColorMask       ( true,                                      "Enable color masking for shadow map" );
IGLUBool				useLambertian      ( true,                                      "Use Lambertian shading" );
IGLUBool				usePhong           ( false,                                     "Use Phong shading" );
IGLUBool                rotateFlag         ( false,                                     "Rotate the object or not" );
IGLUBool                showText           ( true,                                      "Show the text or not" );

IGLUInt                 smapRes            ( 1024,   IGLURange<int>(2,8192),   1,       "Shadow map resolution" );
IGLUInt                 displayMode        ( 0,      IGLURange<int>(0,1),      1,       "Display Mode" );
IGLUInt                 methodFlag         ( 0,      IGLURange<int>(0,1),   1,          "Method ID" );
IGLUInt                 jointNum           (50, IGLURange<int>(2, 8192), 1,             "Joint num on each hair");
IGLUInt                 lineSegNum         (50, IGLURange<int>(2, 8192), 1,             "Hair number");
IGLUInt                 vertexID           (0, IGLURange<int>(0, 2000), 1,              "Debug: vertexID");

IGLUFloat               hairLineWidth           (0.03, IGLURange<float>(0, 100), 0.01,      "Hair line width");
IGLUFloat               hairLineRasterSize      (30.0,      IGLURange<float>(1, 100), 1.0,  "Hair line raster size");
IGLUFloat               hairPointRasterSize     (12.0, IGLURange<float>(1, 100), 1.0,        "Hair point raster size");
IGLUFloat               lineWidth          ( 0.1,    IGLURange<float>(0, 1),  0.01,       "Rasterize line width");
IGLUFloat               proxyCylinderRadius ( 0.015,    IGLURange<float>(0, 1), 0.001,    "Line width");
IGLUFloat               proxySphereRadius   ( 0.015,    IGLURange<float>(0, 1), 0.001,    "Jiont sphere radius");
IGLUFloat               visualizeLayer     ( 0.0,    IGLURange<float>(0,1),   1.0,      "Omni shadow map layer to display" );
IGLUFloat               eyeNear            ( 1.0,    IGLURange<float>(0,10), 0.1,       "Eye near" );
IGLUFloat               eyeFar             ( 10.0,   IGLURange<float>(10,120), 1.0,     "Eye far" );
IGLUFloat               eyeFOV             ( 38.5,   IGLURange<float>(20,150), 1,       "Viewer field-of-view" );
IGLUFloat               lightMovement      ( 3.5,    IGLURange<float>(0,4), 0.04,       "Light movement" );
IGLUFloat               lightNear          ( 0.10,   IGLURange<float>(0.01,1),  0.01,   "Light near plane" );
IGLUFloat               lightFar           ( 8.3,    IGLURange<float>(5.0,100), 0.1,    "Light far plane" );
IGLUFloat               lightFOV           ( 70,     IGLURange<float>(30,150), 1,       "Shadow map field-of-view" );
IGLUFloat               lightIntense       ( 1.0,    IGLURange<float>(0.0,15.0), 0.1,   "Light intensity" );
IGLUFloat               phongAlpha         ( 2.0,    IGLURange<float>(0.1,15.0), 0.01,  "Phong alpha" );
IGLUFloat               constAmbient       ( 0.2,    IGLURange<float>(0.0,1.0), 0.01,   "Constant ambient" );

IGLUFloat               objX               ( 8.7,    IGLURange<float>(-30,30), 0.1,     "Model X" );
IGLUFloat               objY               ( -0.8,   IGLURange<float>(-30,30), 0.1,     "Model Y" );
IGLUFloat               objZ               ( -3.5,   IGLURange<float>(-30,30), 0.1,     "Model Z" );
IGLUFloat               objS               ( 2.0,    IGLURange<float>(0,10), 0.1,       "Model scale" );

IGLUFloat               rotateAngle1        ( 180.0,    IGLURange<float>(0,360), 1.0,   "Obj1 Rotate Amount" );

IGLUFloat               objX2               ( 3.1,   IGLURange<float>(-20,20), 0.1,     "Box X" );
IGLUFloat               objY2               ( -1.4,  IGLURange<float>(-20,20), 0.1,     "Box Y" );
IGLUFloat               objZ2               ( -6.2,  IGLURange<float>(-20,20), 0.1,     "Box Z" );
IGLUFloat               objS2               ( 1.6,   IGLURange<float>(0,2), 0.1,        "Box scale" );

IGLUFloat               rotateAngle2        ( 180.0,   IGLURange<float>(0,360), 1.0,   "Obj2 Rotate Amount" );

IGLUFloat               objX3               ( 0.0,   IGLURange<float>(-10,10), 0.1,     "Object3 X" );      
IGLUFloat               objY3               ( 8.2,   IGLURange<float>(-10,10), 0.1,     "Object3 Y" );
IGLUFloat               objZ3               ( -5.0,  IGLURange<float>(-10,10), 0.1,     "Object3 Z" );
IGLUFloat               objS3               ( 30.0,  IGLURange<float>(0,30), 0.1,       "Object3 scale" );

IGLUFloat               rotateAngle         ( 0.2,   IGLURange<float>(0.1,10.0), 0.1,   "Rotate angle" );

float *vboData;
GLuint vboHair, vaoHair;
GLuint vboQuad, vaoQuad;
GLuint vboCube, vaoCube;
GLuint ac_buffer    = 0;
//GLuint imageCounter = 0;// image unit
//IGLUTexture2D::Ptr imageCounterPtr = NULL;
//GLuint tbo = 0; // transform feedback

// variables for glMultiDrawArray --- since we are drawing line strips
GLint*     firstElement;
GLsizei*   countElement;

int      quadNum = 6;
GLint*   quadFirstElement;
GLsizei* quadCountElement;

GLuint cubeElementbuffer;
//GLint         prevProg;
//const GLsizei elementNum = 1;   // Number of primitives, i.e.: line segments
//int           numVertices = 3;

//float rotateAngle    = 1.0; // Model rotates in the animation
float accumulateAngle  = 0.0;   // Rotate accumulation
//int methodID = 0;

int  lightFlag  = 0;    // Interact with light or not --- initial is 0, saying no.
int  oldDisplayMode = displayMode;
// For timing
double gatheredTime[2];

IGLUGPUTimer::Ptr gpuTimer = 0;

//The scene geometry consists of two objects: a wavefront obj file and a square 
IGLUOBJReader::Ptr      objCow     = 0;
IGLUOBJReader::Ptr      objCube    = 0;
IGLUOBJReader::Ptr      sphere     = 0;
IGLUOBJReader::Ptr      boxReader  = 0;
IGLUOBJReader::Ptr      bigScene   = 0;

//Shader pool --- ATTENTION: All the "perspective" means traditional shadow map, since that is perspective projection on light space
//These are accessed using the enum below
IGLUShaderProgram::Ptr  shaders[10] = {0};
enum ShaderType{ SM_PERSPECTIVE = 0, SHADE_WIREFRAME_GEO, SHADE_WIREFRAME_NO_GEO, SHADE_WIREFRAME_PROXY_GEO, SHADE_WIREFRAME_PROXY_NO_GEO, SHADE_TUBE_GEO, SHADE_TUBE_NO_GEO, SHADE_QUAD};

uiViewInteraction *myViewer = 0;
uiViewInteraction *myLight  = 0;

//Trackball for interaction
IGLUTrackball::Ptr	   ball = 0;
// Trackball for point light source
IGLUTrackball::Ptr viewBall = 0;

//The frame buffer holding our shadow map
IGLUFramebuffer::Ptr    shadowMapFBO, secondDepthFBO, parabSMapFBO, parabSecondDepthFBO;

// Location of the light and the eye
vec3 lightPos = vec3(2.78,5.25,2.795);
vec3 lightAt  = vec3(2.78, 2.73, 2.795);
vec3 eyePos   = vec3(2.78, 2.73, -8.00);
vec3 eyeAt    = vec3(2.78, 2.73, 2.795);

// Matrices for setting up the view from the eye
IGLUMatrix4x4   eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
IGLUMatrix4x4   eyeView = IGLUMatrix4x4::LookAt( eyePos, eyeAt, vec3::YAxis() );

// Matrices for setting up the view from the light (i.e., the shadow map )
IGLUMatrix4x4   lightView = IGLUMatrix4x4::LookAt( lightPos, lightAt, vec3::ZAxis()+vec3::XAxis() );
IGLUMatrix4x4   lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1.0, lightNear, lightFar ); 

// Matrices for positioning our geometry relative to the world origin
IGLUMatrix4x4   objPosition   = IGLUMatrix4x4::Translate( objX, objY, objZ ) * IGLUMatrix4x4::Scale( objS ) * 
                                IGLUMatrix4x4::Rotate( rotateAngle1, vec3::XAxis() );

IGLUMatrix4x4   objPosition2  = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * IGLUMatrix4x4::Scale( objS2 ) * 
                                IGLUMatrix4x4::Rotate( rotateAngle2, vec3::XAxis() );

IGLUMatrix4x4   bigScenePos   = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );

// Generate some line segments
void GenerateLineSegPrimitives(int jointNum, int lineSegNum, float radius);

// Shade 
void shadeWireFrameGeo();
// Shade
void shadeWireFrameProxyGeo();
// Shade
void shadeHairGeo();

// Update eye projection --- update when modifying eye near and eye far
void UpdateEyeProjection();

// Initialize some parameters for better view --- just for convenience
void sceneInitialize();

// Display current using method
void DisplayMethod();

/************ Render loop ************/
void display ( void )	
{
    //if (displayGPUTime) gpuTimer->Start();
    //glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // Simple object rotation
    if( rotateFlag )
    {
        objPosition2     = objPosition2 * IGLUMatrix4x4::Rotate(rotateAngle, vec3::YAxis());
        accumulateAngle += rotateAngle;
    }

    if (displayMode == 0)
        shadeHairGeo();
    else if (displayMode == 1)
        shadeWireFrameProxyGeo();

    DisplayMethod();

    //exit(0);
}


void shadeWireFrameProxyGeo()
{
    //glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glDisable(GL_BLEND);
    //glEnable(GL_DEPTH_TEST);

//gpuTimer->Tick();

    shaders[SHADE_WIREFRAME_PROXY_GEO]["model"] = objPosition * ball->GetMatrix();
    objCow->Draw(shaders[SHADE_WIREFRAME_PROXY_GEO]);
    //shaders[SHADE_WIREFRAME_PROXY_GEO]["model"] = objPosition2 * ball->GetMatrix();
    //objCube->Draw(shaders[SHADE_WIREFRAME_PROXY_GEO]);

//gatheredTime[0] = gpuTimer->Tick();

}

void shadeHairGeo()
{
    //glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //shaders[SHADE_TUBE_NO_GEO]["znear"] = eyeNear;
    // Draw the geometry, modifying per-geometry shader variables as needed
    //shaders[SHADE_PERSPECTIVE]->Enable();
    //shaders[SHADE_TUBE_NO_GEO]["model"] = objPosition * ball->GetMatrix();

    //glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    //glEnable(GL_LINE_SMOOTH);

    // bind vao
    glBindVertexArray(vaoHair);

    shaders[SHADE_TUBE_GEO]["model"] = objPosition * ball->GetMatrix();
    shaders[SHADE_TUBE_GEO]->Enable();

//gpuTimer->Tick();
    glMultiDrawArrays(GL_LINE_STRIP, firstElement, countElement, lineSegNum);
//gatheredTime[0] = gpuTimer->Tick();

    shaders[SHADE_TUBE_GEO]->Disable();
    glBindVertexArray(0);

}


/************************** Call Back Functions ******************************/

void updateLightPos( IGLUMatrix4x4 light_view )
{
    //shaders[SM_PERSPECTIVE]["view"] = light_view;

    //shaders[SHADE_WIREFRAME_GEO]["lightView"]    = light_view;
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightView"] = light_view;
}

void updateEyePos( IGLUMatrix4x4 eye_view )
{
    shaders[SHADE_WIREFRAME_GEO]["view"]          = eye_view;
    //shaders[SHADE_WIREFRAME_NO_GEO]["view"]       = eye_view;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["view"] = eye_view;
    shaders[SHADE_WIREFRAME_PROXY_GEO]["view"]    = eye_view;
    //shaders[SHADE_TUBE_NO_GEO]["view"]            = eye_view;
    shaders[SHADE_TUBE_GEO]["view"]               = eye_view;
    //shaders[SHADE_QUAD]["view"]                   = eye_view;
}


// Track any updates to the trackball matrix when the mouse moves 
void Motion(int x, int y)
{
    if( lightFlag == 0 )
    {
        myViewer->MouseMotion(x, y);
        updateEyePos( myViewer->getEyeView() );
    }
    else
    {
        myLight->MouseMotion(x, y);
        updateLightPos( myLight->getEyeView() );
    }

}

// When the user clicks/releases in the main window, start/stop tracking with our trackball
void Button(int button, int state, int x, int y  )
{
    if( lightFlag == 0 )
        myViewer->MouseButton(button,state,x, y);
    else
        myLight->MouseButton(button,state,x, y);
}

// keyboard call back
void KeyBoard(unsigned char key, int a, int b)
{   
    switch ( key )
    {
    case 'c':
    case 'C':
        lightFlag = 1-lightFlag; // switch between light and view
        break;
    case 'p':
    case 'P':
        std::cout << "light pos and at:\n";
        myLight->getEyePos().Print();
        myLight->getEyeAt().Print();
        std::cout << "eye pos and at:\n";
        myViewer->getEyePos().Print();
        myViewer->getEyeAt().Print();
        break;
    default: break;
    }

    if( lightFlag == 0 )
    {
        myViewer->KeyBoard(key, a, b);
        updateEyePos( myViewer->getEyeView() );

        //myViewer->MoveCountStart();
        //myViewer->setKey( key );
    }
    else
    {
        myLight->KeyBoard(key, a, b);
        updateLightPos( myLight->getEyeView() );
    }

}

// When the user clicks the "Display Framerate" checkbox, we need to update some
//    window parameters, so this callback function needs to be called
void ToggleFramerate( IGLUVariable *uiBool )
{
	myWin->DisplayFramerate( ((IGLUBool *)uiBool)->GetValue() );
}

// Display method
void DisplayMethod()
{
    char buf[64];
    int  fontHeight = 55;

    if (displayMode == 0)
        sprintf(buf, "Line strips (input) by quadratic mesh.");
    else if (displayMode == 1)
        sprintf(buf, "Tubular style wireframe of triangle mesh (cow) by quadratic mesh.");
    if(showText) IGLUDraw::DrawText( IGLU_FONT_FIXED, 0, winSize.Y()-fontHeight, buf, 1.0 );
}
//
void UpdateLineWidth(IGLUVariable *uiFloat)
{
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = 1.0 - ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_GEO]["lineWidth"]    = 1.0 - ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["lineWidth"] = 1.0 - ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = lineWidth;
    //shaders[SHADE_WIREFRAME_GEO]["lineWidth"] = ((IGLUFloat *)uiFloat)->GetValue();
}

void UpdateProxyCylinderRadius(IGLUVariable *uiFloat)
{
    float radiusSquare = ((IGLUFloat *)uiFloat)->GetValue() * ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = 1.0 - ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["cylinderRadiusSquare"] = radiusSquare;
    shaders[SHADE_WIREFRAME_PROXY_GEO]["cylinderRadiusSquare"]    = radiusSquare;
    //shaders[SHADE_TUBE_NO_GEO]["cylinderRadiusSquare"]            = radiusSquare;
    shaders[SHADE_TUBE_GEO]["cylinderRadiusSquare"]               = radiusSquare;
    shaders[SHADE_TUBE_GEO]["radius"]                             = ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = lineWidth;
    //shaders[SHADE_WIREFRAME_GEO]["lineWidth"] = ((IGLUFloat *)uiFloat)->GetValue();
}

void UpdateProxySphereRadius(IGLUVariable *uiFloat)
{
    float radiusSquare = ((IGLUFloat *)uiFloat)->GetValue() * ((IGLUFloat *)uiFloat)->GetValue();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = 1.0 - ((IGLUFloat *)uiFloat)->GetValue();
    shaders[SHADE_WIREFRAME_PROXY_GEO]["sphereRadiusSquare"] = radiusSquare;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["sphereRadiusSquare"] = radiusSquare;
    //shaders[SHADE_TUBE_NO_GEO]["sphereRadiusSquare"] = radiusSquare;
    shaders[SHADE_TUBE_GEO]["sphereRadiusSquare"] = radiusSquare;
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = lineWidth;
    //shaders[SHADE_WIREFRAME_GEO]["lineWidth"] = ((IGLUFloat *)uiFloat)->GetValue();
}

// When the user resizes the shadow map with the slider, we need to actually resize
//    the frame buffer object, so this callback function needs to be called
//void UpdateShadowMapResolution( IGLUVariable *uiInt )
//{
//	int size = ((IGLUInt *)uiInt)->GetValue();
//
//	shadowMapFBO->Resize( size, size );
//    secondDepthFBO->Resize( size, size );
//
//    parabSMapFBO->Resize(size, size);
//    parabSecondDepthFBO->Resize(size, size);
//}

void UpdateEyeProjection( IGLUVariable *uiFloat )
{

    eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
    myViewer->setEyeProj(eyeProj);

    shaders[SHADE_WIREFRAME_GEO]["proj"]          = myViewer->getEyeProj();
    //shaders[SHADE_WIREFRAME_NO_GEO]["proj"]       = myViewer->getEyeProj();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["proj"] = myViewer->getEyeProj();
    shaders[SHADE_WIREFRAME_PROXY_GEO]["proj"]    = myViewer->getEyeProj();
    //shaders[SHADE_TUBE_NO_GEO]["proj"]            = myViewer->getEyeProj();
    shaders[SHADE_TUBE_GEO]["proj"]               = myViewer->getEyeProj();
    //shaders[SHADE_QUAD]["proj"]              = myViewer->getEyeProj();
}

void UpdateLightProjection( IGLUVariable *uiFloat )
{
	lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1, lightNear, lightFar );
    myLight->setEyeProj(lightProj);
	//shaders[SM_PERSPECTIVE]["proj"]          = lightProj;

    //shaders[SHADE_WIREFRAME_GEO]["lightProj"]  = lightProj;
    //shaders[SHADE_WIREFRAME_GEO]["lightFar"]   = lightFar;
    //shaders[SHADE_WIREFRAME_GEO]["lightNear"]  = lightNear;
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightProj"] = lightProj;
    //shaders[SHADE_WIREFRAME_GEO]["lightFar"]   = lightFar;
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightNear"] = lightNear;
}

void UpdateObjPosition( IGLUVariable *uiFloat )
{
    objPosition = IGLUMatrix4x4::Translate( objX, objY, objZ ) * 
                  IGLUMatrix4x4::Scale( objS ) * 
                  IGLUMatrix4x4::Rotate( rotateAngle1, vec3::YAxis() );

    //shaders[SHADE_WIREFRAME_GEO]["model"]          = objPosition * ball->GetMatrix();
    //shaders[SHADE_WIREFRAME_NO_GEO]["model"]       = objPosition * ball->GetMatrix();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["model"] = objPosition * ball->GetMatrix();
}

void UpdateObjPosition2( IGLUVariable *uiFloat )
{
    objPosition2 = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * 
                   IGLUMatrix4x4::Scale( objS2 ) * 
                   IGLUMatrix4x4::Rotate( rotateAngle2, vec3::XAxis() );

    //shaders[SHADE_WIREFRAME_GEO]["model"] = objPosition2 * ball->GetMatrix();
    //shaders[SHADE_WIREFRAME_NO_GEO]["model"] = objPosition2 * ball->GetMatrix();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["model"] = objPosition2 * ball->GetMatrix();
}

//void UpdateObjPosition3( IGLUVariable *uiFloat )
//{
//    bigScenePos  = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );
//}

// Basically we ...
void GenerateLineSegPrimitives(int jointNum, int lineSegNum, float radius)
{

    ////// allocate cpu memory of firstElement & countElement
    firstElement = new GLint[lineSegNum];
    countElement = new GLsizei[lineSegNum];

    int verticesNum = 0; // The vertices offset for each line strip

    // this loop gets the number of vertices and set values for firstElement & countElement simultaneously        
    for (int i = 0; i < lineSegNum; ++i)
    {
        firstElement[i] = verticesNum;
        verticesNum    += jointNum;
        countElement[i] = jointNum;
    }

    int dataSize = verticesNum * 3; // 3 means every vertex has 3 values, in the future we need this adaptive
    vboData = new float[dataSize];

    float phiDelta   = 0.5*3.1415926585897 / (float)jointNum;
    float thetaDelta = 2.0*3.1415926585897 / (float)lineSegNum;

    float meanX = 0.0; float meanY = 0.0; float meanZ = 0.0;

    ////// Generate the curve geometry --- the circle radius is 1.0
    int dataIndex = 0;
    for (int theta = 0; theta < lineSegNum; ++theta)
    for (int phi = 0; phi < jointNum; ++phi)
    {
        float x = radius * sin(phi*phiDelta) * cos(theta*thetaDelta);
        float y = radius * sin(phi*phiDelta) * sin(theta*thetaDelta);
        float z = radius * cos(phi*phiDelta);

        vboData[dataIndex++] = x;
        vboData[dataIndex++] = y;
        vboData[dataIndex++] = z;

        meanX += x;  meanY += y; meanZ += z;
    }
    meanX /= (float)verticesNum;  meanY /= (float)verticesNum; meanZ /= (float)verticesNum;
    // Move the model to the origin
    for (int i = 0; i < dataIndex; i += 3)
    {
        vboData[i]     -= meanX;
        vboData[i + 1] -= meanY;
        vboData[i + 2] -= meanZ;
    }

    ////// Upload al the data to GPU. Generate VAO and VBO.
    // Todo: convert all these to "glDrawIndirect" to reduce CPU latency.
    glGenVertexArrays(1, &vaoHair);
    glBindVertexArray(vaoHair);

    //if (!glIsBuffer(vboHair))
    glGenBuffers(1, &vboHair);

    glBindBuffer(GL_ARRAY_BUFFER, vboHair);
    //glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(float)*9, vboData, GL_STATIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(float)* dataSize, vboData, GL_STATIC_DRAW);

    // Edit the first input stream in the shader --- we only have one vbo or input stream here. So we enable index 0.
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0, 	                                // vertex coordinate attribute --- layout(location = 0)
        3,                	                // number of elements per vertex, here (x,y,z)
        GL_FLOAT,          	                // the type of each element
        GL_FALSE,          	                // no normalization
        0,                                  // non zero only vbo stores interleaved data like: vvvnnnttvvvnnntt...
        0                                   // offset of first element
        );

    glBindVertexArray(0);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (0 != vboData) delete[] vboData;
}


// Basically we create a cube in the form of quads
void GenerateQuad()
{
    ////// allocate cpu memory of firstElement & countElement
    // Note VBO only stores vertices (and maybe other attributs like normal). It is element arrays that stores the indices to those vertices.
    quadFirstElement = new GLint[1];
    quadCountElement = new GLsizei[1];
    int verticesNum = 4;
    quadFirstElement[0] = 0;
    quadCountElement[0] = verticesNum;

    // 3 means every vertex has 3 values (x,y,z)
    int dataSize = verticesNum * 3;
    vboData = new float[dataSize];
    int dataIndex = 0;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = 0.0;

    vboData[dataIndex++] = 0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = 0.0;

    vboData[dataIndex++] = 0.5;
    vboData[dataIndex++] = 0.5;
    vboData[dataIndex++] = 0.0;

    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = 0.5;
    vboData[dataIndex++] = 0.0;

    ////// Upload al the data to GPU. Generate VAO and VBO.
    // Todo: convert all these to "glDrawIndirect" to reduce CPU latency.
    // Basically, "glDrawIndirect" allows us to set some drawing flags without explicitly passing them from CPU side
    glGenVertexArrays(1, &vaoQuad);
    glBindVertexArray(vaoQuad);

    glGenBuffers(1, &vboQuad);

    glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
    //glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(float)*9, vboData, GL_STATIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(float)* dataSize, vboData, GL_STATIC_DRAW);

    // Edit the first input stream in the shader --- we only have one vbo or input stream here. So we enable index 0.
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0, 	                                // vertex coordinate attribute --- layout(location = 0)
        3,                	                // number of elements per vertex, here (x,y,z)
        GL_FLOAT,          	                // the type of each element
        GL_FALSE,          	                // no normalization
        0,                                  // non zero only vbo stores interleaved data like: vvvnnnttvvvnnntt...
        0                                   // offset of first element
        );

    glBindVertexArray(0);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (0 != vboData) delete[] vboData;
}

// Basically we create a cube in the form of quads
void GenerateQuadCube()
{
    ////// allocate cpu memory of firstElement & countElement
    int verticesNum = 8;
    // 3 means every vertex has 3 values (x,y,z)
    int dataSize = verticesNum * 3;
    vboData = new float[dataSize];
    int dataIndex = 0;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = -0.5;

    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = -0.5;

    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] = -0.5;

    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] = -0.5;

    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] =  0.5;

    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] =  0.5;

    vboData[dataIndex++] = -0.5;
    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] =  0.5;

    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] =  0.5;
    vboData[dataIndex++] =  0.5;

    // The indices
    //unsigned int indices[24] = 
    //{
    //    0,1,2,3,
    //    5,1,2,6,
    //    6,2,3,7,
    //    7,3,0,4,
    //    4,0,1,5,
    //    4,5,6,7
    //};
    unsigned int indices[36] =
    {
        2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1
    };

    ////// Upload al the data to GPU. Generate VAO and VBO.
    // Todo: convert all these to "glDrawIndirect" to reduce CPU latency.
    // Basically, "glDrawIndirect" allows us to set some drawing flags without explicitly passing them from CPU side
    glGenVertexArrays(1, &vaoCube);
    glBindVertexArray(vaoCube);

    // Generate indices array for glDrawElements
    glGenBuffers(1, &cubeElementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeElementbuffer);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, 8 * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 36 * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

    // VBO
    glGenBuffers(1, &vboCube);
    glBindBuffer(GL_ARRAY_BUFFER, vboCube);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(float)* dataSize, vboData, GL_STATIC_DRAW);

    // Edit the first input stream in the shader --- we only have one vbo or input stream here. So we enable index 0.
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0, 	                                // vertex coordinate attribute --- layout(location = 0)
        3,                	                // number of elements per vertex, here (x,y,z)
        GL_FLOAT,          	                // the type of each element
        GL_FALSE,          	                // no normalization
        0,                                  // non zero only vbo stores interleaved data like: vvvnnnttvvvnnntt...
        0                                   // offset of first element
        );

    glBindVertexArray(0);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (0 != vboData) delete[] vboData;
}

/************************** Initialize Some Parameters ******************************/
// Initialize some parameters for better view --- just for convenience
void sceneInitialize()
{
    eyeNear = 1.0;
    eyeFar  = 20.0;

    lightNear = 1.0;
    lightFar  = 30.0;

    objX = 2.8;
    objY = 2.2;
    objZ = 0.0;
    objS = 2.0;

    objX2 = 1.5;
    objY2 = 2.2;
    objZ2 = 0.0;
    objS2 = 0.8;

    eyePos = vec3(2.78, 2.73, -8.00);
    eyeAt  = vec3(2.78, 2.73, 2.795);
    eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
    eyeView = IGLUMatrix4x4::LookAt( eyePos, eyeAt, vec3::YAxis() );

    lightPos  = eyePos;
    lightAt   = eyeAt;
    lightView = IGLUMatrix4x4::LookAt( lightPos, lightAt, vec3::ZAxis()+vec3::XAxis() );
    lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1.0, lightNear, lightFar ); 

    // Matrices for positioning our geometry relative to the world origin
    objPosition   = IGLUMatrix4x4::Translate( objX, objY, objZ ) * IGLUMatrix4x4::Scale( objS ) * IGLUMatrix4x4::Rotate( rotateAngle1, vec3::YAxis() );
    objPosition2  = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * IGLUMatrix4x4::Scale( objS2 ) * IGLUMatrix4x4::Rotate( rotateAngle2, vec3::YAxis() );
    bigScenePos   = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );

    ///// For hair drawing
    GenerateLineSegPrimitives(jointNum, lineSegNum, 1.0);

    ///// For quad debug drawing
    GenerateQuad();

    ///// For Cube debug drawing
    GenerateQuadCube();

    //// For atomic buffer
    //unsigned int init_ac_data = 0;
    //glGenBuffers(1, &ac_buffer);
    //glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, ac_buffer);
    // In our case, we don't really need to initiate the counter value. As long as it increases and wrap when overflow, we are good.
    //glBufferData(GL_ATOMIC_COUNTER_BUFFER, sizeof(GLuint), 0, GL_DYNAMIC_DRAW);
    //glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    ///// For image store and load
    //imageCounterPtr = new IGLUTexture2D(0, 1, 1, false, IGLU_TEXTURE_DEFAULT, true);
    
    // bind texture to image units
    //imageCounterPtr->SetImageAccess(IGLU_IMAGE_READ_WRITE);
    //imageCounterPtr->BindToImageUnit(imageCounter);
    //glGenTextures(1, &imageCounter);
    //glBindTexture(GL_TEXTURE_1D, imageCounter);
    //glTexImage1D(GL_TEXTURE_1D, 0, GL_R32I, 1, 0, GL_RED_INTEGER, GL_INT, 0);
    // Bind the texture to image unit 0
    //glBindImageTexture(0, imageCounter, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32I);

    // For transform feedback
    //glGenBuffers(1, &tbo);
    //glBindBuffer(GL_ARRAY_BUFFER, tbo);
    ////glBufferData(GL_ARRAY_BUFFER, sizeof(float)*8*3, 0, GL_STATIC_READ);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(float)*12, 0, GL_STATIC_READ);
    //glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, tbo);
}


// Code that initializes our OpenGL state.  This is guaranteed (by IGLUWindow) to be 
//    called after the OpenGL context & all extensions have been initialized
void OpenGLInitialization( void ){

    // Initialize some scene parameters
    sceneInitialize();

    // Initiate timer
    gpuTimer = new IGLUGPUTimer();
    gpuTimer->Start();
    //Load obj files
    //objCow    = new IGLUOBJReader( "../../../CommonSampleFiles/models/buddha.obj", IGLU_OBJ_UNITIZE );     // lucy100k  cube_V24578_F49152 bunny
    objCow = new IGLUOBJReader("models/cow.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE);     // buddha lucy100k  cube_V24578_F49152 bunny xyzrgb_statuette_1M
    //objCube = new IGLUOBJReader("E:/Graphics/models/cube_6_planes_new.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE);        // lucy100k  cube_V24578_F49152 bunny  cube_6planes bunny_quad
    //sphere    = new IGLUOBJReader( "E:/Graphics/models/simpleSphere.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE); // lucy100k  cube_V24578_F49152 bunny 
    //boxReader = new IGLUOBJReader( "E:/Graphics/models/testBox.obj" );    // testBox_21165_40960 testBox
    //bigScene  = new IGLUOBJReader( "../../../CommonSampleFiles/Models/scene/sponza.obj", IGLU_OBJ_UNITIZE );

    //IGLUOBJMaterialReader::FinalizeMaterialsForRendering();
    // Note: If there is no texture, the following line will crash the program
    //IGLUOBJMaterialReader::s_matlTexArray->SetTextureParameters(IGLU_REPEAT_S | IGLU_REPEAT_T);

    // Create virtual trackballs for free views
    ball     = new IGLUTrackball( myWin->w(), myWin->h() );
    viewBall = new IGLUTrackball( myWin->w(), myWin->h() );

    // Create a free view
    myViewer = new uiViewInteraction(eyePos, eyeAt, eyeProj, ball, 0.6, 0.6);

    // Create a flexible perspective light source
    myLight = new uiViewInteraction(lightPos, lightAt, lightProj, viewBall, 0.6, 0.6);

    // Create a shadow map.  Make sure to use nearest neighbor on the z-buffer/shadow map.  It doesn't
    //    make sense to linearly interpolate (which is the default IGLU behavior for textures)
    //shadowMapFBO = IGLUFramebuffer::Create(GL_R8, smapRes, smapRes, true, false);
    //shadowMapFBO[IGLU_COLOR0].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );
    //shadowMapFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );

    // Load a shader that creates a shadow map
    //shaders[SM_PERSPECTIVE] = new IGLUShaderProgram( "shaders/perspectiveSM/sm_perspective.vert.glsl", "shaders/perspectiveSM/sm_perspective.frag.glsl" );
    //shaders[SM_PERSPECTIVE]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST ); 
    //shaders[SM_PERSPECTIVE]["proj"] = myLight->getEyeProj();
    //shaders[SM_PERSPECTIVE]["view"] = myLight->getEyeView();
    
    // Load a shader shade with traditional shadow map (with geometry shader)
    shaders[SHADE_WIREFRAME_GEO] = new IGLUShaderProgram("shaders/wireFrame/shade_wireFrame.vert.glsl", "shaders/wireFrame/shade_wireFrame.geom.glsl", "shaders/wireFrame/shade_wireFrame.frag.glsl");
    shaders[SHADE_WIREFRAME_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    //shaders[SHADE_WIREFRAME_GEO]["shadowTex"]   = shadowMapFBO[IGLU_DEPTH];
    //shaders[SHADE_WIREFRAME_GEO]["shadowTex2"]  = secondDepthFBO[IGLU_DEPTH];
    //shaders[SHADE_WIREFRAME_GEO]["normalTex"]   = shadowMapFBO[IGLU_COLOR0];
    shaders[SHADE_WIREFRAME_GEO]["view"]        = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_WIREFRAME_GEO]["proj"]        = myViewer->getEyeProj(); //eyeProj;
    shaders[SHADE_WIREFRAME_GEO]["winSize"]     = (float)winSize.X();
    shaders[SHADE_WIREFRAME_GEO]["lineWidth"] = 1.0 - lineWidth;

    // Load a shader shade wireframe (without geometry shader)
    //shaders[SHADE_WIREFRAME_NO_GEO] = new IGLUShaderProgram("shaders/wireFrame/shade_wireFrame_noGeo.vert.glsl", "shaders/wireFrame/shade_wireFrame_noGeo.frag.glsl");
    //shaders[SHADE_WIREFRAME_NO_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    //shaders[SHADE_WIREFRAME_NO_GEO]["shadowTex"] = shadowMapFBO[IGLU_DEPTH];
    //shaderSHADE_PERSPECTIVE_NOGEO]["shadowTex2"]  = secondDepthFBO[IGLU_DEPTH];
    //shaderSHADE_PERSPECTIVE_NOGEO]["normalTex"]   = shadowMapFBO[IGLU_COLOR0];
    //shaders[SHADE_WIREFRAME_NO_GEO]["view"]    = myViewer->getEyeView(); //eyeView;
    //shaders[SHADE_WIREFRAME_NO_GEO]["proj"]    = myViewer->getEyeProj(); //eyeProj;
    //shaders[SHADE_WIREFRAME_NO_GEO]["winSize"] = (float)winSize.X();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lineWidth"] = 1.0 - lineWidth;
    //shaders[SHADE_WIREFRAME_NO_GEO]["vertexID"] = vertexID;
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightProj"] = myLight->getEyeProj();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightView"] = myLight->getEyeView();
    //shaders[SHADE_WIREFRAME_NO_GEO]["lightNear"] = lightNear;
    //shaderSHADE_PERSPECTIVE_NOGEO]["lightFar"]    = lightFar;
    //shaders[SHADE_WIREFRAME_NO_GEO]["viewBound"] = lightNear * tan(lightFOV*float(3.141592653589793238462643) / 360.0f);
    //shaders[SHADE_WIREFRAME_NO_GEO]["matlInfoTex"] = IGLUOBJMaterialReader::s_matlCoefBuf;
    //shaders[SHADE_WIREFRAME_NO_GEO]["textureArray"] = IGLUOBJMaterialReader::s_matlTexArray; 

    // Load a shader shade wireframe with proxy (with geometry shader)
    shaders[SHADE_WIREFRAME_PROXY_GEO] = new IGLUShaderProgram("shaders/wireFrame/shade_wireFrame_proxy_Geo.vert.glsl", "shaders/wireFrame/shade_wireFrame_proxy_Geo.geom.glsl", "shaders/wireFrame/shade_wireFrame_proxy_Geo.frag.glsl");
    shaders[SHADE_WIREFRAME_PROXY_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    shaders[SHADE_WIREFRAME_PROXY_GEO]["view"] = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_WIREFRAME_PROXY_GEO]["proj"] = myViewer->getEyeProj(); //eyeProj;
    shaders[SHADE_WIREFRAME_PROXY_GEO]["cylinderRadiusSquare"] = proxyCylinderRadius*proxyCylinderRadius;
    shaders[SHADE_WIREFRAME_PROXY_GEO]["sphereRadiusSquare"] = proxySphereRadius*proxySphereRadius;

    // Load a shader shade wireframe with proxy (without geometry shader)
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO] = new IGLUShaderProgram("shaders/wireFrame/shade_wireFrame_proxy_noGeo.vert.glsl", "shaders/wireFrame/shade_wireFrame_proxy_noGeo.frag.glsl");
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["view"]      = myViewer->getEyeView(); //eyeView;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["proj"]      = myViewer->getEyeProj(); //eyeProj;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["winSize"]   = (float)winSize.X();
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["cylinderRadiusSquare"] = proxyCylinderRadius*proxyCylinderRadius;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["sphereRadiusSquare"] = proxySphereRadius*proxySphereRadius;
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["viewBound"] = eyeNear * tan(eyeFOV*float(3.141592653589793238462643) / 360.0f);
    //shaders[SHADE_WIREFRAME_PROXY_NO_GEO]["znear"]     = eyeNear;

    // Load a shader to shade tubular structure with proxy (with geometry shader)
    shaders[SHADE_TUBE_GEO] = new IGLUShaderProgram("shaders/geoProxy/shade_tube_Geo.vert.glsl", "shaders/geoProxy/shade_tube_Geo.geom.glsl", "shaders/geoProxy/shade_tube_Geo.frag.glsl");
    shaders[SHADE_TUBE_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    shaders[SHADE_TUBE_GEO]["view"]    = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_TUBE_GEO]["proj"]    = myViewer->getEyeProj(); //eyeProj;
    shaders[SHADE_TUBE_GEO]["winSize"] = (float)winSize.X();
    shaders[SHADE_TUBE_GEO]["znear"]   = eyeNear;
    shaders[SHADE_TUBE_GEO]["viewBound"] = eyeNear * tan(eyeFOV*float(3.141592653589793238462643) / 360.0f);
    shaders[SHADE_TUBE_GEO]["cylinderRadiusSquare"] = proxyCylinderRadius*proxyCylinderRadius;
    shaders[SHADE_TUBE_GEO]["sphereRadiusSquare"] = proxySphereRadius*proxySphereRadius;
    shaders[SHADE_TUBE_GEO]["radius"] = proxyCylinderRadius;

    // Load a shader shade wireframe with proxy (without geometry shader)
    //shaders[SHADE_TUBE_NO_GEO] = new IGLUShaderProgram("shaders/geoProxy/shade_tube_NoGeo.vert.glsl", "shaders/geoProxy/shade_tube_NoGeo.frag.glsl");
    //shaders[SHADE_TUBE_NO_GEO]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    //shaders[SHADE_TUBE_NO_GEO]["view"] = myViewer->getEyeView(); //eyeView;
    //shaders[SHADE_TUBE_NO_GEO]["proj"] = myViewer->getEyeProj(); //eyeProj;
    //shaders[SHADE_TUBE_NO_GEO]["winSize"]   = (float)winSize.X();
    //shaders[SHADE_TUBE_NO_GEO]["znear"]     = eyeNear;
    //shaders[SHADE_TUBE_NO_GEO]["viewBound"] = eyeNear * tan(eyeFOV*float(3.141592653589793238462643) / 360.0f);
    ////shaders[SHADE_TUBE_NO_GEO]["shapeFlag"] = 1;
    //shaders[SHADE_TUBE_NO_GEO]["cylinderRadiusSquare"] = proxyCylinderRadius;
    //shaders[SHADE_TUBE_NO_GEO]["sphereRadiusSquare"] = proxySphereRadius;

    // Load a shader to test quad
    //shaders[SHADE_QUAD] = new IGLUShaderProgram("shaders/wireFrame/shade_wireFrame_quad_noGeo.vert.glsl", "shaders/wireFrame/shade_wireFrame_quad_noGeo.frag.glsl");
    //shaders[SHADE_QUAD]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    //shaders[SHADE_QUAD]["view"] = myViewer->getEyeView(); //eyeView;
    //shaders[SHADE_QUAD]["proj"] = myViewer->getEyeProj(); //eyeProj;
    //const GLchar* feedbackVaryings[] = { "outValue" };
    //glTransformFeedbackVaryings(shaders[SHADE_QUAD]->GetProgramID(), 1, feedbackVaryings, GL_INTERLEAVED_ATTRIBS);
    //shaders[SHADE_QUAD]["sphereRadiusSquare"] = proxySphereRadius;
    //glUseProgram(shaders[SHADE_QUAD]->GetProgramID());
    // We only use image unit 0
    //glUniform1i(glGetUniformLocation(shaders[SHADE_QUAD]->GetProgramID(), "imageCounter"), 0);
    

    ///////// Set up opengl flags
    //glClearColor(0.8, 0.8, 0.8, 1.0);
    glClearColor(0.0, 0.0, 0.0, 0.0);
}


int main(int argc, char** argv)
{
	// Create our main window
	myWin = new IGLUWindow( winSize.X(), winSize.Y(), "Quadratic Mesh (tubular geomery proxy)" );
	myWin->SetWindowProperties( IGLU_WINDOW_NO_RESIZE |	
								IGLU_WINDOW_DOUBLE |
								IGLU_WINDOW_REDRAW_ON_IDLE |
                                IGLU_WINDOW_DEPTH |
								IGLU_WINDOW_W_FRAMERATE ); 
	myWin->SetDisplayCallback( display );  
	myWin->SetIdleCallback( IGLUWindow::NullIdle );
	myWin->SetPreprocessOnGLInit( OpenGLInitialization );
	myWin->SetActiveMotionCallback(Motion);
	myWin->SetMouseButtonCallback(Button);
    myWin->SetKeyboardCallback( KeyBoard );
	myWin->CreateWindow( argc, argv );

	// Create our widget window & add our widgets to it
	uiWin = new IGLUWidgetWindow( 300, 550, "UI Widget Window" );
    uiWin->AddWidget( &displayMode );
    //uiWin->AddButton( &displayGPUTime );
    //uiWin->AddWidget( &displayGPUTime );
    //uiWin->AddWidget( &printLightInfo );
	//uiWin->AddWidget( &zBias );
    //uiWin->AddWidget( &normBiasBound );
    //uiWin->AddWidget( &upperBoundRD );
    //uiWin->AddWidget( &depthBiasRatio );
    //uiWin->AddWidget( &depthBiasConst );
    //uiWin->AddWidget( &visualizeLayer );
    //uiWin->AddWidget( &vertexID);
    //uiWin->AddWidget( &hairPointRasterSize);
    //uiWin->AddWidget( &hairLineRasterSize );
    //uiWin->AddWidget( &lineWidth,      new IGLUVariableCallback(UpdateLineWidth));
    uiWin->AddWidget(&proxyCylinderRadius, new IGLUVariableCallback(UpdateProxyCylinderRadius));
    uiWin->AddWidget(&proxySphereRadius, new IGLUVariableCallback(UpdateProxySphereRadius));
	//uiWin->AddWidget( &smapRes,     new IGLUVariableCallback( UpdateShadowMapResolution ) );
	//uiWin->AddWidgetSpacer();
    //uiWin->AddWidget( &eyeNear,     new IGLUVariableCallback( UpdateEyeProjection ) );
    //uiWin->AddWidget( &eyeFar,      new IGLUVariableCallback( UpdateEyeProjection ) );
    //uiWin->AddWidget( &lightMovement, new IGLUVariableCallback( UpdateLightPosition ) );
    //uiWin->AddWidget( &lightNear,     new IGLUVariableCallback( UpdateLightProjection ) );
    //uiWin->AddWidget( &lightFar,      new IGLUVariableCallback( UpdateLightProjection ) );
    //uiWin->AddWidget( &lightFOV,      new IGLUVariableCallback( UpdateLightProjection ) );
    //uiWin->AddWidget( &lightIntense );
    //uiWin->AddWidgetSpacer();
    //uiWin->AddWidget( &omniNear );
    //uiWin->AddWidget( &omniFar );
    //uiWin->AddWidget( &omniLightX );
    //uiWin->AddWidget( &omniLightY );
    //uiWin->AddWidget( &omniLightZ );
	//uiWin->AddWidgetSpacer();
    uiWin->AddWidget( &objX,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objY,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objZ,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objS,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &rotateAngle1,  new IGLUVariableCallback( UpdateObjPosition ) );
    //uiWin->AddWidget( &objX2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidget( &objY2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidget( &objZ2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidget( &objS2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidget( &rotateAngle2,  new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidgetSpacer();
    //uiWin->AddWidget( &objX3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objY3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objZ3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objS3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &visLightView );
    //uiWin->AddWidget( &centerFlag );
    //uiWin->AddWidget( &adaptiveFlag );
    //uiWin->AddWidget( &useBiasBound );
    //uiWin->AddWidget( &realBoundFlag );
    uiWin->AddWidget( &rotateFlag );
    uiWin->AddWidget( &rotateAngle );
    //uiWin->AddWidget( &shadowFlag );
	//uiWin->AddWidget( &useLambertian );
    //uiWin->AddWidget( &constAmbient );
    //uiWin->AddWidget( &usePhong );
    //uiWin->AddWidget( &phongAlpha );
    //uiWin->AddWidget( &methodFlag );
	//uiWin->AddWidget( &useColorMask );
    uiWin->AddWidget( &showText);
	//uiWin->AddWidgetSpacer();
	uiWin->AddWidget( &displayFPS,    new IGLUVariableCallback( ToggleFramerate ) );
	myWin->SetWidgetWindow( uiWin );

    // Make the widget show by default
    uiWin->show();
	// Start running our IGLU OpenGL program!
	IGLUWindow::Run();
	return 0;
}
