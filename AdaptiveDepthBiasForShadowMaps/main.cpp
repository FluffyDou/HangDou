#include <GL/glew.h>
#include <stdio.h>

// All headers are automatically included from "iglu.h"
#include <GL/glut.h>
#include "iglu.h"

#include "header/uiViewInteraction.h"

using namespace iglu;

ivec2 winSize = ivec2(1024, 1024);

// These are handles to the two windows, should you need to modify their
//   behavior after window creation.
IGLUWindow::Ptr         myWin = 0;
IGLUWidgetWindow::Ptr   uiWin = 0;

// All the variables in the control panel
IGLUBool                printLightInfo     ( false,                                     "Print light info" );
IGLUBool                displayGPUTime     ( false,                                     "Display my timing" );
IGLUBool                displayFPS         ( true,                                      "Display framerate" );
IGLUBool				useColorMask       ( true,                                      "Enable color masking for shadow map" );
IGLUBool				useLambertian      ( true,                                      "Use Lambertian shading" );
IGLUBool				usePhong           ( false,                                     "Use Phong shading" );
IGLUBool				shadowFlag         ( true,                                      "Display shadow" );
IGLUBool				centerFlag         ( false,                                     "Show grid center" );
IGLUBool				adaptiveFlag       ( true,                                      "Use adaptive bias" );
IGLUBool				useBiasBound       ( false,                                     "Use bias bound" );
IGLUBool				realBoundFlag      ( true,                                      "Real depth value bound" );
IGLUBool				visLightView       ( false,                                     "Display light view" );
IGLUBool                rotateFlag         ( false,                                     "Rotate the object or not" );
IGLUBool                showText           ( true,                                      "Show the text or not" );

IGLUInt                 smapRes            ( 1024,   IGLURange<int>(2,8192),   1,       "Shadow map resolution" );
IGLUInt                 displayMode        ( 0,      IGLURange<int>(0,3),     1,        "Display Mode" );
IGLUInt                 methodFlag         ( 0,   IGLURange<int>(0,1),   1,             "Method ID" );

IGLUFloat               visualizeLayer     ( 0.0,    IGLURange<float>(0,1),   1.0,      "Omni shadow map layer to display" );
// This is the constant part for adaptive depth bias
IGLUFloat               constEpsilon       ( 0.003,  IGLURange<float>(0,10), 0.0001,    "Constant depth bias based on scene scale" );
// This is the simple normalized depth bias
IGLUFloat               zBias              ( -0.0001,IGLURange<float>(-0.2,0), 0.0001,  "Simple normalized constant bias" );
IGLUFloat               normBiasBound      ( 0.001,  IGLURange<float>(0,0.02), 0.0001,  "Simple normalized bias lower bound" );
IGLUFloat               eyeNear            ( 0.1,    IGLURange<float>(0,10), 0.1,       "Eye near" );
IGLUFloat               eyeFar             ( 15.0,   IGLURange<float>(10,120), 1.0,     "Eye far" );
IGLUFloat               eyeFOV             ( 38.5,   IGLURange<float>(20,150), 1,       "Viewer field-of-view" );
IGLUFloat               lightMovement      ( 3.5,    IGLURange<float>(0,4), 0.04,       "Light movement" );
IGLUFloat               lightNear          ( 0.10,   IGLURange<float>(0.01,1),  0.01,   "Light near plane" );
IGLUFloat               lightFar           ( 8.3,    IGLURange<float>(5.0,100), 0.1,     "Light far plane" );
IGLUFloat               lightFOV           ( 110,    IGLURange<float>(60,150), 1,       "Shadow map field-of-view" );
IGLUFloat               lightIntense       ( 1.0,    IGLURange<float>(0.0,15.0), 0.1,   "Light intensity" );
IGLUFloat               phongAlpha         ( 2.0,    IGLURange<float>(0.1,15.0), 0.01,  "Phong alpha" );
IGLUFloat               constAmbient       ( 0.2,    IGLURange<float>(0.0,1.0), 0.01,   "Constant ambient" );

IGLUFloat               omniNear           ( 0.10,   IGLURange<float>(0.01,2),  0.01,   "Omni near plane" );
IGLUFloat               omniFar            ( 8.3,   IGLURange<float>(2.0,70), 0.1,      "Omni far plane" );
IGLUFloat               omniLightX         ( 3.02,   IGLURange<float>(-20,20), 0.01,    "Omni x" );
IGLUFloat               omniLightY         ( 3.20,   IGLURange<float>(-20,20), 0.01,    "Omni y" ); 
IGLUFloat               omniLightZ         ( 0.19,  IGLURange<float>(-20,20), 0.01,     "Omni z" );

IGLUFloat               objX               ( 8.7,    IGLURange<float>(-30,30), 0.1,     "Dragon X" );
IGLUFloat               objY               ( -0.8,   IGLURange<float>(-30,30), 0.1,     "Dragon Y" );
IGLUFloat               objZ               ( -3.5,   IGLURange<float>(-30,30), 0.1,     "Dragon Z" );
IGLUFloat               objS               ( 1.6,    IGLURange<float>(0,10), 0.1,       "Dragon scale" );

IGLUFloat               rotateX            ( 0.0,    IGLURange<float>(0,360), 1.0,      "Dragon Rotate Amount" );

IGLUFloat               objX2               ( 3.1,   IGLURange<float>(-10,10), 0.1,     "Box X" );
IGLUFloat               objY2               ( -1.4,  IGLURange<float>(-10,10), 0.1,     "Box Y" );
IGLUFloat               objZ2               ( -6.2,  IGLURange<float>(-10,10), 0.1,     "Box Z" );
IGLUFloat               objS2               ( 1.0,   IGLURange<float>(0,10), 0.1,       "Box scale" );

IGLUFloat               rotateX2            ( 0.0,   IGLURange<float>(0,360), 1.0,      "Box Rotate Amount" );

IGLUFloat               objX3               ( 0.0,   IGLURange<float>(-10,10), 0.1,     "Object3 X" );      
IGLUFloat               objY3               ( 8.2,   IGLURange<float>(-10,10), 0.1,     "Object3 Y" );
IGLUFloat               objZ3               ( -5.0,  IGLURange<float>(-10,10), 0.1,     "Object3 Z" );
IGLUFloat               objS3               ( 30.0,  IGLURange<float>(0,30), 0.1,       "Object3 scale" );

IGLUFloat               rotateAngle         ( 0.2,   IGLURange<float>(0.1,10.0), 0.1,   "Rotate angle" );

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
IGLUOBJReader::Ptr      objCube2   = 0;
IGLUOBJReader::Ptr      sphere     = 0;
IGLUOBJReader::Ptr      boxReader  = 0;
IGLUOBJReader::Ptr      bigScene   = 0;

//Shader pool --- ATTENTION: All the "perspective" means traditional shadow map, since that is perspective projection on light space
//These are accessed using the enum below
IGLUShaderProgram::Ptr  shaders[10] = {0};
enum ShaderType{ SM_PERSPECTIVE = 0, SHADE_PERSPECTIVE, SM_OMNI, SHADE_OMNI, TEX_DISPLAY_OMNI, SM_SECONDDEPTH_PERS, SM_OMNI_SECONDDEPTH };

uiViewInteraction *myViewer = 0;
uiViewInteraction *myLight  = 0;

//Trackball for interaction
IGLUTrackball::Ptr	   ball = 0;
// Trackball for point light source
IGLUTrackball::Ptr viewBall = 0;

//The frame buffer holding our shadow map
IGLUFramebuffer::Ptr    shadowMapFBO, secondDepthFBO, parabSMapFBO, parabSecondDepthFBO;

// Location of the light and the eye
vec3            lightPos = vec3(2.78,5.25,2.795);
vec3            lightAt  = vec3(2.78, 2.73, 2.795);
vec3            eyePos   = vec3(2.78, 2.73, -8.00);
vec3            eyeAt    = vec3(2.78, 2.73, 2.795);

// Matrices for setting up the view from the eye
IGLUMatrix4x4   eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
IGLUMatrix4x4   eyeView = IGLUMatrix4x4::LookAt( eyePos, eyeAt, vec3::YAxis() );

// Matrices for setting up the view from the light (i.e., the shadow map )
IGLUMatrix4x4   lightView = IGLUMatrix4x4::LookAt( lightPos, lightAt, vec3::ZAxis()+vec3::XAxis() );
IGLUMatrix4x4   lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1.0, lightNear, lightFar ); 

// Matrices for positioning our geometry relative to the world origin
IGLUMatrix4x4   objPosition   = IGLUMatrix4x4::Translate( objX, objY, objZ ) * IGLUMatrix4x4::Scale( objS ) * 
                                IGLUMatrix4x4::Rotate( rotateX, vec3::XAxis() );

IGLUMatrix4x4   objPosition2  = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * IGLUMatrix4x4::Scale( objS2 ) * 
                                IGLUMatrix4x4::Rotate( rotateX2, vec3::XAxis() );

IGLUMatrix4x4   bigScenePos   = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );


// Generate perspective shadow map
void SMPerspecitve();
// Omni Shadow Map
void SMOmni();
// Shade with omni light source
void ShadeOmni();
// Shade with perspective light source
void shadePersPecitve();
// Update eye projection --- update when modifying eye near and eye far
void UpdateEyeProjection();

// Initialize some parameters for better view --- just for convenience
void sceneInitialize();

// Display current using method
void DisplayMethod();

/************ Render loop ************/
void display ( void )	
{
    //glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // Simple object rotation
    if( rotateFlag )
    {
        objPosition2     = objPosition2 * IGLUMatrix4x4::Rotate(rotateAngle, vec3::YAxis());
        accumulateAngle += rotateAngle;
    }

    // Render scene with traditional shadow map
    if( displayMode == 0 )
    {
        //glEnable(GL_CULL_FACE);
        //glCullFace(GL_BACK);

        if(oldDisplayMode != displayMode)
            constEpsilon = 0.0059;

        // Generate shadow map
        SMPerspecitve();

        // Shade the scene
        shadePersPecitve();

        //glDisable(GL_CULL_FACE);

    }
    else if( displayMode == 1 ) // Visualize the shadow map
    {	
        // Generate shadow map
        SMPerspecitve();

        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        IGLUDraw::Fullscreen( shadowMapFBO[IGLU_DEPTH] );
        //IGLUDraw::Fullscreen( shadowMapFBO[IGLU_COLOR0] );
    }
    else if( displayMode == 2 ) // Shade with Omni Light
    {	
        //glEnable(GL_CULL_FACE);
        //glCullFace(GL_BACK);        

        if(oldDisplayMode != displayMode)
            constEpsilon = 0.0574;

        // Create Omni shadow map
        SMOmni();

        // Shade the scene with Omni light source
        ShadeOmni();

        //glDisable(GL_CULL_FACE);
    }
    else if( displayMode == 3 ) // Visualize Omni shadow map
    {	
        SMOmni();

        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        shaders[TEX_DISPLAY_OMNI]["layer"] = visualizeLayer;
        IGLUDraw::Fullscreen( shaders[TEX_DISPLAY_OMNI], parabSMapFBO[IGLU_DEPTH], "inputTex" );
    }
    else
    {
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    }

    // Keep the initial parameter when start new display mode
    oldDisplayMode = displayMode;

    // Display the method we choose
    DisplayMethod();

    //exit(0);
}


// Generate traditional shadow map
void SMPerspecitve()
{
    //Render into the shadow map.  This can be done with or without a color mask.
    //   When no color channels are enabled for rendering, creating the shadow map
    //   is significantly faster (usually 2-4x), so color masking is usually a good
    //   idea for when rendering only to a depth buffer.
    shadowMapFBO->Bind();
    //shadowMapFBO->Clear();

    if (useColorMask)
    {
        glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
        //glDisable(GL_BLEND);
    }

    shadowMapFBO->Clear();

    // Draw our objects in the correct position
    shaders[SM_PERSPECTIVE]["model"] = objPosition * ball->GetMatrix();
    objCow->Draw( shaders[SM_PERSPECTIVE] );

    shaders[SM_PERSPECTIVE]["model"] = objPosition2 * ball->GetMatrix();
    objCube2->Draw( shaders[SM_PERSPECTIVE] );

    shaders[SM_PERSPECTIVE]["model"] = IGLUMatrix4x4::Identity();
    boxReader->Draw( shaders[SM_PERSPECTIVE] );

    //shaders[SM_PERSPECTIVE]["model"] = bigScenePos;
    //bigScene->Draw( shaders[SM_PERSPECTIVE] );

    // Finish rendering our shadow map
    if (useColorMask) 
    {
        glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
        //glEnable(GL_BLEND);
    }

    shadowMapFBO->Unbind();
}


// Omni Shadow Map
void SMOmni()
{
    parabSMapFBO->Bind();

    if (useColorMask)
    {
        glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
        //glDisable(GL_BLEND);
    }

    parabSMapFBO->Clear();
    
    shaders[SM_OMNI]["nearDist"]   = omniNear;
    shaders[SM_OMNI]["farDist"]    = omniFar;
    shaders[SM_OMNI]["wsLightPos"] = vec4(omniLightX, omniLightY, omniLightZ, 1.0);

    shaders[SM_OMNI]["model"] = objPosition * ball->GetMatrix();
    objCow->Draw( shaders[SM_OMNI] );

    shaders[SM_OMNI]["model"] = objPosition2 * ball->GetMatrix();
    objCube->Draw( shaders[SM_OMNI] ); 

    //shaders[SM_OMNI]["model"] = bigScenePos;
    //bigScene->Draw( shaders[SM_OMNI] );

    // Finish rendering our shadow map
    if (useColorMask) 
    {
        glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
        //glEnable(GL_BLEND);
    }

    parabSMapFBO->Unbind();
}


// Shade with omni light source
void ShadeOmni()
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    vec4 esLightPos = myViewer->getEyeView() * vec4( omniLightX, omniLightY, omniLightZ, 1.0f );

    shaders[SHADE_OMNI]["useLambertian"] = useLambertian ? 1.0f : 0.0f;
    shaders[SHADE_OMNI]["usePhong"]      = usePhong ? 1.0f : 0.0f;
    shaders[SHADE_OMNI]["adaptiveFlag"]  = adaptiveFlag;

    //shaders[SHADE_OMNI]["methodFlag"]    = methodFlag;
    shaders[SHADE_OMNI]["phongAlpha"]    = phongAlpha;
    shaders[SHADE_OMNI]["lightIntense"]  = lightIntense;
    shaders[SHADE_OMNI]["constAmbient"]  = constAmbient;
    //shaders[SHADE_OMNI]["shadowFlag"]    = shadowFlag;
    //shaders[SHADE_OMNI]["normBiasBound"] = normBiasBound;
    //shaders[SHADE_OMNI]["realBoundFlag"] = realBoundFlag;
    //shaders[SHADE_OMNI]["upperBoundRD"] = upperBoundRD;
    shaders[SHADE_OMNI]["constantBias"] = constEpsilon;
    shaders[SHADE_OMNI]["smBufferRes"]  = float(smapRes);
    //shaders[SHADE_OMNI]["centerFlag"]   = centerFlag;
    shaders[SHADE_OMNI]["zBias"]        = zBias;
    shaders[SHADE_OMNI]["nearDist"]     = omniNear;
    shaders[SHADE_OMNI]["farDist"]      = omniFar;
    shaders[SHADE_OMNI]["wsLightPos"]   = vec4(omniLightX, omniLightY, omniLightZ, 1.0);
    shaders[SHADE_OMNI]["esLightPos"]   = esLightPos;

    shaders[SHADE_OMNI]["view"]         = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_OMNI]["proj"]         = myViewer->getEyeProj(); //eyeProj;

    // Visualize the light
    shaders[SHADE_OMNI]["model"] = IGLUMatrix4x4::Translate(vec3(omniLightX, omniLightY, omniLightZ)) * IGLUMatrix4x4::Scale(0.1);
    sphere->Draw( shaders[SHADE_OMNI] );
    // Draw object 1
    shaders[SHADE_OMNI]["model"] = objPosition * ball->GetMatrix();
    objCow->Draw( shaders[SHADE_OMNI] );
    // Draw object 2
    shaders[SHADE_OMNI]["model"] = objPosition2 * ball->GetMatrix();
    objCube->Draw( shaders[SHADE_OMNI] );

    // Draw the box
    shaders[SHADE_OMNI]["model"] = IGLUMatrix4x4::Identity();
    boxReader->Draw( shaders[SHADE_OMNI] );

    // Draw the big scene
    //shaders[SHADE_OMNI]["model"] = bigScenePos;
    //bigScene->Draw( shaders[SHADE_OMNI] );


}


// Shade with traditional shadow map
void shadePersPecitve()
{

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // Compute the eye-space position of the light this frame
    vec4 esLightPos = myViewer->getEyeView() * vec4( myLight->getEyePos(), 1.0f );

    // Setup some values constant over all the geometry in the frame
    shaders[SHADE_PERSPECTIVE]["useLambertian"] = useLambertian ? 1.0f : 0.0f;
    shaders[SHADE_PERSPECTIVE]["usePhong"]      = usePhong ? 1.0f : 0.0f;
    shaders[SHADE_PERSPECTIVE]["esLightPos"]    = esLightPos;

    //shaders[SHADE_PERSPECTIVE]["viewFrus"]    = eyeNear * tan( eyeFOV*float( 3.141592653589793238462643 )/360.0f );
    shaders[SHADE_PERSPECTIVE]["viewBound"]   = lightNear * tan( lightFOV*float( 3.141592653589793238462643 )/360.0f );
    //shaders[SHADE_PERSPECTIVE]["methodFlag"]    = methodFlag;
    shaders[SHADE_PERSPECTIVE]["constAmbient"]  = constAmbient;
    shaders[SHADE_PERSPECTIVE]["phongAlpha"]    = phongAlpha;
    shaders[SHADE_PERSPECTIVE]["smBufferRes"]   = float(smapRes);
    shaders[SHADE_PERSPECTIVE]["lightView"]     = myLight->getEyeView();
    shaders[SHADE_PERSPECTIVE]["lightProj"]     = myLight->getEyeProj();
    shaders[SHADE_PERSPECTIVE]["zBias"]         = zBias;
    //shaders[SHADE_PERSPECTIVE]["upperBoundRD"]  = upperBoundRD;
    shaders[SHADE_PERSPECTIVE]["constantBias"]  = constEpsilon;
    //shaders[SHADE_PERSPECTIVE]["shadowFlag"]    = shadowFlag;
    //shaders[SHADE_PERSPECTIVE]["centerFlag"]    = centerFlag;
    shaders[SHADE_PERSPECTIVE]["adaptiveFlag"]  = adaptiveFlag;
    //shaders[SHADE_PERSPECTIVE]["useBiasBound"]  = useBiasBound;
    shaders[SHADE_PERSPECTIVE]["lightIntense"]  = lightIntense;
    //shaders[SHADE_PERSPECTIVE]["normBiasBound"] = normBiasBound;
    //shaders[SHADE_PERSPECTIVE]["realBoundFlag"] = realBoundFlag;

    // Change between light view and camera view --- for debug
    if( visLightView )
    {
        shaders[SHADE_PERSPECTIVE]["view"]        = myLight->getEyeView(); //eyeView;
        shaders[SHADE_PERSPECTIVE]["proj"]        = myLight->getEyeProj(); //eyeProj;
    }
    else
    {
        shaders[SHADE_PERSPECTIVE]["view"]        = myViewer->getEyeView(); //eyeView;
        shaders[SHADE_PERSPECTIVE]["proj"]        = myViewer->getEyeProj(); //eyeProj;
    }

    // Draw the geometry, modifying per-geometry shader variables as needed
    //shaders[SHADE_PERSPECTIVE]->Enable();
    shaders[SHADE_PERSPECTIVE]["model"] = objPosition * ball->GetMatrix();
    objCow->Draw( shaders[SHADE_PERSPECTIVE] ); 
    shaders[SHADE_PERSPECTIVE]["model"] = objPosition2 * ball->GetMatrix();
    objCube2->Draw( shaders[SHADE_PERSPECTIVE] );

    // Visualize the light
    shaders[SHADE_PERSPECTIVE]["model"] = IGLUMatrix4x4::Translate(myLight->getEyePos()) * IGLUMatrix4x4::Scale(0.1);
    sphere->Draw( shaders[SHADE_PERSPECTIVE] );

    shaders[SHADE_PERSPECTIVE]["model"] = IGLUMatrix4x4::Identity();
    boxReader->Draw( shaders[SHADE_PERSPECTIVE] );

    //shaders[SHADE_PERSPECTIVE]["model"] = bigScenePos;
    //bigScene->Draw( shaders[SHADE_PERSPECTIVE] );

    //shaders[SHADE_PERSPECTIVE]->Disable();

}

/************************** Call Back Functions ******************************/

void updateLightPos( IGLUMatrix4x4 light_view )
{
    shaders[SM_PERSPECTIVE]["view"]      = light_view;
    shaders[SM_SECONDDEPTH_PERS]["view"] = light_view;

    shaders[SHADE_PERSPECTIVE]["lightView"] = light_view;    
}

void updateEyePos( IGLUMatrix4x4 eye_view )
{
    shaders[SHADE_PERSPECTIVE]["view"] = eye_view;
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
        myLight->getEyePos().Print();
        myLight->getEyeAt().Print();
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
    if(0 == displayMode)
    {
        if( adaptiveFlag == true )
            sprintf( buf, "Traditional Shadow map (Adaptive Depth Bias)" );
        else if( adaptiveFlag == false )
            sprintf( buf, "Traditional Shadow map (Constant Depth Bias)" );
    }
    else if(1 == displayMode)
    {
        sprintf( buf, "Traditional Shadow map" );
    }
    else if(2 == displayMode)
    {
        if( adaptiveFlag == true )
            sprintf( buf, "Omni Shadow map (Adaptive Depth Bias)" );
        else if( adaptiveFlag == false )
            sprintf( buf, "Omni Shadow map (Constant Depth Bias)" );
    }
    else if(3 == displayMode)
    {
        sprintf( buf, "Omni Shadow map" );
    }

    if(showText)
        IGLUDraw::DrawText( IGLU_FONT_FIXED, 0, winSize.Y()-fontHeight, buf, 1.0 );
}

// When the user resizes the shadow map with the slider, we need to actually resize
//    the frame buffer object, so this callback function needs to be called
void UpdateShadowMapResolution( IGLUVariable *uiInt )
{
	int size = ((IGLUInt *)uiInt)->GetValue();

	shadowMapFBO->Resize( size, size );
    secondDepthFBO->Resize( size, size );

    parabSMapFBO->Resize(size, size);
    parabSecondDepthFBO->Resize(size, size);
}

void UpdateEyeProjection( IGLUVariable *uiFloat )
{

    eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
    myViewer->setEyeProj(eyeProj);

    shaders[SHADE_PERSPECTIVE]["proj"]  = myViewer->getEyeProj();
    shaders[SHADE_OMNI]["proj"]         = myViewer->getEyeProj();
    
}


// Move the light via a slider
void UpdateLightPosition( IGLUVariable *uiFloat )
{
	float pos = ((IGLUFloat *)uiFloat)->GetValue();
	if (pos >=0 && pos<1)
		lightPos = vec3(5.560,(1.0f-pos)*5.35+0.1,2.795);
	else if (pos >= 1 && pos<2)
		lightPos = vec3((2.0-pos)*5.460+0.1, 0.1, 2.795);
	else if (pos >= 2 && pos<3)
		lightPos = vec3( 0.1, (pos-2.0)*5.35+0.1, 2.795);
	else if (pos >= 3 && pos<4)
		lightPos = vec3( (pos-3.0)*5.460+0.1, 5.45, 2.795);

	//lightView = IGLUMatrix4x4::LookAt( lightPos, vec3(2.78, 2.73, 2.795), vec3::ZAxis()+vec3::XAxis() );
	shaders[SM_PERSPECTIVE]["view"]         = myLight->getEyeView();
    shaders[SHADE_PERSPECTIVE]["lightView"] = myLight->getEyeView();
}

void UpdateLightProjection( IGLUVariable *uiFloat )
{
	lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1, lightNear, lightFar );
    myLight->setEyeProj(lightProj);
	shaders[SM_PERSPECTIVE]["proj"]          = lightProj;
    shaders[SM_SECONDDEPTH_PERS]["proj"]     = lightProj;

    shaders[SHADE_PERSPECTIVE]["lightProj"]  = lightProj;
    //shaders[SHADE_PERSPECTIVE]["lightFar"]   = lightFar;
    shaders[SHADE_PERSPECTIVE]["lightNear"]  = lightNear;
}

void UpdateObjPosition( IGLUVariable *uiFloat )
{
    objPosition = IGLUMatrix4x4::Translate( objX, objY, objZ ) * 
                  IGLUMatrix4x4::Scale( objS ) * 
                  IGLUMatrix4x4::Rotate( rotateX, vec3::XAxis() );
}

void UpdateObjPosition2( IGLUVariable *uiFloat )
{
    objPosition2 = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * 
                   IGLUMatrix4x4::Scale( objS2 ) * 
                   IGLUMatrix4x4::Rotate( rotateX2, vec3::XAxis() );
}

void UpdateObjPosition3( IGLUVariable *uiFloat )
{
    bigScenePos  = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );
}


/************************** Initialize Some Parameters ******************************/
// Initialize some parameters for better view --- just for convenience
void sceneInitialize()
{
    constEpsilon       = 0.006;   // 0.0001 0.0059

    eyeNear = 1.0;
    eyeFar  = 20.0;    // 85.0   15.0

    lightNear = 1.0;   // 1.0
    lightFar  = 20.0;  // 70.0 20.0

    omniFar    = 25.5;    // 25.5   8.3
    omniLightX = 3.02;    // 5.85   2.78
    omniLightY = 2.7;     // -1.7   3.77
    omniLightZ = 3.21;    // -6.32  2.79

    objX = 1.6;    // 8.7    1.6
    objY = 1.1;    // -0.8   1.1
    objZ = 2.79;   // -3.5  2.795
    objS = 1.0;    // 1.6   1.0

    objX2 = 4.3;   // 3.1     4.3
    objY2 = 1.0;   // -1.4    1.0
    objZ2 = 2.79;  // -6.2    2.795
    objS2 = 1.0;   // 1.0     1.0

    eyeProj = IGLUMatrix4x4::Perspective( eyeFOV, 1.0, eyeNear, eyeFar );
    eyeView = IGLUMatrix4x4::LookAt( eyePos, eyeAt, vec3::YAxis() );

    lightPos = vec3(2.78,5.25,2.795);   // vec3( -1.9402, 5.5244, -5.2119)  vec3(2.78,5.25,2.795);
    lightAt  = vec3(2.78, 2.73, 2.795); // vec3( -1.2378, 4.8229, -5.0918)  vec3(2.78, 2.73, 2.795);

    //lightView = IGLUMatrix4x4::LookAt( lightPos, lightAt, vec3::ZAxis()+vec3::XAxis() );
    //lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1.0, lightNear, lightFar ); 

    // Matrices for positioning our geometry relative to the world origin
    objPosition   = IGLUMatrix4x4::Translate( objX, objY, objZ ) * IGLUMatrix4x4::Scale( objS ) * IGLUMatrix4x4::Rotate( rotateX, vec3::XAxis() );
    objPosition2  = IGLUMatrix4x4::Translate( objX2, objY2, objZ2 ) * IGLUMatrix4x4::Scale( objS2 ) * IGLUMatrix4x4::Rotate( rotateX2, vec3::XAxis() );
    bigScenePos   = IGLUMatrix4x4::Translate( objX3, objY3, objZ3 ) * IGLUMatrix4x4::Scale( objS3 );
}


// Code that initializes our OpenGL state.  This is guaranteed (by IGLUWindow) to be 
//    called after the OpenGL context & all extensions have been initialized
void OpenGLInitialization( void ){

    // Initialize some scene parameters
    sceneInitialize();

    // Initiate timer
    gpuTimer = new IGLUGPUTimer();

    //Load obj files
    objCow       = new IGLUOBJReader( "./models/dragon.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE );     // lucy100k  cube_V24578_F49152 bunny
    //objCow       = new IGLUOBJReader( "../../../../CommonSampleFiles/Models/princeton/hippo.obj", IGLU_OBJ_UNITIZE );     // lucy100k  cube_V24578_F49152 bunny armadillo dragon
    objCube      = new IGLUOBJReader( "./models/cube_V24578_F49152.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE);        // lucy100k  cube_V24578_F49152 bunny  cube_6planes
    objCube2     = new IGLUOBJReader( "./models/cube_6planes.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE);        // lucy100k  cube_V24578_F49152 bunny  cube_6planes
    //objCube      = new IGLUOBJReader( "../../../../CommonSampleFiles/Models/screen/screen.obj", IGLU_OBJ_UNITIZE );        // lucy100k  cube_V24578_F49152 bunny 
    //objCube      = new IGLUOBJReader( "../../../../CommonSampleFiles/Models/dragon.obj", IGLU_OBJ_UNITIZE );        // lucy100k  cube_V24578_F49152 bunny 
    sphere       = new IGLUOBJReader( "./models/simpleSphere.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE); // lucy100k  cube_V24578_F49152 bunny 
    boxReader    = new IGLUOBJReader( "./models/testBox.obj" );    // testBox_21165_40960 testBox
    //bigScene     = new IGLUOBJReader( "../../../CommonSampleFiles/Models/scene/sponza.obj", IGLU_OBJ_UNITIZE );
    //bigScene     = new IGLUOBJReader( "../../../../CommonSampleFiles/Models/Interior/Interior_701.obj", IGLU_OBJ_UNITIZE );
    IGLUOBJMaterialReader::FinalizeMaterialsForRendering();
    //IGLUOBJMaterialReader::s_matlTexArray->SetTextureParameters(IGLU_REPEAT_S | IGLU_REPEAT_T);
    //std::cout<<"r is: "<<eyeNear * tan( eyeFOV*float( 3.141592653589793238462643 )/360.0f )<<std::endl;

    // Create a virtual trackball
    ball = new IGLUTrackball( myWin->w(), myWin->h() );
    viewBall = new IGLUTrackball( myWin->w(), myWin->h() );

    // Create a free view
    myViewer = new uiViewInteraction(eyePos, eyeAt, eyeProj, ball, 0.6, 0.6);

    // Create a flexible perspective light source
    myLight = new uiViewInteraction(lightPos, lightAt, lightProj, viewBall, 0.6, 0.6);

    // Create a shadow map.  Make sure to use nearest neighbor on the z-buffer/shadow map.  It doesn't
    //    make sense to linearly interpolate (which is the default IGLU behavior for textures)
    shadowMapFBO = IGLUFramebuffer::Create(GL_R8, smapRes, smapRes, true, false);
    shadowMapFBO[IGLU_COLOR0].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );
    shadowMapFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );

    // The second nearest depth buffer
    secondDepthFBO = IGLUFramebuffer::Create(GL_R8, smapRes, smapRes, true, false);
    secondDepthFBO[IGLU_COLOR0].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );
    secondDepthFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );

    // FBO for omni shadow map
    parabSMapFBO = IGLUFramebuffer::CreateArray(GL_R8, smapRes, smapRes, 2, true, false);
    parabSMapFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );

    // FBO for omni second depth
    parabSecondDepthFBO = IGLUFramebuffer::CreateArray(GL_R8, smapRes, smapRes, 2, true, false);
    parabSecondDepthFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );

    // Load in a shader to visualize the texture array --- omni shadow map
    // Load a shader for display omni-directional shadow map
    shaders[TEX_DISPLAY_OMNI] = new IGLUShaderProgram( "shaders/paraboloidSM/displayShadowArray.vert.glsl", "shaders/paraboloidSM/displayShadowArray.frag.glsl" );

    // Load in a shader to create a omni shadow map
    shaders[SM_OMNI] = new IGLUShaderProgram( "shaders/paraboloidSM/paraboloidShadowMapTransform.vert.glsl", 
                                              "shaders/paraboloidSM/paraboloidShadowMapTransform.geom.glsl", 
                                              "shaders/paraboloidSM/paraboloidShadowMapTransform.frag.glsl" );

    shaders[SM_OMNI]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST ); 

    // Load in a shader to create omni second depth 
    shaders[SM_OMNI_SECONDDEPTH] = new IGLUShaderProgram( "shaders/paraboloidSM/paraboloidSecondDepth.vert.glsl", 
                                                          "shaders/paraboloidSM/paraboloidSecondDepth.geom.glsl", 
                                                          "shaders/paraboloidSM/paraboloidSecondDepth.frag.glsl" );

    shaders[SM_OMNI_SECONDDEPTH]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST ); 
    shaders[SM_OMNI_SECONDDEPTH]["firstDepth"] = parabSMapFBO[IGLU_DEPTH];

    // Load a shader that creates a shadow map
    shaders[SM_PERSPECTIVE] = new IGLUShaderProgram( "shaders/perspectiveSM/sm_perspective.vert.glsl", "shaders/perspectiveSM/sm_perspective.frag.glsl" );
    shaders[SM_PERSPECTIVE]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST ); 
    shaders[SM_PERSPECTIVE]["proj"] = myLight->getEyeProj();
    shaders[SM_PERSPECTIVE]["view"] = myLight->getEyeView();

    // Load in the texture to generate second nearest depth
    shaders[SM_SECONDDEPTH_PERS] = new IGLUShaderProgram( "shaders/perspectiveSM/sm_perspective_secondDepth.vert.glsl", "shaders/perspectiveSM/sm_perspective_secondDepth.frag.glsl" );
    shaders[SM_SECONDDEPTH_PERS]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST );
    shaders[SM_SECONDDEPTH_PERS]["firstDepth"] = shadowMapFBO[IGLU_DEPTH];
    shaders[SM_SECONDDEPTH_PERS]["proj"]       = myLight->getEyeProj();
    shaders[SM_SECONDDEPTH_PERS]["view"]       = myLight->getEyeView();
    

    // Load a shader that queries a shadow map
    shaders[SHADE_PERSPECTIVE] = new IGLUShaderProgram("shaders/perspectiveSM/shade_perspective.vert.glsl", "shaders/perspectiveSM/shade_perspective.frag.glsl");
    shaders[SHADE_PERSPECTIVE]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    shaders[SHADE_PERSPECTIVE]["shadowTex"]   = shadowMapFBO[IGLU_DEPTH];
    //shaders[SHADE_PERSPECTIVE]["shadowTex2"]  = secondDepthFBO[IGLU_DEPTH];
    //shaders[SHADE_PERSPECTIVE]["normalTex"]   = shadowMapFBO[IGLU_COLOR0];
    shaders[SHADE_PERSPECTIVE]["view"]        = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_PERSPECTIVE]["proj"]        = myViewer->getEyeProj(); //eyeProj;
    shaders[SHADE_PERSPECTIVE]["lightProj"]   = myLight->getEyeProj();
    shaders[SHADE_PERSPECTIVE]["lightView"]   = myLight->getEyeView();
    shaders[SHADE_PERSPECTIVE]["lightNear"]   = lightNear;
    //shaders[SHADE_PERSPECTIVE]["lightFar"]    = lightFar;
    shaders[SHADE_PERSPECTIVE]["viewBound"]   = lightNear * tan( lightFOV*float( 3.141592653589793238462643 )/360.0f );
    shaders[SHADE_PERSPECTIVE]["matlInfoTex"] = IGLUOBJMaterialReader::s_matlCoefBuf;
    //shaders[SHADE_PERSPECTIVE]["textureArray"] = IGLUOBJMaterialReader::s_matlTexArray;

    // Load a shader to render scene with media using omni-direction light    
    shaders[SHADE_OMNI] = new IGLUShaderProgram("shaders/paraboloidSM/shade_paraboloid.vert.glsl", "shaders/paraboloidSM/shade_paraboloid.frag.glsl");
    shaders[SHADE_OMNI]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
    shaders[SHADE_OMNI]["shadowTex"]    = parabSMapFBO[IGLU_DEPTH];
    //shaders[SHADE_OMNI]["secondDepth"]  = parabSecondDepthFBO[IGLU_DEPTH];
    shaders[SHADE_OMNI]["view"]         = myViewer->getEyeView(); //eyeView;
    shaders[SHADE_OMNI]["proj"]         = myViewer->getEyeProj(); //eyeProj;
    shaders[SHADE_OMNI]["matlInfoTex"]  = IGLUOBJMaterialReader::s_matlCoefBuf;
    //shaders[SHADE_OMNI]["textureArray"] = IGLUOBJMaterialReader::s_matlTexArray;    
}


int main(int argc, char** argv)
{
	// Create our main window
	myWin = new IGLUWindow( winSize.X(), winSize.Y(), "Simple IGLU Example:  Create and Use a Shadow Map" );
	myWin->SetWindowProperties( IGLU_WINDOW_NO_RESIZE |	
								IGLU_WINDOW_DOUBLE |
								IGLU_WINDOW_REDRAW_ON_IDLE |
								IGLU_WINDOW_W_FRAMERATE ); 
	myWin->SetDisplayCallback( display );  
	myWin->SetIdleCallback( IGLUWindow::NullIdle );
	myWin->SetPreprocessOnGLInit( OpenGLInitialization );
	myWin->SetActiveMotionCallback(Motion);
	myWin->SetMouseButtonCallback(Button);
    myWin->SetKeyboardCallback( KeyBoard );
	myWin->CreateWindow( argc, argv );

	// Create our widget window & add our widgets to it
	uiWin = new IGLUWidgetWindow( 300, 700, "UI Widget Window" );
    uiWin->AddWidget( &displayMode );
    //uiWin->AddWidget( &displayGPUTime );
    //uiWin->AddWidget( &printLightInfo );
	uiWin->AddWidget( &zBias );
    //uiWin->AddWidget( &normBiasBound );
    //uiWin->AddWidget( &upperBoundRD );
    uiWin->AddWidget( &constEpsilon );
    uiWin->AddWidget( &visualizeLayer );
	uiWin->AddWidget( &smapRes,       new IGLUVariableCallback( UpdateShadowMapResolution ) );
	uiWin->AddWidgetSpacer();
    uiWin->AddWidget( &eyeNear,     new IGLUVariableCallback( UpdateEyeProjection ) );
    uiWin->AddWidget( &eyeFar,      new IGLUVariableCallback( UpdateEyeProjection ) );
    //uiWin->AddWidget( &lightMovement, new IGLUVariableCallback( UpdateLightPosition ) );
    uiWin->AddWidget( &lightNear,     new IGLUVariableCallback( UpdateLightProjection ) );
    uiWin->AddWidget( &lightFar,      new IGLUVariableCallback( UpdateLightProjection ) );
    uiWin->AddWidget( &lightFOV,      new IGLUVariableCallback( UpdateLightProjection ) );
    uiWin->AddWidget( &lightIntense );
    uiWin->AddWidgetSpacer();
    uiWin->AddWidget( &omniNear );
    uiWin->AddWidget( &omniFar );
    uiWin->AddWidget( &omniLightX );
    uiWin->AddWidget( &omniLightY );
    uiWin->AddWidget( &omniLightZ );
	uiWin->AddWidgetSpacer();
    uiWin->AddWidget( &objX,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objY,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objZ,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objS,     new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &rotateX,  new IGLUVariableCallback( UpdateObjPosition ) );
    uiWin->AddWidget( &objX2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    uiWin->AddWidget( &objY2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    uiWin->AddWidget( &objZ2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    uiWin->AddWidget( &objS2,    new IGLUVariableCallback( UpdateObjPosition2 ) );
    uiWin->AddWidget( &rotateX2,  new IGLUVariableCallback( UpdateObjPosition2 ) );
    //uiWin->AddWidgetSpacer();
    //uiWin->AddWidget( &objX3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objY3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objZ3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &objS3,    new IGLUVariableCallback( UpdateObjPosition3 ) );
    //uiWin->AddWidget( &visLightView );
    //uiWin->AddWidget( &centerFlag );
    uiWin->AddWidget( &adaptiveFlag );
    //uiWin->AddWidget( &useBiasBound );
    //uiWin->AddWidget( &realBoundFlag );
    uiWin->AddWidget( &rotateFlag );
    uiWin->AddWidget( &rotateAngle );
    //uiWin->AddWidget( &shadowFlag );
	uiWin->AddWidget( &useLambertian );
    uiWin->AddWidget( &constAmbient );
    uiWin->AddWidget( &usePhong );
    uiWin->AddWidget( &phongAlpha );
    //uiWin->AddWidget( &methodFlag );
	uiWin->AddWidget( &useColorMask );
    uiWin->AddWidget( &showText);
	uiWin->AddWidgetSpacer();
	uiWin->AddWidget( &displayFPS,    new IGLUVariableCallback( ToggleFramerate ) );
	myWin->SetWidgetWindow( uiWin );

	// Start running our IGLU OpenGL program!
	IGLUWindow::Run();
	return 0;
}
