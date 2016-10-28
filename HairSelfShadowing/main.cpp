#include <GL/glew.h>
#include <stdio.h>

// All headers are automatically included from "iglu.h"
#include <GL/glut.h>
#include "iglu.h"
using namespace iglu;

// These are handles to the two windows, should you need to modify their
//   behavior after window creation.
IGLUWindow::Ptr         myWin = 0;
IGLUWidgetWindow::Ptr   uiWin = 0;

IGLUInt                 shadowFlag         ( 0,      IGLURange<int>(0,1),  1,            "For draw scene: Occupancy SM/ Normal SM " );
IGLUFloat               fragAlpha          ( 0.999,  IGLURange<float>(0.0, 1.0),  0.001, "For draw scene: Alpha value per fragment" );
IGLUFloat               intenseBias        ( 32.0,   IGLURange<float>(0.0, 96.0), 1.0,   "For draw scene&tex: Brightness bias(denominator)" );
IGLUInt                 layerT             ( 0,      IGLURange<int>(0,2),  1,            "For draw tex: Frag num/ near surf/ far surf" );
IGLUInt                 layerS             ( 0,      IGLURange<int>(0,1),  1,            "For draw scene: FragNum trans/ alpha trans" );
IGLUBool                displayFPS         ( true,                                       "Display framerate" );
IGLUInt                 displayMode        ( 0,      IGLURange<int>(0,1),  1,            "Draw scene / draw textures" );
IGLUBool				useColorMask       ( false,                                      "Enable color masking for shadow map" );
IGLUBool				useLambertian      ( false,                                      "Use Lambertian shading" );
IGLUFloat               zBias              ( -0.005, IGLURange<float>(-0.02,0), 0.0002,  "Shadow map bias" );
IGLUInt                 smapRes            ( 1024,   IGLURange<int>(128,4096),  1,       "Shadow map resolution" );
IGLUFloat               lightMovement      ( 3.5,    IGLURange<float>(0,4), 0.04,        "Light movement" );
IGLUFloat               lightNear          ( 0.1,    IGLURange<float>(0.1,3),  0.1,      "Light near plane" );
IGLUFloat               lightFar           ( 9.3,    IGLURange<float>(5.0,10), 0.1,      "Light far plane" );
IGLUFloat               lightFOV           ( 110,    IGLURange<float>(60,150), 1,        "Shadow map field-of-view" );

//The scene geometry consists of two objects: a wavefront obj file and a square 
IGLUOBJReader::Ptr      objCow     = 0;
IGLUOBJReader::Ptr      boxReader  = 0;

//We use two shaders--one to create the shadow map, another to query the shadow map
//These are accessed using the enum below
IGLUShaderProgram::Ptr  shaders[3]     = {0, 0, 0};
IGLUShaderProgram::Ptr  transShader[2] = {0, 0};

enum ShaderType{ CREATE = 0, QUERY = 1, DISPLAY = 2 };

//Trackball for interaction
IGLUTrackball::Ptr		ball = 0;

//The frame buffer holding our shadow map
IGLUFramebuffer::Ptr    duelDepthFBO;
IGLUFramebuffer::Ptr    occupancyFBO;

// Location of the light and the eye
vec3            lightPos = vec3(2.78,5.35,2.795);
vec3            eyePos   = vec3(0.0,0.0,2.0);

// Matrices for setting up the view from the eye
IGLUMatrix4x4   eyeProj = IGLUMatrix4x4::Perspective( 38.5, 1.0, 5, 20 );
IGLUMatrix4x4   eyeView = IGLUMatrix4x4::LookAt( vec3(2.78, 2.73, -8.00), 
												vec3(2.78, 2.73, 2.795), 
												vec3::YAxis() );

// Matrices for setting up the view from the light (i.e., the shadow map
IGLUMatrix4x4   lightView = IGLUMatrix4x4::LookAt( lightPos, vec3(2.78, 2.73, 2.795), vec3::ZAxis()+vec3::XAxis() );
IGLUMatrix4x4   lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1, lightNear, lightFar ); 

// Matrices for positioning our geometry relative to the world origin
IGLUMatrix4x4   objPosition   = IGLUMatrix4x4::Translate( 2.78, 2.2, 2.795 ) * IGLUMatrix4x4::Scale( 2.0f ) * 
                                IGLUMatrix4x4::Rotate( -45, vec3::YAxis() );
IGLUMatrix4x4   planePosition = IGLUMatrix4x4::Translate( -0.45*vec3::YAxis() ) * IGLUMatrix4x4::Scale( 2.0f ) *
                                IGLUMatrix4x4::Rotate( 90, vec3::XAxis() );


void display ( void )	
{	
    /****************************************** Get the buffers we needed ******************************************/
	//Render into the shadow map.  This can be done with or without a color mask.
	//   When no color channels are enabled for rendering, creating the shadow map
	//   is significantly faster (usually 2-4x), so color masking is usually a good
	//   idea for when rendering only to a depth buffer.

    /////////////////// First step: get the near depth and the far depth
	duelDepthFBO->Bind();
	    //if (useColorMask) 
	    	//glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );

        //////////// One pass version with depth test disabled (min & max blending for color & alpha channel)
        // we use min for color channel and max for alpha channel. So we clear the frame buffer color below
        glClearColor(0, 0, 1, 0);
	    duelDepthFBO->Clear();               

		// Draw our objects in the correct position
		shaders[CREATE]->Enable();

            // set up the blending function to obtain duel depth
            glBlendEquationSeparate(GL_MIN, GL_MAX);

            // disable the depth test to get duel depth
            glDisable(GL_DEPTH_TEST);
            //glDepthFunc(GL_ALWAYS);

			shaders[CREATE]["model"] = objPosition * ball->GetMatrix();
			objCow->Draw( shaders[CREATE] ); 

			//shaders[CREATE]["model"] = IGLUMatrix4x4::Identity();
			//boxReader->Draw( shaders[CREATE] );
        
            // enable the depth test back
            glEnable(GL_DEPTH_TEST);
            //glDepthFunc(GL_LESS);
		shaders[CREATE]->Disable();

        //////////// Two pass version with depth test enabled (min & max blending for color & alpha channel)
        // to be implemented

        // Finish rendering our shadow map
        //if (useColorMask) 
            //glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
    duelDepthFBO->Unbind();

    /////////////////// Second step: record the transmittance function
    occupancyFBO->Bind();

        glClearColor(0, 0, 0, 0);
        occupancyFBO->Clear();

        transShader[CREATE]->Enable();

            // draw the fragment occupancy into the second attached texture
            // disable the depth test & blend to get the transmittance function
            glDisable(GL_DEPTH_TEST);            
            glEnable(GL_COLOR_LOGIC_OP); // this operation will automatically disable GL_BLEND
            glLogicOp(GL_OR);

            //glActiveTexture(GL_TEXTURE0);
            //glBindTexture(GL_TEXTURE_2D, duelDepthFBO[IGLU_COLOR0].GetTextureID());
            transShader[CREATE]["model"] = objPosition * ball->GetMatrix();            
            objCow->Draw( transShader[CREATE] ); 

            //transShader[CREATE]["model"] = IGLUMatrix4x4::Identity();
            //boxReader->Draw( shaders[CREATE] );

            // enable the depth test & blend back
            glEnable(GL_DEPTH_TEST);            
            glDisable(GL_COLOR_LOGIC_OP);
        transShader[CREATE]->Disable();

    occupancyFBO->Unbind();


    /****************************************** Draw something on the screen ******************************************/    

    // enable the blend again
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glBlendEquation(GL_FUNC_ADD);

    //glClearDepthf(0);
    glClearColor(0, 0, 0, 0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );    

	// Did the user ask to just draw the shadow map?
	//if( displayMode == 1 )
		//IGLUDraw::Fullscreen( duelDepthFBO[IGLU_DEPTH] );
    if( displayMode == 1 )
    {
        shaders[DISPLAY]["layer"]       = layerT;
        shaders[DISPLAY]["intenseBias"] = intenseBias;
        IGLUDraw::Fullscreen( shaders[DISPLAY], occupancyFBO[IGLU_COLOR], "depthTex" );
    }
	else if( displayMode == 0 )// Draw a full scene models with shadows
	{	
		// Compute the eye-space position of the light this frame
		vec4 esLightPos = eyeView * vec4( lightPos.X(), lightPos.Y(), lightPos.Z(), 1.0f );

		// Setup some values constant over all the geometry in the frame
		shaders[QUERY]["useLambertian"] = useLambertian ? 1.0f : 0.0f;
		shaders[QUERY]["esLightPos"]    = esLightPos;
        shaders[QUERY]["zBias"]         = zBias;
        shaders[QUERY]["fragAlpha"]     = fragAlpha;
        shaders[QUERY]["shadowFlag"]    = shadowFlag;
        shaders[QUERY]["intenseBias"]   = intenseBias;
        shaders[QUERY]["layer"]         = layerS;

		// Draw the geometry, modifying per-geometry shader variables as needed
		shaders[QUERY]->Enable();
            //glDepthFunc(GL_GEQUAL);
			shaders[QUERY]["model"] = objPosition * ball->GetMatrix();
			objCow->Draw( shaders[QUERY] ); 

			shaders[QUERY]["model"] = IGLUMatrix4x4::Identity();
			boxReader->Draw( shaders[QUERY] );
            //glDepthFunc(GL_LESS);
		shaders[QUERY]->Disable();
	}

    //exit(0);
}


// Code that initializes our OpenGL state.  This is guaranteed (by IGLUWindow) to be 
//    called after the OpenGL context & all extensions have been initialized
void OpenGLInitialization( void ){
	
	//Load obj files
	objCow       = new IGLUOBJReader( "models/hairball.obj", IGLU_OBJ_UNITIZE | IGLU_OBJ_COMPACT_STORAGE); // heptoroid hairball cow
	boxReader    = new IGLUOBJReader( "models/testBox.obj" );
	IGLUOBJMaterialReader::FinalizeMaterialsForRendering();

	// Create a virtual trackball
	ball         = new IGLUTrackball( myWin->w(), myWin->h() );
	
    // The texture storing duel depth ---- z_near & z_far
	duelDepthFBO = IGLUFramebuffer::Create(GL_RGBA16F, smapRes, smapRes, true, false);
	duelDepthFBO[IGLU_DEPTH].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST );
    duelDepthFBO[IGLU_COLOR].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST | IGLU_CLAMP_TO_EDGE_S | IGLU_CLAMP_TO_EDGE_T );    
    
    // The texture storing the occupancy map
    occupancyFBO = IGLUFramebuffer::Create(GL_RGBA32UI, smapRes, smapRes, false, false);
    occupancyFBO[IGLU_COLOR].SetTextureParameters( IGLU_MAG_NEAREST | IGLU_MIN_NEAREST | IGLU_CLAMP_TO_EDGE_S | IGLU_CLAMP_TO_EDGE_T ); 

	//Load a shader that creates a shadow map
	shaders[CREATE] = new IGLUShaderProgram( "shaders/duelDepthCreate.vert.glsl", "shaders/duelDepthCreate.frag.glsl" );
	shaders[CREATE]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST | IGLU_GLSL_BLEND); 
	shaders[CREATE]["proj"] = lightProj;
	shaders[CREATE]["view"] = lightView;

    //Load a shader that display the duel shadow map
    shaders[DISPLAY] = new IGLUShaderProgram( "shaders/displaySM.vert.glsl", "shaders/displaySM.frag.glsl" );
    shaders[DISPLAY]->SetProgramEnables( IGLU_GLSL_DEPTH_TEST ); 
    shaders[DISPLAY]["layer"]       = layerT;
    shaders[DISPLAY]["intenseBias"] = intenseBias;

	// Load a shader that queries a shadow map
	shaders[QUERY] = new IGLUShaderProgram("shaders/query-shadowmap.vert.glsl", "shaders/query-shadowmap.frag.glsl");
	shaders[QUERY]->SetProgramEnables(IGLU_GLSL_DEPTH_TEST);
	shaders[QUERY]["shadowTex"]   = occupancyFBO[IGLU_COLOR];
    shaders[QUERY]["view"]        = eyeView;
	shaders[QUERY]["proj"]        = eyeProj;
    shaders[QUERY]["lightView"]   = lightView;
    shaders[QUERY]["lightProj"]   = lightProj;
    shaders[QUERY]["fragAlpha"]   = fragAlpha;
    shaders[QUERY]["shadowFlag"]  = shadowFlag;
    shaders[QUERY]["zBias"]       = zBias;
    shaders[QUERY]["intenseBias"] = intenseBias;
    shaders[QUERY]["znear"]       = lightNear;
    shaders[QUERY]["zfar"]        = lightFar;
    shaders[QUERY]["layer"]       = layerS;
	shaders[QUERY]["matlInfoTex"] = IGLUOBJMaterialReader::s_matlCoefBuf;
    shaders[QUERY]["matlTexArray"] = IGLUOBJMaterialReader::s_matlTexArray;

    // Load a shader that record the opacity bit wisely ---- way to store the transmittance function
    transShader[CREATE] = new IGLUShaderProgram( "shaders/transCreate.vert.glsl", "shaders/transCreate.frag.glsl" ); 
    transShader[CREATE]->SetProgramEnables(IGLU_GLSL_BLEND); 
    transShader[CREATE]["proj"]        = lightProj;
    transShader[CREATE]["view"]        = lightView;
    transShader[CREATE]["znear"]       = lightNear;
    transShader[CREATE]["zfar"]        = lightFar;
    transShader[CREATE]["screenRes"]   = vec2(smapRes);
    transShader[CREATE]["depthRange"]  = duelDepthFBO[IGLU_COLOR];
}


// Track any updates to the trackball matrix when the mouse moves 
void Motion(int x, int y)
{
	ball->UpdateOnMotion(x, y);
}

// When the user clicks/releases in the main window, start/stop tracking with our trackball
void Button(int button, int state, int x, int y  )
{
	if(IGLU_EVENT_DOWN == state)
		ball->SetOnClick( x, y );
	else
		ball->Release();
}

// When the user clicks the "Display Framerate" checkbox, we need to update some
//    window parameters, so this callback function needs to be called
void ToggleFramerate( void *uiBool )
{
	myWin->DisplayFramerate( ((IGLUBool *)uiBool)->GetValue() );
}

// When the user resizes the shadow map with the slider, we need to actually resize
//    the framebuffer object, so this callback function needs to be called
void UpdateShadowMapResolution( void *uiInt )
{
	int size = ((IGLUInt *)uiInt)->GetValue();
	duelDepthFBO->Resize( size, size );
}

// Move the light via a slider
void UpdateLightPosition( void *uiFloat )
{
	float pos = ((IGLUFloat *)uiFloat)->GetValue();
	if (pos >=0 && pos<1)
		lightPos = vec3(5.560,                (1.0f-pos)*5.35+0.1, 2.795);
	else if (pos >= 1 && pos<2)
		lightPos = vec3((2.0-pos)*5.460+0.1,  0.1,                 2.795);
	else if (pos >= 2 && pos<3)
		lightPos = vec3( 0.1,                 (pos-2.0)*5.35+0.1,  2.795);
	else if (pos >= 3 && pos<4)
		lightPos = vec3( (pos-3.0)*5.460+0.1, 5.35,                2.795);

	lightView = IGLUMatrix4x4::LookAt( lightPos, vec3(2.78, 2.73, 2.795), vec3::ZAxis()+vec3::XAxis() );
	shaders[CREATE]["view"]     = lightView;    
    transShader[CREATE]["view"] = lightView;
    shaders[QUERY]["lightView"] = lightView;
}

void UpdateLightProjection( void *uiFloat )
{
	lightProj = IGLUMatrix4x4::Perspective( lightFOV, 1, lightNear, lightFar ); 
	shaders[CREATE]["proj"]     = lightProj;
    transShader[CREATE]["proj"] = lightProj;
    shaders[QUERY]["lightProj"] = lightProj;
}


int main(int argc, char** argv)
{
	// Create our main window
	myWin = new IGLUWindow( 768, 768, "Simple IGLU Example:  Create and Use a Shadow Map" );
	myWin->SetWindowProperties( IGLU_WINDOW_NO_RESIZE |	
								IGLU_WINDOW_DOUBLE |
								IGLU_WINDOW_REDRAW_ON_IDLE |
								IGLU_WINDOW_W_FRAMERATE ); 
	myWin->SetDisplayCallback( display );  
	myWin->SetIdleCallback( IGLUWindow::NullIdle );
	myWin->SetPreprocessOnGLInit( OpenGLInitialization );
	myWin->SetActiveMotionCallback(Motion);
	myWin->SetMouseButtonCallback(Button);
	myWin->CreateWindow( argc, argv );     
	
	// Create our widget window & add our widgets to it
	uiWin = new IGLUWidgetWindow( 350, 720, "UI Widget Window" );
	uiWin->AddWidget( &zBias );
	uiWin->AddWidget( &smapRes );
	uiWin->AddWidgetSpacer();
	uiWin->AddWidget( &lightMovement );
	uiWin->AddWidget( &lightNear );
	uiWin->AddWidget( &lightFar );
	uiWin->AddWidget( &lightFOV );        
    uiWin->AddWidgetSpacer();
    uiWin->AddWidget( &displayMode );
    uiWin->AddWidget( &shadowFlag );
    uiWin->AddWidget( &layerS );    
    uiWin->AddWidget( &fragAlpha );
    uiWin->AddWidget( &layerT );        
    uiWin->AddWidget( &intenseBias );
	uiWin->AddWidgetSpacer();
	uiWin->AddWidget( &useLambertian );
	uiWin->AddWidget( &useColorMask );	
	uiWin->AddWidgetSpacer();
	uiWin->AddWidget( &displayFPS );
	myWin->SetWidgetWindow( uiWin );
	
	// Some of our widget/variables need callbacks to affect OpenGL state
	//    when they change, so setup the callback functions here.
	displayFPS.SetCallback( ToggleFramerate );
	smapRes.SetCallback( UpdateShadowMapResolution );
	lightMovement.SetCallback( UpdateLightPosition );
	lightNear.SetCallback( UpdateLightProjection );
	lightFar.SetCallback( UpdateLightProjection );
	lightFOV.SetCallback( UpdateLightProjection );

	// Start running our IGLU OpenGL program!
	IGLUWindow::Run();
	return 0;
}
