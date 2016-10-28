/***************************************************************************/
/* raytracer.cpp                                                           */
/* -------------                                                           */
/*                                                                         */
/* The main routine for the MC rendering ray tracing framework.            */
/*                                                                         */
/* Hang Dou (03/25/2012)                                                   */
/***************************************************************************/


/* Include all the necessary #include files from a single common #include  */
#include "raytracer.h"


/* This is where the program starts executing */
int main(int argc, char** argv)
{
    // create a frame rate record object
    fRate = new FrameRate( 15 );

    // create the scene with geometries and a camera
	myScene = SceneSetup();

    /************* Initialize OpenGL window **************/
    glutInit(&argc, argv);
    glutInitDisplayMode(  GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA  );
    glutInitWindowSize( myScene->GetWidth(), myScene->GetHeight() );
    glutInitWindowPosition( 50, 10 );
    glutCreateWindow("MC Ray Tracing");

    //////////// initialize GLEW ///////////
    GLenum err = glewInit();
    if (GLEW_OK != err) 
        printf("Error: %s\n", glewGetErrorString(err));
    printf("Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    // print out the OpenGL version
    const char* version = (const char*)glGetString(GL_VERSION);
    printf("OpenGL Version is£º %s\n\n", version);
    const GLubyte* pShaderVersion = glGetString(GL_SHADING_LANGUAGE_VERSION);
    printf("GLSL Version is£º %s\n\n", pShaderVersion);

    // create a FBO with its attached render texture
    //renderBuffer = RenderBufInit(myScene->GetWidth(), myScene->GetHeight() );

    // create a buffer on CPU to store results from ray tracing --- 4 channels: RGBA
    std::cout<<"(+) Begin to create sample vector pool.....\n";
    myImage              = new Image(myScene->GetWidth(), myScene->GetHeight(), 4);
    myImage->GetBuffer() = new float[myImage->GetXRes() * myImage->GetYRes() * myImage->GetChanNum()];
    myImage->vecPool     = new SampleVector[myScene->GetWidth()*myScene->GetHeight()*SPP];
    std::cout<<"    Create successfully!\n";

    std::cout<<"(+) Begin to create average sample vector.....\n";
    featureImage              = new Image(myScene->GetWidth(), myScene->GetHeight(), 4);
    featureImage->GetBuffer() = new float[featureImage->GetXRes() * featureImage->GetYRes()];
    featureImage->vecPool     = new SampleVector[myScene->GetWidth()*myScene->GetHeight()];
    std::cout<<"    Create successfully!\n";

    // load in the shader for render the texture
    textureShader = TexShaderSetup("glsl/texVShader.glsl", "glsl/texFShader.glsl");
    textureShader->InitAttribute( myImage->GetTexture() );
    textureShader->updateUniformTex( myImage->GetTexture() );

    /************ register OpenGL call back functions ***************/
    glutDisplayFunc( myDisplay );
    glutReshapeFunc( myReshape );
    glutIdleFunc( myIdle );
    glutKeyboardFunc( myKeys );
    //glutMouseFunc( myMouseButtonCallback );
    //glutMotionFunc( myMouseMotionCallback );

    // Open the window and start running!
    glutMainLoop();
	return 0;
}


// Display function
void myDisplay()
{
    fRate->StartFrame();

    // trace the scene and store the image into the frame buffer texture
    if( sampleCount < SPP)
    {
        if( sampleCount == SPP-1)   
            std::cout<<"Now pass "<< SPP <<" completed."<< std::endl;

        // shoot rays to trace for every pixel
        ShootRays( *myImage, myScene, sampleCount );        
        sampleCount++;

        // update the texture to draw in frame buffer
        myImage->updateTex();
    }
    else if( stageFlag1 == 0 ) // after tracing SPP times, do computation based on the sample vectors
    {
        // after tracing enough samples, we generate the filter and display features
        std::cout<<"(+) Now averaging sample vectors....\n";

        AverageSampleVectors( *featureImage, *myImage );

        std::cout<<"    Done.\n";
        stageFlag1 = 1;
        renderFlag = SAMPLECOLOR;
    }
        

    // choose which texture to render --- each texture carries a feature
    if( renderFlag == SCREENPOS )
    {
        UpdateCanvus( *featureImage, SCREENPOS );
    }
    else if( renderFlag == RANPARA )
    {
        UpdateCanvus( *featureImage, RANPARA );
    }
    else if( renderFlag == WORLDSPACECOORD )
    {
        UpdateCanvus( *featureImage, WORLDSPACECOORD );
    }
    else if( renderFlag == SURFACENORMAL )
    {
        UpdateCanvus( *featureImage, SURFACENORMAL );
    }
    else if( renderFlag == TEXVALUE )
    {
        UpdateCanvus( *featureImage, TEXVALUE );
    }
    else if( renderFlag == SAMPLECOLOR )
    {
        UpdateCanvus( *featureImage, SAMPLECOLOR );
    }

    // display the scene
    textureShader->Enable();

        glClear(GL_COLOR_BUFFER_BIT |  GL_DEPTH_BUFFER_BIT);
        if( sampleCount < SPP || renderFlag == RENDERRESULT )
        {
            textureShader->updateUniformTex( myImage->GetTexture() );
            textureShader->RenderTexture( myImage->GetTexture() );
        }
        else
        {
            textureShader->updateUniformTex( featureImage->GetTexture() );
            textureShader->RenderTexture( featureImage->GetTexture() );
        }

    textureShader->Disable();

    // Display a frame rate timer, based on the speed recorded for this frame
    DisplayTimer( fRate->EndFrame() );

    glFlush();
    glutSwapBuffers();
}


/* This function takes in an image and a scene, traces a ray through every */
/*    pixel in the scene, and stores the color in the image.               */
void ShootRays( Image &image, Scene *scene, int samplePass )
{
    // For every pixel in the image --- y is row, x is column
    for (int y=0; y < image.GetYRes(); y++)
        for (int x=0; x < image.GetXRes(); x++)
        {
            // Create a ray through this pixel --- origin, direction, init bounce time
            Ray ray( scene->camera->GetEye(), 
                            scene->camera->GenerateRay( (float)x, (float)y ), 0 );

            // before trace the ray, format the sample vector --- kinda expensive
            //scene->X->Format(scene->bounceLimit);

            /*************** record the screen position of the certain sample *************/
            scene->X->screenPos[0] = vec2( (float)x, (float)y );
            
            // trace that ray into the scene
            Color pixelColor = scene->TraceRay( ray );

            /*************** record the sample color of the certain sample *************/
            (scene->X->sampleColor)[0] = pixelColor.ColorVec();

            /***************** store the sample vector ****************************/
            image.vecPool[x + y*image.GetYRes() + image.GetXRes()*image.GetYRes()*samplePass] = *(scene->X);

            // 0+a/1, (a+b)/2, (((a+b )/2) *2 + c) / 3, ......
            float pathScale = samplePass != 0 ? 1.0/(float)(samplePass+1): 1.0 ;

            // accumulate the color and divide the result by the number of path --- store as FLOAT
            image.Buffer(x, y, 0) = (image.Buffer(x, y, 0)*(float)samplePass + pixelColor.Red())   * pathScale;
            image.Buffer(x, y, 1) = (image.Buffer(x, y, 1)*(float)samplePass + pixelColor.Green()) * pathScale;
            image.Buffer(x, y, 2) = (image.Buffer(x, y, 2)*(float)samplePass + pixelColor.Blue())  * pathScale;
            image.Buffer(x, y, 3) = 1.0;

            // tone mapping
            pixelColor = Color( image.Buffer(x, y, 0), image.Buffer(x, y, 1), image.Buffer(x, y, 2) );
            float scales = 1.0 / ( 1.0 + pixelColor.Luminance() );
            pixelColor = pixelColor * scales;
            // Gama correction
            pixelColor = Color( powf( pixelColor.Red(), 1.0f/2.2), 
                                    powf( pixelColor.Green(), 1.0f/2.2), 
                                        powf( pixelColor.Blue(), 1.0f/2.2) );

            // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
            image(x, y, 0) = unsigned char( MIN(pixelColor.Red(),   1.0) * 255 ); // R
            image(x, y, 1) = unsigned char( MIN(pixelColor.Green(), 1.0) * 255 ); // G
            image(x, y, 2) = unsigned char( MIN(pixelColor.Blue(),  1.0) * 255 ); // B
            image(x, y, 3) = 255;                                                 // A

        }        
}


// update the texture of Image for rendering 
void UpdateCanvus( Image &image, int vectorFlag )
{
    if( vectorFlag == SCREENPOS ) // render the screen position
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {
                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].screenPos[0].Y() / (float)SCREENHEIGHT, 1.0) * 255 );   // R
                image(x, y, 1) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].screenPos[0].X() / (float)SCREENWIDTH,  1.0) * 255 );   // G
                image(x, y, 2) = 0;                                                                                                              // B
                image(x, y, 3) = 255;                                                                                                            // A             
            }
        // put the render result onto GPU --- glTexture
        image.updateTex();
    }
    else if( vectorFlag == RANPARA )
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {
                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].ranPara[0].X(), 1.0) * 255 );   // R
                image(x, y, 1) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].ranPara[0].Y(), 1.0) * 255 );   // G
                image(x, y, 2) = 0;                                                                                      // B
                image(x, y, 3) = 255;                                                                                    // A             
            }
            // put the render result onto GPU --- glTexture
        image.updateTex(NEAREST);
    }
    else if( vectorFlag == WORLDSPACECOORD )
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {
                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].worldSpaceCoord[0].X() / 1000.0, 1.0) * 255 );   // R
                image(x, y, 1) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].worldSpaceCoord[0].Y() / 1000.0, 1.0) * 255 );   // G
                image(x, y, 2) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].worldSpaceCoord[0].Z() / 1000.0, 1.0) * 255 );   // B                                                                                          // B
                image(x, y, 3) = 255;                                                                                                     // A             
            }
            // put the render result onto GPU --- glTexture
            image.updateTex(LINEAR);
    }
    else if( vectorFlag == SURFACENORMAL )
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {
                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(abs(image.vecPool[x + y*image.GetXRes()].surfaceNorm[0].X()), 1.0) * 255 );   // R
                image(x, y, 1) = unsigned char( MIN(abs(image.vecPool[x + y*image.GetXRes()].surfaceNorm[0].Y()), 1.0) * 255 );   // G
                image(x, y, 2) = unsigned char( MIN(abs(image.vecPool[x + y*image.GetXRes()].surfaceNorm[0].Z()), 1.0) * 255 );   // B                                                                                          // B
                image(x, y, 3) = 255;                                                                                        // A             
            }
            // put the render result onto GPU --- glTexture
            image.updateTex(LINEAR);
    }
    else if( vectorFlag == TEXVALUE )
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {
                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].texValue[0].X(), 1.0) * 255 );   // R
                image(x, y, 1) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].texValue[0].Y(), 1.0) * 255 );   // G
                image(x, y, 2) = unsigned char( MIN(image.vecPool[x + y*image.GetXRes()].texValue[0].Z(), 1.0) * 255 );   // B                                                                                          // B
                image(x, y, 3) = 255;                                                                                     // A             
            }
            // put the render result onto GPU --- glTexture
            image.updateTex(LINEAR);
    }
    else if( vectorFlag == SAMPLECOLOR )
    {
        // For every pixel in the image
        for (int y=0; y < image.GetYRes(); y++)
            for (int x=0; x < image.GetXRes(); x++)
            {

                // tone mapping
                Color pixelColor = Color( image.vecPool[x + y*image.GetXRes()].sampleColor[0].X(), 
                                             image.vecPool[x + y*image.GetXRes()].sampleColor[0].Y(),
                                                image.vecPool[x + y*image.GetXRes()].sampleColor[0].Z() );

                float scales = 1.0 / ( 1.0 + pixelColor.Luminance() );
                pixelColor = pixelColor * scales;
                // Gama correction
                pixelColor = Color( powf( pixelColor.Red(), 1.0f/2.2), 
                                        powf( pixelColor.Green(), 1.0f/2.2), 
                                            powf( pixelColor.Blue(), 1.0f/2.2) );

                // store the color into an unsigned char array which is bound to GL render texture --- render as BYTE
                image(x, y, 0) = unsigned char( MIN(pixelColor.Red(),   1.0) * 255 ); // R
                image(x, y, 1) = unsigned char( MIN(pixelColor.Green(), 1.0) * 255 ); // G
                image(x, y, 2) = unsigned char( MIN(pixelColor.Blue(),  1.0) * 255 ); // B
                image(x, y, 3) = 255;                                                 // A             
            }
            // put the render result onto GPU --- glTexture
            image.updateTex(LINEAR);
    }
}


// average the sample vector values to render them out
void AverageSampleVectors(  Image &destination, Image &source )
{
    // For every pixel in the image
    for (int y=0; y < source.GetYRes(); y++)
        for (int x=0; x < source.GetXRes(); x++)
        {
            for(int s = 0; s < SPP; s++)
            {
                // accumulate the sample vector value
                destination.vecPool[x + y*source.GetXRes()].screenPos[0]        =   source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].screenPos[0];
                destination.vecPool[x + y*source.GetXRes()].sampleColor[0]     +=   source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].sampleColor[0];
                destination.vecPool[x + y*source.GetXRes()].texValue[0]        +=   source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].texValue[0];

                destination.vecPool[x + y*source.GetXRes()].ranPara[0]         += ( source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].ranPara[0] + 
                                                                                    source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].ranPara[1] ) * 0.5;
                destination.vecPool[x + y*source.GetXRes()].worldSpaceCoord[0] += ( source.vecPool[x + y*source.GetXRes() +
                                                                                        s*source.GetXRes()*source.GetYRes()].worldSpaceCoord[0] + 
                                                                                    source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].worldSpaceCoord[1] ) * 0.5;
                destination.vecPool[x + y*source.GetXRes()].surfaceNorm[0]     += ( source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].surfaceNorm[0] +
                                                                                    source.vecPool[x + y*source.GetXRes() + 
                                                                                        s*source.GetXRes()*source.GetYRes()].surfaceNorm[1] ) * 0.5;
            }

            // average the sample vector value
            destination.vecPool[x + y*source.GetXRes()].ranPara[0]         = destination.vecPool[x + y*source.GetXRes()].ranPara[0]         * (1.0/(float)SPP);
            destination.vecPool[x + y*source.GetXRes()].worldSpaceCoord[0] = destination.vecPool[x + y*source.GetXRes()].worldSpaceCoord[0] * (1.0/(float)SPP);
            destination.vecPool[x + y*source.GetXRes()].surfaceNorm[0]     = destination.vecPool[x + y*source.GetXRes()].surfaceNorm[0]     * (1.0/(float)SPP);
            destination.vecPool[x + y*source.GetXRes()].texValue[0]        = destination.vecPool[x + y*source.GetXRes()].texValue[0]        * (1.0/(float)SPP);
            destination.vecPool[x + y*source.GetXRes()].sampleColor[0]     = destination.vecPool[x + y*source.GetXRes()].sampleColor[0]     * (1.0/(float)SPP);

        }
}
// initialize the Render Texture for frame buffer object
//frameBuf* RenderBufInit(int screenWidth, int screenHeigth)
//{
//    std::cout << "Creating frame buffer to render in....\n";
//
//    frameBuf* renderBuf = new frameBuf(screenWidth, screenHeigth);
//
//    std::cout<< "Create successfully!\n";
//
//    return renderBuf;
//}


// initialize Shader
texShader* TexShaderSetup(char* vshader, char* fshader)
{

    std::cout << "Loading in texture shader....\n";

    texShader* shader = new texShader(vshader, fshader);

    std::cout<< "Load in successfully!\n";

    return shader;
}


/* Create a scene to render! */
Scene *SceneSetup( void )
{
    std::cout << "Setting up the scene...\n";

    // Create the main scene data structure
    Scene *scn = new Scene();

    //scn->light = new Light(vec3(4.0, 4.0, 10.0), vec3(1.0, 1.0, 1.0));
    //scn->light = new Light( vec3(343.0, 548.8, 332.0), vec3(1.0, 1.0, 1.0) );
    scn->light = new Light( vec3(343.0, 400.0, 132.0), vec3(18.4f, 15.6f, 8.0f) );

    // set ambient light for the scene
    scn->sceneAmbient = vec3(0.0, 0.0, 0.0);

    // set bounce limit for the scene --- 2 bounce for now
    scn->bounceLimit = 2;

    // Set the scene's background color
    //scn->backgroundColor = Color( 0.462f, 0.725f, 0.0f ); 
    scn->backgroundColor = Color( 0.0f, 0.0f, 0.0f );

    // Create a camera for the scene
    scn->camera = new Camera( LOOKEYE, LOOKAT, LOOKUP, FOV, SCREENHEIGHT, SCREENWIDTH );

    // Initiate the scene sample vector --- used to record all the features getting by the tracing ray
    // need to get format whenever before start tracing
    scn->X = new SampleVector();

    // Create a group to contain all the scene geometry, tell the scene about it
    Group *grp = new Group();
    scn->geometry = grp;

    // Create a material type for our geometry
    //Material *mtl = new AmbientOcclusionMaterial( 10 );
    Material *mt0        = new ConstantColorMaterial( Color(0.0, 0.0, 0.8) );
    Material *mtCenter   = new LambertianMaterial( Color(0.76, 0.75, 0.5) );
    Material *mtLeft     = new LambertianMaterial( Color(0.63, 0.06, 0.04) );
    Material *mtRight    = new LambertianMaterial( Color(0.15, 0.48, 0.09) );
    Material *mtCenterMC = new MCLambertianMaterial( Color(0.76, 0.75, 0.5) );
    Material *mtLeftMC   = new MCLambertianMaterial( Color(0.63, 0.06, 0.04) );
    Material *mtRightMC  = new MCLambertianMaterial( Color(0.15, 0.48, 0.09) );

    Material *mtLight  = new LightMaterial( Color(18.4f, 15.6f, 8.0f) );

    Material* mtGroup[5];
    mtGroup[0] = mtCenterMC;
    mtGroup[1] = mtCenterMC;
    mtGroup[2] = mtLeftMC;
    mtGroup[3] = mtRightMC;
    mtGroup[4] = mtLight;

    // Create a built in Cornell box
    createCornellBox(grp, mtGroup);

    // Read in an obj file
    readInObjFile(grp, mtGroup);


    std::cout << "Setup successfully!\n";
    return scn;
}

// Read in an obj file
void readInObjFile(Group* grp, Material** mt)
{
    GLMmodel *m = glmReadOBJ( "../models/testBox.obj" );
    glmUnitize( m );             // Scale & translate cow to [-1..1]^3
    glmFacetNormals( m );        // Compute facet normals
    glmVertexNormals( m, 90.0 ); // Compute vertex normals

    Primitive *prm = 0;

    for (unsigned int i=0; i < m->numtriangles; i++)
    {
        // Get the first triangle vertex
        int index0 = m->triangles[i].vindices[0];
        vec3 vert0( &m->vertices[3*index0] );

        // Get the second triangle vertex
        int index1 = m->triangles[i].vindices[1];
        vec3 vert1( &m->vertices[3*index1] );

        // Get the third triangle vertex
        int index2 = m->triangles[i].vindices[2];
        vec3 vert2( &m->vertices[3*index2] );

        // Do something with the three triangle vertices!
        prm = new Triangle( mt[2], 50.0*vert0 + vec3(200,150,150), 50.0*vert1 + vec3(200,150,150), 50.0*vert2 + vec3(200,150,150));
        grp->Add( prm );
    }
}


// create a Cornell box
void createCornellBox(Group* grp, Material** mt)
{
    Primitive *prm = 0;
    //prm = new Sphere( mtl, vec3( 0.0, 0.0, 0.0 ), 5.0 );
    //grp->Add( prm );
    prm = new Sphere( mt[1], vec3( 376.0, 276.0, 476.0 ), 90.0 );
    grp->Add( prm );

    prm = new Sphere( mt[0], vec3( 176.0, 176.0, 300.0 ), 60.0 );
    grp->Add( prm );

    // area light
    prm = new Triangle( mt[4], vec3( 343.0, 548.0, 227.0 ), vec3( 213.0, 548.0, 332.0 ), vec3( 343.0, 548.0, 332.0 ) );
    grp->Add( prm );

    prm = new Triangle( mt[4], vec3( 343.0, 548.0, 227.0 ), vec3( 213.0, 548.0, 227.0 ), vec3( 213.0, 548.0, 332.0 ) );
    grp->Add( prm );

    // floor
    prm = new Triangle( mt[1], vec3( 0.0, 0.0, 0.0 ), vec3( 552.8, 0.0, 0.0 ), vec3( 0.0, 0.0, 559.2 ) );
    grp->Add( prm );

    prm = new Triangle( mt[1], vec3( 0.0, 0.0, 559.2 ), vec3( 552.8, 0.0, 0.0 ), vec3( 549.6, 0.0, 559.2 ) );
    grp->Add( prm );

    // ceiling
    prm = new Triangle( mt[1], vec3( 556.0, 548.8, 559.2 ), vec3( 556.0, 548.8, 0.0 ), vec3( 0.0, 548.8, 559.2 ));
    grp->Add( prm );

    prm = new Triangle( mt[1], vec3( 0.0, 548.8, 559.2 ), vec3( 556.0, 548.8, 0.0 ), vec3( 0.0, 548.8, 0.0 ) );
    grp->Add( prm );

    // back wall
    prm = new Triangle( mt[1], vec3( 549.6, 0.0, 559.2 ), vec3( 0.0, 548.8, 559.2 ), vec3( 0.0, 0.0, 559.2 ) );
    grp->Add( prm );

    prm = new Triangle( mt[1], vec3( 549.6, 0.0, 559.2 ), vec3( 556.0, 548.8, 559.2 ), vec3( 0.0, 548.8, 559.2 ) );
    grp->Add( prm );

    // right wall
    prm = new Triangle( mt[3], vec3( 0.0, 0.0, 559.2 ), vec3( 0.0, 548.8, 0.0 ), vec3( 0.0, 0.0, 0.0 ) );
    grp->Add( prm );

    prm = new Triangle( mt[3], vec3( 0.0, 0.0, 559.2 ), vec3( 0.0, 548.8, 559.2 ), vec3( 0.0, 548.8, 0.0 ) );
    grp->Add( prm );

    // left wall
    prm = new Triangle( mt[2], vec3( 552.8, 0.0, 0.0 ), vec3( 556.0, 548.8, 559.2 ), vec3( 549.6, 0.0, 559.2 ) );
    grp->Add( prm );

    prm = new Triangle( mt[2], vec3( 552.8, 0.0, 0.0 ), vec3( 556.0, 548.8, 0.0 ), vec3( 556.0, 548.8, 559.2 ) );
    grp->Add( prm );
}

/***************GL other Call Back Functions***************************/


// reshape function
void myReshape(int w, int h)
{
    //myScene->camera->SetScreenRatio(w, h);

    glViewport( 0, 0, (GLsizei)w, (GLsizei)h );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( 0, (GLsizei)w, 0, (GLsizei)h );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glutPostRedisplay();
}


void myKeys(unsigned char key, int x, int y)
{

    switch ( key )
    {
        case 'W':
        case 'w':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3( 0.0, 10.0, 0.0) );
            break;
        case 'S':
        case 's':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3( 0.0, -10.0, 0.0) );
            break;
        case 'A':
        case 'a':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3( 10.0, 0.0, 0.0) );
            break;
        case 'D':
        case 'd':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3(-10.0, 0.0, 0.0) );
            break;
        case 'Z':
        case 'z':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3( 0.0, 0.0, 10.0) );
            break;
        case 'X':
        case 'x':
            myScene->light->SetPosition( myScene->light->GetPosition() + vec3( 0.0, 0.0, -10.0) );
            break;
        case '1':
            renderFlag = SCREENPOS;
            break;
        case '2':
            renderFlag = RANPARA;
            break;
        case '3':
            renderFlag = WORLDSPACECOORD;
            break;
        case '4':
            renderFlag = SURFACENORMAL;
            break;
        case '5':
            renderFlag = TEXVALUE;
            break;
        case '6':
            renderFlag = SAMPLECOLOR;
            break;
        case '7':
            renderFlag = RENDERRESULT;
            break;
        case 'Q':
        case 'q':
        case  27:
            exit(0);
            break;
        default:
            break;
    }

}


// Idle function
void myIdle()
{
    glutPostRedisplay();
};


// Displays a timer on the lower left corner of the screen
void DisplayTimer( float fps )
{
    char buf[1024];
    glColor3f(0.9, 0.9, 0.9);
    glRasterPos2i(3,10);
    sprintf( buf, "%.2f fps", fps );
    int len = (int) strlen(buf);
    for(int i=0; i<len; i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buf[i]);
}

