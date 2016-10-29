/*************************************************************************
** uiViewInteraction.cpp                                                **
** --------------------------                                           **
**                                                                      **
** This is a class for view interaction                                 **
**                                                                      **
** Hang Dou (04/06/2012)                                                **
*************************************************************************/

#include "uiViewInteraction.h"

void uiViewInteraction::MouseMotion( int x, int y )
{
    m_Ball->UpdateOnMotion(x, y);

    if( 0 != m_EyeRotateFlag )
    {
        if(x - m_OldPos.X() > m_MotionEpsilon)
            rotateEyeHorizontal( -m_RotateAmount );
        else if(m_OldPos.X() - x > m_MotionEpsilon)
            rotateEyeHorizontal(  m_RotateAmount );

        if(y - m_OldPos.Y() > m_MotionEpsilon)
            rotateEyeVertical( -m_RotateAmount );
        else if(m_OldPos.Y() -y > m_MotionEpsilon)
            rotateEyeVertical(  m_RotateAmount );
    }

    m_OldPos = vec2(x,y);
    //m_Motion = m_OldPos;

}

// When the user clicks/releases in the main window, start/stop tracking with our trackball
void uiViewInteraction::MouseButton( int button, int state, int x, int y )
{
    // rotate the track ball
    if(IGLU_EVENT_DOWN == state && IGLU_EVENT_LEFT_BUTTON == button)
    {
        m_Ball->SetOnClick( x, y );
    }
    else if(IGLU_EVENT_UP == state && IGLU_EVENT_LEFT_BUTTON == button)
    {
        m_Ball->Release();
    }

    // rotate the eye view
    if(IGLU_EVENT_DOWN == state && IGLU_EVENT_RIGHT_BUTTON == button)
    {
        m_EyeRotateFlag = 1;
    }
    if(IGLU_EVENT_UP == state && IGLU_EVENT_RIGHT_BUTTON == button)
    {
        m_EyeRotateFlag = 0;
    }

    m_MouseState  = state;
    m_MouseButton = button;

    m_OldPos = vec2(x,y);
    //m_Motion = m_OldPos;
}


// Move the eye in the scene
void uiViewInteraction::updateEyePos(float forwardAndBack, float leftAndRight)
{
    vec3 eyeForward = (m_eyeAt - m_eyePos).vNormalize();

    vec3 eyeUp = abs( eyeForward.Dot(vec3::YAxis()) ) < 0.999                        ? 
                    (( eyeForward.Cross(vec3::YAxis()) ).Cross(eyeForward)).vNormalize() :
                    (( eyeForward.Cross(vec3::XAxis()) ).Cross(eyeForward)).vNormalize() ;   

    vec3 eyeRight = (eyeForward.Cross(eyeUp)).vNormalize();

    // move eye along up and forward
    vec3 moveAmout = forwardAndBack*0.2*eyeForward + leftAndRight*0.2*eyeRight;
    m_eyePos = m_eyePos + moveAmout;
    m_eyeAt  = m_eyeAt  + moveAmout;
    // update eyeView
    m_eyeView = IGLUMatrix4x4::LookAt( m_eyePos, m_eyeAt, eyeUp );

}


// rotate the view direction for the eye horizontally
void uiViewInteraction::rotateEyeHorizontal(float angle)
{

    // compute the eye basis
    vec3 eyeForward = (m_eyeAt - m_eyePos).vNormalize();

    vec3 eyeUp = abs( eyeForward.Dot(vec3::YAxis()) ) < 0.999                         ? 
                    (( eyeForward.Cross(vec3::YAxis()) ).Cross(eyeForward)).vNormalize() :
                    (( eyeForward.Cross(vec3::XAxis()) ).Cross(eyeForward)).vNormalize() ;   

    // new forward direction for eye
    eyeForward = ( IGLUMatrix4x4::Rotate(angle, eyeUp) * vec4(eyeForward, 0.0) ).xyz();

    // update look at for eye
    m_eyeAt = m_eyePos + eyeForward;

    //// New forward
    //m_eyeForward = ( ( IGLUMatrix4x4::Rotate(angle, vec3::YAxis()) * vec4(m_eyeForward, 0.0) ).xyz() ).vNormalize();
    //// New At
    //m_eyeAt      = m_eyePos + m_eyeForward;
    //// New Right
    //m_eyeRight   = ( m_eyeForward.Cross( m_eyeUp ) ).vNormalize();
    ////// New Up
    ////m_eyeUp = abs(  m_eyeForward.Dot(vec3::YAxis()) ) < 0.999                         ? 
    ////             (( m_eyeForward.Cross(vec3::YAxis()) ).Cross(m_eyeForward)).vNormalize() :
    ////             (( m_eyeForward.Cross(vec3::XAxis()) ).Cross(m_eyeForward)).vNormalize() ; 

    // update eyeView matrix
    m_eyeView = IGLUMatrix4x4::LookAt( m_eyePos, m_eyeAt, eyeUp );
}


// rotate the view direction for the eye vertically
void uiViewInteraction::rotateEyeVertical( float angle )
{
    // compute the eye basis
    vec3 eyeForward = (m_eyeAt - m_eyePos).vNormalize();

    vec3 eyeUp = abs( eyeForward.Dot(vec3::YAxis()) ) < 0.999                            ? 
                    (( eyeForward.Cross(vec3::YAxis()) ).Cross(eyeForward)).vNormalize() :
                    (( eyeForward.Cross(vec3::XAxis()) ).Cross(eyeForward)).vNormalize() ;   
    
    vec3 eyeRight = (eyeForward.Cross(eyeUp)).vNormalize();

    // new forward direction for eye
    eyeForward = ( IGLUMatrix4x4::Rotate(angle, eyeRight) * vec4(eyeForward, 0.0) ).xyz();

    // update look at for eye
    m_eyeAt = m_eyePos + eyeForward;
    // new eye up
    eyeUp = abs(  eyeForward.Dot(vec3::YAxis()) ) < 0.999                         ? 
               (( eyeForward.Cross(vec3::YAxis()) ).Cross(eyeForward)).vNormalize() :
               (( eyeForward.Cross(vec3::XAxis()) ).Cross(eyeForward)).vNormalize() ; 

    //// New forward
    //m_eyeForward = ( ( IGLUMatrix4x4::Rotate(angle, m_eyeRight) * vec4(m_eyeForward, 0.0) ).xyz() ).vNormalize();
    //// New At
    //m_eyeAt = m_eyePos + m_eyeForward;
    //// New Up
    ////m_eyeUp = abs(  m_eyeForward.Dot(vec3::YAxis()) ) < 0.999                         ? 
    ////             (( m_eyeForward.Cross(vec3::YAxis()) ).Cross(m_eyeForward)).vNormalize() :
    ////             (( m_eyeForward.Cross(vec3::XAxis()) ).Cross(m_eyeForward)).vNormalize() ; 

    //m_eyeUp = ( m_eyeRight.Cross( m_eyeForward ) ).vNormalize();

    // update eyeView matrix
    m_eyeView = IGLUMatrix4x4::LookAt( m_eyePos, m_eyeAt, eyeUp );
}


void uiViewInteraction::setEyeProj( IGLUMatrix4x4 eyeProj )
{
    m_eyeProj = eyeProj;
}


void uiViewInteraction::KeyBoard(unsigned char key, int a, int b)
{
    switch ( key )
    {
    case 'a':
    case 'A':
        updateEyePos(  0.0, -m_MoveAmount );
        break;
    case 'd':
    case 'D':
        updateEyePos(  0.0,  m_MoveAmount );
        break;
    case 'w':
    case 'W':
        updateEyePos(  m_MoveAmount,  0.0 );
        break;
    case 's':
    case 'S':
        updateEyePos( -m_MoveAmount,  0.0 );
        break;
    case 'f':
    case 'F':
        rotateEyeHorizontal(  m_RotateAmount );
        break;
    case 'h':
    case 'H':
        rotateEyeHorizontal( -m_RotateAmount );
        break;
    case 't':
    case 'T':
        rotateEyeVertical(  m_RotateAmount );
        break;
    case 'g':
    case 'G':
        rotateEyeVertical( -m_RotateAmount );
        break;
    default: break;
    }

}