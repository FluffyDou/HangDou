/*************************************************************************
** uiViewInteraction.h                                                  **
** --------------------------                                           **
**                                                                      **
** This is a class for view interaction                                 **
**                                                                      **
** Hang Dou (04/06/2012)                                                **
*************************************************************************/

#ifndef MY_UI_VIEW_INTERATION_H
#define MY_UI_VIEW_INTERATION_H

#include "iglu.h"
#include <gl/glut.h>

using namespace iglu;

class uiViewInteraction
{

private:

    IGLUTrackball::Ptr m_Ball;
    IGLUTrackball::Ptr m_Ball2; // ball2 is designed for track ball move

    int   m_EyeRotateFlag, m_MouseState, m_MouseButton;

    float m_MotionEpsilon;
    float m_RotateAmount, m_MoveAmount;
    float m_eyeNear, m_eyeFar, m_eyeFov;

    vec2  m_OldPos; // screen position for mouse interaction
    //vec2 m_Motion; // For motion tract

    vec3  m_eyePos, m_eyeAt;
    vec3 m_eyeUp, m_eyeRight, m_eyeForward;

    IGLUMatrix4x4 m_eyeView, m_eyeProj;

    void updateEyePos(float forwardAndBack, float leftAndRight);
    void rotateEyeHorizontal(float angle);
    void rotateEyeVertical(float angle);

public:

    uiViewInteraction(): m_MotionEpsilon(0.001), m_Ball(0), m_EyeRotateFlag(0), m_RotateAmount(0.4), m_MoveAmount(5.0), m_MouseState(0), m_MouseButton(0), 
                         m_OldPos(vec2(0,0))
                         {}

    uiViewInteraction( vec3 eyePos, vec3 eyeAt, IGLUMatrix4x4 eyeProj, IGLUTrackball::Ptr ball ): m_MotionEpsilon(0.001), m_Ball(0), 
        m_EyeRotateFlag(0), m_RotateAmount(0.4), m_MoveAmount(5.0), m_MouseState(0), m_MouseButton(0), m_OldPos(vec2(0,0))
    {
        m_eyeForward = (eyeAt - eyePos).vNormalize();

        m_eyeUp      = abs(  m_eyeForward.Dot(vec3::YAxis()) ) < 0.999                             ? 
                          (( m_eyeForward.Cross(vec3::YAxis()) ).Cross(m_eyeForward)).vNormalize() :
                          (( m_eyeForward.Cross(vec3::XAxis()) ).Cross(m_eyeForward)).vNormalize() ;

        m_eyeRight = (m_eyeForward.Cross(m_eyeUp)).vNormalize();

        m_eyePos   = eyePos;
        m_eyeAt    = m_eyePos + (eyeAt - m_eyePos).vNormalize();

        m_eyeView = IGLUMatrix4x4::LookAt(m_eyePos, m_eyeAt, m_eyeUp);
        m_eyeProj = eyeProj;
        m_Ball    = ball;
    }

    uiViewInteraction( vec3 eyePos, vec3 eyeAt, IGLUMatrix4x4 eyeProj, IGLUTrackball::Ptr ball, float rotateAmount, float moveAmount ): m_MotionEpsilon(0.001), m_Ball(0), 
        m_EyeRotateFlag(0), m_RotateAmount(rotateAmount), m_MoveAmount(moveAmount), m_MouseState(0), m_MouseButton(0), m_OldPos(vec2(0,0))
    {
        m_eyeForward = (eyeAt - eyePos).vNormalize();

        m_eyeUp      = abs(  m_eyeForward.Dot(vec3::YAxis()) ) < 0.999                             ? 
            (( m_eyeForward.Cross(vec3::YAxis()) ).Cross(m_eyeForward)).vNormalize() :
        (( m_eyeForward.Cross(vec3::XAxis()) ).Cross(m_eyeForward)).vNormalize() ;

        m_eyeRight = (m_eyeForward.Cross(m_eyeUp)).vNormalize();

        m_eyePos   = eyePos;
        m_eyeAt    = m_eyePos + (eyeAt - m_eyePos).vNormalize();

        m_eyeView = IGLUMatrix4x4::LookAt(m_eyePos, m_eyeAt, m_eyeUp);
        m_eyeProj = eyeProj;
        m_Ball    = ball;
    }

    ~uiViewInteraction();
    // Call back function
    void MouseMotion( int x, int y );
    void MouseButton( int button, int state, int x, int y );
    void KeyBoard(unsigned char key, int a, int b);
    // Getter
    IGLUMatrix4x4 getEyeView()        {return m_eyeView;}
    IGLUMatrix4x4 getEyeProj()        {return m_eyeProj;}
    IGLUTrackball::Ptr getTrackBall() {return m_Ball;}
    vec3 getEyePos()                  {return m_eyePos;}
    vec3 getEyeAt()                   {return m_eyeAt;}
    vec3 getForward()                 {return (m_eyeAt - m_eyePos).vNormalize();}
    // Setter
    void setEyeProj( IGLUMatrix4x4 eyeProj );

};


#endif