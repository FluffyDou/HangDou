/*********************************************************************/
// IPFP.h
// This is a class implementing IPFP. Basically it takes in an initial
// solution, initial score, data number (in our matching format case) 
// and score matrix.

// ATTENTION! IPFP relies on Hungarian which is not applicable for double
// line direction case in a straight forward way.

// Hang 11/17/2014
/*********************************************************************/

#pragma once

#include <iostream>
#include "Eigen\Dense"
#include "hungarian.h"

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class IPFP
{
public:

    IPFP(void);

    IPFP(MatrixXd matrix, VectorXd result, double score, int maxStepNum, int moduleLen, int dataLen);

    ~IPFP(void);

    /*** Setters ***/
    void Reset(MatrixXd matrix, VectorXd result, double score, int maxStepNum, int moduleLen, int dataLen);

    void SetScoreMatrix(MatrixXd matrix)   { m_scoreMatrix = matrix; }

    void SetInitialResult(VectorXd result) { m_initialIndicator = result; }

    void SetInitialScore(double score)     { m_initialScore = score; }

    void SetMaxStepNum(int maxStepNum)     { m_maxStepNum = maxStepNum; }

    void SetModuleLen(int moduleLen)       { m_moduleLen = moduleLen; }

    void SetDataLen(int dataLen)           { m_dataLen = dataLen; }

    /*** Getters ***/
    VectorXd GetRefinedIndicator() { return m_refinedIndicator; }

    VectorXi GetRefinedResult() { return m_RefinedResult; }

    double GetRefinedScore()    { return m_RefinedScore; }

    /*** Methods ***/
    //
    void GenIndicatorThroughHungarian();
    //
    void DoIPFP();
    // Initialize member variables --- Also initiate matrix and vectors. ToDo.
    void Clear();

private:

    MatrixXd m_scoreMatrix;
    VectorXd m_initialIndicator; // Indicator vector (binary or continuous)
    double   m_initialScore;
    int      m_maxStepNum;
    int      m_moduleLen;
    int      m_dataLen;

    VectorXd m_refinedIndicator;
    VectorXi m_RefinedResult;   // In our format --- ei[fi]
    double   m_RefinedScore;
};

