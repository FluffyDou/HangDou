/*********************************************************************/
// ScoreMatrix2DPointCloud.h
// This the class holding the score matrix for 2D point cloud.
// This class is very similar to that for single line direction. The
// only difference is we only care about midpoint distance now.

// Hang 11/19/2014
/*********************************************************************/

#pragma once

//#include <omp.h>
#include <iostream>
#include <array>
#include "helper.h"
#include "ScoreMatrix.h"
#include "PointCloud2D.h"

class ScoreMatrix2DPointCloud :public ScoreMatrix
{
public:

    typedef std::array<int,4>       ivector4;
    typedef std::array<int,2>       ivector2;
    typedef std::array<double, 7>   dvector7;
    typedef std::array<Vector3f,2>  lineSeg;
    typedef Vector3f                point3D;

    ScoreMatrix2DPointCloud(void);

    ScoreMatrix2DPointCloud(PointCloud2D* moduleSet, PointCloud2D* dataSet);

    ~ScoreMatrix2DPointCloud(void);

    /*** Setters ***/
    // Set the module set
    void SetModuleSet( PointCloud2D *moduleSet ) { m_moduleSet = moduleSet; Reset(); }
    // Set the data set
    void SetDataSet( PointCloud2D *dataSet )     { m_dataSet = dataSet; Reset(); }
    // Set both module and data set
    void SetModuleAndDataSet(PointCloud2D *moduleSet, PointCloud2D *dataSet) { m_moduleSet = moduleSet; m_dataSet = dataSet; Reset(); }
    // Set parameters to compute the score. --- change parameter list to double. ToDo.
    void SetParameters( double wf, double tdm, double tdd, double td, double ta, double tt, double tl )
    {
        // Record the thresholds
        m_parameterList[0] = wf; m_parameterList[1] = tdm; m_parameterList[2] = tdd; m_parameterList[3] = td;
        m_parameterList[4] = ta; m_parameterList[5] = tt;  m_parameterList[6] = tl;
        // Pre-compute the corresponding sigma values
        m_sigmaList[0] = ComputeSigma(wf); m_sigmaList[1] = ComputeSigma(tdm); m_sigmaList[2] = ComputeSigma(tdd);
        m_sigmaList[3] = ComputeSigma(td); m_sigmaList[4] = ComputeSigma(ta);  m_sigmaList[5] = ComputeSigma(tt);
        m_sigmaList[6] = ComputeSigma(tl);
    }

    // Set up real thresholds --- not used for computing sigma
    void SetThresholds( double wf, double tdm, double tdd, double td, double ta, double tt, double tl )
    {
        m_parameterList[0] = wf; m_parameterList[1] = tdm; m_parameterList[2] = tdd; m_parameterList[3] = td;
        m_parameterList[4] = ta; m_parameterList[5] = tt;  m_parameterList[6] = tl;
    }

    // Set the match truth --- to compute the correct rate
    // The match truth vector should has the same length as input "dataHelix"
    void SetMatchTruth(std::vector<int> matchTruth) { m_matchTruth = matchTruth; }

    // Set the score function
    void SetScoreFun(int flag) { m_scoreFunFlag = flag; }

    // Construct the score matrix --- add parameter list as std::array. ToDo.
    // Call this after module and data read in as well as matrix initialization
    void GenScoreMatrix();

    /*** Getters ***/
    // Get Parameter list
    const dvector7& GetParaList() const { return m_parameterList; }
    // Get match truth
    const std::vector<int> & GetMatchTruth() const { return m_matchTruth; }
    // Get score function we use
    const int GetScoreFun() const { return m_scoreFunFlag; }

    //const double GetMatrixSum() const { return m_matrixSum; }

private:

    // Module set and Data set. ATTENTION! We don't actually
    // hold the sample set in this class. They are interfaces
    // for us to access the sample set when computing entries
    // in score matrix. DO NOT delete them in Reset. Just set
    // The pointers to 0.

    // The size of the score matrix should be decided by module
    // set and data set.
    PointCloud2D *m_moduleSet;
    PointCloud2D *m_dataSet;

    // Threshold parameter vector: distance weight flag, module bounding box, data bounding box, midpoint distance, alpha&beta, theta, helix length
    // wf[0],tdm[1],tdd[2],td[3],ta[4],tt[5],tl[6]
    dvector7 m_parameterList;
    dvector7 m_sigmaList;
    // Those entries' scores lower than this value will be set to 0.
    double m_epsilon;

    // Flag to choose score function --- Gaussian by default
    int m_scoreFunFlag;

    // To hold the matching truth. To compute the accuracy rate.
    std::vector<int> m_matchTruth;

    // Matrix low triangle plus diagonal sum --- for debug
    //double m_matrixSum;

private:    

    // This function clear or release all the existing member variables
    void Clear();

    // Reset the matrix and its size after reading into new module, data or module and data.
    // Basically, delete the old matrix, create a new one and update dimension information.
    // This should follow "SetModuleSet" and "SetDataSet" --- ToDo
    void Reset();

    // Given matrix entry index (x,y), return the corresponding element
    // indices (maybe other information), i,j and k,l.
    ivector4 ComputeIndices(int x, int y);

    // Compute the flag for each index:
    // 0 means original line direction 1 means invert line direction
    //ivector2 ComputeFlags(int x, int y);

    // Fill in the entry of the score matrix
    double ComputeEntry(int x, int y);

    // The score function to use in ComputeEntry. It takes the distance
    // of two assignments and a kernel parameter and return the score
    // evaluated. 
    double ScoreFun(double diff, double sigma, int flag);

    // Score function using Gaussian Kernel
    double ScoreFunGau(double diff, double sigma);
    // Score function using Epanechnikov Kernel
    double ScoreFunEpan(double diff, double sigma);
    // Score function using Quadratic
    double ScoreQuadratic(double diff, double sigma);
    // Evaluate the sigma for different kernel
    double ComputeSigma(double threshold);
};

