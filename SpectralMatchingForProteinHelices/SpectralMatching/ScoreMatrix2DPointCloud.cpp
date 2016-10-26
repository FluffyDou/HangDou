/*********************************************************************/
// ScoreMatrix2DPointCloud.cpp
// This the class holding the score matrix for 2D point cloud.
// This class is very similar to that for single line direction. The
// only difference is we only care about midpoint distance now.

// Hang 11/19/2014
/*********************************************************************/

#include "ScoreMatrix2DPointCloud.h"


ScoreMatrix2DPointCloud::ScoreMatrix2DPointCloud(void)
{
    Clear();
}


ScoreMatrix2DPointCloud::ScoreMatrix2DPointCloud(PointCloud2D *moduleSet, PointCloud2D *dataSet)
{
    SetModuleAndDataSet(moduleSet, dataSet);
}


ScoreMatrix2DPointCloud::~ScoreMatrix2DPointCloud(void)
{
    // ATTENTION! Do not delete m_moduleSet and m_dataSet.
    // This class only access to them rather than owning them.
}



// This function clear or release all the existing member variables
void ScoreMatrix2DPointCloud::Clear()
{
    //if( 0 != m_scoreM )
    //    delete m_scoreM;

    m_scoreM.resize(0,0);

    ivector2 dim = {0, 0};
    m_matrixDim  = dim;

    m_moduleSet = 0;
    m_dataSet   = 0;

    m_epsilon = 0.000001;

    m_scoreFunFlag = 0;

    //m_matrixSum = 0.0;
}

// Construct the score matrix
void ScoreMatrix2DPointCloud::GenScoreMatrix()
{
    // For each entry, compute the score.
    // Since the matrix is symmetric, only compute the lower triangle.
    // Use OpenMP to do this with multi-thread. ToDo.
    // Use OpenAL to do this in parallel. ToDo.

//    double tempSum = 0.0;
    //m_matrixSum = 0.0;
    for(int i = 0; i < m_dim; ++i)
        for(int j = 0; j < i || j == i; ++j)
        {
            m_scoreM(i, j) = ComputeEntry(j, i);
//            tempSum += m_scoreM(i, j);
        }

//        std::cout<<"Matrix lower triangle Sum: "<<tempSum<<"\n";
}


// Fill in the entry of the score matrix --- To Modify: no flags
// ATTENTION!!! Now the matrix is for TPS point model in SM paper.
double ScoreMatrix2DPointCloud::ComputeEntry(int x, int y)
{
    // fi,ej; fk,el
    // Compute the indices: i[0],j[1],k[2],l[3]
    ivector4 indices = ComputeIndices(x, y);
    int i = indices[0];
    int j = indices[1];
    int k = indices[2];
    int l = indices[3];

    // Compute score for diagonal elements: i == k, j == l
    // Conflict features excluded

    if( i != k && j != l)
    {
        // Points from module and data set
        const Vector2f & dataPoint_i   = m_dataSet->GetPoint(i);
        const Vector2f & modulePoint_j = m_moduleSet->GetPoint(j);
        const Vector2f & dataPoint_k   = m_dataSet->GetPoint(k);
        const Vector2f & modulePoint_l = m_moduleSet->GetPoint(l);

        /***** Midpoint distance difference *****/
        double mLen  = (modulePoint_j - modulePoint_l).norm();
        double dLen  = (dataPoint_i - dataPoint_k).norm();

        //double diffD = abs(mLen - dLen);
        //if( (m_parameterList[3] > 0 || m_parameterList[3] == 0) && (diffD > m_parameterList[3] || diffD == m_parameterList[3]) ) 
        //    return 0;

        //double finalScore = ScoreFun(diffD, m_sigmaList[3], m_scoreFunFlag);

        double tempEpsilon = 0.000001;
        double diffD       = abs((mLen - dLen)/(dLen+tempEpsilon));
        if( (m_parameterList[3] > 0 || m_parameterList[3] == 0) && (diffD > m_parameterList[3] || diffD == m_parameterList[3]) ) 
            return 0;

        double tempGama = 0.5;
        double finalScore = ScoreFun(diffD, m_sigmaList[3], m_scoreFunFlag) * (1.0 - tempGama);

        // .... 
        Vector2f dvik = dataPoint_i - dataPoint_k;
        Vector2f dvjl = modulePoint_j - modulePoint_l;

        double diffA = VecAng(dvik, dvjl);
        // We use parameter 4 to control the angle threshold
        if( (m_parameterList[4] > 0 || m_parameterList[4] == 0) && (diffA > m_parameterList[4] || diffA == m_parameterList[4]) ) 
            return 0;

        finalScore += ScoreFun(diffA, m_sigmaList[4], m_scoreFunFlag) * tempGama;

        return finalScore > m_epsilon ? finalScore : 0;
    }

    // Else, return 0 (conflict entry)
    return 0;
}


// Given matrix entry index (x,y), return the corresponding element
// indices (maybe other information), i,j and k,l.
ScoreMatrix2DPointCloud::ivector4 ScoreMatrix2DPointCloud::ComputeIndices(int x, int y)
{

    int i     = y / m_moduleNum;
    int j     = y % m_moduleNum;

    int k     = x / m_moduleNum;
    int l     = x % m_moduleNum;

    ivector4 indices = {i,j,k,l};
    return indices;
}

// Return the score evaluated by picked score function
// ATTENTION! Really should change to use function pointer! ToDo.
double ScoreMatrix2DPointCloud::ScoreFun(double diff, double sigma, int flag)
{
    // Choose the score function to use in ComputeEntry.
    switch( flag )
    {
    case 0: 
        return ScoreFunGau(diff, sigma);
        break;
    case 1: 
        return ScoreFunEpan(diff, sigma);
        break;
    case 2: 
        return ScoreQuadratic(diff, sigma);
        break;
    default:
        std::cout<<"Invalid input flag ....";
        break;
    }

    return 0;
}


// Score function using Gaussian Kernel
double ScoreMatrix2DPointCloud::ScoreFunGau(double diff, double sigma)
{
    //std::cout<<diff<<", "<<-(diff*diff)<<", "<<(2.0*sigma*sigma)<<", "<<-(diff*diff)/(2.0*sigma*sigma)<<"\n";
    return exp( -(diff*diff)/(2.0*sigma*sigma) );
}
// Score function using Epanechnikov Kernel
double ScoreMatrix2DPointCloud::ScoreFunEpan(double diff, double sigma)
{
    // We assume: threshold == 3 * sigma
    return 1.0 - (diff*diff)/(9.0*sigma*sigma);
}
// Score function using unnormalized Epanechnikov Kernel
double ScoreMatrix2DPointCloud::ScoreQuadratic(double diff, double sigma)
{
    // We assume: threshold == 3 * sigma
    return 4.5 - (diff*diff)/(2.0*sigma*sigma);
}

// Compute the sigma value based on threshold
double ScoreMatrix2DPointCloud::ComputeSigma(double threshold)
{
    return threshold / 3.0;
}

void ScoreMatrix2DPointCloud::Reset()
{
    // Update matrix dimension
    m_moduleNum = m_moduleSet->GetPointNum();
    m_dataNum   = m_dataSet->GetPointNum();
    m_dim       = m_dataNum * m_moduleNum;

    // Update the matrix
    //m_scoreM = new MatrixXd(m_dim, m_dim);
    m_scoreM.resize(m_dim, m_dim);

    // Set score function back to Gaussian
    m_scoreFunFlag = 0;

    //m_matrixSum = 0.0;
}