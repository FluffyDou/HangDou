/*********************************************************************/
// ScoreMatrixDoubleLineDir.cpp
// This the class holding the score matrix for helix line segment, with 
// single line directions.

// Hang 11/14/2014
/*********************************************************************/

#include "ScoreMatrixSingleLineDir.h"


ScoreMatrixSingleLineDir::ScoreMatrixSingleLineDir(void)
{
    Clear();
}


ScoreMatrixSingleLineDir::ScoreMatrixSingleLineDir(HelixGeometry *moduleSet, HelixGeometry *dataSet)
{
    SetModuleAndDataSet(moduleSet, dataSet);
}


ScoreMatrixSingleLineDir::~ScoreMatrixSingleLineDir(void)
{
    // ATTENTION! Do not delete m_moduleSet and m_dataSet.
    // This class only access to them rather than owning them.
}



// This function clear or release all the existing member variables
void ScoreMatrixSingleLineDir::Clear()
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
void ScoreMatrixSingleLineDir::GenScoreMatrix()
{
    // For each entry, compute the score.
    // Since the matrix is symmetric, only compute the lower triangle.
    // Use OpenMP to do this with multi-thread. ToDo.
    // Use OpenAL to do this in parallel. ToDo.

    //double tempSum = 0.0;
    //m_matrixSum = 0.0;
    for(int i = 0; i < m_dim; ++i)
        for(int j = 0; j < i || j == i; ++j)
        {
            m_scoreM(i, j) = ComputeEntry(j, i);
            //tempSum += m_scoreM(i, j);
        }

    //std::cout<<"Matrix lower triangle Sum: "<<tempSum<<"\n";
}


// Construct the score matrix for Orient Single Line Direction
void ScoreMatrixSingleLineDir::GenScoreMatrixOrientLine()
{
    for(int i = 0; i < m_dim; ++i)
        for(int j = 0; j < i || j == i; ++j)
        {
            m_scoreM(i, j) = ComputeEntryOrientLine(j, i);
            //tempSum += m_scoreM(i, j);
        }
}


// Fill in the entry of the score matrix --- To Modify: no flags
double ScoreMatrixSingleLineDir::ComputeEntry(int x, int y)
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

    if(x == y)
    {
        // Length difference between two lines in one assignment
        const LineSegment & dataLineSeg   = m_dataSet->GetLineSegment(i);
        const LineSegment & moduleLineSeg = m_moduleSet->GetLineSegment(j);

        double diffL = (double) abs( dataLineSeg.GetLen() - moduleLineSeg.GetLen() );
        // Check if past threshold
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffL > m_parameterList[6] ) 
            return 0;

        double finalScore = ScoreFun(diffL, m_sigmaList[6], m_scoreFunFlag);

        return finalScore > m_epsilon ? finalScore : 0;

    }
    // Score for off diagonal elements
    // Conflict features excluded
    else if( i != k && j != l)
    {
        // Line segments from module and data set
        const LineSegment & dataLine_i   = m_dataSet->GetLineSegment(i);
        const LineSegment & moduleLine_j = m_moduleSet->GetLineSegment(j);
        const LineSegment & dataLine_k   = m_dataSet->GetLineSegment(k);
        const LineSegment & moduleLine_l = m_moduleSet->GetLineSegment(l);

        // Mid points for fi, ej ; fk, el. All can be precomputed. ToDo.
        const Vector3f & mMidPj = moduleLine_j.GetMidPoint();
        const Vector3f & mMidPl = moduleLine_l.GetMidPoint();
        const Vector3f & dMidPi = dataLine_i.GetMidPoint();
        const Vector3f & dMidPk = dataLine_k.GetMidPoint();

        // Line vectors in module set. mvj & mvl are precomputed.
        // For module set, the line direction is fixed : first point -> second point
        const Vector3f mvjl = (mMidPl - mMidPj).normalized(); // mid line vector
        const Vector3f mvj  = moduleLine_j.GetLineDir();
        const Vector3f mvl  = moduleLine_l.GetLineDir();

        // There are two line directions for data set
        Vector3f dvik = (dMidPk - dMidPi).normalized(); // mid line vector
        Vector3f dvi1  = dataLine_i.GetLineDir();
        Vector3f dvk1  = dataLine_k.GetLineDir();
        Vector3f dvi2  = -dvi1;
        Vector3f dvk2  = -dvk1;

        /***** Midpoint distance difference *****/
        double mLen  = EuDis(mMidPj, mMidPl);
        double dLen  = EuDis(dMidPi, dMidPk);
        double diffD = abs(mLen - dLen);
        if( (m_parameterList[3] > 0 || m_parameterList[3] == 0) && diffD > m_parameterList[3] ) 
            return 0;

        // Distance weight threshold.
        if( (mLen > m_parameterList[1] || dLen > m_parameterList[2]) && 0 != m_parameterList[0] && (m_parameterList[1] > 0 || m_parameterList[1] == 0) && (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
            return 0;

        /***** Helix length difference *****/
        // All the helix length can be precomputed. ToDo.
        // Difference between i, j
        double iLen = dataLine_i.GetLen();
        double jLen = moduleLine_j.GetLen();
        double diffLIJ = abs(iLen - jLen);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffLIJ > m_parameterList[6] ) 
            return 0;

        // Difference between k, l
        double kLen = dataLine_k.GetLen();
        double lLen = moduleLine_l.GetLen();
        double diffLKL = abs(kLen - lLen);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffLKL > (double)m_parameterList[6] ) 
            return 0;

        /***** Alpha, Beta and Theta for Module side *****/
        // Alpha
        double mAlpha = VecAng(mvj, mvjl);
        // Beta
        double mBeta = VecAng(mvl, -mvjl);
        // Theta for module set
        Vector3f tempV1 = mvj - mvjl * mvj.dot(mvjl);
        Vector3f tempV2 = mvl - mvjl * mvl.dot(mvjl);
        double mTheta   = VecAng(tempV1, tempV2);

        /***** Alpha, Beta and Theta for data side *****/
        // Loop to choose the direction with minimum angle difference.
        double diffA = 0.0; double diffB = 0.0; double diffT = 0.0;

        // If all the angle differences pass the threshold.
        bool passThreshFlag = true;
        //  4 paris in all: i1, k1; i1, k2; i2, k1; i2, k2
        std::vector< std::pair<Vector3f, Vector3f> > ikVecPair(4);
        ikVecPair[0].first = dvi1; ikVecPair[0].second = dvk1;
        ikVecPair[1].first = dvi1; ikVecPair[1].second = dvk2;
        ikVecPair[2].first = dvi2; ikVecPair[2].second = dvk1;
        ikVecPair[3].first = dvi2; ikVecPair[3].second = dvk2;

        // The largest angle sum should be less than or equal to 3.
        double minAngleSum = 100000.0;
        // Loop over the 4 pairs to pick out smallest angle difference
        for(auto it = ikVecPair.begin(); it != ikVecPair.end(); ++it)
        {
            // data alpha
            double dAlpha    = VecAng(it->first, dvik);
            double tempDiffA = abs(mAlpha - dAlpha);
            if( (m_parameterList[4] > 0 || m_parameterList[4] == 0) && tempDiffA > m_parameterList[4] ) 
                continue;
            // data beta
            double dBeta     = VecAng(it->second, -dvik);
            double tempDiffB = abs(mBeta - dBeta);
            if( (m_parameterList[4] > 0 || m_parameterList[4] == 0) && tempDiffB > m_parameterList[4] ) 
                continue;
            // data theta
            tempV1 = it->first - dvik * it->first.dot(dvik);
            tempV2 = it->second - dvik * it->second.dot(dvik);
            double dTheta    = VecAng(tempV1, tempV2);
            double tempDiffT = abs(mTheta - dTheta);
            if( (m_parameterList[5] > 0 || m_parameterList[5] == 0) && tempDiffT > m_parameterList[5] )
                continue;

            // For the current pair, all the angle differences are within the threshold
            passThreshFlag = false;
            // Record the angle difference with smallest angle difference sum
            double tempSum = tempDiffA + tempDiffB + tempDiffT;
            if( tempSum < minAngleSum )
            {
                minAngleSum = tempSum;
                diffA = tempDiffA; diffB = tempDiffB; diffT = tempDiffT;
            }
        }

        // None of the pair combination is within the threshold
        if( passThreshFlag )
            return 0;

        // Multiply the feature scores together
        double finalScore = 1.0;
        if(  (m_parameterList[4] > 0 || m_parameterList[4] == 0) )
            finalScore *= ScoreFun(diffA, m_sigmaList[4], m_scoreFunFlag) * ScoreFun(diffB, m_sigmaList[4], m_scoreFunFlag);
        if( (m_parameterList[5] > 0 || m_parameterList[5] == 0) )
            finalScore *= ScoreFun(diffT, m_sigmaList[5], m_scoreFunFlag);
        if( (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
            finalScore *= ScoreFun(diffD, m_sigmaList[3], 0);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == m_scoreFunFlag) )
            finalScore *= ScoreFun(diffLIJ, m_sigmaList[6], m_scoreFunFlag) * ScoreFun(diffLKL, m_sigmaList[6], m_scoreFunFlag);

        // Distance weight
        if( m_parameterList[0] != 0 )
        {
            if( (m_parameterList[1] > 0 || m_parameterList[1] == 0) )
                finalScore *= ScoreFun(mLen, m_sigmaList[1], m_scoreFunFlag);
            if( (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
                finalScore *= ScoreFun(dLen, m_sigmaList[2], m_scoreFunFlag);
        }

        return finalScore > m_epsilon ? finalScore : 0;
    }

    // Else, return 0 (conflict entry)
    return 0;
}


// Fill in the entry of the score matrix for Orient Single Line Direction --- Direction: first point to second point
double ScoreMatrixSingleLineDir::ComputeEntryOrientLine(int x, int y)
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

    if(x == y)
    {
        // Length difference between two lines in one assignment
        const LineSegment & dataLineSeg   = m_dataSet->GetLineSegment(i);
        const LineSegment & moduleLineSeg = m_moduleSet->GetLineSegment(j);

        double diffL = (double) abs( dataLineSeg.GetLen() - moduleLineSeg.GetLen() );
        // Check if past threshold
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffL > m_parameterList[6] ) 
            return 0;

        double finalScore = ScoreFun(diffL, m_sigmaList[6], m_scoreFunFlag);

        return finalScore > m_epsilon ? finalScore : 0;

    }
    // Score for off diagonal elements
    // Conflict features excluded
    else if( i != k && j != l)
    {
        // Line segments from module and data set
        const LineSegment & dataLine_i   = m_dataSet->GetLineSegment(i);
        const LineSegment & moduleLine_j = m_moduleSet->GetLineSegment(j);
        const LineSegment & dataLine_k   = m_dataSet->GetLineSegment(k);
        const LineSegment & moduleLine_l = m_moduleSet->GetLineSegment(l);

        // Mid points for fi, ej ; fk, el. All can be precomputed. ToDo.
        const Vector3f & mMidPj = moduleLine_j.GetMidPoint();
        const Vector3f & mMidPl = moduleLine_l.GetMidPoint();
        const Vector3f & dMidPi = dataLine_i.GetMidPoint();
        const Vector3f & dMidPk = dataLine_k.GetMidPoint();

        // Line vectors in module set. mvj & mvl can be precomputed. ToDo.
        // For module set, the line direction is fixed : first point -> second point
        const Vector3f mvjl = (mMidPl - mMidPj).normalized(); // mid line vector
        const Vector3f mvj  = moduleLine_j.GetLineDir();
        const Vector3f mvl  = moduleLine_l.GetLineDir();

        // Line vectors in data set. dvi & dvk can be precomputed. ToDo.
        Vector3f dvik = (dMidPk - dMidPi).normalized(); // mid line vector
        Vector3f dvi  = dataLine_i.GetLineDir();
        Vector3f dvk  = dataLine_k.GetLineDir();

        /***** Midpoint distance difference *****/
        double mLen  = EuDis(mMidPj, mMidPl);
        double dLen  = EuDis(dMidPi, dMidPk);
        double diffD = abs(mLen - dLen);
        if( (m_parameterList[3] > 0 || m_parameterList[3] == 0) && diffD > m_parameterList[3] ) 
            return 0;

        // Distance weight threshold.
        if( (mLen > m_parameterList[1] || dLen > m_parameterList[2]) && 0 != m_parameterList[0] && (m_parameterList[1] > 0 || m_parameterList[1] == 0) && (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
            return 0;

        /***** Helix length difference *****/
        // All the helix length are precomputed.
        // Difference between i, j
        double iLen = dataLine_i.GetLen();
        double jLen = moduleLine_j.GetLen();
        double diffLIJ = abs(iLen - jLen);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffLIJ > m_parameterList[6] ) 
            return 0;

        // Difference between k, l
        double kLen = dataLine_k.GetLen();
        double lLen = moduleLine_l.GetLen();
        double diffLKL = abs(kLen - lLen);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == 0) && diffLKL > (double)m_parameterList[6] ) 
            return 0;

        /***** Alpha difference *****/
        // Alpha : j & mid line, i & mid line
        double mAlpha = VecAng(mvj, mvjl);
        double dAlpha = VecAng(dvi, dvik);
        double diffA  = abs(mAlpha - dAlpha);
        if( (m_parameterList[4] > 0 || m_parameterList[4] == 0) && diffA > m_parameterList[4] ) 
            return 0;

        /***** Beta difference *****/
        // Beta : l & mid line, k & mid line
        double mBeta = VecAng(mvl, -mvjl);
        double dBeta = VecAng(dvk, -dvik);
        double diffB = abs(mBeta - dBeta);
        if( (m_parameterList[4] > 0 || m_parameterList[4] == 0) && diffB > m_parameterList[4] ) 
            return 0;

        /***** Theta difference *****/
        // Gorgon projects the two line vector to the middle line vector.
        // Then "connect" the project point with line end point.
        // Then it uses the vector to compute the theta.

        // Theta for module set
        Vector3f tempV1 = mvj - mvjl * mvj.dot(mvjl);
        Vector3f tempV2 = mvl - mvjl * mvl.dot(mvjl);
        double mTheta   = VecAng(tempV1, tempV2);
        // Theta for data set
        tempV1 = dvi - dvik * dvi.dot(dvik);
        tempV2 = dvk - dvik * dvk.dot(dvik);
        double dTheta = VecAng(tempV1, tempV2);
        double diffT  = abs(mTheta - dTheta);
        if( (m_parameterList[5] > 0 || m_parameterList[5] == 0) && diffT > m_parameterList[5] )
            return 0;

        /***** Compute the final score *****/

        // Multiply the feature scores together
        double finalScore = 1.0;
        if(  (m_parameterList[4] > 0 || m_parameterList[4] == 0) )
            finalScore *= ScoreFun(diffA, m_sigmaList[4], m_scoreFunFlag) * ScoreFun(diffB, m_sigmaList[4], m_scoreFunFlag);
        if( (m_parameterList[5] > 0 || m_parameterList[5] == 0) )
            finalScore *= ScoreFun(diffT, m_sigmaList[5], m_scoreFunFlag);
        if( (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
            finalScore *= ScoreFun(diffD, m_sigmaList[3], 0);
        if( (m_parameterList[6] > 0 || m_parameterList[6] == m_scoreFunFlag) )
            finalScore *= ScoreFun(diffLIJ, m_sigmaList[6], m_scoreFunFlag) * ScoreFun(diffLKL, m_sigmaList[6], m_scoreFunFlag);

        // Distance weight
        if( m_parameterList[0] != 0 )
        {
            if( (m_parameterList[1] > 0 || m_parameterList[1] == 0) )
                finalScore *= ScoreFun(mLen, m_sigmaList[1], m_scoreFunFlag);
            if( (m_parameterList[2] > 0 || m_parameterList[2] == 0) )
                finalScore *= ScoreFun(dLen, m_sigmaList[2], m_scoreFunFlag);
        }

        return finalScore > m_epsilon ? finalScore : 0;
    }

    // Else, return 0 (conflict entry)
    return 0;
}


// Given matrix entry index (x,y), return the corresponding element
// indices (maybe other information), i,j and k,l.
ScoreMatrixSingleLineDir::ivector4 ScoreMatrixSingleLineDir::ComputeIndices(int x, int y)
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
double ScoreMatrixSingleLineDir::ScoreFun(double diff, double sigma, int flag)
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
    default:
        std::cout<<"Invalid input flag ....";
        break;
    }

    return 0;
}


// Score function using Gaussian Kernel
double ScoreMatrixSingleLineDir::ScoreFunGau(double diff, double sigma)
{
    //std::cout<<diff<<", "<<-(diff*diff)<<", "<<(2.0*sigma*sigma)<<", "<<-(diff*diff)/(2.0*sigma*sigma)<<"\n";
    return exp( -(diff*diff)/(2.0*sigma*sigma) );
}
// Score function using Epanechnikov Kernel
double ScoreMatrixSingleLineDir::ScoreFunEpan(double diff, double sigma)
{
    // We assume: threshold == 3 * sigma
    return 1.0 - (diff*diff)/(9.0*sigma*sigma);
}
// Compute the sigma value based on threshold
double ScoreMatrixSingleLineDir::ComputeSigma(double threshold)
{

    return threshold / 3.0;
}

void ScoreMatrixSingleLineDir::Reset()
{
    // Update matrix dimension
    m_moduleNum = m_moduleSet->GetHelixNum();
    m_dataNum   = m_dataSet->GetHelixNum();
    m_dim       = m_dataNum * m_moduleNum;

    // Update the matrix
    //m_scoreM = new MatrixXd(m_dim, m_dim);
    m_scoreM.resize(m_dim, m_dim);

    // Set score function back to Gaussian
    m_scoreFunFlag = 0;

    //m_matrixSum = 0.0;
}