/*********************************************************************/
// IPFP.cpp
// This is a class implementing IPFP. Basically it takes in an initial
// solution, initial score, data number (in our matching format case) 
// and score matrix.

// ATTENTION! IPFP relies on Hungarian which is not applicable for double
// line direction case in a straight forward way.

// Hang 11/17/2014
/*********************************************************************/

#include "IPFP.h"

IPFP::IPFP(void)
{
    Clear();
}

IPFP::IPFP(MatrixXd matrix, VectorXd result, double score, int maxStepNum, int moduleLen, int dataLen)
{
    SetScoreMatrix(matrix);

    SetInitialResult(result);

    SetInitialScore(score);

    SetMaxStepNum(maxStepNum);

    SetModuleLen(moduleLen);

    SetDataLen(dataLen);
}

IPFP::~IPFP(void)
{
}

void IPFP::Reset(MatrixXd matrix, VectorXd result, double score, int maxStepNum, int moduleLen, int dataLen)
{
    SetScoreMatrix(matrix);

    SetInitialResult(result);

    SetInitialScore(score);

    SetMaxStepNum(maxStepNum);

    SetModuleLen(moduleLen);

    SetDataLen(dataLen);
}

void IPFP::Clear()
{
    m_initialScore = 0.0;
    m_maxStepNum   = 0;
    m_RefinedScore = 0.0;
    m_moduleLen    = 0;
    m_dataLen      = 0;
}

// Add comments later. ToDo.
//
void IPFP::DoIPFP()
{
    VectorXd xOld  = m_initialIndicator;
    VectorXd xNew  = xOld;
    VectorXd xStar = xOld;

    double sStar = m_initialScore;

    for(int i = 0; i < m_maxStepNum; ++i)
    {
        /****** Hungarian step to obtain bNew******/
        // Initialize the Hungarian solver
        hungarian_problem_t p;
        // Set up Hungarian cost matrix
        double** m = array_to_matrix(((-1.0)*m_scoreMatrix*xOld).array().data(), m_dataLen, m_moduleLen);
        // Initialize the gungarian_problem using the cost matrix
        hungarian_init(&p, m , m_dataLen, m_moduleLen, HUNGARIAN_MODE_MINIMIZE_COST);
        // Solve the Hungarian
        hungarian_solve(&p);        

        // Convert the output from Hungarian solver into VectorXd format
        VectorXd bNew = VectorXd::Zero(m_initialIndicator.size());
        for(int i = 0; i < m_dataLen; ++i)
        {
            // Search for each row (f1_e1, f1_e2, ... ,f1_en) in p
            for(int j = 0; j < m_moduleLen; ++j)
                bNew[i*m_moduleLen + j] = (double)p.assignment[i][j];
        }

        //printf("bNew sum: %f\n", bNew.array().sum());

        /* free used memory */
        hungarian_free(&p);
        if( NULL != m)
        {
            for(int h=0; h<m_dataLen; h++)
                if(NULL != m[h])
                    free(m[h]);
            free(m);
            m = NULL;
        }


        /****** Other steps ******/
        double c = xOld.transpose() * m_scoreMatrix * (bNew - xOld);
        double d = (bNew - xOld).transpose() * m_scoreMatrix * (bNew - xOld);
        if(d > 0 || 0 == d)
        {
            xNew = bNew;
        }
        else
        {
            double r = std::min(-c/d, 1.0);
            xNew = xOld + r*(bNew - xOld);
        }

        double tempScore = bNew.transpose() * m_scoreMatrix * bNew;
        if(tempScore > sStar || tempScore == sStar)
        {
            sStar = tempScore;
            xStar = bNew;
        }

        if((xNew - xOld).squaredNorm() == 0)
            break;

        xOld = xNew;
    }

    m_RefinedScore     = sStar;
    m_refinedIndicator = xStar;

    // Convert the output from Hungarian solver into our format: ej[fi]
    m_RefinedResult = VectorXi::Zero(m_dataLen);
    for(int i = 0; i < m_dataLen; ++i)
    {
        // Search for each row (f1_e1, f1_e2, ... ,f1_en) in p
        for(int j = 0; j < m_moduleLen; ++j)
            if( 0 != xStar[i*m_moduleLen + j] )
                m_RefinedResult[i] = j+1;
    }
}