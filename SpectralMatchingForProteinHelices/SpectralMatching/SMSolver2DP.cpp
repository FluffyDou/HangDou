/*********************************************************************/
// SMSolver2DP.cpp
// This is a class holding solvers for Spectral Matching for 2D point
// cloud, including solving for eigenvectors, binarizing the result,
// computing indicator vector and computing matrix score sum (x'Mx)

// Hang 11/21/2014
/*********************************************************************/

#include "SMSolver2DP.h"


SMSolver2DP::SMSolver2DP(void)
{
    m_matrixSize    = 0;
    m_FunctionScore = 0.0;
}


SMSolver2DP::~SMSolver2DP(void)
{
}

// Compute the score sum for Score Function
void SMSolver2DP::ComputeFunScore( const ScoreMatrix2DPointCloud & scoreMatrix )
{
    const int mLen = scoreMatrix.GetModuleNum();
    const int dLen = scoreMatrix.GetDataNum();

    // Compute indicator vector based on matching truth
    VectorXd indicatorVec = VectorXd::Zero( scoreMatrix.GetDim() );

    for(int i = 0; i < dLen; ++i)
    {
        if(0 != m_matchResult[i])
            indicatorVec[ mLen*i + m_matchResult[i]-1 ] = 1.0;
    }

    // Compute x'Mx
    MatrixXd tempMatrix = scoreMatrix.GetMatrix();
    tempMatrix.triangularView<Eigen::Upper>() = tempMatrix.transpose();

    m_FunctionScore  = indicatorVec.transpose() * tempMatrix * indicatorVec;
    // Return the normalized function score: (x'Mx)/(x'x)
    //m_FunctionScore /=  indicatorVec.dot(indicatorVec);
}

// Solve for the principle eigenvector
void SMSolver2DP::SolveForEigenvectors( const MatrixXd & matrix )
{
    m_solver.compute(matrix);
    // Record the matrix size
    m_matrixSize = matrix.cols();
    // Record the principle eigenvector
    m_prinEigenVec = m_solver.eigenvectors().col(m_matrixSize-1);
}

// Do Simple Greedy on double line direction score matrix
void SMSolver2DP::SimpleGreedyPick( const ScoreMatrix2DPointCloud & scoreMatrix )
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    //m_matchFlag   = m_matchResult;

    // Keep a remain index array, denoting remaining rows and columns
    // Initially, all the rows remain.
    std::vector<bool> remainIndex(matrix.cols(), true);
    // This is for matrix shrink
    std::vector<bool> remainIndexOld = remainIndex;

    // Array to hold the padding back eigenvector
    std::vector<double> paddedEigenvec(matrix.cols(), 0);

    // Prune all the zero rows --- only need to be done once
    PruneMatrix(matrix, remainIndex, remainIndexOld);

    //Solve the matrix (remain matrix) for principle eigenvector
    VectorXd prunedVec = m_solver.compute(matrix).eigenvectors().col(matrix.cols()-1);
    // Get absolute value of the principle eigenvector
    prunedVec = prunedVec.cwiseAbs();

    //Pad the eigenvector back to original length
    PadEigenvec(prunedVec, paddedEigenvec, remainIndex);

    // Make a pair list of the eigenvector so that after sorting
    // We can still find each element's original index
    std::vector<std::pair<double, int>> finalEigenvec(paddedEigenvec.size());
    for(int i = 0; i < (int)finalEigenvec.size(); ++i)
    {
        // First value and Second index
        finalEigenvec[i].first  = paddedEigenvec[i];
        finalEigenvec[i].second = i;
    }

    // Sort the principle eigenvector
    std::sort(finalEigenvec.begin(), finalEigenvec.end(), IndexCompare);

    // In case dataLen is larger than moduleLen
    int resultLen = std::min(moduleLen, dataLen);
    int pickedNum = 0;
    // Loop to pick out the largest element and exclude conflicts
    for(int i = 0; i < matrix.cols(); ++i)
    {
        // If the largest remaining element is zero, we can stop
        if(0 == finalEigenvec[i].first || pickedNum == resultLen)
            break;

        // This one has been excluded
        if(0 == paddedEigenvec[finalEigenvec[i].second])
            continue;

        // The remaining largest element
        int acceptIndex = finalEigenvec[i].second;

        // Append the picked assignment to the result
        std::array<int, 2> indices = ObtainFeatureIndices(acceptIndex, moduleLen);
        m_matchResult[indices[0]] = indices[1]+1;
        //int tempFlag = ObtainFeatureFlag(acceptIndex, moduleLen);
        //m_matchFlag[indices[0]] = tempFlag;

        // Exclude conflict assignments
        // Exclude all the fi and fi'
        int start = moduleLen*indices[0];
        int end   = start + moduleLen;
        for( int h = start; h < end; ++h )
        {
            paddedEigenvec[h] = 0;                
        }

        // Exclude all the ej
        start = 0;
        end   = dataLen;
        for( int h = start; h < end; ++h )
        {
            paddedEigenvec[ h*moduleLen + indices[1] ] = 0.0;
        }

        ++pickedNum;
    }

    //// Output results
    //for(auto it = m_matchResult.begin(); it != m_matchResult.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
    //// Output Flags
    //for(auto it = m_matchFlag.begin(); it != m_matchFlag.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
}

// Do Iterative Pick on double line direction score matrix
void SMSolver2DP::IterativePicking(const ScoreMatrix2DPointCloud & scoreMatrix)
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    //m_matchFlag   = m_matchResult;

    // Keep a remain index array, denoting remaining rows and columns
    // Initially, all the rows remain.
    std::vector<bool> remainIndex(matrix.cols(), true);
    // This is for matrix shrink
    std::vector<bool> remainIndexOld = remainIndex;

    // Array to hold the padding back eigenvector
    std::vector<double> paddedEigenvec(matrix.cols(), 0);

    //std::cout<<"matrix size: "<<matrix.cols()<<"\n";

    // Prune all the zero rows --- only need to be done once
    PruneMatrix(matrix, remainIndex, remainIndexOld);    
    
    // While the matrix size (remain matrix) is still larger than 0
    while( matrix.cols() > 0 )
    {
        //std::cout<<"zero rows:\n";
        //for(int i = 0; i < (int)remainIndex.size(); ++i)
        //{
        //    if(remainIndex[i] == false)
        //        std::cout<<i+1<<", ";
        //}
        //std::cout<<"\n";

        //Solve the matrix (remain matrix) for principle eigenvector
        VectorXd prunedVec = m_solver.compute(matrix).eigenvectors().col(matrix.cols()-1);
        // Get absolute value of the principle eigenvector
        prunedVec = prunedVec.cwiseAbs();

        //std::cout<<"matrix size: "<<matrix.cols()<<"\n";
        //std::cout<<"eigenvector sum: "<<prunedVec.sum()<<"\n";

        //Pad the eigenvector back to original length
        PadEigenvec(prunedVec, paddedEigenvec, remainIndex);

        //double sum = 0.0;
        //double templargest = -1;
        //int    tempIndex = 0;
        //for(int i = 0; i < (int)paddedEigenvec.size(); ++i)
        //{
        //    sum += paddedEigenvec[i];
        //    if(paddedEigenvec[i] > templargest)
        //    {
        //        templargest = paddedEigenvec[i];
        //        tempIndex = i;
        //    }
        //}
        //std::cout<<"Sum after padding: "<<sum<<"\n";
        //std::cout<<"largest index: "<<tempIndex+1<<"\n";

        // Find or pick the largest element (its index) of the padded eigenvector
        // NOTE: So far we don't take equal value into consideration.
        int largestEleIndex = std::distance(paddedEigenvec.begin(), std::max_element(paddedEigenvec.begin(), paddedEigenvec.end()));
        //std::cout<<"(+) picked index: "<<largestEleIndex+1<<"\n";




        //// Pick the largest in Hungarian choice
        ///*******************************************************************/
        //for(auto it = paddedEigenvec.begin(); it != paddedEigenvec.end(); ++it)
        //{
        //    *it = -*it;
        //}
        //// Initialize the hungarian solver
        //hungarian_problem_t p;
        //// Set up hungarian cost matrix
        //double** m = array_to_matrix(paddedEigenvec.data(), dataLen, moduleLen);
        //// Initialize the gungarian_problem using the cost matrix
        //hungarian_init(&p, m , dataLen, moduleLen, HUNGARIAN_MODE_MINIMIZE_COST);
        //// Solve the hungarian
        //hungarian_solve(&p);

        //// Find the largest assignment in result from Hungarian
        //double tempLargest = 0.0;
        //int    tempLargestIndex = 0;
        //// Convert the output from Hungarian solver into our format: ej[fi]
        //for(int i = 0; i < dataLen; ++i)
        //{
        //    // Search for each row in p
        //    for(int j = 0; j < moduleLen; ++j) 
        //        if( 0 != p.assignment[i][j] )
        //        {
        //            // Find the corresponding value in padded eigenvector
        //            int    tempIndex = moduleLen*i +j;
        //            double tempVal   = -paddedEigenvec[tempIndex];
        //            if(tempVal > tempLargest)
        //            {
        //                tempLargest      = tempVal;
        //                tempLargestIndex = tempIndex;
        //            }
        //            break;
        //        }
        //}

        //int largestEleIndex = tempLargestIndex;

        ///* free used memory */
        //hungarian_free(&p);
        //if(NULL != m)
        //{
        //    for(int i=0; i<dataLen; i++)
        //        if(NULL != m[i])
        //            free(m[i]);

        //    free(m);
        //    m = NULL;
        //}

        /*******************************************************************/




        // Append the picked assignment to the result
        std::array<int, 2> tempResult = ObtainFeatureIndices(largestEleIndex, moduleLen);
        m_matchResult[tempResult[0]] = tempResult[1]+1;
        //int tempFlag = ObtainFeatureFlag(largestEleIndex, moduleLen);
        //m_matchFlag[tempResult[0]] = tempFlag;

        // Update the remaining rows --- add zero rows
        ExcludeAssigns(largestEleIndex, remainIndex, scoreMatrix.GetModuleNum(), scoreMatrix.GetDataNum());

        // Shrink the matrix
        ShrinkMatrix(largestEleIndex, matrix, remainIndex, remainIndexOld);

        // To locate rows (columns) index, for matrix shrink
        remainIndexOld = remainIndex;
    }

    //// Output results
    //for(auto it = m_matchResult.begin(); it != m_matchResult.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
    //// Output Flags
    //for(auto it = m_matchFlag.begin(); it != m_matchFlag.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
}

// Do Simple Greedy on double line direction score matrix
void SMSolver2DP::HungarianPick( const ScoreMatrix2DPointCloud & scoreMatrix )
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    //m_matchFlag   = m_matchResult;

    // Keep a remain index array, denoting remaining rows and columns
    // Initially, all the rows remain.
    std::vector<bool> remainIndex(matrix.cols(), true);
    // This is for matrix shrink
    std::vector<bool> remainIndexOld = remainIndex;

    // Array to hold the padding back eigenvector
    std::vector<double> paddedEigenvec(matrix.cols(), 0);

    // Prune all the zero rows --- only need to be done once
    PruneMatrix(matrix, remainIndex, remainIndexOld);

    //Solve the matrix (remain matrix) for principle eigenvector
    VectorXd prunedVec = m_solver.compute(matrix).eigenvectors().col(matrix.cols()-1);
    // Get absolute value of the principle eigenvector
    prunedVec = prunedVec.cwiseAbs();

    //Pad the eigenvector back to original length
    PadEigenvec(prunedVec, paddedEigenvec, remainIndex);

    /**** Use Munkres method to binarize ****/
    // Wrap that into a class. ToDo.
    // Change all positive eigenvector to all negative
    for(auto it = paddedEigenvec.begin(); it != paddedEigenvec.end(); ++it)
    {
        *it = -*it;
    }

    // Initialize the hungarian solver
    hungarian_problem_t p;
    // Set up hungarian cost matrix
    double** m = array_to_matrix(paddedEigenvec.data(), dataLen, moduleLen);
    // Initialize the gungarian_problem using the cost matrix
    hungarian_init(&p, m , dataLen, moduleLen, HUNGARIAN_MODE_MINIMIZE_COST);
    // Solve the hungarian
    hungarian_solve(&p);

    // Convert the output from Hungarian solver into our format: ej[fi]
    for(int i = 0; i < dataLen; ++i)
    {
        // Search for each row in p
        for(int j = 0; j < moduleLen; ++j) 
            if( 0 != p.assignment[i][j] )
            {
                // The matching index starts from 1 --- the same as Mathematica
                m_matchResult[i] = j+1;
                break;
            }
    }

    /* free used memory */
    hungarian_free(&p);
    if(NULL != m)
    {
        for(int i=0; i<dataLen; i++)
            if(NULL != m[i])
                free(m[i]);

        free(m);
        m = NULL;
    }


    //// Output results
    //for(auto it = m_matchResult.begin(); it != m_matchResult.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
    //// Output Flags
    //for(auto it = m_matchFlag.begin(); it != m_matchFlag.end(); ++it)
    //    std::cout<<*it<<",";
    //std::cout<<"\n";
}

// Prune the zero columns and rows.
// ATTENTION! only the lower triangle matrix is valid
void SMSolver2DP::PruneMatrix(MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld)
{    
    /////// Record the zero rows
    for( int i = 0; i < (int)matrix.rows(); ++i )
    {
        int j = 0;
        // Check along the row till diagonal
        while(j < i)
        {
            if( 0 != matrix.row(i)[j])
            {
                j = -1;
                break;
            }
            ++j;
        }
        // If not all zero go and check the next row (column)
        if(j < 0)
            continue;

        // Check along the column till matrix edge
        while( j < (int)matrix.rows() )
        {
            if( 0 != matrix.col(i)[j])
            {
                j = -1;
                break;
            }
            ++j;
        }
        // If not all zero go and check the next row (column)
        if(j < 0)
            continue;

        remainIndex[i] = false;
    }

    /////// Shrink the matrix
    // checking the "remain index" array. If 1 than move to top, else keep looking
    int offset = 0;
    int newMatrixDim = 0;
    for(int i = 0, j = 0; i < (int)remainIndex.size(); ++i)
    {
        // This row (column) remains, move it to the top
        if(true == remainIndex[i])
        {
            if(offset > 0)
            {
                matrix.row(j-offset) = matrix.row(j);
                matrix.col(j-offset) = matrix.col(j);
            }

            ++newMatrixDim;
        }

        // This row (column) is excluded in this iteration.
        if(true == remainIndexOld[i])
        {
            if(false == remainIndex[i])
                ++offset;
            ++j;
        }

        // We have finished checking all the rows (columns) for the remaining matrix
        if(j == matrix.cols())
            break;
    }

    // Trunk the matrix
    matrix.conservativeResize(newMatrixDim, newMatrixDim);

    // Update old_remain_index
    remainIndexOld = remainIndex;
}

// Pad the eigenvector back to original length
void SMSolver2DP::PadEigenvec(const VectorXd & prunedVec, std::vector<double> & paddedEigenvec, const std::vector<bool> & remainIndex)
{
    // Clear existing padded Eigenvector --- set to all 0
    //std::fill(paddedEigenvec.begin(), paddedEigenvec.end(), 0);

    int j = 0;
    for(int i = 0; i < (int)paddedEigenvec.size(); ++i)
    {
        if(true == remainIndex[i])
        {
            paddedEigenvec[i] = prunedVec[j];
            ++j;
        }
        else
            paddedEigenvec[i] = 0.0;

        //if( j == (int)prunedVec.size() )
        //    break;
    }

    //std::cout<<"size of prunedVec: "<<prunedVec.size()<<"\n";
    //std::cout<<"num of index: "<<remainIndex.size()<<"\n";
    //std::cout<<"num of remain index: "<<j<<"\n";
}

// Pick out conflict assignments. Set the corresponding element in remaining array to be 0.
// Here we need the scoreMatrix class mainly to make use of its indices computation function.
// This is not a good design. Redesign it later. ToDo.
void SMSolver2DP::ExcludeAssigns(int pickedIndex, std::vector<bool> & remainIndex, int moduleLen, int dataLen)
{
    // Pick out i,j (fi,ej)
    std::array<int, 2> indices = ObtainFeatureIndices(pickedIndex, moduleLen);

    // Exclude all the fi and fi'
    int start = moduleLen*indices[0];
    int end   = start + moduleLen;
    for( int i = start; i < end; ++i )
    {
        remainIndex[i] = false;                
    }

    // Exclude all the ej
    start = 0;
    end   = dataLen;
    for( int i = start; i < end; ++i )
    {
        remainIndex[ i*moduleLen + indices[1] ] = false;
    }
}

// Add the contribution of the picked row back to the matrix and shrink the matrix
// by copying the remaining rows and columns to top left and resize the matrix
void SMSolver2DP::ShrinkMatrix(int largestEleIndex, MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld)
{   
    // Locate the row index of the picked assignment for shrinked matrix
    int picked = std::count( remainIndexOld.begin(), remainIndexOld.begin()+largestEleIndex, true );

    // Add the contribution back --- ATTENTION! upper triangle of the matrix is undefined
    for(int i = 0; i < picked; ++i)
        matrix.diagonal()[i] += 2.0*matrix.row(picked)[i];
    for(int i = picked; i < matrix.cols(); ++i)
        matrix.diagonal()[i] += 2.0*matrix.col(picked)[i];

    // Shrink the matrix
    // checking the "remain index" array. If 1 than move to top, else keep looking
    int offset = 0;
    int newMatrixDim = 0;
    for(int i = 0, j = 0; i < (int)remainIndex.size(); ++i)
    {
        // This row (column) remains, move it to the top
        if(true == remainIndex[i])
        {
            if(offset > 0)
            {
                matrix.row(j-offset) = matrix.row(j);
                matrix.col(j-offset) = matrix.col(j);
            }

            ++newMatrixDim;
        }

        // This row (column) is excluded in this iteration.
        if(true == remainIndexOld[i])
        {
            if(false == remainIndex[i])
                ++offset;
            ++j;
        }

        // We have finished checking all the rows (columns) for the remaining matrix
        if(j == matrix.cols())
            break;
    }

    // Trunk the matrix
    matrix.conservativeResize(newMatrixDim, newMatrixDim);

    // Only set rows and columns to 0 rather than removing them
    // Compare the performance with 

    ////// So far, we can not get right results in this way. It is slower anyway.
    //VectorXd zeros = VectorXd::Zero(matrix.cols());

    //for(int i = 0, j = 0; i < (int)remainIndex.size(); ++i)
    //{
    //    if( true == remainIndexOld[i] && false == remainIndex[i])
    //    {
    //        matrix.row(i) = zeros;
    //        matrix.col(i) = zeros;
    //    }
    //}

}

// Return the full matrix --- the member matrix is only half valid
MatrixXd SMSolver2DP::GetFullMatrix(const ScoreMatrix2DPointCloud & scoreMatrix)
{
    MatrixXd tempMatrix = scoreMatrix.GetMatrix();
    tempMatrix.triangularView<Eigen::Upper>() = tempMatrix.transpose();

    return tempMatrix;
}

// Return the indicator Vector
VectorXd SMSolver2DP::GetIndicatorVec(const ScoreMatrix2DPointCloud & scoreMatrix)
{
    const int mLen = scoreMatrix.GetModuleNum();
    const int dLen = scoreMatrix.GetDataNum();

    // Compute indicator vector based on matching truth
    VectorXd indicatorVec = VectorXd::Zero( scoreMatrix.GetDim() );

    for(int i = 0; i < dLen; ++i)
    {
        if(0 != m_matchResult[i])
            indicatorVec[ mLen*i + m_matchResult[i]-1 ] = 1.0;
    }

    return indicatorVec;
}