/*********************************************************************/
// SMSolverDLD.cpp
// This is a class holding solvers for Spectral Matching for double line
// direction, including solving for eigenvectors, binarizing the result,
// computing indicator vector and computing matrix score sum (x'Mx)

// Hang 10/14/2014
/*********************************************************************/

#include "SMSolverDLD.h"


SMSolverDLD::SMSolverDLD(void)
{
    m_matrixSize    = 0;
    m_FunctionScore = 0.0;
}


SMSolverDLD::~SMSolverDLD(void)
{
}

// Compute the score sum for Score Function
void SMSolverDLD::ComputeFunScore( const ScoreMatrixDoubleLineDir & scoreMatrix )
{
    const int mLen = scoreMatrix.GetModuleNum();
    const int dLen = scoreMatrix.GetDataNum();

    // Compute indicator vector based on matching truth
    VectorXd indicatorVec = VectorXd::Zero( scoreMatrix.GetDim() );

    for(int i = 0; i < dLen; ++i)
    {
        if(0 != m_matchResult[i])
            indicatorVec[ 2*mLen*i + m_matchResult[i]-1 + m_matchFlag[i]*mLen ] = 1.0;
    }

    // Compute x'Mx
    MatrixXd tempMatrix = scoreMatrix.GetMatrix();
    tempMatrix.triangularView<Eigen::Upper>() = tempMatrix.transpose();

    m_FunctionScore  = indicatorVec.transpose() * tempMatrix * indicatorVec;
    // Return the normalized function score: (x'Mx)/(x'x)
    //m_FunctionScore /=  indicatorVec.dot(indicatorVec);
}

// Solve for the principle eigenvector
void SMSolverDLD::SolveForEigenvectors( const MatrixXd & matrix )
{
    m_solver.compute(matrix);
    // Record the matrix size
    m_matrixSize = matrix.cols();
    // Record the principle eigenvector
    m_prinEigenVec = m_solver.eigenvectors().col(m_matrixSize-1);
}


// Do Simple Greedy on double line direction score matrix
void SMSolverDLD::SimpleGreedyPick( const ScoreMatrixDoubleLineDir & scoreMatrix )
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    m_matchFlag   = m_matchResult;

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
    // Note: a better way to do this greedy is by hash map
    int resultLen = moduleLen < dataLen ? moduleLen : dataLen;
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
        // The matching index starts from 1 --- the same as Mathematica
        m_matchResult[indices[0]] = indices[1]+1;
        int tempFlag = ObtainFeatureFlag(acceptIndex, moduleLen);
        m_matchFlag[indices[0]] = tempFlag;

        // Exclude conflict assignments
        // Exclude all the fi and fi'
        int start = 2*moduleLen*indices[0];
        int end   = start + 2*moduleLen;
        for( int h = start; h < end; ++h )
        {
            paddedEigenvec[h] = 0;                
        }

        // Exclude all the ej
        start = 0;
        end   = 2*dataLen;
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
void SMSolverDLD::IterativePicking(const ScoreMatrixDoubleLineDir & scoreMatrix)
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    m_matchFlag   = m_matchResult;

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
        // NOTE: So far we don't take equal value into consideration (the sort is not guaranteed to be stationary).
        int largestEleIndex = std::distance(paddedEigenvec.begin(), std::max_element(paddedEigenvec.begin(), paddedEigenvec.end()));
        //std::cout<<"(+) picked index: "<<largestEleIndex+1<<"\n";

        // Append the picked assignment to the result
        std::array<int, 2> tempResult = ObtainFeatureIndices(largestEleIndex, moduleLen);
        m_matchResult[tempResult[0]] = tempResult[1]+1;
        int tempFlag = ObtainFeatureFlag(largestEleIndex, moduleLen);
        m_matchFlag[tempResult[0]] = tempFlag;

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

// Try picking more than one assignment each iteration
void SMSolverDLD::IterativePickingNewStop(const ScoreMatrixDoubleLineDir & scoreMatrix, double eigenRatio)
{
    // Make a copy of the score matrix since we need to prune it. 
    MatrixXd matrix = scoreMatrix.GetMatrix();

    // Initialize the result
    const int moduleLen = scoreMatrix.GetModuleNum();
    const int dataLen   = scoreMatrix.GetDataNum();

    // The match result is the same in Mathematica --- index starts from 1
    m_matchResult = std::vector<int>(dataLen, 0);
    m_matchFlag   = m_matchResult;

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

    //int iterationRound = 0;
    // While the matrix size (remain matrix) is still larger than 0

    //// Initialize the matlab library
    //if( !LeadingEigenvecInitialize() )
    //{
    //    std::cerr << "Could not initialize the library properly."<< std::endl;
    //    return;
    //}

    while( matrix.cols() > 0 )
    {
        //++iterationRound;

        // Solve the matrix (remain matrix) for principle eigenvector
        VectorXd prunedVec = m_solver.compute(matrix).eigenvectors().col(matrix.cols()-1);
        // Get absolute value of the principle eigenvector
        prunedVec = prunedVec.cwiseAbs();

        /***********************************************************/
        /*************** Matlab and MKL Hybrid *********************/
        /***********************************************************/

        //VectorXd prunedVec;
        //if( matrix.cols() < 400 ) // Use Eigen / mkl to solve
        //{
        //    // Solve the matrix (remain matrix) for principle eigenvector
        //    prunedVec = m_solver.compute(matrix).eigenvectors().col(matrix.cols()-1);
        //    // Get absolute value of the principle eigenvector
        //    prunedVec = prunedVec.cwiseAbs();
        //}
        //else // Use Matlab kernel to solve for leading eigenvector
        //{
        //    // Resize the eigenvector
        //    prunedVec.resize(matrix.cols());

        //    // Input matrix
        //    mwArray matrix_mat(matrix.cols(), matrix.cols(), mxDOUBLE_CLASS, mxREAL);

        //    // Output
        //    mwArray out;

        //    matrix_mat.SetData(matrix.data(), matrix.cols()*matrix.cols());

        //    // Matlab kernel computation --- 1 means one output
        //    LeadingEigenVec(1, out, matrix_mat);

        //    out.GetData(prunedVec.data(), matrix.cols());
        //}

        /***********************************************************/
        /********************** Matlab Only ************************/
        /***********************************************************/

        //// Resize the eigenvector
        //VectorXd prunedVec(matrix.cols());

        //// Input matrix
        //mwArray matrix_mat(matrix.cols(), matrix.cols(), mxDOUBLE_CLASS, mxREAL);

        //// Output
        //mwArray out;

        //matrix_mat.SetData(matrix.data(), matrix.cols()*matrix.cols());

        //// Matlab kernel computation --- 1 means one output
        //LeadingEigenVec(1, out, matrix_mat);

        //out.GetData(prunedVec.data(), matrix.cols());

        /***********************************************************/
        /***********************************************************/
        /***********************************************************/

        //std::cout<<"matrix size: "<<matrix.cols()<<"\n";
        //std::cout<<"eigenvector sum: "<<prunedVec.sum()<<"\n";

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

        // Keep picking until we meet stop criteria
        // In case dataLen is larger than moduleLen
        // NOTE: a better way to do this greedy is by hash map
        int resultLen = moduleLen < dataLen ? moduleLen : dataLen;
        int pickedNum = 0;
        int assignNumThreshold = resultLen;
        double eigenLowerBound = eigenRatio * finalEigenvec[0].first;

        for(int i = 0; i < (int)finalEigenvec.size() && pickedNum < assignNumThreshold; ++i)
        {
            // If the largest remaining element is zero, we can stop
            if(0 == finalEigenvec[i].first || pickedNum == resultLen || eigenLowerBound > finalEigenvec[i].first)
                break;

            ///////// Old hybrid picking
            //// This one has been excluded --- the next largest is in conflict with accepted assignments
            //if(0 == paddedEigenvec[finalEigenvec[i].second])
            //    break;

            /////// New hybrid picking. 
            // This one has been excluded --- the next largest is in conflict with accepted assignments
            if(0 == paddedEigenvec[finalEigenvec[i].second])
                continue;

            // The remaining largest element
            int acceptIndex = finalEigenvec[i].second;

            // Append the picked assignment to the result
            std::array<int, 2> indices = ObtainFeatureIndices(acceptIndex, moduleLen);
            //m_matchResult[indices[0]] = indices[1]+1;

            // The matching index starts from 1 --- the same as Mathematica
            m_matchResult[indices[0]] = indices[1]+1;
            int tempFlag = ObtainFeatureFlag(acceptIndex, moduleLen);
            m_matchFlag[indices[0]] = tempFlag;

            // Update the remaining rows --- add zero rows
            ExcludeAssigns(acceptIndex, remainIndex, scoreMatrix.GetModuleNum(), scoreMatrix.GetDataNum());

            // Shrink the matrix
            ShrinkMatrix(acceptIndex, matrix, remainIndex, remainIndexOld);

            // To locate rows (columns) index, for matrix shrink
            remainIndexOld = remainIndex;

            // Exclude conflict assignments
            // Exclude all the fi and fi'
            int start = 2*moduleLen*indices[0];
            int end   = start + 2*moduleLen;
            for( int h = start; h < end; ++h )
            {
                paddedEigenvec[h] = 0;                
            }

            // Exclude all the ej
            start = 0;
            end   = 2*dataLen;
            for( int h = start; h < end; ++h )
            {
                paddedEigenvec[ h*moduleLen + indices[1] ] = 0.0;
            }

            ++pickedNum;
        }
    }
    
    //// End of matlab kernel
    //LeadingEigenvecTerminate();

    //// Output information
    //std::cout<<"Actual iteration rounds: "<<iterationRound<<"/"<<std::min(moduleLen, dataLen)<<"\n";

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
void SMSolverDLD::PruneMatrix(MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld)
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

//
void SMSolverDLD::PadEigenvec(const VectorXd & prunedVec, std::vector<double> & paddedEigenvec, const std::vector<bool> & remainIndex)
{
    // Clear existing padded Eigenvector --- set to all 0
    // Not necessary
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
void SMSolverDLD::ExcludeAssigns(int pickedIndex, std::vector<bool> & remainIndex, int moduleLen, int dataLen)
{
    // Pick out i,j (fi,ej)
    std::array<int, 2> indices = ObtainFeatureIndices(pickedIndex, moduleLen);

    // Exclude all the fi and fi'
    int start = 2*moduleLen*indices[0];
    int end   = start + 2*moduleLen;
    for( int i = start; i < end; ++i )
    {
        remainIndex[i] = false;                
    }

    // Exclude all the ej
    start = 0;
    end   = 2*dataLen;
    for( int i = start; i < end; ++i )
    {
        remainIndex[ i*moduleLen + indices[1] ] = false;
    }
}

// Add the contribution of the picked row back to the matrix and shrink the matrix
// by copying the remaining rows and columns to top left and resize the matrix
void SMSolverDLD::ShrinkMatrix(int largestEleIndex, MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld)
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
MatrixXd SMSolverDLD::GetFullMatrix(const ScoreMatrixDoubleLineDir & scoreMatrix)
{
    MatrixXd tempMatrix = scoreMatrix.GetMatrix();
    tempMatrix.triangularView<Eigen::Upper>() = tempMatrix.transpose();

    return tempMatrix;
}

// Return the indicator Vector
VectorXd SMSolverDLD::GetIndicatorVec(const ScoreMatrixDoubleLineDir & scoreMatrix)
{
    const int mLen = scoreMatrix.GetModuleNum();
    const int dLen = scoreMatrix.GetDataNum();

    // Compute indicator vector based on matching truth
    VectorXd indicatorVec = VectorXd::Zero( scoreMatrix.GetDim() );

    for(int i = 0; i < dLen; ++i)
    {
        if(0 != m_matchResult[i])
            indicatorVec[ 2*mLen*i + m_matchResult[i]-1 + m_matchFlag[i]*mLen ] = 1.0;
    }

    return indicatorVec;
}