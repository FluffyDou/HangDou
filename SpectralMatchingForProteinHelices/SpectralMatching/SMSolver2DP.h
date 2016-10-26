/*********************************************************************/
// SMSolver2DP.h
// This is a class holding solvers for Spectral Matching for 2D point
// cloud, including solving for eigenvectors, binarizing the result,
// computing indicator vector and computing matrix score sum (x'Mx)

// Hang 11/21/2014
/*********************************************************************/

#pragma once

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <utility>
#include <array>
#include <vector>
#include <algorithm>
#include "Eigen\Dense"
//#include "Eigen\Eigenvalues"

#include "ScoreMatrix2DPointCloud.h"
#include "hungarian.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;

class SMSolver2DP
{
public:

    SMSolver2DP(void);

    ~SMSolver2DP(void);

    /*** Setters ***/

    // Solve for the principle eigenvector
    void SolveForEigenvectors( const MatrixXd & matrix );

    // Use iterative picking to binarize the result
    // Add an parameter so that we can do one set to many sets matching. ToDo.
    void IterativePicking( const ScoreMatrix2DPointCloud & scoreMatrix );

    // Simple greedy method to binarize the result
    void SimpleGreedyPick( const ScoreMatrix2DPointCloud & scoreMatrix );

    // Hungarian method to binarize the result. ToDo.
    void HungarianPick( const ScoreMatrix2DPointCloud & scoreMatrix );

    // Compute the score sum for Score Function
    void ComputeFunScore( const ScoreMatrix2DPointCloud & scoreMatrix );

    /*** Getters ***/
    // Return principle eigenvector
    const VectorXd& GetPrinEigenVec() const { return m_prinEigenVec; }
    // Return the solver to retrieve other info
    const SelfAdjointEigenSolver<MatrixXd> & GetSolver() const { return m_solver; }

    // Return match result and flag
    const std::vector<int> & GetMatchResult() const { return m_matchResult; }
    //const std::vector<int> & GetMatchFlag()   const { return m_matchFlag; }

    // Return function score
    const double GetFunScore() const { return m_FunctionScore; }

    // Return the full matrix --- the member matrix is only half valid
    MatrixXd GetFullMatrix(const ScoreMatrix2DPointCloud & scoreMatrix);
    // Return the indicator Vector
    VectorXd GetIndicatorVec(const ScoreMatrix2DPointCloud & scoreMatrix);

private:

    // Prune all the zero rows and columns
    // Return an array containing relation between old and new indices
    void PruneMatrix(MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld);

    // Pad the eigenvector back to original length
    void PadEigenvec(const VectorXd & prunedVec, std::vector<double> & paddedEigenvec, const std::vector<bool> & remainIndex);

    // Pick out conflict assignments. Set the corresponding element in remaining array to be 0.
    void ExcludeAssigns(int pickedIndex, std::vector<bool> & remainIndex, int moduleLen, int dataLen);

    // Obtain feature indices (i, j) from the index in the indicator vector
    // Each element in the indicator vector represents a pair of assignment i,j.
    std::array<int, 2> ObtainFeatureIndices(int index, int moduleLen)
    {
        int i     = index / moduleLen;
        int j     = index % moduleLen;

        std::array<int, 2> indices = {i,j};
        return indices;
    }

    // Obtain the line direction flag of "fi" from the index in the indicator vector
    //int ObtainFeatureFlag(int index, int moduleLen)
    //{
    //    int iFlag = index % (2*moduleLen) / moduleLen;
    //    return iFlag;
    //}

    static bool IndexCompare(const std::pair<double, int> & a, const std::pair<double, int> & b) {return a.first > b.first;}

    // Add the contribution of the picked row back to the matrix and shrink the matrix
    // by copying the remaining rows and columns to top left and resize the matrix
    void ShrinkMatrix(int largestEleIndex, MatrixXd & matrix, std::vector<bool> & remainIndex, std::vector<bool> & remainIndexOld);

private:
    // Eigen solver
    SelfAdjointEigenSolver<MatrixXd> m_solver;

    // Principle eigenvector --- for simple binary
    VectorXd m_prinEigenVec;

    // Matrix size
    int m_matrixSize;

    // Matching result --- starts from 1. 0 means no matching
    std::vector<int> m_matchResult;

    // Matching result flag --- 0 and 1
    //std::vector<int> m_matchFlag;

    // Score function score
    double m_FunctionScore;
};

