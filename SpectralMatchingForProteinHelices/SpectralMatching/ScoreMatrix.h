/*********************************************************************/
// This is a class holding a score matrix. This is an abstract class.
// It only holds the matrix. Functions will be in its inherited class.
// The score matrix represents assignment pairs of two sample sets,
// module and data. With SM, we want to find the mapping from elements
// in data to elements in module.
// The score matrix for SM is always symmetric. Each entry represents
// the score (mutual agreement) of two nodes. Each node represents one
// pair of match in module and data.

// Hang 10/10/2014
/*********************************************************************/

#pragma once

#include <array>
#include "Eigen\Dense"

using Eigen::MatrixXd;

class ScoreMatrix
{
public:

    typedef std::array<int,4> ivector4;
    typedef std::array<int,2> ivector2;

    ScoreMatrix(void);

    virtual ~ScoreMatrix(void);

    /*** Setters ***/
    // Construct the score matrix
    virtual void GenScoreMatrix() = 0;

    /*** Getters ***/
    // Return the reference of score matrix
    const MatrixXd& GetMatrix() const { return m_scoreM; }

    // Return the matrix dimension
    const int GetDim() const { return m_dim; }

    // Return the matrix dimension
    const int GetModuleNum() const { return m_moduleNum; }

    // Return the matrix dimension
    const int GetDataNum() const { return m_dataNum; }

protected:

    // Fill in the entry of the score matrix
    virtual double ComputeEntry(int x, int y) = 0;

    // Given matrix entry index (x,y), return the corresponding element
    // indices (maybe other information), i,j and k,l.
    virtual ivector4 ComputeIndices(int x, int y) = 0;

    // The score function to use in ComputeEntry. It takes the distance
    // of two assignments and a kernel parameter and return the score
    // evaluated. The flag is used to choose among kernel types in the
    // subclasses.
    virtual double ScoreFun(double diff, double sigma, int flag) = 0;

protected:

    // Score matrix
    MatrixXd m_scoreM;
    // Matrix Dimension
    ivector2 m_matrixDim;

    // Dimension of the matrix --- row# and column#
    int m_dim;
    int m_moduleNum;
    int m_dataNum;
};

