
// The score matrix represents assignment pairs of two sample sets,
// module and data. With SM, we want to find the mapping from elements
// in data to elements in module.
// The score matrix for SM is always symmetric. Each entry represents
// the score (mutual agreement) of two nodes. Each node represents one
// pair of match in module and data.

// Hang 10/10/2014
/*********************************************************************/

#include "ScoreMatrix.h"


ScoreMatrix::ScoreMatrix(void)
{
}


ScoreMatrix::~ScoreMatrix(void)
{
    // Should be unnecessary
    //if( 0 != m_scoreM )
    //    delete m_scoreM;
}
