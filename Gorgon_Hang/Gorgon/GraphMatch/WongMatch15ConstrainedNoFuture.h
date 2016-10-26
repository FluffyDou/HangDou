// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Graph Matching algorithm based on the following paper:
//				  Wong et al., An algorithm for graph optimal monomorphism, 
//				  IEEE Transactions on Systems, Man, and Cybernetics, Vol 20, No 3, May/June 1990, pp628-636
//				  Modified to support missing helixes and constrained nodes
//				  Contains no future function



#ifndef WONGMATCH15CONSTRAINEDNOFUTURE_H
#define WONGMATCH15CONSTRAINEDNOFUTURE_H

#include "StandardGraph.h"
#include "LinkedNode.h"
#include "NodeList.h"
#include "GlobalConstants.h"
#include "PathGenerator.h"
#include "Structures.h"
#include <ctime>
//#include <SkeletonMaker/PriorityQueue.h>
#include <Foundation/GorgonPriorityQueue.h>
#include "SSECorrespondenceResult.h"

#include <fstream>
#include <sstream>

#include "Eigen/Sparse"
#include <iomanip>

#include "mat.h"
#include "LeadingEigenVecForSparseM.h"

using namespace wustl_mm::Foundation;
namespace wustl_mm {
	namespace GraphMatch {

		class WongMatch15ConstrainedNoFuture{
		public:
			StandardGraph * patternGraph;
			StandardGraph * baseGraph;
		public:
			WongMatch15ConstrainedNoFuture(StandardGraph * patternGraph, StandardGraph * baseGraph);
			WongMatch15ConstrainedNoFuture(StandardGraph * patternGraph, StandardGraph * baseGraph, int missingHelixCount, int missingSheetCount);
			~WongMatch15ConstrainedNoFuture();
			int RunMatching(clock_t startTime);
			SSECorrespondenceResult GetResult(int rank);
			void SaveResults();


		private:
		#ifdef VERBOSE
			clock_t timeInGetA;
			clock_t timeInGetB;
			clock_t timeInQueue;
		#endif
			LinkedNode * currentNode;
			//PriorityQueue<LinkedNode, double> * queue;
			GorgonPriorityQueue<double, LinkedNode *> * queue;
			vector<LinkedNodeStub*> usedNodes;
			vector<SSECorrespondenceResult> solutions;
			int missingHelixCount;
			int missingSheetCount;
			int expandCount;
			int foundCount;
			int longestMatch;
			PathGenerator * pathGenerator;

		private:
			void Init(StandardGraph * patternGraph, StandardGraph * baseGraph);
			double GetC(int p, int qp);
			double GetC(int j, int p, int qj, int qp);
			double GetC(int p, int qp, LinkedNodeStub * currentNode);
			double GetCost(int d, int m, int qj, int qp, bool debugMsg);
			double GetPenaltyCost(int d, int m, bool debugMsg);
			double GetCPrime(int a, int b, int c, int d);
			double GetK(int p, int qp);
			double GetKPrime(int i, int q);
			double GetF();
			void PopBestNode(); // Gets the best (first) node from the active nodes list.
			bool ExpandNode(LinkedNodeStub * currentStub);  // Expands all the children of the current node.
			void ComputeSolutionCost(int solution[], bool extraMessages);
			void AnalyzeResults(int results[][MAX_NODES], int groundTruth[]);
			void NormalizeGraphs();
			void NormalizeSheets();
			unsigned long long EncodeNode(unsigned long long bitmap, int node);
			void PrintNodeConcise(LinkedNode * node, int rank, bool endOfLine, bool printCostBreakdown);
            int bestMatches[RESULT_COUNT][MAX_NODES];

            typedef struct MatchPair{
                int p;
                int q;
                MatchPair(){}
                MatchPair(int pi, int qi) : p(pi), q(qi){}

            } MatchPair;

            // Use spectral matching to replace the searching method
            void   AddAugmentEdges(StandardGraph * const patternGraph, const int & missHelicesInVolumeMap);
            void   SpectralMatch(StandardGraph * const patternGraph, StandardGraph * const baseGraph, std::vector<MatchPair> & simgpleGreedyResult, std::vector<std::vector<MatchPair>> & orderGreedyResults, Eigen::SparseMatrix<double, Eigen::ColMajor> & affinityMatrix);
            void   BuildAffinityScoreMatrix(StandardGraph * const patternGraph, StandardGraph * const baseGraph, std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, const double & tempSigma);
            void   SolveSMSimpleGreedyByMatlab(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, std::vector<MatchPair> & matchResult, const int & seqNodeNum, const int & volNodeNum);            
            void   SolveSMOrderGreedyByMatlab(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, std::vector< std::vector<MatchPair> > & matchResults, const int & seqNodeNum, const int & volNodeNum, const int & t1);
            void   SolveSMIteratively(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue);
            double GetAffinityScore(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const int & pi, const int & pj, const int & qi, const int & qj);
            double GetEdgePairCost(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const int & pi, const int & pj, const int & qi, const int & qj);
            double GaussianKernel(double const & value, double const & sigma) { return exp(-value*value / sigma); }

            void   PrintMatchResult(const std::vector<MatchPair> & matchResult, int seqNodeNum, double cost, double score, int acc, int trueAcc, int rank);
            void   ObtainBestNodeMatch(const std::vector<std::vector<MatchPair>> & results, const std::vector<MatchPair> & trueResult, std::vector<MatchPair> & bestResult, int & bestAcc, int & bestRank);
            void   ObtainBestOrderMatch(const std::vector<std::vector<MatchPair>> & results, const std::vector<MatchPair> & trueResult, std::vector<MatchPair> & bestResult, int & bestAcc, int & bestRank);

            double ComputeMatchCost(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const std::vector<MatchPair> & matchResult);
            double ComputeMatchScore(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const std::vector<MatchPair> & matchResult, const Eigen::SparseMatrix<double, Eigen::ColMajor> & affinityMatrix);

            static bool IndexCompareStdPairDIDescend(const std::pair<double, int> & a, const std::pair<double, int> & b) { return a.first > b.first; }
            static bool IndexCompareMatchPairIIAscend(const MatchPair & a, const MatchPair & b) { return a.p < b.p; }
		};

		WongMatch15ConstrainedNoFuture::WongMatch15ConstrainedNoFuture(StandardGraph * patternGraph, StandardGraph * baseGraph) {
			Init(patternGraph, baseGraph);
		}

		WongMatch15ConstrainedNoFuture::WongMatch15ConstrainedNoFuture(StandardGraph * patternGraph, StandardGraph * baseGraph, int missingHelixCount, int missingSheetCount) {
			Init(patternGraph, baseGraph);
			this->missingHelixCount = missingHelixCount;
			this->missingSheetCount = missingSheetCount;
		}

		WongMatch15ConstrainedNoFuture::~WongMatch15ConstrainedNoFuture() {
			for(unsigned int i = 0; i < usedNodes.size(); i++) {
				delete usedNodes[i];
			}
			usedNodes.clear();
			solutions.clear();
			//int queueSize = queue->getLength();
			//double tempKey;
			//LinkedNode * tempNode;
			//for(int i = 0; i < queueSize; i++) {		
			//	queue->remove(tempNode, tempKey);
			//	delete tempNode;
			//}
			LinkedNode * tempNode;
			while(!queue->IsEmpty()) {
				tempNode = queue->PopFirst();
				delete tempNode;
			}

			delete queue;
			
			delete pathGenerator;
		}

        // Deprecated
        // Suppose the graph only have helix nodes for now.
        // The graph index starts at 0, even nodes are helix start and odd nodes are helix end.
        void WongMatch15ConstrainedNoFuture::AddAugmentEdges(StandardGraph * const patternGraph, const int & missHelicesInVolumeMap)
        {
            for (int i = 1; i < patternGraph->GetNodeCount(); i += 2)
                for (int j = 1; j <= missHelicesInVolumeMap; ++j)
                {
                    // A helix consists of two nodes. Thus the skip step is 2
                    int otherEndIndex = i + 2*j+1;
                    if (otherEndIndex >= patternGraph->GetNodeCount())
                        break;
                    
                    patternGraph->adjacencyMatrix[i][otherEndIndex][0] = GRAPHNODE_HELIX;
                    patternGraph->adjacencyMatrix[otherEndIndex][i][0] = GRAPHNODE_HELIX;

                    double tempDistance = 0.0;
                    // Add the number of skipped helices between node i and otherEndIndex
                    for (int h = i; h < otherEndIndex; ++h)
                        tempDistance += patternGraph->adjacencyMatrix[h][h+1][1];

                    patternGraph->adjacencyMatrix[i][otherEndIndex][1] = tempDistance;
                    patternGraph->adjacencyMatrix[otherEndIndex][i][1] = tempDistance;
                }
        }

        // GS (sequence graph) node: p; GC (volume graph) node: q
        // The affinity score matrix has the row (column): 
        // (p1 q1), (p1, q2), ... ,(p1, qm), (p2 q1), (p2, q2), ... , (p2, qm), ...
        // So far we use matlab for sparse eigen decomposition
        // Note: patternGraph is sequence graph and baseGraph is volume graph
        void WongMatch15ConstrainedNoFuture::BuildAffinityScoreMatrix(StandardGraph * const patternGraph, StandardGraph * const baseGraph, std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, const double & tempSigma)
        {
            const double tmpEpsilon = 0.000001;
            //const double tempSigma  = 40.0;
            const double nodeScore  = 1.0;

            //int matrixDim = patternGraph->GetNodeCount() * baseGraph->GetNodeCount();

            /////// Todo: do not loop on dimension, loop on existing edges instead. (IMPORTANT)
            /////// Todo: put the diagonal elements in order in the end and record the start and end indices for iterative picking.
            /////// Todo: remove all zero rows and columns and record the padding indices
            //// For each matrix row
            //for (int i = 0; i < matrixDim; ++i)
            //    // For each matrix column
            //    for (int j = 0; j < matrixDim; ++j)
            //    {
            //        int pi = i / baseGraph->GetNodeCount(), pj = j / baseGraph->GetNodeCount();
            //        int qi = i % baseGraph->GetNodeCount(), qj = j % baseGraph->GetNodeCount();

            //        // Assign the affinity score for this two matches: score = e^(edge cost + node cost + volume graph missing penalty)
            //        double tempScore = GetAffinityScore(patternGraph, baseGraph, pi, pj, qi, qj);
            //        // Valid score (usually larger than 0)
            //        if (tempScore > tmpEpsilon)
            //        {
            //            rowIndex.push_back((double)i);
            //            colIndex.push_back((double)j);
            //            matValue.push_back(tempScore);
            //        }
            //    }

            struct graphEdges
            {
                int qi;
                int qj;
                double length;
                bool euclideanFlag;

                graphEdges() : qi(0), qj(0), length(0.0), euclideanFlag(false) {}
                graphEdges(int i_i, int i_j, double i_l, bool i_f) : qi(i_i), qj(i_j), length(i_l), euclideanFlag(i_f) {}
            };

            // Prepare two lists for loop edges and helix edges in volume graph. 
            // The edges are non-duplicated and index ascending: (1,2) not (2,1).
            std::vector<graphEdges> helixEdges, loopEdges;

            for (int i = 0; i < baseGraph->GetNodeCount(); ++i)
            for (int j = 0; j < i; ++j)
            {
                if ((int)(baseGraph->adjacencyMatrix[i][j][0]+0.01) == GRAPHEDGE_HELIX)
                    helixEdges.push_back(graphEdges(j, i, baseGraph->adjacencyMatrix[i][j][1], false));
                else if ((int)(baseGraph->adjacencyMatrix[i][j][0] + 0.01) == GRAPHEDGE_LOOP)
                    loopEdges.push_back(graphEdges(j, i, baseGraph->adjacencyMatrix[i][j][1], false));
                else if ((int)(baseGraph->adjacencyMatrix[i][j][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)
                    loopEdges.push_back(graphEdges(j, i, baseGraph->adjacencyMatrix[i][j][1], true));
            }


            for (int pi = 0; pi < patternGraph->GetNodeCount(); ++pi)
            {
                if (pi % 2 == 0) // Helix edge
                {
                    // The other end (node) of the helix
                    int pj = pi + 1;

                    // Compute all the affinity score matrx entries involving the current "loop" edge in sequence graph and all "loop" edges in volume graph
                    for (int k = 0; k < helixEdges.size(); ++k)
                    {
                        // The entry index in the sparse affinity score matrix
                        double i1 = pi*baseGraph->GetNodeCount() + helixEdges[k].qi;
                        double j1 = pj*baseGraph->GetNodeCount() + helixEdges[k].qj;
                        double i2 = pi*baseGraph->GetNodeCount() + helixEdges[k].qj;
                        double j2 = pj*baseGraph->GetNodeCount() + helixEdges[k].qi;
                        // Compute the entry score --- GetCost function take-in index starts from 1

                        double cost = GetCost(pi + 1, pj - pi, helixEdges[k].qi + 1, helixEdges[k].qj + 1, false);

                        // Should be a useless protection
                        if (cost < 0) continue;

                        cost += GetPenaltyCost(pi, pj - pi, false);
                        
                        // No sheet in concern. Node affinity should be always zero
                        // cost += GetC(pi+1, helixEdges[k].qi+1) + GetC(pj+1, helixEdges[k].qj+1);
                        double score = GaussianKernel(cost, tempSigma);

                        if (score > tmpEpsilon)
                        {
                            // For edge pair (pi, qi) and (pj, qj)
                            rowIndex.push_back(i1);
                            colIndex.push_back(j1);
                            matValue.push_back(score);
                            // We are building a sysmetric matrix
                            rowIndex.push_back(j1);
                            colIndex.push_back(i1);
                            matValue.push_back(score);

                            // For edge pair (pi, qj) and (pj, qi)
                            rowIndex.push_back(i2);
                            colIndex.push_back(j2);
                            matValue.push_back(score);
                            // We are building a sysmetric matrix
                            rowIndex.push_back(j2);
                            colIndex.push_back(i2);
                            matValue.push_back(score);
                        }
                    }
                }
                else // Loop edge
                {
                    for (int skippedHelixNum = 0; skippedHelixNum <= missingHelixCount; ++skippedHelixNum)                    
                    {
                        // A helix consists of two nodes. Thus the skip step is 2
                        int pj = pi + 2 * skippedHelixNum + 1;
                        if (pj >= patternGraph->GetNodeCount())
                            break;

                        // Compute all the affinity score matrx entries involving the current "loop" edge in sequence graph and all "loop" edges in volume graph
                        for (int k = 0; k < loopEdges.size(); ++k)
                        {
                            // The entry index in the sparse affinity score matrix
                            double i1 = pi*baseGraph->GetNodeCount() + loopEdges[k].qi;
                            double j1 = pj*baseGraph->GetNodeCount() + loopEdges[k].qj;
                            double i2 = pi*baseGraph->GetNodeCount() + loopEdges[k].qj;
                            double j2 = pj*baseGraph->GetNodeCount() + loopEdges[k].qi;

                            double cost = GetCost(pi + 1, pj - pi, loopEdges[k].qi + 1, loopEdges[k].qj + 1, false);

                            // Should be a useless protection
                            if (cost < 0) continue;

                            cost += GetPenaltyCost(pi, pj - pi, false);

                            // No sheet in concern. Node affinity should be always zero
                            // cost += GetC(pi+1, helixEdges[k].qi+1) + GetC(pj+1, helixEdges[k].qj+1);
                            double score = GaussianKernel(cost, tempSigma);

                            if (score > tmpEpsilon)
                            {
                                // For edge pair (pi, qi) and (pj, qj)
                                rowIndex.push_back(i1);
                                colIndex.push_back(j1);
                                matValue.push_back(score);
                                // We are building a sysmetric matrix
                                rowIndex.push_back(j1);
                                colIndex.push_back(i1);
                                matValue.push_back(score);

                                // For edge pair (pi, qj) and (pj, qi)
                                rowIndex.push_back(i2);
                                colIndex.push_back(j2);
                                matValue.push_back(score);
                                // We are building a sysmetric matrix
                                rowIndex.push_back(j2);
                                colIndex.push_back(i2);
                                matValue.push_back(score);
                            }
                        }
                    }
                }
            }

/***********************
            ///////////// This is Hang's implementation of getting the graph cost according to the paper /////////////////////
            ///////////// but the resulting cost "sometimes" differ a bit from Gorgon's original code    /////////////////////
            // Loop over "loop edge" pairs --- start from index 1
            for (int pi = 1; pi < patternGraph->GetNodeCount(); pi += 2)
            for (int skippedHelixNum = 0; skippedHelixNum <= missingHelixCount; ++skippedHelixNum)
            {
                // A helix consists of two nodes. Thus the skip step is 2
                int pj = pi + 2*skippedHelixNum + 1;
                if (pj >= patternGraph->GetNodeCount())
                    break;

                double sequenceEdgeLength    = 0.0;
                double missingHelicesPenalty = 0.0;
                // Compute the edge length in sequence graph
                for (int h = pi; h < pj; ++h)
                {
                    sequenceEdgeLength += patternGraph->adjacencyMatrix[h][h + 1][1];

                    // Compute the missing helices penalty
                    //if (h % 2 == 0)
                    //    missingHelicesPenalty += patternGraph->adjacencyMatrix[h][h + 1][1];
                }

                missingHelicesPenalty = GetPenaltyCost(pi, pj-pi, false);

                // Compute all the affinity score matrx entries involving the current "loop" edge in sequence graph and all "loop" edges in volume graph
                for (int k = 0; k < loopEdges.size(); ++k)
                {
                    double weight = LOOP_WEIGHT_COEFFICIENT;
                    weight = loopEdges[k].euclideanFlag ? weight * EUCLIDEAN_LOOP_PENALTY : weight;
                    double score = GaussianKernel(weight*abs(sequenceEdgeLength - loopEdges[k].length) + missingHelicesPenalty, tempSigma);

                    if (score > tmpEpsilon)
                    {
                        // The entry index in the sparse affinity score matrix
                        double i1 = pi*baseGraph->GetNodeCount() + loopEdges[k].qi;
                        double j1 = pj*baseGraph->GetNodeCount() + loopEdges[k].qj;
                        double i2 = pi*baseGraph->GetNodeCount() + loopEdges[k].qj;
                        double j2 = pj*baseGraph->GetNodeCount() + loopEdges[k].qi;

                        // For edge pair (pi, qi) and (pj, qj)
                        rowIndex.push_back(i1);
                        colIndex.push_back(j1);
                        matValue.push_back(score);
                        // We are building a sysmetric matrix
                        rowIndex.push_back(j1);
                        colIndex.push_back(i1);
                        matValue.push_back(score);

                        // For edge pair (pi, qj) and (pj, qi)
                        rowIndex.push_back(i2);
                        colIndex.push_back(j2);
                        matValue.push_back(score);
                        // We are building a sysmetric matrix
                        rowIndex.push_back(j2);
                        colIndex.push_back(i2);
                        matValue.push_back(score);
                    }
                }
            }

            // Loop over "helix edge" pairs
            for (int pi = 0; pi < patternGraph->GetNodeCount(); pi += 2)
            {
                // A helix consists of two nodes. Thus the skip step is 2
                int pj = pi + 1;

                // Compute the edge (helix) length in sequence graph
                double sequenceEdgeLength = patternGraph->adjacencyMatrix[pi][pj][1];

                // Compute all the affinity score matrx entries involving the current "loop" edge in sequence graph and all "loop" edges in volume graph
                for (int k = 0; k < helixEdges.size(); ++k)
                {
                    double weight = HELIX_WEIGHT_COEFFICIENT;
                    double score = GaussianKernel(weight*(sequenceEdgeLength - helixEdges[k].length), tempSigma);

                    if (score > tmpEpsilon)
                    {
                        // The entry index in the sparse affinity score matrix
                        double i1 = pi*baseGraph->GetNodeCount() + helixEdges[k].qi;
                        double j1 = pj*baseGraph->GetNodeCount() + helixEdges[k].qj;
                        double i2 = pi*baseGraph->GetNodeCount() + helixEdges[k].qj;
                        double j2 = pj*baseGraph->GetNodeCount() + helixEdges[k].qi;

                        // For edge pair (pi, qi) and (pj, qj)
                        rowIndex.push_back(i1);
                        colIndex.push_back(j1);
                        matValue.push_back(score);
                        // We are building a sysmetric matrix
                        rowIndex.push_back(j1);
                        colIndex.push_back(i1);
                        matValue.push_back(score);

                        // For edge pair (pi, qj) and (pj, qi)
                        rowIndex.push_back(i2);
                        colIndex.push_back(j2);
                        matValue.push_back(score);
                        // We are building a sysmetric matrix
                        rowIndex.push_back(j2);
                        colIndex.push_back(i2);
                        matValue.push_back(score);
                    }
                }
            }
****************/

            // For the current node-to-node correspondence, the diagonal of the affinity score matrix is all 1,
            // which is good since the matrix will usually be full rank.
            // ToDo: Replace the implementation with the insert function
            // Loop over node pairs
            for (int i = 0; i < patternGraph->GetNodeCount()*baseGraph->GetNodeCount(); ++i)
            {
                rowIndex.push_back((double)i);
                colIndex.push_back((double)i);
                matValue.push_back(nodeScore);
            }

            // NOTE: matlab can handle all zero rows and columns for eigen decomposition.
            // However, we might want to prune all zero rows and columns in the future
            // For this affinity matrix, the diagonal is always all 1 and thus not all zero rows and columns.
        }

        // Deprecated
        // The affinity scores are normalized to [0,1]
        double WongMatch15ConstrainedNoFuture::GetAffinityScore(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const int & pi, const int & pj, const int & qi, const int & qj)
        {
            const double constNodeScore = 1.0;

            // Affinity matrix diagonal --- only node score
            if (pi == pj && qi == qj)
                return constNodeScore;

            // Both edges in sequence and volume graph exit
            if (pi != pj && qi != qj && patternGraph->EdgeExists(pi, pj) && baseGraph->EdgeExists(qi, qj) && int(patternGraph->adjacencyMatrix[pi][qi][0] + 0.01) == int(baseGraph->adjacencyMatrix[pj][qj][0] + 0.01))
            {
                return GetEdgePairCost(patternGraph, baseGraph, pi, pj, qi, qj);
            }

            // Conflicting match or existing invalid edges for the two matches 
            return 0.0;
        }

        // Deprecated
        double WongMatch15ConstrainedNoFuture::GetEdgePairCost(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const int & pi, const int & pj, const int & qi, const int & qj)
        {
            return 0.0;
        }

        //
        void WongMatch15ConstrainedNoFuture::SolveSMSimpleGreedyByMatlab(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, std::vector<MatchPair> & matchResult, const int & seqNodeNum, const int & volNodeNum)
        {
            // Enable Matlab kernel
            if (!LeadingEigenVecForSparseMInitialize())
            {
                std::cerr << "Error calling LeadingEigenVecForSparseMInitialize..." << std::endl;
                return;
            }

clock_t matLabstartTime = clock();

            mwArray out;

            mwArray i(1, rowIndex.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray j(1, colIndex.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray s(1, matValue.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray n(1, 1, mxDOUBLE_CLASS, mxREAL);

            // Pass in row indices to Matlab
            i.SetData(rowIndex.data(), rowIndex.size());
            // Pass in column indices to Matlab
            j.SetData(colIndex.data(), colIndex.size());
            // Pass in matrix entry values to Matlab
            s.SetData(matValue.data(), matValue.size());

            // Pass in matrix dimension to Matlab
            //double dim = mLaplacianMatrix.rows();
            double dim = seqNodeNum*volNodeNum;
            n.SetData(&dim, 1);

            ////// Matlab kernel computation
            LeadingEigenVecForSparseM(1, out, i, j, s, n);

            // Grab the data out from Matlab
            //Eigen::VectorXd currEigenVec(currentRowsCols.size());
            std::vector<double> confidenceVector(dim);
            out.GetData(confidenceVector.data(), confidenceVector.size());

clock_t matlabElapesed = clock() - matLabstartTime;
std::cout << " Matrix sparsity: " << (float)matValue.size() / (float)(dim*dim) << endl;
std::cout << " Matlab uses (excluding matlab init) " << matlabElapesed / 1000 << "sec" << matlabElapesed % 1000 << "mm" << endl;

            std::cout << "Continuous solution from Matlab: \n";
            std::cout << std::fixed << std::setprecision(10);
            for (auto & it : confidenceVector)
            {
                std::cout << it << ", ";
            }
            std::cout << "\n";

            // Disable Matlab kernel
            LeadingEigenVecForSparseMTerminate();

            // Binarize the confidenceVector (continuous solution) --- sequence index start from 0
            // Very simple. Each time, pick the correspondence with highest score and the correspondence in the same helix. Reject conflicting ones and keep going.
            // Maybe set up a threshold to check if both correspondences are good enough. If not, record the score and keep checking.
            // Make a pair list of the eigenvector so that after sorting
            // We can still find each element's original index

clock_t greedyBinarizeStartTime = clock();

            // Container to store sorted eigenvector elements
            std::vector<std::pair<double, int>> finalEigenvec(confidenceVector.size());
            for (int i = 0; i < (int)finalEigenvec.size(); ++i)
            {
                // First value and Second index
                finalEigenvec[i].first = confidenceVector[i];
                finalEigenvec[i].second = i;
            }

            // Sort the principle eigenvector
            std::sort(finalEigenvec.begin(), finalEigenvec.end(), IndexCompareStdPairDIDescend);

            // In case dataLen is larger than moduleLen
            // Note: a better way to do this greedy is by hash map
            int resultLen = seqNodeNum < volNodeNum ? seqNodeNum : volNodeNum;
            int pickedNum = 0;
            // Loop to pick out the largest element and exclude conflicts
            for (int i = 0; i < dim; ++i)
            {
                // If the largest remaining element is zero, we can stop
                // finalEigenvec[i].first is always larger than or equal to zero
                if (0.0 == finalEigenvec[i].first || pickedNum == resultLen)
                    break;

                // This one has been excluded
                if (0.0 == confidenceVector[finalEigenvec[i].second])
                    continue;

                // The remaining largest element
                int acceptIndex = finalEigenvec[i].second;

                // Append the picked assignment to the result. Note we accept two correspondences (for one helix) at once.
                int pi = acceptIndex / volNodeNum;
                int qi = acceptIndex % volNodeNum;
                int pj = pi%2 == 0 ? pi+1 : pi-1;  
                int qj = qi%2 == 0 ? qi+1 : qi-1;

                // Check if this helix match is allowed
                if (GetCost(min(pi,pj)+1, abs(pj-pi), qi+1, qj+1, false) < 0)
                {
                    // Exclude both end nodes of this helix
                    confidenceVector[acceptIndex] = 0.0;
                    confidenceVector[pj*volNodeNum + qj] = 0.0;
                    continue;
                }

                matchResult.push_back(MatchPair(pi, qi));
                matchResult.push_back(MatchPair(pj, qj));

                // Exclude conflict assignments for pi and pj chunk
                for (int h = pi*volNodeNum; h < (pi + 1)*volNodeNum; ++h)
                {
                    confidenceVector[h] = 0.0;
                }
                for (int h = pj*volNodeNum; h < (pj + 1)*volNodeNum; ++h)
                {
                    confidenceVector[h] = 0.0;
                }

                // Exclude conflict assignments for qi and qj slots
                for (int h = qi; h < dim; h += volNodeNum)
                {
                    confidenceVector[h] = 0.0;
                }
                for (int h = qj; h < dim; h += volNodeNum)
                {
                    confidenceVector[h] = 0.0;
                }
                pickedNum += 2;
            }

clock_t binarizeElapesed = clock() - greedyBinarizeStartTime;
std::cout << "Binarization time: " << binarizeElapesed / 1000 << "sec" << binarizeElapesed % 1000 << "mm" << endl;

            // Print out the match result
            std::cout << "Simple binarization result: \n";
            std::sort(matchResult.begin(), matchResult.end(), IndexCompareMatchPairIIAscend);
            int tempI = 0;
            for (int j = 0; j < seqNodeNum; ++j)
            {
                int pi = matchResult[tempI].p;
                int qi = matchResult[tempI].q;
                if (pi != j || tempI == matchResult.size())
                {
                    std::cout << "-1 ";
                    continue;
                }
                std::cout << qi+1 << " ";
                ++tempI;
            }

            // Compute the match score
            std::cout << std::fixed << std::setprecision(10);
            std::vector<std::vector<int>> invalidConnection;
            double tempCost = 0.0;
            for (int i = 0; i < matchResult.size(); ++i)
            {
                // p represents sequence node, q represents volume node
                int pi = matchResult[i].p, qi = matchResult[i].q;

                if (i + 1 < matchResult.size())
                {
                    bool invalidFlag = true;
                    int piPlus1 = matchResult[i + 1].p, qiPlus1 = matchResult[i + 1].q;
                    if ((int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_HELIX || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)
                    {
                        double currCost = GetCost(pi + 1, piPlus1 - pi, qi + 1, qiPlus1 + 1, false);
                        invalidFlag = false;
                        if (currCost < 0)
                        {
                            invalidFlag = true;
                        }
                        else
                        {
                            tempCost += currCost;
                            tempCost += GetPenaltyCost(pi, piPlus1 - pi, false);
                        }
                    }

                    // The connection is invalid, record it
                    if (invalidFlag)
                    {
                        tempCost += 1000.0;

                        std::vector<int> tempPair;
                        tempPair.push_back(qi+1);
                        tempPair.push_back(qiPlus1+1);
                        invalidConnection.push_back(tempPair);
                    }
                }
                //if (GetC(pi + 1, qi + 1) > 0)
                //    std::cout << "Oops, GetC does return non-zero value: " << GetC(pi + 1, qi + 1) << "\n";

                tempCost += GetC(pi+1, qi+1);
            }
            // Add the penalty of missing the first helix
            if (!matchResult.empty() && matchResult.front().p != 0)
            {
                tempCost += GetPenaltyCost(0, matchResult.front().p, false);
            }
            // Add the penalty of missing the last helix
            if (!matchResult.empty() && matchResult.back().p != patternGraph->GetNodeCount() - 1)
            {
                tempCost += GetPenaltyCost(matchResult.back().p, patternGraph->GetNodeCount() - matchResult.back().p, false);
            }
            
            std::cout << tempCost << "\n";
            if (!invalidConnection.empty())
            {
                std::cout << "Invalid edge connections:";
                for (int i = 0; i < invalidConnection.size(); ++i)
                    std::cout << "(" << invalidConnection[i][0] << "," << invalidConnection[i][1] << ") ";
                std::cout << "\n";
            }
        }

        // 
        void WongMatch15ConstrainedNoFuture::SolveSMOrderGreedyByMatlab(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue, std::vector< std::vector<MatchPair> > & matchResults, const int & seqNodeNum, const int & volNodeNum, const int & t1)
        {
            // Enable Matlab kernel
            if (!LeadingEigenVecForSparseMInitialize())
            {
                std::cerr << "Error calling LeadingEigenVecForSparseMInitialize..." << std::endl;
                return;
            }

clock_t matLabstartTime = clock();

            mwArray out;

            mwArray i(1, rowIndex.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray j(1, colIndex.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray s(1, matValue.size(), mxDOUBLE_CLASS, mxREAL);
            mwArray n(1, 1, mxDOUBLE_CLASS, mxREAL);

            // Pass in row indices to Matlab
            i.SetData(rowIndex.data(), rowIndex.size());
            // Pass in column indices to Matlab
            j.SetData(colIndex.data(), colIndex.size());
            // Pass in matrix entry values to Matlab
            s.SetData(matValue.data(), matValue.size());

            // Pass in matrix dimension to Matlab
            //double dim = mLaplacianMatrix.rows();
            double dim = seqNodeNum*volNodeNum;
            n.SetData(&dim, 1);

            ////// Matlab kernel computation
            LeadingEigenVecForSparseM(1, out, i, j, s, n);

            // Grab the data out from Matlab
            //Eigen::VectorXd currEigenVec(currentRowsCols.size());
            std::vector<double> confidenceVector(dim);
            out.GetData(confidenceVector.data(), confidenceVector.size());

clock_t matlabElapesed = clock() - matLabstartTime;
std::cout << " Matrix sparsity: " << (float)matValue.size() / (float)(dim*dim) << endl;
std::cout << " Matlab uses (excluding matlab init) " << matlabElapesed / 1000 << "sec" << matlabElapesed % 1000 << "mm" << endl;

            //std::cout << "Continuous solution from Matlab: \n";
            //std::cout << std::fixed << std::setprecision(10);
            //for (auto & it : confidenceVector)
            //{
            //    std::cout << it << ", ";
            //}
            //std::cout << "\n";

            // Disable Matlab kernel
            LeadingEigenVecForSparseMTerminate();

            // Binarize the confidenceVector (continuous solution) --- sequence index start from 0
            // Very simple. Each time, pick the correspondence with highest score and the correspondence in the same helix. Reject conflicting ones and keep going.
            // Maybe set up a threshold to check if both correspondences are good enough. If not, record the score and keep checking.
            // Make a pair list of the eigenvector so that after sorting
            // We can still find each element's original index

clock_t greedyBinarizeStartTime = clock();

            int mostConfidentIndex = -1;
            double mostConfidentScore = t1;
            for (int i = 0; i < confidenceVector.size(); ++i)
            {
                if (confidenceVector[i] > mostConfidentScore)
                {
                    mostConfidentIndex = i;
                    mostConfidentScore = confidenceVector[i];
                }
            }

            // In case vol node number is larger than seq node number
            // Note: maybe a better way to do this greedy is by real hash map
            int resultLen = seqNodeNum < volNodeNum ? seqNodeNum : volNodeNum;

            // Find the highest score sequence node and push all its possible matches into root stack
            if (-1 == mostConfidentIndex)
            {
                std::cout << "Error: highest score of eigenvector is smaller than the threshold...\n";
                return;
            }

            std::vector<int> resultRootStack;
            int startIndex = (mostConfidentIndex / volNodeNum) * volNodeNum;
            int endIndex = startIndex + volNodeNum;
            // Push all the potential matches into root stack
            for (int i = startIndex; i < endIndex; ++i)
            {
                if (confidenceVector[i] > t1)
                    resultRootStack.push_back(i);
            }

            // Match result list
            //std::vector<std::vector<MatchPair>> matchResults;

            while (!resultRootStack.empty())
            {
                int currStartIndex = resultRootStack.back();
                resultRootStack.pop_back();

                int  pickedNum = 0;
                bool stopFlag  = false;

                std::vector<double> currConfidentVec = confidenceVector;

                // Append the picked assignment to the result. Note we accept two correspondences (for one helix) at once.
                int pi = currStartIndex / volNodeNum;
                int qi = currStartIndex % volNodeNum;
                int pj = pi % 2 == 0 ? pi + 1 : pi - 1;
                int qj = qi % 2 == 0 ? qi + 1 : qi - 1;

                // Check if this helix match is allowed
                if (GetCost(min(pi, pj) + 1, abs(pj - pi), qi + 1, qj + 1, false) < 0)
                {
                    currConfidentVec[currStartIndex] = 0.0;
                    currConfidentVec[pj*volNodeNum + qj] = 0.0;
                    continue;
                }

                std::vector<MatchPair> currMatchResult;
                // Accept current index
                currMatchResult.push_back(MatchPair(pi, qi));
                currMatchResult.push_back(MatchPair(pj, qj));
                pickedNum += 2;
                if (pickedNum == resultLen)
                {
                    matchResults.push_back(currMatchResult);
                    continue;
                }

                // Exclude conflict assignments for pi and pj chunk
                for (int h = pi*volNodeNum; h < (pi + 1)*volNodeNum; ++h)
                {
                    currConfidentVec[h] = 0.0;
                }
                for (int h = pj*volNodeNum; h < (pj + 1)*volNodeNum; ++h)
                {
                    currConfidentVec[h] = 0.0;
                }
                // Exclude conflict assignments for qi and qj slots
                for (int h = qi; h < dim; h += volNodeNum)
                {
                    currConfidentVec[h] = 0.0;
                }
                for (int h = qj; h < dim; h += volNodeNum)
                {
                    currConfidentVec[h] = 0.0;
                }

                // Forward picking
                int fStart      = pi > pj ? pi : pj;
                int lasVolNode  = pi > pj ? qi : qj;
                int lastSeqNode = fStart;
                // For each sequence node --- i means node in the sequence graph
                for (int i = fStart + 1; i < seqNodeNum; i += 2)
                {
                    int startIndex = i*volNodeNum;
                    int tempMaxIndex = -1;
                    double tempMaxScore = t1;
                    // For each possible match of that sequence node
                    for (int j = startIndex; j < startIndex+volNodeNum; ++j)
                    {
                        int currVolNode = j % volNodeNum;
                        int currSeqNode = j / volNodeNum;
                        // We find the match with the highest score and a valid path in the volume connecting it with its previous node
                        //if (currConfidentVec[j] > tempMaxScore && ((int)(baseGraph->adjacencyMatrix[lasVolNode][currVolNode][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[lasVolNode][currVolNode][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN) && (GetCost(lastSeqNode+1,currSeqNode-lastSeqNode,lasVolNode+1,currVolNode+1,false)>=0))
                        if (currConfidentVec[j] > tempMaxScore && GetCost(lastSeqNode + 1, currSeqNode - lastSeqNode, lasVolNode + 1, currVolNode + 1, false) >= 0)
                        {
                            tempMaxScore = currConfidentVec[j];
                            tempMaxIndex = j;
                        }
                        // Visited this match, exclude it
                        currConfidentVec[j] = 0.0;
                    }

                    // We do find a correspondence for successive sequence node
                    if (tempMaxIndex != -1)
                    {
                        // Accept the match
                        int pi = tempMaxIndex / volNodeNum;
                        int qi = tempMaxIndex % volNodeNum;
                        int pj = pi % 2 == 0 ? pi + 1 : pi - 1;
                        int qj = qi % 2 == 0 ? qi + 1 : qi - 1;

                        // Oops, invalid match due to helix length constraint
                        if (GetCost(min(pi, pj) + 1, abs(pi - pj), qi + 1, qj + 1, false) < 0)
                        {
                            currConfidentVec[pj*volNodeNum + qj] = 0.0;
                            continue;
                        }

                        currMatchResult.push_back(MatchPair(pi, qi));
                        currMatchResult.push_back(MatchPair(pj, qj));
                        pickedNum += 2;
                        if (pickedNum == resultLen)
                        {
                            stopFlag = true;
                            break;
                        }
                        // Exclude conflict matches
                        for (int h = qi; h < dim; h += volNodeNum)
                        {
                            currConfidentVec[h] = 0.0;
                        }
                        for (int h = qj; h < dim; h += volNodeNum)
                        {
                            currConfidentVec[h] = 0.0;
                        }

                        // Update last node for path checking
                        lasVolNode = pi > pj ? qi : qj;
                        lastSeqNode = pi > pj ? pi : pj;
                    }
                } // end for

                // We have found enough matches
                if (stopFlag)
                {
                    matchResults.push_back(currMatchResult);
                    continue;
                }

                // Backward picking
                int bStart  = pi < pj ? pi : pj;
                lasVolNode  = pi < pj ? qi : qj;
                lastSeqNode = bStart;
                for (int i = bStart - 1; i > 0; i -= 2)
                {
                    int startIndex = i*volNodeNum;
                    int tempMaxIndex = -1;
                    double tempMaxScore = t1;
                    // For each possible match of that sequence node
                    for (int j = startIndex; j < startIndex + volNodeNum; ++j)
                    {
                        int currVolNode = j%volNodeNum;
                        int currSeqNode = j/volNodeNum;
                        // We find the match with the highest score and a valid path in the volume connecting it with its previous node
                        //if (currConfidentVec[j] > tempMaxScore && ((int)(baseGraph->adjacencyMatrix[lasVolNode][currVolNode][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[lasVolNode][currVolNode][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN) && (GetCost(currSeqNode+1, lastSeqNode-currSeqNode, lasVolNode+1,currVolNode+1,false) >= 0))
                        if (currConfidentVec[j] > tempMaxScore && GetCost(currSeqNode + 1, lastSeqNode - currSeqNode, lasVolNode + 1, currVolNode + 1, false) >= 0)
                        {
                            tempMaxScore = currConfidentVec[j];
                            tempMaxIndex = j;
                        }
                        // Visited this match, exclude it
                        currConfidentVec[j] = 0.0;
                    }

                    // We do find the correspondence for successive sequence node
                    if (tempMaxIndex != -1)
                    {
                        // Accept the match
                        int pi = tempMaxIndex / volNodeNum;
                        int qi = tempMaxIndex % volNodeNum;
                        int pj = pi % 2 == 0 ? pi + 1 : pi - 1;
                        int qj = qi % 2 == 0 ? qi + 1 : qi - 1;

                        // Oops, invalid match due to helix length constraint
                        if (GetCost(min(pi, pj) + 1, abs(pi - pj), qi + 1, qj + 1, false) < 0)
                        {
                            currConfidentVec[pj*volNodeNum + qj] = 0.0;
                            continue;
                        }

                        currMatchResult.push_back(MatchPair(pi, qi));
                        currMatchResult.push_back(MatchPair(pj, qj));
                        pickedNum += 2;
                        if (pickedNum == resultLen)
                            break;

                        // Exclude conflict matches
                        for (int h = qi; h < dim; h += volNodeNum)
                        {
                            currConfidentVec[h] = 0.0;
                        }
                        for (int h = qj; h < dim; h += volNodeNum)
                        {
                            currConfidentVec[h] = 0.0;
                        }

                        // Update last node for path checking
                        lasVolNode  = pi < pj ? qi : qj;
                        lastSeqNode = pi < pj ? pi : pj;
                    }
                } // end for

                matchResults.push_back(currMatchResult);
            }// end while

clock_t binarizeElapesed = clock() - greedyBinarizeStartTime;
std::cout << "Order binarization time: " << binarizeElapesed / 1000 << "sec" << binarizeElapesed % 1000 << "mm" << endl;

            std::cout << "Order binarization results:\n";
            for (int i = 0; i < matchResults.size(); ++i)
            {
                std::vector<MatchPair> currResult = matchResults[i];
                if (currResult.empty()) continue;

                // Print out the match result
                std::sort(currResult.begin(), currResult.end(), IndexCompareMatchPairIIAscend);
                int tempI = 0;
                for (int j = 0; j < seqNodeNum; ++j)
                {
                    int pi = currResult[tempI].p;
                    int qi = currResult[tempI].q;
                    if (pi != j || tempI == currResult.size())
                    {
                        std::cout << "-1 ";
                        continue;
                    }
                    std::cout << qi+1 << " ";
                    ++tempI;
                }

                // Compute the match cost
                double tempCost = 0.0;                
                for (int i = 0; i < currResult.size(); ++i)
                {
                    // p represents sequence node, q represents volume node
                    int pi = currResult[i].p, qi = currResult[i].q;

                    if (i + 1 < currResult.size())
                    {
                        bool invalidFlag = true;
                        int piPlus1 = currResult[i + 1].p, qiPlus1 = currResult[i + 1].q;
                        if ((int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_HELIX || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)
                        {
                            double currCost = GetCost(pi + 1, piPlus1 - pi, qi + 1, qiPlus1 + 1, false);
                            invalidFlag = false;
                            if (currCost < 0)
                            {
                                invalidFlag = true;
                            }
                            else
                            {
                                tempCost += currCost;
                                tempCost += GetPenaltyCost(pi, piPlus1 - pi, false);
                            }
                        }

                        // The connection is invalid, record it
                        if (invalidFlag)
                        {
                            tempCost += 1000.0;
                        }
                    }
                    tempCost += GetC(pi + 1, qi + 1);
                }
                // Add the penalty of missing the first helix
                if (!currResult.empty() && currResult.front().p != 0)
                {
                    tempCost += GetPenaltyCost(0, currResult.front().p, false);
                }
                // Add the penalty of missing the last helix
                if (!currResult.empty() && currResult.back().p != patternGraph->GetNodeCount() - 1)
                {
                    tempCost += GetPenaltyCost(currResult.back().p, patternGraph->GetNodeCount() - currResult.back().p, false);
                }

                // Print out the match cost
                std::cout << tempCost << "\n";
            }

        }

        void WongMatch15ConstrainedNoFuture::SolveSMIteratively(std::vector<double> & rowIndex, std::vector<double> & colIndex, std::vector<double> & matValue)
        {}

        // GS (sequence graph) node: p; GC (volume graph) node: q
        // So far we use matlab for sparse eigen decomposition
        // Note: patternGraph is sequence graph and baseGraph is volume graph
        void WongMatch15ConstrainedNoFuture::SpectralMatch(StandardGraph * patternGraph, StandardGraph * baseGraph, std::vector<MatchPair> & simgpleGreedyResult, std::vector<std::vector<MatchPair>> & orderGreedyResults, Eigen::SparseMatrix<double, Eigen::ColMajor> & affinityMatrix)
        {
            /***** Add augment edges in the pattern graph based on the possible missing helices in the volume graph *****/
            // Deprecated
            //AddAugmentEdges(patternGraph, missingHelixCount);

            /***** Build the sparse affinity score matrix *****/
            std::vector<double> rowIndex, colIndex, matValue;

            const double tempSigma = 800.0; // 100 800
clock_t buildAffinityGraphStartTime = clock();
            BuildAffinityScoreMatrix(patternGraph, baseGraph, rowIndex, colIndex, matValue, tempSigma);
clock_t buildAffinityGraphElapse = clock() - buildAffinityGraphStartTime;
std::cout << " Time to build affinity matrix: " << buildAffinityGraphElapse / 1000 << "sec" << buildAffinityGraphElapse % 1000 << "mm" << endl;

            // Store the affinity matrix for score computing in the future
            std::vector< Eigen::Triplet<double> > tripletList;
            for (int i = 0; i < matValue.size(); ++i)
                tripletList.push_back(Eigen::Triplet<double>(rowIndex[i], colIndex[i], matValue[i]));

            affinityMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
            affinityMatrix.makeCompressed();

            //Eigen::SparseMatrix<double, Eigen::ColMajor> mat(patternGraph->GetNodeCount() * baseGraph->GetNodeCount(), patternGraph->GetNodeCount() * baseGraph->GetNodeCount());
            //for (int i = 0; i < rowIndex.size(); ++i)
            //{
            //    mat.insert(rowIndex[i], colIndex[i]) = matValue[i];
            //}
            //std::cout << "Affinity score matrix:\n" << mat << "\n";

            //for (auto & t : rowIndex)
            //{
            //    std::cout << t << ",";
            //}
            //std::cout << "\n";
            //for (auto & t : colIndex)
            //{
            //    std::cout << t << ",";
            //}
            //std::cout << "\n";
            //for (auto & t : matValue)
            //{
            //    std::cout << t << ",";
            //}
            //std::cout << "\n";

            /***** Solve for continuous solution *****/
            SolveSMSimpleGreedyByMatlab(rowIndex, colIndex, matValue, simgpleGreedyResult, patternGraph->GetNodeCount(), baseGraph->GetNodeCount());

            const double tempEpsilon = 0.000001;
            SolveSMOrderGreedyByMatlab(rowIndex, colIndex, matValue, orderGreedyResults, patternGraph->GetNodeCount(), baseGraph->GetNodeCount(), tempEpsilon);

            //////////////// Cost of ground truth --- by reading in a file named "match_truth.txt" //////////////
            //std::string inputGroundTruthFile = "match_truth.txt";
            //std::ifstream truthFile(inputGroundTruthFile);
            //std::string   line;

            //if (!truthFile.is_open())
            //    std::cout << "Unable to open match_truth.txt\n";
            //else
            //{
            //    // Only suppose to have one line
            //    std::getline(truthFile, line);
            //}
            //std::istringstream ss(line);
            //std::string token;
            //std::vector<MatchPair> trueResult;
            //int index = 0;
            //while (std::getline(ss, token, ' '))
            //{
            //    if (token.size() > 0)
            //    {
            //        int tempMatch = std::stoi(token);
            //        if (tempMatch > 0)
            //            trueResult.push_back(MatchPair(index, tempMatch-1));
            //    }
            //    ++index;
            //}

            //std::cout << std::fixed << std::setprecision(10);
            //std::vector<std::vector<int>> invalidConnection;
            //double tempCost = 0.0;            
            //for (int i = 0; i < trueResult.size(); ++i)
            //{
            //    // p represents sequence node, q represents volume node
            //    int pi = trueResult[i].p, qi = trueResult[i].q;

            //    if (i + 1 < trueResult.size())
            //    {
            //        bool invalidFlag = true;
            //        int piPlus1 = trueResult[i + 1].p, qiPlus1 = trueResult[i + 1].q;
            //        if ((int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_HELIX || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)
            //        {
            //            double currCost = GetCost(pi + 1, piPlus1 - pi, qi + 1, qiPlus1 + 1, false);
            //            invalidFlag = false;
            //            if (currCost < 0)
            //            {
            //                invalidFlag = true;
            //            }
            //            else
            //            {
            //                tempCost += currCost;
            //                tempCost += GetPenaltyCost(pi, piPlus1 - pi, false);
            //            }
            //        }

            //        // The connection is invalid, record it
            //        if (invalidFlag)
            //        {
            //            tempCost += 1000.0;

            //            std::vector<int> tempPair;
            //            tempPair.push_back(qi+1);
            //            tempPair.push_back(qiPlus1+1);
            //            invalidConnection.push_back(tempPair);
            //        }
            //    }
            //    if (GetC(pi + 1, qi + 1) > 0)
            //        std::cout << "Oops, GetC does return non-zero value: " << GetC(pi + 1, qi + 1) << "\n";

            //    tempCost += GetC(pi+1, qi+1);
            //}
            //// Add the penalty of missing the first helix
            //if (!trueResult.empty() && trueResult.front().p != 0)
            //{
            //    tempCost += GetPenaltyCost(0, trueResult.front().p, false);
            //}
            //// Add the penalty of missing the last helix
            //if (!trueResult.empty() && trueResult.back().p != patternGraph->GetNodeCount() - 1)
            //{
            //    tempCost += GetPenaltyCost(trueResult.back().p, patternGraph->GetNodeCount() - trueResult.back().p, false);
            //}
            //
            //std::cout << "Ground truth cost: " << tempCost << "\n";

            //if (!invalidConnection.empty())
            //{
            //    std::cout << "Invalid edge connections:";
            //    for (int i = 0; i < invalidConnection.size(); ++i)
            //        std::cout << "(" << invalidConnection[i][0] << "," << invalidConnection[i][1] << ") ";
            //    std::cout << "\n";
            //}
            


            /***** Feed the matching solution to Gorgon *****/
            // (The solutions may contains multiple matchings with different scores in the future)

            //std::cout << "Grid search for a good sigma:\n";
            //for (int i = 1; i < 200; ++i)
            //{
            //    std::vector<double> rowIndex, colIndex, matValue;
            //    std::vector<MatchPair> matchResult;

            //    const double tempSigma = 5.0;
            //    BuildAffinityScoreMatrix(patternGraph, baseGraph, rowIndex, colIndex, matValue, (double)i * tempSigma);
            //    //std::cout << "With sigma: " << (double)i * tempSigma << "\n";
            //    SolveSMSimpleGreedyByMatlab(rowIndex, colIndex, matValue, matchResult, patternGraph->GetNodeCount(), baseGraph->GetNodeCount());
            //}
            //std::cout << "Done with grid search:\n";

        }

        // Note: when counting correct "order", we count edges connecting the correct node order
        // E.g., 1 2 3 4 5 6 7 8 --> 1 (2 3) (4 5) (6 7) 8 : check if we have (2,3) (4,5) (6,7)
        void WongMatch15ConstrainedNoFuture::ObtainBestOrderMatch(const std::vector<std::vector<MatchPair>> & results, const std::vector<MatchPair> & trueResult, std::vector<MatchPair> & bestResult, int & bestAcc, int & bestRank)
        {
            for (int i = 0; i < results.size(); ++i)
            {
                int tempAcc = 0;
                // Convert the result into the compact version (no -1)
                std::vector<MatchPair> tempResult = results[i];
                std::sort(tempResult.begin(), tempResult.end(), IndexCompareMatchPairIIAscend);
                // Compute accuracy
                for (int h = 1; h < tempResult.size()-1; h += 2)
                for (int g = 1; g < trueResult.size()-1; g += 2)
                {
                    if (tempResult[h].q == trueResult[g].q && tempResult[h+1].q == trueResult[g+1].q)
                    {
                        ++tempAcc;
                        continue;
                    }
                }
                // Update best result
                if (tempAcc > bestAcc)
                {
                    bestAcc = tempAcc;
                    bestRank = i + 1;
                    bestResult = tempResult;
                }
            }
        }

        void WongMatch15ConstrainedNoFuture::ObtainBestNodeMatch(const std::vector<std::vector<MatchPair>> & results, const std::vector<MatchPair> & trueResult, std::vector<MatchPair> & bestResult, int & bestAcc, int & bestRank)
        {
            for (int i = 0; i < results.size(); ++i)
            {
                int tempAcc = 0;
                // Convert the result into the compact version (no -1)
                std::vector<MatchPair> tempResult = results[i];
                std::sort(tempResult.begin(), tempResult.end(), IndexCompareMatchPairIIAscend);
                // Compute accuracy
                for (int h = 0; h < tempResult.size(); ++h)
                for (int g = 0; g < trueResult.size(); ++g)
                {
                    if (tempResult[h].p == trueResult[g].p && tempResult[h].q == trueResult[g].q)
                    {
                        ++tempAcc;
                        continue;
                    }
                }
                // Update best result
                if (tempAcc > bestAcc)
                {
                    bestAcc = tempAcc;
                    bestRank = i + 1;
                    bestResult = tempResult;
                }
            }
        }

        void WongMatch15ConstrainedNoFuture::PrintMatchResult(const std::vector<MatchPair> & matchResult, int seqNodeNum, double cost, double score, int acc, int trueAcc, int rank)
        {
            int tempI = 0;
            for (int j = 0; j < seqNodeNum; ++j)
            {
                if (tempI >= matchResult.size())
                {
                        std::cout << "-1\t";
                        continue;
                }

                int pi = matchResult[tempI].p;
                int qi = matchResult[tempI].q;
                if (pi != j || tempI == matchResult.size())
                {
                    std::cout << "-1\t";
                    continue;
                }
                std::cout << qi + 1 << "\t";
                ++tempI;
            }

            std::cout << cost << "\t" << score << "\t" << acc << "\\" << trueAcc << "\t" << rank <<"\n";
        }

        double WongMatch15ConstrainedNoFuture::ComputeMatchScore(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const std::vector<MatchPair> & matchResult, const Eigen::SparseMatrix<double, Eigen::ColMajor> & affinityMatrix)
        {
            double tempScore = 0.0;
            int n = patternGraph->GetNodeCount();
            int m = baseGraph->GetNodeCount();
            // Build indicator vector
            Eigen::VectorXd indicatorVec = Eigen::VectorXd::Zero(n*m);
            for (int i = 0; i < matchResult.size(); ++i)
                indicatorVec[m*matchResult[i].p + matchResult[i].q] = 1.0;

            return indicatorVec.transpose() * affinityMatrix * indicatorVec;
        }

        double WongMatch15ConstrainedNoFuture::ComputeMatchCost(StandardGraph * const patternGraph, StandardGraph * const baseGraph, const std::vector<MatchPair> & matchResult)
        {
            std::vector<std::vector<int>> invalidConnection;
            double tempCost = 0.0;
            for (int i = 0; i < matchResult.size(); ++i)
            {
                // p represents sequence node, q represents volume node
                int pi = matchResult[i].p, qi = matchResult[i].q;

                if (i + 1 < matchResult.size())
                {
                    bool invalidFlag = true;
                    int piPlus1 = matchResult[i + 1].p, qiPlus1 = matchResult[i + 1].q;
                    if ((int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_HELIX || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP || (int)(baseGraph->adjacencyMatrix[qi][qiPlus1][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)
                    {
                        double currCost = GetCost(pi + 1, piPlus1 - pi, qi + 1, qiPlus1 + 1, false);
                        invalidFlag = false;
                        if (currCost < 0)
                        {
                            invalidFlag = true;
                        }
                        else
                        {
                            tempCost += currCost;
                            tempCost += GetPenaltyCost(pi, piPlus1 - pi, false);
                        }
                    }

                    // The connection is invalid, record it
                    if (invalidFlag)
                    {
                        tempCost += 1000.0;

                        std::vector<int> tempPair;
                        tempPair.push_back(qi + 1);
                        tempPair.push_back(qiPlus1 + 1);
                        invalidConnection.push_back(tempPair);
                    }
                }
                if (GetC(pi + 1, qi + 1) > 0)
                    std::cout << "Oops, GetC does return non-zero value: " << GetC(pi + 1, qi + 1) << "\n";

                tempCost += GetC(pi + 1, qi + 1);
            }
            // Add the penalty of missing the first helix
            if (!matchResult.empty() && matchResult.front().p != 0)
            {
                tempCost += GetPenaltyCost(0, matchResult.front().p, false);
            }
            // Add the penalty of missing the last helix
            if (!matchResult.empty() && matchResult.back().p != patternGraph->GetNodeCount() - 1)
            {
                tempCost += GetPenaltyCost(matchResult.back().p, patternGraph->GetNodeCount() - matchResult.back().p, false);
            }

            //std::cout << "cost: " << tempCost << "\n";
            if (!invalidConnection.empty())
            {
                std::cout << "Invalid edge connections:";
                for (int i = 0; i < invalidConnection.size(); ++i)
                    std::cout << "(" << invalidConnection[i][0] << "," << invalidConnection[i][1] << ") ";
                std::cout << "\n";
            }            

            return tempCost;
        }

		void WongMatch15ConstrainedNoFuture::Init(StandardGraph * patternGraph, StandardGraph * baseGraph) {	

#ifdef VERBOSE
			cout << "Initializing search" << endl;
			cout << "Base graph has " << baseGraph->GetHelixCount() << " helices and " << baseGraph->GetSheetCount() << " sheets." << endl;
			cout << "Pattern graph has " << patternGraph->GetHelixCount() << " helices and " << patternGraph->GetSheetCount() << " sheets." << endl;
#endif // VERBOSE
			usedNodes.clear();
#ifdef VERBOSE
			cout << "Creating priority queue" << endl;
#endif // VERBOSE
			//queue = new PriorityQueue<LinkedNode, double> (PRIORITYQUEUESIZE);
			queue = new GorgonPriorityQueue<double, LinkedNode *>(false);
#ifdef VERBOSE
			cout << "Loading pattern graph" << endl;
#endif // VERBOSE
			this->patternGraph = patternGraph;
#ifdef VERBOSE
			cout << "Loading base graph" << endl;
#endif // VERBOSE
			this->baseGraph = baseGraph;
			expandCount = 0;

#ifdef VERBOSE
			cout << "Finding the number of missing helices and sheets" << endl;
#endif // VERBOSE

            //cout << "init missingHelixCount " << missingHelixCount << "\n";

            // Deprecated
			//missingHelixCount = (patternGraph->GetNodeCount() - baseGraph->GetNodeCount()) / 2;
			//if(missingHelixCount < 0)  {
			//	missingHelixCount = 0;
			//}

			if(missingSheetCount < 0)  {
				missingSheetCount = 0;
			}
			
			// TEMPORARY CODE!
			//EUCLIDEAN_VOXEL_TO_PDB_RATIO *= 1.5;
			//EUCLIDEAN_VOXEL_TO_PDB_RATIO = 3.0;
			//EUCLIDEAN_VOXEL_TO_PDB_RATIO = 1.5 * LOOP_C_ALPHA_TO_ANGSTROMS; // this is what i think it should be.
			// END TEMPORARY CODE

			// new method of counting missing sheets and helices
			// count helix nodes in base graph
			int baseHelixNodes = 0;
			int baseSheetNodes = 0;
			int patternHelixNodes = 0;
			int patternSheetNodes = 0;            

            //NOTE!! Currently we don't take sheet into consideration
			for (int i = 0; i < baseGraph->GetNodeCount(); i++) {
				//cout << "base graph node " << i << " has type " << (int)(baseGraph->adjacencyMatrix[i][i][0]) << endl;
				switch((int)(baseGraph->adjacencyMatrix[i][i][0] + 0.01)) {
					case(GRAPHNODE_HELIX) : 
						baseHelixNodes++;
						break;
					case(GRAPHNODE_SHEET):
						baseSheetNodes++;
						break;
					default:
						break;
				}
			}

            // Get the missing helix number in the volume graph
            int volHelicesNum = baseGraph->GetHelixCount();
            for (int i = 0; i < baseGraph->skeletonHelixes.size(); ++i)
            {
                if (baseGraph->skeletonHelixes[i]->cornerCells.size() < 2)
                {
                    cout << "non-occupying-voxel helicx id (baseGraph->skeletonHelixes[i]): " << i << "\n";
                    --volHelicesNum;
                }
            }

            // May not need this protection
            if(volHelicesNum < 0)
                volHelicesNum = 0;

            cout << "non-occupying-voxel helicx number: " << baseGraph->GetHelixCount() - volHelicesNum << "\n";

#ifdef VERBOSE
			cout << "base graph has " << baseHelixNodes << " helix nodes and " << baseSheetNodes << " sheet nodes." << endl;
#endif // VERBOSE

            // Try to change so that those non-voxel occupied helices does not count in volume graph
			for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
				//cout << "pattern graph node " << i << " has type " << (int)(patternGraph->adjacencyMatrix[i][i][0]) << endl;
				switch((int)(patternGraph->adjacencyMatrix[i][i][0] + 0.01)) {
					case(GRAPHNODE_HELIX) : 
						patternHelixNodes++;
						break;
					case(GRAPHNODE_SHEET):
						patternSheetNodes++;
						break;
					default:
						break;
				}
			}
#ifdef VERBOSE
			cout << "pattern graph has " << patternHelixNodes << " helix nodes and " << patternSheetNodes << " sheet nodes." << endl;
#endif // VERBOSE

			//missingHelixCount = (patternHelixNodes - baseHelixNodes) / 2;
            missingHelixCount = (patternHelixNodes - 2*volHelicesNum) / 2;
            
            // Note: Do not check the missing helices box in Gorgon now. Will result in non-correspondence.
            //if (MISSING_HELIX_COUNT != -1)
            //    missingHelixCount = MISSING_HELIX_COUNT;

            if (missingHelixCount < 0)  
            {
                cout << "Error setting missing helix to be negative. Now setting it to zero.\n";
                missingHelixCount = 0;
            }

			// allow all strands to be missing
			//missingSheetCount = patternSheetNodes - baseSheetNodes;
			missingSheetCount = patternSheetNodes;

#ifdef VERBOSE
			cout << "missing helix count is " << missingHelixCount << ", missing sheet count is " << missingSheetCount << endl;
			cout << "missing helix penalty is " << MISSING_HELIX_PENALTY << ", missing sheet penalty is " << MISSING_SHEET_PENALTY << endl;
#endif // VERBOSE

			if(!PERFORMANCE_COMPARISON_MODE) {
				NormalizeGraphs();
				NormalizeSheets();
			}
			foundCount = 0;
			longestMatch = 0;
#ifdef VERBOSE
			timeInGetA = 0;
			timeInGetB = 0;
			timeInQueue = 0;
#endif

			// create and set up a new node to start the search
			currentNode = new LinkedNode();
			for(int j = 1; j <= patternGraph->nodeCount; j++) {
				LinkedNode::AddNodeToBitmap(currentNode->m1Bitmap, j);
			}
			for(int j = 1; j <= baseGraph->nodeCount; j++) {
				LinkedNode::AddNodeToBitmap(currentNode->m2Bitmap, j);
			}
			//queue->add(currentNode, currentNode->cost);
			queue->Add(currentNode->cost, currentNode);
			pathGenerator = new PathGenerator(baseGraph);
		}



		// searches for correspondences between the pattern graph and base graph.
		int WongMatch15ConstrainedNoFuture::RunMatching(clock_t startTime) {

            //// My test ////
            cout << "Aha! sequence graph matrix:\n";
            patternGraph->PrintGraph();

            cout << "Aha! volume graph matrix:\n";
            baseGraph->PrintGraph();

            std::vector<MatchPair> simgpleGreedyResult;
            std::vector< std::vector<MatchPair> > orderGreedyResults;
            std::vector< std::vector<MatchPair> > gorgonResults;

            int affinitMatirxDim = patternGraph->GetNodeCount() * baseGraph->GetNodeCount();
            Eigen::SparseMatrix<double, Eigen::ColMajor> affinityMatrix(affinitMatirxDim, affinitMatirxDim);

            SpectralMatch(patternGraph, baseGraph, simgpleGreedyResult, orderGreedyResults, affinityMatrix);

            /////////////// Analyze spectral result ///////////////////

            // Read in the ground truth
            std::string   inputGroundTruthFile = "match_truth.txt";
            std::ifstream truthFile(inputGroundTruthFile);
            std::string   line;
            if (!truthFile.is_open())
                std::cout << "Unable to open match_truth.txt\n";
            else
            {
                // Only suppose to have one line
                std::getline(truthFile, line);
            }
            std::istringstream ss(line);
            std::string token;
            std::vector<MatchPair> trueResult;
            int index = 0;
            while (std::getline(ss, token, ' '))
            {
                if (token.size() > 0)
                {
                    int tempMatch = std::stoi(token);
                    if (tempMatch > 0)
                        trueResult.push_back(MatchPair(index, tempMatch - 1));
                }
                ++index;
            }

            //// For order-greedy result
            int orderGreedyBestMatchAcc = -1;
            int orderGreedyBestMatchRank = -1;
            std::vector<MatchPair> orderGreedyBestMatchResult;
            ObtainBestNodeMatch(orderGreedyResults, trueResult, orderGreedyBestMatchResult, orderGreedyBestMatchAcc, orderGreedyBestMatchRank);

            //// For simple greedy result
            int simpleGreedyBestMatchAcc = 0;
            int simpleGreedyBestMatchRank = 1;
            // Convert the result into the compact version (no -1)
            std::sort(simgpleGreedyResult.begin(), simgpleGreedyResult.end(), IndexCompareMatchPairIIAscend);
            // Compute accuracy
            for (int h = 0; h < simgpleGreedyResult.size(); ++h)
            for (int g = 0; g < trueResult.size(); ++g)
            {
                if (simgpleGreedyResult[h].p == trueResult[g].p && simgpleGreedyResult[h].q == trueResult[g].q)
                {
                    ++simpleGreedyBestMatchAcc;
                    continue;
                }
            }
            int simpleGreedyBestOrderAcc = 0;
            int simpleGreedyBestOrderRank = 1;
            for (int h = 1; h < simgpleGreedyResult.size() - 1; h += 2)
            for (int g = 1; g < trueResult.size() - 1; g += 2)
            {
                if (simgpleGreedyResult[h].q == trueResult[g].q && simgpleGreedyResult[h + 1].q == trueResult[g + 1].q)
                {
                    ++simpleGreedyBestOrderAcc;
                    continue;
                }
            }

            // Compute the cost for ground truth, order-greedy and simple-greedy
            std::cout << "Compute truthMatch cost...\n";
            double truthCost = ComputeMatchCost(patternGraph, baseGraph, trueResult);
            std::cout << "Compute orderGreedyBestMatchCost cost...\n";
            double orderGreedyBestMatchCost = ComputeMatchCost(patternGraph, baseGraph, orderGreedyBestMatchResult);
            std::cout << "Compute simpleGreedyMatchCost cost...\n";
            double simpleGreedyMatchCost = ComputeMatchCost(patternGraph, baseGraph, simgpleGreedyResult);

            // Compute the score for ground truth, order-greedy, simple-greedy and gorgon-highest-cost result
            double truthScore = ComputeMatchScore(patternGraph, baseGraph, trueResult, affinityMatrix);
            double orderGreedyBestMatchScore = ComputeMatchScore(patternGraph, baseGraph, orderGreedyBestMatchResult, affinityMatrix);
            double simpleGreedyMatchScore = ComputeMatchScore(patternGraph, baseGraph, simgpleGreedyResult, affinityMatrix);

            // Display best correct-match-node results
            std::cout << "Results analysis:\n" << "Truth cost and score: " << truthCost << ", " << truthScore << "\n";
            std::cout << "Best correct-match-node results:\n";
            std::cout << "Truth, simple-geedy, order-greedy:\n";
            //std::cout << std::fixed << std::setprecision(10);
            PrintMatchResult(trueResult, patternGraph->GetNodeCount(), truthCost, truthScore, trueResult.size(), trueResult.size(), 1);
            PrintMatchResult(simgpleGreedyResult, patternGraph->GetNodeCount(), simpleGreedyMatchCost, simpleGreedyMatchScore, simpleGreedyBestMatchAcc, trueResult.size(), simpleGreedyBestMatchRank);
            PrintMatchResult(orderGreedyBestMatchResult, patternGraph->GetNodeCount(), orderGreedyBestMatchCost, orderGreedyBestMatchScore, orderGreedyBestMatchAcc, trueResult.size(), orderGreedyBestMatchRank);

            // Make sure we print in the console (before exit)
            flushall();

#ifdef VERBOSE
			cout << "Starting to search for correspondences." << endl;
			DisplayConstants();
#endif // VERBOSE
			bool continueLoop = true;
			clock_t finishTime;

clock_t myRecordTime = clock();

			// repeat the following loop until all results are found
			while(continueLoop)
			{
				PopBestNode();		
				if(currentNode == NULL) {
					break;
				}

				// if currentNode contains a complete sequence match, add it to the solutions list
				if(currentNode->depth == patternGraph->nodeCount) {
					finishTime = clock();
					foundCount++;
					//currentNode->PrintNodeConcise(foundCount, false);
					printf("\n");
					PrintNodeConcise(currentNode,foundCount, false, false);
					//printf(": (%d expanded) (%f seconds) (%fkB Memory) (%d queue size) (%d parent size)\n", expandCount, (double) (finishTime - startTime) / (double) CLOCKS_PER_SEC, (queue->getLength() * sizeof(LinkedNode) + usedNodes.size() * sizeof(LinkedNodeStub)) / 1024.0, queue->getLength(), (int)usedNodes.size());
					//printf(": (%d expanded) (%f seconds) (%d parent size)\n", expandCount, (double) (finishTime - startTime) / (double) CLOCKS_PER_SEC, (int)usedNodes.size());
                    //std::cout << "Expanded nodes so far: " << expandCount << endl;

#ifdef _WIN32
					flushall();
#endif
					int numHelices = baseGraph->GetHelixCount();
					solutions.push_back(SSECorrespondenceResult(currentNode, numHelices));

#ifdef MAKE_FINAL_MRC
					char fileName[80];
					sprintf(fileName, "Solution%d.mrc", foundCount);
					pathGenerator->GenerateGraph(currentNode, fileName);
#endif
				// otherwise, expand currentNode and adds its children to usedNodes
				} else {
					LinkedNodeStub * currentStub = new LinkedNodeStub(currentNode);
					if(ExpandNode(currentStub)) {
						usedNodes.push_back(currentStub);
					} else {
						delete currentStub;
					}
					delete currentNode;
				}		
				// continue until desired number of results are found
				continueLoop = (foundCount < RESULT_COUNT) && (!queue->IsEmpty());
			}

			//Cleaning up memory
			for(unsigned int i = 0; i < usedNodes.size(); i++) {
				delete usedNodes[i];
			}
			usedNodes.clear();

			//int queueSize = queue->getLength();
			//double tempKey;
			//LinkedNode * tempNode;
			//for(int i = 0; i < queueSize; i++) {		
			//	queue->remove(tempNode, tempKey);
			//	delete tempNode;
			//}
			LinkedNode * tempNode;
			while(!queue->IsEmpty()) {
				tempNode = queue->PopFirst();
				delete tempNode;
			}

clock_t searchElapesed = clock() - myRecordTime;
std::cout << "\nSearching time (my record): " << searchElapesed / 1000 << "sec" << searchElapesed % 1000 << "mm" << endl;
std::cout << "Expanded nodes in all: " << expandCount << endl;

#ifdef VERBOSE
			cout << "Finished the correspondence search. Found " << foundCount << " results." << endl;
#endif // VERBOSE

			/*vector<int> groundTruth;
			for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
				groundTruth[i]=SOLUTION[i];
			}
			ComputeSolutionCost(groundTruth);*/
			//cout << "The ground truth solution is" << endl;
			//cout << "**      ";
			//for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
			//	cout.width(2);
			//	cout << SOLUTION[i] << " ";
			//}
			////ComputeSolutionCost(SOLUTION, false);
			//ComputeSolutionCost(SOLUTION, false);
			//AnalyzeResults(bestMatches, SOLUTION);
			//cout << endl;


            /////// Analyze results from different methods /////////

            // Grab gorgon results from "bestMatches" array
            for (int i = 0; i < foundCount; ++i)
            {
                std::vector<MatchPair> tempResult;
                for (int j = 0; j < patternGraph->GetNodeCount(); ++j)
                {
                    if (bestMatches[i][j] > 0)
                        tempResult.push_back(MatchPair(j, bestMatches[i][j] - 1));
                }
                gorgonResults.push_back(tempResult);
            }

            int gorgonBestMatchAcc  = -1;
            int gorgonBestMatchRank = -1;
            int gorgonBestCostAcc   = -1;
            int gorgonBestCostRank  = -1;
            std::vector<MatchPair> gorgonBestMatchResult;
            std::vector< std::vector<MatchPair> > tempResultContainer;
            tempResultContainer.push_back(gorgonResults[0]);
            std::vector<MatchPair> gorgonBestCostResult;
            ObtainBestNodeMatch(gorgonResults, trueResult, gorgonBestMatchResult, gorgonBestMatchAcc, gorgonBestMatchRank);
            ObtainBestNodeMatch(tempResultContainer, trueResult, gorgonBestCostResult, gorgonBestCostAcc, gorgonBestCostRank);
            std::cout << "Compute gorgonBestMatch cost...\n";
            double gorgonBestMatchCost = ComputeMatchCost(patternGraph, baseGraph, gorgonBestMatchResult);
            std::cout << "Compute gorgonBestCost cost...\n";
            double gorgonBestCostCost = ComputeMatchCost(patternGraph, baseGraph, gorgonBestCostResult);
            std::cout << "Compute gorgonBestMatch score...\n";
            double gorgonBestMatchScore = ComputeMatchScore(patternGraph, baseGraph, gorgonBestMatchResult, affinityMatrix);
            std::cout << "Compute gorgonBestCost score...\n";
            double gorgonBestCostScore = ComputeMatchScore(patternGraph, baseGraph, gorgonBestCostResult, affinityMatrix);
            std::cout << "Truth, gorgon-best-match, gorgon-best-cost:\n";
            PrintMatchResult(trueResult, patternGraph->GetNodeCount(), truthCost, truthScore, trueResult.size(), trueResult.size(), 1);
            PrintMatchResult(gorgonBestMatchResult, patternGraph->GetNodeCount(), gorgonBestMatchCost, gorgonBestMatchScore, gorgonBestMatchAcc, trueResult.size(), gorgonBestMatchRank);
            PrintMatchResult(gorgonBestCostResult, patternGraph->GetNodeCount(), gorgonBestCostCost, gorgonBestCostScore, gorgonBestCostAcc, trueResult.size(), gorgonBestCostRank);
/*
            //////// A: Find the best correct-match-node result from Gorgon and order-greedy --- record accuracy, Gorgon rank and order-greedy dummy rank

            //// For gorgon result
            int gorgonBestMatchAcc  = -1;
            int gorgonBestMatchRank = -1;
            std::vector<MatchPair> gorgonBestMatchResult;
            ObtainBestNodeMatch(gorgonResults, trueResult, gorgonBestMatchResult, gorgonBestMatchAcc, gorgonBestMatchRank);

            //// For order-greedy result
            int orderGreedyBestMatchAcc  = -1;
            int orderGreedyBestMatchRank = -1;
            std::vector<MatchPair> orderGreedyBestMatchResult;
            ObtainBestNodeMatch(orderGreedyResults, trueResult, orderGreedyBestMatchResult, orderGreedyBestMatchAcc, orderGreedyBestMatchRank);

            //// For simple greedy result
            int simpleGreedyBestMatchAcc  = 0;
            int simpleGreedyBestMatchRank = 1;
            // Convert the result into the compact version (no -1)
            std::sort(simgpleGreedyResult.begin(), simgpleGreedyResult.end(), IndexCompareMatchPairIIAscend);
            // Compute accuracy
            for (int h = 0; h < simgpleGreedyResult.size(); ++h)
            for (int g = 0; g < trueResult.size(); ++g)
            {
                if (simgpleGreedyResult[h].p == trueResult[g].p && simgpleGreedyResult[h].q == trueResult[g].q)
                {
                    ++simpleGreedyBestMatchAcc;
                    continue;
                }
            }
            int simpleGreedyBestOrderAcc  = 0;
            int simpleGreedyBestOrderRank = 1;
            for (int h = 1; h < simgpleGreedyResult.size()-1; h += 2)
            for (int g = 1; g < trueResult.size()-1; g += 2)
            {
                if (simgpleGreedyResult[h].q == trueResult[g].q && simgpleGreedyResult[h+1].q == trueResult[g+1].q)
                {
                    ++simpleGreedyBestOrderAcc;
                    continue;
                }
            }

            //////// B: Find the best correct-order result from Gorgon and order-greedy --- record accuracy, Gorgon rank and order-greedy dummy rank
            //// For gorgon result
            int gorgonBestOrderAcc = -1;
            int gorgonBestOrderRank = -1;
            std::vector<MatchPair> gorgonBestOrderResult;
            ObtainBestOrderMatch(gorgonResults, trueResult, gorgonBestOrderResult, gorgonBestOrderAcc, gorgonBestOrderRank);

            //// For order-greedy result
            int orderGreedyBestOrderAcc = -1;
            int orderGreedyBestOrderRank = -1;
            std::vector<MatchPair> orderGreedyBestOrderResult;
            ObtainBestOrderMatch(orderGreedyResults, trueResult, orderGreedyBestOrderResult, orderGreedyBestOrderAcc, orderGreedyBestOrderRank);

            // Compute the cost for ground truth, A, B and simple-greedy
            std::cout << "Compute truthMatch cost...\n";
            double truthCost = ComputeMatchCost(patternGraph, baseGraph, trueResult);
            std::cout << "Compute gorgonBestMatch cost...\n";
            double gorgonBestMatchCost = ComputeMatchCost(patternGraph, baseGraph, gorgonBestMatchResult);
            std::cout << "Compute orderGreedyBestMatchCost cost...\n";
            double orderGreedyBestMatchCost = ComputeMatchCost(patternGraph, baseGraph, orderGreedyBestMatchResult);
            std::cout << "Compute simpleGreedyMatchCost cost...\n";
            double simpleGreedyMatchCost = ComputeMatchCost(patternGraph, baseGraph, simgpleGreedyResult);

            double simpleGreedyOrderCost = simpleGreedyMatchCost;
            std::cout << "Compute gorgonBestOrderCost...\n";
            double gorgonBestOrderCost = ComputeMatchCost(patternGraph, baseGraph, gorgonBestOrderResult);
            std::cout << "Compute orderGreedyBestOrderCost...\n";
            double orderGreedyBestOrderCost = ComputeMatchCost(patternGraph, baseGraph, orderGreedyBestOrderResult);

            // Compute the score for ground truth, A, B and simple-greedy
            double truthScore = ComputeMatchScore(patternGraph, baseGraph, trueResult, affinityMatrix);
            double gorgonBestMatchScore = ComputeMatchScore(patternGraph, baseGraph, gorgonBestMatchResult, affinityMatrix);
            double orderGreedyBestMatchScore = ComputeMatchScore(patternGraph, baseGraph, orderGreedyBestMatchResult, affinityMatrix);
            double simpleGreedyMatchScore = ComputeMatchScore(patternGraph, baseGraph, simgpleGreedyResult, affinityMatrix);

            double simpleGreedyOrderScore = simpleGreedyMatchScore;
            double gorgonBestOrderScore = ComputeMatchScore(patternGraph, baseGraph, gorgonBestOrderResult, affinityMatrix);
            double orderGreedyBestOrderScore = ComputeMatchScore(patternGraph, baseGraph, orderGreedyBestOrderResult, affinityMatrix);

            // Display best correct-match-node results
            std::cout << "Results analysis:\n" << "Truth cost and score: " << truthCost << ", " << truthScore << "\n";
            std::cout << "Best correct-match-node results:\n";
            std::cout << "Truth, Gorgon, simple-geedy, order-greedy:\n";
            //std::cout << std::fixed << std::setprecision(10);
            PrintMatchResult(trueResult, patternGraph->GetNodeCount(), truthCost, truthScore, trueResult.size(), trueResult.size(), 1);
            PrintMatchResult(gorgonBestMatchResult, patternGraph->GetNodeCount(), gorgonBestMatchCost, gorgonBestMatchScore, gorgonBestMatchAcc, trueResult.size(), gorgonBestMatchRank);            
            PrintMatchResult(simgpleGreedyResult, patternGraph->GetNodeCount(), simpleGreedyMatchCost, simpleGreedyMatchScore, simpleGreedyBestMatchAcc, trueResult.size(), simpleGreedyBestMatchRank);
            PrintMatchResult(orderGreedyBestMatchResult, patternGraph->GetNodeCount(), orderGreedyBestMatchCost, orderGreedyBestMatchScore, orderGreedyBestMatchAcc, trueResult.size(), orderGreedyBestMatchRank);

            // Display best correct-order results
            std::cout << "Best correct-order results:\n";
            std::cout << "Gorgon, simple-geedy, order-greedy:\n";
            // Note now we count edges with correct connectivity
            PrintMatchResult(gorgonBestOrderResult, patternGraph->GetNodeCount(), gorgonBestOrderCost, gorgonBestOrderScore, gorgonBestOrderAcc, trueResult.size()/2 - 1, gorgonBestOrderRank);            
            PrintMatchResult(simgpleGreedyResult, patternGraph->GetNodeCount(), simpleGreedyOrderCost, simpleGreedyOrderScore, simpleGreedyBestOrderAcc, trueResult.size()/2 - 1, simpleGreedyBestOrderRank);
            PrintMatchResult(orderGreedyBestOrderResult, patternGraph->GetNodeCount(), orderGreedyBestOrderCost, orderGreedyBestOrderScore, orderGreedyBestOrderAcc, trueResult.size() / 2 - 1, orderGreedyBestOrderRank);

*/

            // Make sure we print in the console (before exit)
            flushall();
            
			return foundCount;
		}

		// returns one of the results of a correspondence search
		SSECorrespondenceResult WongMatch15ConstrainedNoFuture::GetResult(int rank) {
			//if(rank <= (int)solutions.size() && (rank >= 1)) {
				return solutions[rank-1];
			//} else {
			//	return NULL;
			//}
		}

		// prints correspondence search results
		void WongMatch15ConstrainedNoFuture::SaveResults(){
			//printf("\t");
			//for(int i = 0; i < currentNode->n1Top; i++) {
			//	printf("%d ", currentNode->n2[i]);
			//}
			//printf(" - %f : (%d expanded)\n", currentNode->cost, expandCount);
			//printf("%d\t\t", expandCount);	
#ifdef VERBOSE
			printf("Time taken in GetA %f\n", timeInGetA / (double)CLOCKS_PER_SEC);
			printf("Time taken in GetB %f\n", timeInGetB / (double)CLOCKS_PER_SEC);
			printf("Time taken in Queue %f\n", timeInQueue / (double)CLOCKS_PER_SEC);
#endif

		}

		// returns the cost of matching node p in the pattern graph to node qp in the base graph
		// this method does not include any cost for matching strands to sheets.
		double WongMatch15ConstrainedNoFuture::GetC(int p, int qp) {
			//return GetC(p, p, qp, qp);
			double cost = GetC(p, p, qp, qp);
			
			// if sheet-to-strand match, compute the cost of the match based on the unused sheet capacity and the strand length
			if( (int)(patternGraph->adjacencyMatrix[p-1][p-1][0] + 0.01) == GRAPHNODE_SHEET && 
				(int)(baseGraph->adjacencyMatrix[qp-1][qp-1][0] + 0.01) == GRAPHNODE_SHEET ) 
            {
				cost = 0;
				//cout << " ... original cost = " << cost << endl;
				cost = abs(patternGraph->nodeWeights[p-1] - baseGraph->nodeWeights[qp-1]);
				//cout << " ... cost to match strand " << p << " to sheet " << qp << " is " << cost << endl;
			}

			return cost;

		}

		// returns the cost of matching node p in the pattern graph to node qp in the base graph
		// this method includes the cost of matching strands to sheets.
		double WongMatch15ConstrainedNoFuture::GetC(int p, int qp, LinkedNodeStub * currentNode) {
			double cost = GetC(p, p, qp, qp);
			
			// if sheet-to-strand match, compute the cost of the match based on the unused sheet capacity and the strand length
			if( (int)(patternGraph->adjacencyMatrix[p-1][p-1][0] + 0.01) == GRAPHNODE_SHEET && 
				(int)(baseGraph->adjacencyMatrix[qp-1][qp-1][0] + 0.01) == GRAPHNODE_SHEET ) {
				cost = 0;
				//cout << " ... original cost = " << cost << endl;
				cost = abs(patternGraph->nodeWeights[p-1] - baseGraph->nodeWeights[qp-1]);
				//cout << " ... cost to match strand " << p << " to sheet " << qp << " is " << cost << endl;
			}
			return cost;

		}

		// returns the cost of matching edge j,p in the pattern graph to edge qj,qp in the base graph.
		// when j=p and qj=qp, returns the cost of matching node j to node qj.
		// three possible cases:
		//   j == p and qj == qp -- node match cost
		//   j != p and qj == qp -- edge match cost, special case where same sheet revisited by two consecutive nodes
		//   j != p and qj != qp -- edge match cost
		// note: only the first case is ever used, as all calls to this method have j=p and qj=qp.
		double WongMatch15ConstrainedNoFuture::GetC(int j, int p, int qj, int qp) {

			double jpCost;
			double qjqpCost;
			double typeCost = 0;

			// if no edge exists between j and p
			if(patternGraph->adjacencyMatrix[j-1][p-1][1] == MAXINT) {
				jpCost = 1000;
			} else {
				// if edge exists or if j == p
				jpCost = patternGraph->adjacencyMatrix[j-1][p-1][1];
			}

			// if no edge exists between qj and qp
			if(baseGraph->adjacencyMatrix[qj-1][qp-1][1] == MAXINT) {
				qjqpCost = 1000;
			} else {
				// if edge exists or if qj == qp
				qjqpCost = baseGraph->adjacencyMatrix[qj-1][qp-1][1];
			}

			// if edge types or node types do not match. cost set here may be recomputed in next block.
			if(patternGraph->adjacencyMatrix[j-1][p-1][0] != baseGraph->adjacencyMatrix[qj-1][qp-1][0]) {
				typeCost = 1000;
			}

			// case where a sheet is revisited by two consecutive strands
			if( (int)(patternGraph->adjacencyMatrix[j-1][p-1][0] + 0.01) == GRAPHEDGE_LOOP && 
				(int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] + 0.01) == GRAPHNODE_SHEET) {
				//cout << "----> revisiting sheet " << qj-1 << ". type cost 0, total cost " << fabs(jpCost - qjqpCost) << endl;
				typeCost = 0;
			}
			
			//cout << "cost of pattern " << patternGraph->adjacencyMatrix[j-1][p-1][0] << " and base " << baseGraph->adjacencyMatrix[qj-1][qp-1][0] << " is " << fabs(jpCost - qjqpCost) + typeCost << endl;
			return fabs(jpCost - qjqpCost) + typeCost;
		}

		// returns the cost of matching a loop from the pattern graph (with missing helices) to a loop in the base graph.
		// d is the start match in the pattern graph
		// m is the number of missing helices or sheets in the pattern graph
		// qj is the start node in the base graph
		// qp is the end node in the base graph
		double WongMatch15ConstrainedNoFuture::GetCost(int d, int m, int qj, int qp, bool debugMsg) {
			// TODO: Fix patthernLength and baseLength for sheet-to-sheet case.
			double patternLength = 0;
			double baseLength;

			// Count number of skipped helices and sheets
			int skippedHelices = 0;
			int skippedSheets = 0;

			// Adding the length of the skipped helixes
			for(int i = 1; i < m; i++) {
				patternLength += patternGraph->adjacencyMatrix[d+i-1][d+i-1][1];
				if (patternGraph->adjacencyMatrix[d+i-1][d+i-1][0] == GRAPHNODE_SHEET) {
					patternLength += patternGraph->nodeWeights[d+i-1];
					skippedSheets++;
#ifdef VERBOSE
					if (debugMsg) { cout << "  -- found strand " << d+i << ", adding " << patternGraph->nodeWeights[d+i-1] << " to patternLength" << endl; }
#endif // VERBOSE
				} else {
					skippedHelices++;
				}
			}

			// Adding the length of the edges
			// Check that first and last edges are both loops, a requirement for all edges in the pattern graph
			bool firstIsLoop = false;
			bool lastIsLoop = false;
			for(int i = 0; i < m; i++) {
				lastIsLoop = ((int)(patternGraph->adjacencyMatrix[d+i-1][d+i][0] + 0.01) == GRAPHEDGE_LOOP) ;
				if(i==0) {
					firstIsLoop = lastIsLoop;
				}
				patternLength += patternGraph->adjacencyMatrix[d+i-1][d+i][1];	
				
				/*
				// rescale helices if they are skipped
				if ((int)(patternGraph->adjacencyMatrix[d+i-1][d+i][0] + 0.01) == GRAPHEDGE_HELIX) {
					patternLength += (HELIX_C_ALPHA_TO_ANGSTROMS / LOOP_C_ALPHA_TO_ANGSTROMS)*(patternGraph->adjacencyMatrix[d+i-1][d+i][1]);	
					//patternLength += (LOOP_C_ALPHA_TO_ANGSTROMS / HELIX_C_ALPHA_TO_ANGSTROMS)*(patternGraph->adjacencyMatrix[d+i-1][d+i][1]);	
				} else {
					patternLength += patternGraph->adjacencyMatrix[d+i-1][d+i][1];	
				}
				*/
			}
			// TODO: Fix, has bug. But getting closer.
			skippedHelices = skippedHelices / 2;
			//cout << "d=" << d << ", m=" << m << ", skipH=" << skippedHelices << ", skipS=" << skippedSheets << endl;

			bool euclideanEstimate = false;
			double weight = 1.0;

			// if edge begins with an unmatched node in the base graph
			if(qj == -1) { // special handling for missing helixes at the ends
				baseLength = 0;
				switch((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01)) {
					case(GRAPHEDGE_HELIX) : 
						weight = HELIX_WEIGHT_COEFFICIENT;
						break;
					case(GRAPHEDGE_LOOP):
						weight = LOOP_WEIGHT_COEFFICIENT;
						break;
					case(GRAPHEDGE_SHEET):
						weight = SHEET_WEIGHT_COEFFICIENT;
						break;
				}
			}
			else {
				assert(baseGraph->EdgeExists(qj-1,qp-1));		
				baseLength = baseGraph->adjacencyMatrix[qj-1][qp-1][1];
				euclideanEstimate = ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] + 0.01) == GRAPHEDGE_LOOP_EUCLIDEAN);
				switch((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] + 0.01)) {
					case(GRAPHEDGE_HELIX) : 
						weight = HELIX_WEIGHT_COEFFICIENT;
						break;
					case(GRAPHEDGE_LOOP):
						weight = LOOP_WEIGHT_COEFFICIENT;
						break;
					case(GRAPHNODE_SHEET): // two strands in a row match to the same sheet
						//cout << "---> sheet to sheet case. parameters: " << d << "," << m << "," << qj << "," << qp << endl;
						weight = LOOP_WEIGHT_COEFFICIENT;
						break;
				}
				weight = euclideanEstimate? weight * EUCLIDEAN_LOOP_PENALTY: weight;
			}

			if(m == 1) { // not a skip edge
				if( (qj != -1) && // first node in pattern graph is matched
					!( ((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01) == (int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] + 0.01)) ) && // types don't match exactly
					!( ((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01) == GRAPHEDGE_LOOP) && ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) == GRAPHEDGE_LOOP_EUCLIDEAN)) && // not a loop-Euclidianloop match
					!( ((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01) == GRAPHEDGE_LOOP) && ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) == GRAPHNODE_SHEET)) ) 	{ // not a loop-sheet match
                    //cout << "Negative edge cost value due to AAA\n";
					return -1;
				}		
#ifdef VERBOSE
				if (debugMsg) { cout << "  -- euclidean dist = " << baseGraph->euclideanMatrix[qj-1][qp-1] << ", patternLength = " << patternLength << ", loop fudge factor = " << EUCLIDEAN_VOXEL_TO_PDB_RATIO / LOOP_C_ALPHA_TO_ANGSTROMS << ", helix fudge factor = " << EUCLIDEAN_VOXEL_TO_PDB_RATIO / HELIX_C_ALPHA_TO_ANGSTROMS << endl; }
				if (debugMsg) { cout << "  -- scalar ratio required = " << baseGraph->euclideanMatrix[qj-1][qp-1] / patternLength << ", additive headroom = " << baseGraph->euclideanMatrix[qj-1][qp-1] - patternLength << endl; }
#endif // VERBOSE

				if((qj != -1) && ((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01) == GRAPHEDGE_HELIX) && (baseGraph->euclideanMatrix[qj-1][qp-1] > (patternLength * EUCLIDEAN_VOXEL_TO_PDB_RATIO / HELIX_C_ALPHA_TO_ANGSTROMS )) ){					
#ifdef VERBOSE
					if (debugMsg) { cout << "  -- -- -- NOT ALLOWED (HELIX) -- -- -- " << endl; }
#endif // VERBOSE
                    //cout << "Negative edge cost value due to BBB\n";                    
					return -1;
				} else 
				if((qj != -1) && ((int)(patternGraph->adjacencyMatrix[d-1][d][0] + 0.01) == GRAPHEDGE_LOOP)) {
					if (((int)(patternGraph->adjacencyMatrix[d-1][d-1][0] + 0.01) == GRAPHNODE_SHEET || (int)(patternGraph->adjacencyMatrix[d][d][0] + 0.01) == GRAPHNODE_SHEET) && (baseGraph->euclideanMatrix[qj-1][qp-1] > (patternLength * 1.0 * EUCLIDEAN_VOXEL_TO_PDB_RATIO / LOOP_C_ALPHA_TO_ANGSTROMS )) ){					
#ifdef VERBOSE
						if (debugMsg) { cout << "  -- -- -- NOT ALLOWED (LOOP WITH STRAND) -- -- -- " << endl; }		
#endif // VERBOSE
                        //cout << "Negative edge cost value due to CCC\n";
						return -1;
					}
					if (((int)(patternGraph->adjacencyMatrix[d-1][d-1][0] + 0.01) != GRAPHNODE_SHEET && (int)(patternGraph->adjacencyMatrix[d][d][0] + 0.01) != GRAPHNODE_SHEET) && (baseGraph->euclideanMatrix[qj-1][qp-1] > (patternLength * EUCLIDEAN_VOXEL_TO_PDB_RATIO / LOOP_C_ALPHA_TO_ANGSTROMS )) ){					
#ifdef VERBOSE
						if (debugMsg) { cout << "  -- -- -- NOT ALLOWED (LOOP) -- -- -- " << endl; }		
#endif // VERBOSE
                        //cout << "Negative edge cost value due to DDD\n";
						return -1;
					}
				}
			} else { // a skip edge
				// not sure if these checks really help or if they just waste time
				if( !(firstIsLoop && lastIsLoop) || // pattern graph edge doesn't start and end with loops OR 
					( ( ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) != GRAPHEDGE_LOOP)) &&			// (base graph edge not a loop AND
					( ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) != GRAPHEDGE_LOOP_EUCLIDEAN)) &&	// base graph edge not a Euclidian loop AND
					( ((int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) != GRAPHNODE_SHEET)) ) ) {		// base graph edge not a sheet)
                    //cout << "Negative edge cost value due to A\n";
                    //cout << "Negative edge cost value due to EEE\n";
					return -1;
				}
				// check here to sum up the parts of the skip edge and compare to the euclidian distance, if it's a euclidian edge in the base graph
				if((qj != -1) && (baseGraph->euclideanMatrix[qj-1][qp-1] > (patternLength * EUCLIDEAN_VOXEL_TO_PDB_RATIO ))){
                    //cout << "Negative edge cost value due to B\n";
                    //cout << "Negative edge cost value due to FFF\n";
					return -1;
				}
			}


			switch(COST_FUNCTION)
			{
			
			//case(1):
			//	if ( (int)(baseGraph->adjacencyMatrix[qj-1][qp-1][0] +0.01) == GRAPHEDGE_HELIX ) {
			//		return weight * fabs(patternLength - baseLength);
			//	} else {
			//		// TODO: Remove this arbitrary weight!
			//		return 50.0*weight * fabs(patternLength - baseLength) / (patternLength + baseLength);
			//	}
			//	break;
			case(1):
				return weight * fabs(patternLength - baseLength);
				break;
			case(2):
				return weight * fabs(patternLength - baseLength) / (patternLength + baseLength);
				break;
			case(3):
				return weight * pow((patternLength - baseLength),2);
				break;
			}
			// this line should be unreachable
			return 0;
		}

		double WongMatch15ConstrainedNoFuture::GetF() {
			return currentNode->costGStar;	
		}

		void WongMatch15ConstrainedNoFuture::PopBestNode(){
		#ifdef VERBOSE
			clock_t start = clock();
		#endif
			double cost;
			//queue->remove(currentNode, cost);
			queue->PopFirst(cost, currentNode);
		#ifdef VERBOSE
			timeInQueue += clock() - start;
		#endif
		}


		// add in penalties for skipped helices and sheets
		// m is the number of nodes involved in the match. m=1 is no skipped helices or sheets.
		double WongMatch15ConstrainedNoFuture::GetPenaltyCost(int d, int m, bool debugMsg) {
			//if (d==0) {cout << "d=" << d << ", m=" << m << endl; debugMsg=true;}
			double cost = 0.0;
			int lastPatternNode = patternGraph->GetNodeCount() - 1;
			bool startAtBeginning = ( d == 0 );
			bool finishAtEnd = ( d + m-1 == lastPatternNode );
			bool pastFirst = true;
			bool firstHelixFound = false;
			for(int k = d; k < d + m-1; k++) {
				//if (debugMsg) {cout << "+++++++++++++++++++inside loop k=" << k << endl;}
				//cout << "  GetPenaltyCost(" << d << "," << m << "). k=" << k << ". " << endl;
				// add penalties for all skipped helices
				if((int)(patternGraph->adjacencyMatrix[k][k+1][0] + 0.01) == GRAPHEDGE_HELIX) {
					//cout << "    GetPenaltyCost(" << d << "," << m << "). SKIP HELIX. k=" << k << endl;
					cost += MISSING_HELIX_PENALTY;
					cost += patternGraph->adjacencyMatrix[k][k+1][1] * MISSING_HELIX_PENALTY_SCALED;
#ifdef VERBOSE
					if (debugMsg) { cout << "  -- adding missing helix penalties: fixed=" << MISSING_HELIX_PENALTY << ", scaled=" << patternGraph->nodeWeights[k] * MISSING_HELIX_PENALTY_SCALED << endl; }
#endif // VERBOSE

					if (startAtBeginning && !firstHelixFound) {
						cost += START_END_MISSING_HELIX_PENALTY;
#ifdef VERBOSE
						if (debugMsg) { cout << "  -- adding start_end_miss_helix_pen" << endl; }
#endif // VERBOSE
					}
					if (finishAtEnd && !firstHelixFound) {
						cost += START_END_MISSING_HELIX_PENALTY;
#ifdef VERBOSE
						if (debugMsg) { cout << "  -- adding start_end_miss_helix_pen" << endl; }
#endif // VERBOSE
					}
					firstHelixFound = true;
				}
				// add penalties for skipped strands, unless the strand falls at the beginning of the sequence and is not the first node
				else if( (startAtBeginning || pastFirst) && ((int)(patternGraph->adjacencyMatrix[k][k][0] + 0.01) == GRAPHNODE_SHEET) ) {
					cost += MISSING_SHEET_PENALTY;
					cost += patternGraph->nodeWeights[k] * MISSING_SHEET_PENALTY_SCALED;
#ifdef VERBOSE
					if (debugMsg) { cout << "  -- adding missing sheet penalties: fixed=" << MISSING_SHEET_PENALTY << ", scaled=" << patternGraph->nodeWeights[k] * MISSING_SHEET_PENALTY_SCALED << endl; }
#endif // VERBOSE
				}
				//if (startAtBeginning && debugMsg) { cout << "STARTATBEGIN" << endl;}
				pastFirst = true;
			}

			if (finishAtEnd && patternGraph->adjacencyMatrix[lastPatternNode-1][lastPatternNode-1][0] + 0.01 == GRAPHNODE_SHEET){
					cost += MISSING_SHEET_PENALTY;
					cost += patternGraph->nodeWeights[lastPatternNode-1] * MISSING_SHEET_PENALTY_SCALED;
			}
			//if (debugMsg) { cout << "  -- returning cost=" << cost << endl; }
			return cost;
		}


		// expand a node.
		// checks the base graph for edges between this node and every other.
		// if an edge is found, match the pattern graph to that edge and add the match to the queue.
		// also match edges that include skip edges in the pattern graph
		// costs of matches are determined by the GetC method
		bool WongMatch15ConstrainedNoFuture::ExpandNode(LinkedNodeStub * currentStub) {
			bool expanded = false;
			expandCount++;

			LinkedNode * temp;
			double edgeCost;
#ifdef VERBOSE
			if(longestMatch < currentNode->depth) {
				longestMatch = currentNode->depth;
				//printf(" %d elements matched! (%f kB Memory Used)\n", longestMatch, (queue->getLength() * sizeof(LinkedNode) + usedNodes.size() * sizeof(LinkedNodeStub)) / 1024.0);
				printf(" %d elements matched!\n", longestMatch);
				//cout << "WongMatch15ConstrainedNoFuture::ExpandNode: " << longestMatch << " elements expanded (" << ((queue->getLength() * sizeof(LinkedNode) + usedNodes.size() * sizeof(LinkedNodeStub)) / 1024.0) << " kB memory used)" << endl;
			}
#endif //VERBOSE
			
			int currentM1Top = patternGraph->nodeCount - currentNode->depth; // remaining unmatched nodes in sequence
			bool notConstrained;
			//currentNode->PrintNodeConcise(0, true, true);

			// Expanding nodes with a real terminal node
			// for every node i in baseGraph
			for(int i = 1; i <= baseGraph->nodeCount; i++) {
				// if: 
				//   currentNode is at level 0 of tree
				//   or 
				//   i is in the currentNode bitmap, and there is an edge in baseGraph between currentNode and node i
				if((currentNode->depth == 0) || 
					(LinkedNode::IsNodeInBitmap(currentNode->m2Bitmap, i) && (baseGraph->EdgeExists(currentNode->n2Node-1, i-1)))) {						
					int skippedHelixNodes = 0;
					int skippedSheetNodes = 0;

                    //cout << "Well, (enter first level) currentNode->depth: " << (int)currentNode->depth << "\n";
                    if (currentM1Top <= 0)
                        break;

					for (int j = 0; (j <= currentM1Top) && (skippedHelixNodes + currentNode->missingHelixNodesUsed <= missingHelixCount * 2) && (skippedSheetNodes + currentNode->missingSheetNodesUsed <= missingSheetCount); ) {
						// i is the node from baseGraph being matched to currentNode
						// j is the number of missing helix or sheet nodes from patternGraph to be skipped for this match

						notConstrained = true;

						for(int k = currentNode->n1Node + 1; k <= currentNode->n1Node + j; k++) {
							notConstrained = notConstrained && IsNodeAssignmentAllowed(k, -1);
						}
						notConstrained = notConstrained && IsNodeAssignmentAllowed(currentNode->n1Node + j + 1, i);

						// used later for case where first helix is missing
						bool firstMissing = false;

						if(notConstrained) {	
							// store the current node as temp
							temp = currentNode; 
							// create new current node. i is the index of the new node(?), j is the number of skipped nodes.


							// check whether i is a revisitable node (a sheet)
							bool revisitable = ( (int)(baseGraph->adjacencyMatrix[i-1][i-1][0] + 0.01) == GRAPHNODE_SHEET );

							if(((temp->depth == 0) && (j > 0)) || 
								((patternGraph->nodeCount - currentNode->depth == 0) && (currentNode->n2Node == -1))) {
									if (skippedHelixNodes == 0 && patternGraph->adjacencyMatrix[0][0][0] == GRAPHNODE_HELIX) {
										skippedHelixNodes = 1;
#ifdef VERBOSE
										cout << "node skipped. adding one to skippedHelixNodes. result is " << skippedHelixNodes << endl;
#endif // VERBOSE
									}
							}

							// generate a current node, marking it as revisitable or not depending on result from test
							// the constructor marches forward along the sequence, skipping j nodes
							currentNode = new LinkedNode(currentNode, currentStub, i, skippedHelixNodes, skippedSheetNodes, revisitable);

							currentNode->costGStar = 0;

							// if previous node was at top of tree and it was skipped
							if(((temp->depth == 0) && (j > 0)) || 
								((patternGraph->nodeCount - currentNode->depth == 0) && (currentNode->n2Node == -1))) {
								if (skippedHelixNodes > 0) {
									//cout << "first helix is missing." << endl;
									firstMissing = true;
								}
							}	

							// if previous node was at top of tree
							if(temp->depth == 0) {
								edgeCost = 0;
							} else {								
								edgeCost = GetCost(temp->n1Node, j+1, temp->n2Node, currentNode->n2Node, false);
							}
							
                            //cout << "edgeCost: " << edgeCost << "\n";
                            //cout << "next node depth: " << (int)currentNode->depth << "\n";
                            //cout << "Ah: (second level) current node depth: " << (int)temp->depth << "\n";
                            //cout << "Ah: next node depth: " << (int)currentNode->depth << "\n";

							// if this is an allowed match:
							if(edgeCost >= 0) {

                                //cout << "Ah: current node depth: " << (int)temp->depth << "\n";
                                //cout << "Ah: indeed generate next node depth: " << (int)currentNode->depth << "\n";

								//worked! currentNode->costGStar += temp->costGStar + edgeCost + GetC(currentNode->n1Node, currentNode->n2Node, currentNode);
								currentNode->costGStar += temp->costGStar + edgeCost + GetC(currentNode->n1Node, currentNode->n2Node);

								// add costs for skipped helices and sheets
								currentNode->costGStar += GetPenaltyCost(temp->n1Node, j+1, false);
								
								currentNode->cost = GetF();			
								//currentNode->PrintNodeConcise(-1, true, true);
								//queue->add(currentNode, currentNode->cost);
								queue->Add(currentNode->cost, currentNode);
								expanded = true;
							} else { // not an allowed match
								delete currentNode;
							}
							// set currentNode pointer back to the current node
							currentNode = temp;	
						}
						// if this node is a helix, increment j by one more to prepare for the next iteration
						switch ( (int)(patternGraph->adjacencyMatrix[currentNode->n1Node + j][currentNode->n1Node + j][0] + 0.01)) {
							case GRAPHNODE_HELIX:
								skippedHelixNodes += 2;
								j+=2;
								break;
							default:
								skippedSheetNodes += 1;
								j++;
								break;
						}
					} // end for j
				}
			}


			// Expanding nodes with a dummy terminal node:
			// Count the number of sheets and helices remaining in the helix, and then try to make
			// a long skip edge to pass over all of them.
			
			// count the number of helix and sheet nodes that are not yet matched
			int remainingSheetNodes = 0;
			int remainingHelixNodes = 0;
			for (int l = currentNode->n1Node + 1; l <= patternGraph->nodeCount; l++) {
				if (patternGraph->adjacencyMatrix[l-1][l-1][0] == GRAPHNODE_HELIX) {
					remainingHelixNodes++;
				} else if (patternGraph->adjacencyMatrix[l-1][l-1][0] == GRAPHNODE_SHEET) {
					remainingSheetNodes++;
				}
			}

			// if possible, create an edge to jump to the end of the sequence
			if(2*missingHelixCount - currentNode->missingHelixNodesUsed >= remainingHelixNodes && 2*missingSheetCount - currentNode->missingSheetNodesUsed >= remainingSheetNodes) {
				notConstrained = true;
				for(int k = currentNode->n1Node + 1; k <= patternGraph->nodeCount; k++) {
					notConstrained = notConstrained && IsNodeAssignmentAllowed(k, -1);
				}

				if(notConstrained) {
					temp = currentNode;
					currentNode = new LinkedNode(temp);
					currentNode->depth = (char)patternGraph->nodeCount;
					currentNode->costGStar = temp->costGStar;
					currentNode->costGStar += GetPenaltyCost(temp->n1Node, remainingHelixNodes + remainingSheetNodes, false);
					currentNode->cost = currentNode->costGStar;
					//queue->add(currentNode, currentNode->cost);
					queue->Add(currentNode->cost, currentNode);
					currentNode = temp;
				}
			}
			return expanded;
		}

		void WongMatch15ConstrainedNoFuture::AnalyzeResults(int results[][MAX_NODES], int groundTruth[]) {

			int numNodes = patternGraph->GetNodeCount();
			int nh=0, ns=0;
			// count the number of helices and sheets
			for (int i=0; i < numNodes; i++) {
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_HELIX) {
					nh++;
				}
				else if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
					ns++;
				}
			}

			// compute the % of helices and sheets correctly predicted by each result
			int nCorrectHelices[RESULT_COUNT];
			int nCorrectHelicesFlipped[RESULT_COUNT];
			int nCorrectSheets[RESULT_COUNT];

			for (int res = 0; res < RESULT_COUNT; res++) {
				// Check if all helices or all sheets were correctly matched
				for (int i=0; i < numNodes; i++) {
					if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_HELIX) {
						if (results[res][i]==SOLUTION[i]) {
							nCorrectHelices[res]++;
						}
						if (results[res][i]==SOLUTION[i] || results[res][i]==SOLUTION[max(0,i-1)] || results[res][i]==SOLUTION[min(numNodes-1,i+1)]) {
							nCorrectHelicesFlipped[res]++;
						}
					}
					if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
						if (results[res][i]==SOLUTION[i]) {
							nCorrectSheets[res]++;
						}
					}
				}
			}

			cout << endl;
			//for (int i = 0; i < foundCount; i++){
			//	cout << "Result " << i << " has " << nCorrectHelices[i]/2 << " correct helices (" << nCorrectHelicesFlipped[i]/2 << " if direction is ignored) and " << nCorrectSheets[i] << " correct strands." << endl;
			//}
			//for (int i = 0; i < foundCount; i++){
			//	cout << "Result " << i << " has " << (double)nCorrectHelices[i]/(nh) << " correct helices (" << (double)nCorrectHelicesFlipped[i]/nh << " if direction is ignored) and " << (double)nCorrectSheets[i]/ns << " correct strands." << endl;
			//}

			// compute the average % of correct helices and sheets over all results
			

			// count up how many votes each node receives
			int votes[MAX_NODES][MAX_NODES];
			for(int i = 0; i < MAX_NODES; i++) {
				for(int j = 0; j < MAX_NODES; j++) {
					votes[i][j] = 0;
				}
			}
			for (int i=0; i<foundCount; i++) {
				for (int j=0; j<numNodes; j++) {
					int thisVote = results[i][j];
					thisVote = max(thisVote,0); // store votes for -1 at location 0, which is otherwise unused
					//cout << "casting vote " << thisVote << " at votes[" << j << "][" << thisVote << "]" << endl;
					votes[j][thisVote]++;
				}
			}

			/*
			cout << "votes: " << endl;
			for (int i = 0; i < MAX_NODES; i++) {
				for (int j = 0; j < numNodes; j++) {
					cout.width(2);
					cout << votes[j][i];
				}
				cout << endl;
			}*/
			
			// report the % of votes recieved by the ground truth
			int tots=0, toth=0, tothf=0;
			for (int i = 0; i < foundCount; i++) {
				tots += nCorrectSheets[i];
				toth += nCorrectHelices[i];
				tothf += nCorrectHelicesFlipped[i];
			}
			cout << "Average helix predictions accuracy is " << (double)toth/(foundCount*nh) << " (" << (double)tothf/(foundCount*nh) << " if direction is ignored). Average sheet accuracy is " << (double)tots/(foundCount*ns) << endl;


			// compute # of results that had perfect helix and sheets scores
			int nAllH=0, nAllHf=0, nAllS=0;
			for (int i=0; i < foundCount; i++) {
				if (nCorrectSheets[i]==ns) {nAllS++;}
				if (nCorrectHelices[i]==nh) {nAllH++;}
				if (nCorrectHelicesFlipped[i]==nh) {nAllHf++;}
			}
			cout << "Results with all helices correct: " << nAllH << " (" << nAllHf << " counting flips). Results with all sheets correct: " << nAllS << endl;

			// count # of votes that each ground truth node gets
			int groundTruthVotes[MAX_NODES];
			for (int i = 0; i < numNodes; i++) {
				int truthKey=max(groundTruth[i],0);
				groundTruthVotes[i]=votes[i][truthKey];
			}
			cout << "The ground truth solution gets the following number of votes:" << endl;
			cout << " seq #  ";
			for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
				cout.width(2);
				cout << i+1 << " ";
			}
			cout << endl;
			cout << " truth  ";
			for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
				cout.width(2);
				cout << groundTruth[i] << " ";
			}
			cout << endl;
			cout << " votes  ";
			for (int i = 0; i < patternGraph->GetNodeCount(); i++) {
				cout.width(2);
				cout << groundTruthVotes[i] << " ";
			}
			cout << endl;

			// find which structure gets the most votes for each node
			int maxVotesIndex[MAX_NODES];
			int maxVotesCount[MAX_NODES];
			for (int i = 0; i < numNodes; i++) {
				int maxCount = -1;
				int maxCountIndex = -1;
				for (int j = 0; j < MAX_NODES; j++) {
					if (votes[i][j]>maxCount) {
						maxCount = votes[i][j];
						//cout << "maxcount[" << i << "][" << j << "]=" << maxCount << endl;
						maxCountIndex=j;
					}
				}
				maxVotesIndex[i]=maxCountIndex;
				maxVotesCount[i]=maxCount;
			}
			cout << "The max vote solution gets the following number of votes:" << endl;
			cout << " index  ";
			for (int i = 0; i < numNodes; i++) {
				cout.width(2);
				if (maxVotesIndex[i]==0) {
					cout << -1 << " ";
				} else {
					cout << maxVotesIndex[i] << " ";
				}
			}
			cout << endl;
			cout << " votes  ";
			for (int i = 0; i < numNodes; i++) {
				cout.width(2);
				cout << maxVotesCount[i] << " ";
			}
			cout << endl;

			// compute success rate of max vote guess
			double maxVoteSuccess=0.0;
			double maxVoteSuccessHelix=0.0;
			double maxVoteSuccessHelixFlipped=0.0;
			double maxVoteSuccessSheet=0.0;
			for (int i = 0; i < numNodes; i++) {
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_HELIX) {
					if (maxVotesIndex[i]==0 && groundTruth[i]==-1) {
						maxVoteSuccess += 1.0;
						maxVoteSuccessHelix += 1.0;
						maxVoteSuccessHelixFlipped += 1.0;
						//cout << " exact match (missing Helix) detected i=" << i << endl;
					}
					if (maxVotesIndex[i] == groundTruth[i]) {
						maxVoteSuccess += 1.0;
						maxVoteSuccessHelix += 1.0;
						maxVoteSuccessHelixFlipped += 1.0;
						//cout << " exact match detected i=" << i << endl;
					}
					// first test for a flipped helix
					if (i>0 && maxVotesIndex[i-1]==groundTruth[i] && maxVotesIndex[i]==groundTruth[i-1]) {
						// another test to make sure these are really the same helix
						int mini=min(maxVotesIndex[i],maxVotesIndex[i-1]);
						int maxi=max(maxVotesIndex[i],maxVotesIndex[i-1]);
						if (mini%2==1 && maxi==mini+1) {
							maxVoteSuccessHelixFlipped += 1.0;
							//cout << " flip1 detected i=" << i << endl;
						}
					}
					// first test for a flipped helix
					if (i<numNodes-1 && maxVotesIndex[i+1]==groundTruth[i] && maxVotesIndex[i]==groundTruth[i+1]) {
						// another test to make sure these are really the same helix
						int mini=min(maxVotesIndex[i],maxVotesIndex[i+1]);
						int maxi=max(maxVotesIndex[i],maxVotesIndex[i+1]);
						if (mini%2==1 && maxi==mini+1) {
							maxVoteSuccessHelixFlipped += 1.0;
							//cout << " flip2 detected i=" << i << endl;
						}
					}
				} else if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
					if (maxVotesIndex[i] == groundTruth[i]) {
						maxVoteSuccess += 1.0;
						maxVoteSuccessSheet += 1.0;
					}
				}
			}
			maxVoteSuccess /= (double)numNodes;
			maxVoteSuccessHelix /= (double)nh;
			maxVoteSuccessHelixFlipped /= (double)nh;
			maxVoteSuccessSheet /= (double)ns;
			cout << "The max vote guess success rate is " << maxVoteSuccess << ". Helix: " << maxVoteSuccessHelix << " (" << maxVoteSuccessHelixFlipped << " allowing flips). Sheet: " << maxVoteSuccessSheet << "." << endl;

			// compare the # of votes of ground truth with # of votes received by the most-voted 
			double voteRatio[MAX_NODES];
			double averageVoteRatio=0.0;
			double averageVoteRatioHelix=0.0;
			double averageVoteRatioSheet=0.0;
			for (int i = 0; i < numNodes; i++) {
				voteRatio[i]=(double)groundTruthVotes[i]/(double)maxVotesCount[i];
				averageVoteRatio+=voteRatio[i];
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_HELIX) {
					averageVoteRatioHelix+=voteRatio[i];
				}
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
					averageVoteRatioSheet+=voteRatio[i];
				}
			}
			averageVoteRatio /= numNodes;
			averageVoteRatioHelix /= nh;
			averageVoteRatioSheet /= ns;
			cout << " % true ";
			for (int i = 0; i < numNodes; i++) {
				cout.width(2);
				cout << voteRatio[i] << " ";
			}
			cout << endl;
			cout << "Average vote ratio is " << averageVoteRatio << ". Helices: " << averageVoteRatioHelix << ". Sheets: " << averageVoteRatioSheet << endl;
			




			// print the voting rank of the ground truth for each node




		}

		// Compute the cost of the ground truth solution which is submitted by the user.
		void WongMatch15ConstrainedNoFuture::ComputeSolutionCost(int solution[], bool extraMessages) {
			if (extraMessages) {cout << "starting ComputeSolutionCost" << endl;}
			int n1=0, n2=1;
			double edgeCost = 0.0;

			double helixCost = 0.0;
			double helixPenaltyCost = 0.0;
			double loopCost = 0.0;
			double loopPenaltyCost = 0.0;
			double sheetCost = 0.0;
			double sheetPenaltyCost = 0.0;
			double skipPenaltyCost = 0.0;

			double edgePenaltyCost = 0.0;
			double nodeCost = 0.0;
			int skippedHelixNodes = 0;
			int skippedSheetNodes = 0;

			int numNodes = patternGraph->GetNodeCount();

			// iterate over all correspondences, adding each to the previous solution
			while (n2 < numNodes) { 

				// check if first node is skipped helix or sheet
				if (solution[n1] == -1) {
					if (extraMessages) {cout << "skipped node found at " << n1+1 << " with adj matrix value " << patternGraph->adjacencyMatrix[n1][n1][0] << endl;}
					if (patternGraph->adjacencyMatrix[n1][n1][0] == GRAPHNODE_HELIX) {
						skippedHelixNodes++;
					}
					if (patternGraph->adjacencyMatrix[n1][n1][0] == GRAPHNODE_SHEET) {
						skippedSheetNodes++;
					}
				}

				// find the end of the current correspondence
				//cout << "begin while block. n1 = " << n1 << ", n2 = " << n2 << endl;
				//while (solution[n2] == -1 && n2 < numNodes) {
				while (solution[n2] == -1 && n2 < numNodes-1) {
					//cout << "skipped node found at " << n2+1 << " with adj matrix value " << patternGraph->adjacencyMatrix[n2][n2][0] << endl;
					if (patternGraph->adjacencyMatrix[n2][n2][0] == GRAPHNODE_HELIX) {
						skippedHelixNodes++;
					}
					if (patternGraph->adjacencyMatrix[n2][n2][0] == GRAPHNODE_SHEET) {
						skippedSheetNodes++;
					}
					n2++;
				}
	
				//cout << "after advancing n2, n1 = " << n1 << ", n2 = " << n2 << ", skippedHN = " << skippedHelixNodes << ", skippedSN = " << skippedSheetNodes << endl;

				// add edge cost
				if (extraMessages) {cout << "adding (" << n1+1 << "," << n2+1 << "," << solution[n1] << "," << solution[n2] << ")" << endl;}
				double singleEdgeCost = 1000;
				double singleEdgePenaltyCost = 0;

				// if edge exists in base graph, find the cost of this correspondence.
				if (solution[n1] == -1) {
					singleEdgeCost = 0;
					//singleEdgePenaltyCost = GetPenaltyCost(n1+1, n2-n1, extraMessages);
					//singleEdgePenaltyCost = GetPenaltyCost(n1, n2-n1+1, extraMessages);
					singleEdgePenaltyCost = GetPenaltyCost(n1, n2-n1+1, extraMessages);
					//if (extraMessages) {cout << "  GetPenaltyCost("<<n1<<","<<n2-n1<<")="<<singleEdgePenaltyCost<<endl;}
					if (extraMessages) {cout << "  GetPenaltyCost("<<n1<<","<<n2-n1+1<<")="<<singleEdgePenaltyCost<<endl;}
					//if (extraMessages) {cout << "  No edge cost for initial skip edge" << endl;}
					/*if (patternGraph->adjacencyMatrix[n1][n1][0] == GRAPHNODE_SHEET) {
						sheetPenaltyCost += singleEdgePenaltyCost;
						skipPenaltyCost += singleEdgePenaltyCost;
						singleEdgePenaltyCost = MISSING_SHEET_PENALTY;
					}*/
					skipPenaltyCost += singleEdgePenaltyCost;
				} else if (baseGraph->EdgeExists(solution[n1]-1, solution[n2]-1)) {
					if (solution[n2] == -1 && n2==numNodes-1) {
						// last edge is skip edge
						singleEdgeCost = 0; 
					} else {
						singleEdgeCost = GetCost(n1+1, n2-n1, solution[n1], solution[n2], extraMessages);
					}
					singleEdgePenaltyCost = GetPenaltyCost(n1+1, n2-n1, extraMessages);
					//if (patternGraph->adjacencyMatrix[n1][n2][0] == GRAPHEDGE_LOOP) { // misses skip edges
					if (baseGraph->adjacencyMatrix[solution[n1]-1][solution[n2]-1][0] == GRAPHEDGE_HELIX) {
						//cout << "[H" << solution[n1] << "-" << solution[n2] << "]";
						helixCost += singleEdgeCost;
						helixPenaltyCost += singleEdgePenaltyCost;
						skipPenaltyCost += singleEdgePenaltyCost;
					} else {
						loopCost += singleEdgeCost;
						skipPenaltyCost += singleEdgePenaltyCost;
					}
					if (extraMessages) {cout << "  GetCost("<<n1+1<<","<<n2-n1<<","<<solution[n1]<<","<<solution[n2]<<")="<<singleEdgeCost<<endl;}
					if (extraMessages) {cout << "  GetPenaltyCost("<<n1+1<<","<<n2-n1<<")="<<singleEdgePenaltyCost<<endl;}
					if (singleEdgeCost == -1) {cout << "  MATCH FROM NODE " << solution[n1] << " TO NODE " << solution[n2] << " IS NOT ALLOWED (CUTOFF OR TYPE MISMATCH?)" << endl;} 
				} else {
					cout << "  BASE GRAPH DOES NOT HAVE AN EDGE FROM NODE " << solution[n1] << " TO NODE " << solution[n2] << ". THIS SOLUTION NOT POSSIBLE!" << endl;
				}

				// check if first or last helix is unmatched
				if ((n1 == 0 && singleEdgeCost == -1) || (n2 == numNodes && singleEdgeCost == -1)){
					if (extraMessages) {cout << "  first helix or sheet is unmatched. adding penalty of " << singleEdgeCost << endl;}
				} else {
					if (extraMessages) {cout << "  cost of this addition is " << singleEdgeCost << endl;}
				}

				// add node cost
				double singleNodeCost = GetC(n2+1, solution[n2]);
				// if at beginning of sequence, check first node
				if (patternGraph->adjacencyMatrix[n2][n2][0] == GRAPHNODE_SHEET) {
					sheetCost += singleNodeCost;
				}
				if (n1==0 && patternGraph->adjacencyMatrix[n1][n1][0] == GRAPHNODE_SHEET) {
					double firstNodeCost = GetC(n1+1, solution[n1]);
					sheetCost += firstNodeCost;
					singleNodeCost += firstNodeCost;
				}
				if (extraMessages) {cout << "  node cost for nodes " << n2+1 << " and " << solution[n2] << " is " << singleNodeCost << endl;}

				// add the costs from this iteration to the running totals
				edgeCost += singleEdgeCost;
				edgePenaltyCost += singleEdgePenaltyCost;
				nodeCost += singleNodeCost;

				// prepare for next iteration
				n1 = n2;
				n2++;
			}

			// Check if all helices or all sheets were correctly matched
			bool sheetsCorrect=true;
			bool helicesCorrect=true;
			bool helicesCorrectFlipped=true;
			double ratioCorrectHelices = 0;
			double ratioCorrectHelicesFlipped = 0;
			double ratioCorrectSheets = 0;
			int ns = 0, nh=0;
			for (int i=0; i < numNodes; i++) {
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_HELIX) {
					nh++;
					if (solution[i]==SOLUTION[i]) {
						ratioCorrectHelices++;
					} else {
						helicesCorrect=false;
					}
					if (solution[i]==SOLUTION[i] || solution[i]==SOLUTION[max(0,i-1)] || solution[i]==SOLUTION[min(numNodes-1,i+1)]) {
						ratioCorrectHelicesFlipped++;
					} else {
						helicesCorrectFlipped=false;
					}

					/*
					helicesCorrect = helicesCorrect && (solution[i]==SOLUTION[i]);
					helicesCorrectFlipped = helicesCorrectFlipped && 
						(solution[i]==SOLUTION[i] || 
						solution[i]==SOLUTION[max(0,i-1)] ||
						solution[i]==SOLUTION[min(numNodes-1,i+1)]);
					*/
					
				}
				if (patternGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
					ns++;
					if (solution[i]==SOLUTION[i]) {
						ratioCorrectSheets++;
					} else {
						sheetsCorrect=false;
					}
					//sheetsCorrect = sheetsCorrect && (solution[i]==SOLUTION[i]);
				}
			}
			int incorrectHelixCount = (nh - (int)(ratioCorrectHelices+0.1))/2;
			int incorrectStrandCount = ns - (int)(ratioCorrectSheets+0.1);
			if (ns>0) {ratioCorrectSheets /= ns;}
			if (nh>0) {ratioCorrectHelices /= nh;}
			if (nh>0) {ratioCorrectHelicesFlipped /= nh;}

			char sheetChar;
			if (sheetsCorrect) {
				sheetChar='S';
			} else {
				sheetChar='-';
			}

			char helixChar;
			if (helicesCorrect) {
				helixChar='H';
			} else if (helicesCorrectFlipped) {
				helixChar='h';
			} else {
				helixChar='-';
			}

			if (extraMessages) {cout << "total edge cost is " << edgeCost << endl;}
			if (extraMessages) {cout << "total edge penalty cost is " << edgePenaltyCost << endl;}
			if (extraMessages) {cout << "total node cost is " << nodeCost << endl;}

			double helixPenalty = skippedHelixNodes/2 * MISSING_HELIX_PENALTY; 
			double sheetPenalty = skippedSheetNodes * MISSING_SHEET_PENALTY;

			if (extraMessages) {cout << "missing helices: " << skippedHelixNodes/2 << " contribute cost of " << helixPenalty << endl;}
			if (extraMessages) {cout << "missing sheets:  " << skippedSheetNodes << " contribute cost of " << sheetPenalty << endl;}
			if (extraMessages) {cout << "together, missing helices and sheets contribute cost of " << sheetPenalty + helixPenalty << endl;}
			if (extraMessages) {cout << "algorithm thinks missing helices and sheets should contribute cost of " << edgePenaltyCost << endl;}


			double cost = edgeCost + nodeCost + helixPenalty + sheetPenalty;

			if (extraMessages) {cout << "total cost is " << cost << endl;}
			cout << "C=" << loopCost+helixCost+sheetCost+skipPenaltyCost << "(" << helixChar << sheetChar << ")(h" << incorrectHelixCount << ",s" << incorrectStrandCount << ")(" << ratioCorrectHelices << "," << ratioCorrectHelicesFlipped << "," << ratioCorrectSheets << ")(L=" << loopCost << ",H=" << helixCost << ",S=" << sheetCost << ",P=" << skipPenaltyCost << ")";
			if (extraMessages) {cout << endl;}
		}

		void WongMatch15ConstrainedNoFuture::NormalizeGraphs() {
#ifdef VERBOSE
			printf("Normalizing Graphs\n");
#endif // VERBOSE


#ifdef VERBOSE
			printf("\tNormalizing the base graph from Angstroms to amino acids\nNormalized Graph:\n");
#endif // VERBOSE
			for(int i = 0; i < baseGraph->nodeCount; i++) {
				for(int j = 0; j < baseGraph->nodeCount; j++) {
					// base graph
					if(baseGraph->adjacencyMatrix[i][j][1] != MAXINT && baseGraph->adjacencyMatrix[i][j][0] == GRAPHEDGE_HELIX) {
						baseGraph->SetCost(i+1,j+1, baseGraph->adjacencyMatrix[i][j][1] / HELIX_C_ALPHA_TO_ANGSTROMS);
					} else if(baseGraph->adjacencyMatrix[i][j][1] != MAXINT) {
						baseGraph->SetCost(i+1,j+1, baseGraph->adjacencyMatrix[i][j][1] / LOOP_C_ALPHA_TO_ANGSTROMS);
					}
					// euclidean distance matrix
					if(baseGraph->adjacencyMatrix[i][j][0] == GRAPHEDGE_HELIX) {
						baseGraph->euclideanMatrix[i][j] = baseGraph->euclideanMatrix[i][j] / HELIX_C_ALPHA_TO_ANGSTROMS;
					} else {
						baseGraph->euclideanMatrix[i][j] = baseGraph->euclideanMatrix[i][j] / LOOP_C_ALPHA_TO_ANGSTROMS;
					}
				}
			}	


#ifdef VERBOSE
			baseGraph->PrintGraph();
#endif // VERBOSE
		}




		void WongMatch15ConstrainedNoFuture::NormalizeSheets() {
#ifdef VERBOSE
			printf("\tNormalizing the sheet nodes in the base graph based on sheet ratio\nNormalized Graph:\n");
#endif // VERBOSE
			// TODO: Also normalize the sheet capacity here?
			double totalSheetSize = 0;
			double totalStrandLength = 0;

			for(int i = 0; i < (int)baseGraph->skeletonHelixes.size(); i++) {
				if (baseGraph->skeletonHelixes[i]->geometricShapeType == GRAPHEDGE_SHEET) {
					totalSheetSize += (double)baseGraph->skeletonHelixes[i]->length;
#ifdef VERBOSE
					cout << "after sheet " << i << ", total sheet size is now " << totalSheetSize << endl;
#endif // VERBOSE
				}
			}

			for(int i = 0; i < (int)patternGraph->pdbStructures.size(); i++) {
				//cout << "strand " << i << " has type " << patternGraph->pdbStructures[i]->secondaryStructureType << endl;
				if(patternGraph->pdbStructures[i]->secondaryStructureType == GRAPHEDGE_SHEET) {
					totalStrandLength += patternGraph->pdbStructures[i]->GetLengthResidues();
					//totalStrandLength += patternGraph->pdbStructures[i]->GetLengthAngstroms();
#ifdef VERBOSE
					cout << "After adding strand " << i << " with length " << patternGraph->pdbStructures[i]->GetLengthResidues() << ", total strand length is now " << totalStrandLength << endl;
#endif // VERBOSE
				}
			}

			// scale the sheet sizes so that the units are in amino acids
			double ratio = totalStrandLength / totalSheetSize;
#ifdef VERBOSE
			cout << "sheet sizes must be scaled by a factor of " << ratio << endl;
#endif // VERBOSE


#ifdef VERBOSE
			printf("\tNormalizing the base graph sheets from voxels to scaled voxels\nNormalized Graph:\n");
#endif // VERBOSE
			for(int i = 0; i < baseGraph->nodeCount; i++) {
				if(baseGraph->adjacencyMatrix[i][i][1] != MAXINT && baseGraph->adjacencyMatrix[i][i][0] == GRAPHNODE_SHEET) {
					// scale the sheet weight to the # of amino acids
					baseGraph->SetCost(i+1, baseGraph->nodeWeights[i] * ratio);
					// take sqrt for matching algorithm
					baseGraph->SetCost(i+1, sqrt(baseGraph->nodeWeights[i]));
				}
			}	
			
			
#ifdef VERBOSE
			baseGraph->PrintGraph();
#endif // VERBOSE
		}



		unsigned long long WongMatch15ConstrainedNoFuture::EncodeNode(unsigned long long bitmap, int node) {
			if(node == -1)
				return bitmap;

			return (bitmap | ((unsigned long long)1 << node));

		}

		// code copied from LinkedNode::PrintNodeConcise
		// Adding a breakdown of the cost into loops, nodes, and helices
		void WongMatch15ConstrainedNoFuture::PrintNodeConcise(LinkedNode * node, int rank, bool endOfLine, bool printCostBreakdown) {
			bool used[MAX_NODES];
			int n1[MAX_NODES];
			int n2[MAX_NODES];
			int top = 0;
			for(int i = 0; i < MAX_NODES; i++) {
				used[i] = false;
			}

            std::vector<double> costList;
			LinkedNodeStub * currentNode = node;
            
			bool continueLoop = true;
			while(continueLoop) {
				if(currentNode->parentNode == NULL) {
					 break;
				}

				n1[top] = currentNode->n1Node;
				n2[top] = currentNode->n2Node;
				used[(int)currentNode->n1Node] = true;
				top++;
				currentNode = currentNode->parentNode;		
			}

			for(int i = 1; i <= node->depth; i++) {
				if(!used[i]) {
					n1[top] = i;
					n2[top] = -1;
					top++;
				}
			}

			int minIndex;
			int temp;
			for(int i = 0; i < top - 1; i++) {
				minIndex = i;
				for(int j = i+1; j < top; j++) {
					if(n1[minIndex] > n1[j]) {
						minIndex = j;
					}
				}
				temp = n1[minIndex];
				n1[minIndex] = n1[i];
				n1[i] = temp;

				temp = n2[minIndex];
				n2[minIndex] = n2[i];
				n2[i] = temp;
			}
            
//			if(node->IsUserSpecifiedSolution()) {
//				printf("**");
//			} else {
//				printf("  ");
//			}
			

//			if(rank != -1) {
//				printf("%d)", rank);
//			}
//			printf(" ");
			for(int i = 0; i < top; i++) {
				printf("%d ", n2[i]);
			}

			// print the cost of the current solution
			if(INCLUDE_STRANDS) {
				ComputeSolutionCost(n2,false);
			}
			for (int i = 0; i < MAX_NODES; i++){
				bestMatches[rank-1][i]=n2[i];
			}

			if(printCostBreakdown) {
				printf(" - %f = %f + %f", node->cost, node->costGStar, node->cost - node->costGStar);
			} else {
				//printf(" - %f", node->cost);
                printf(" %f", node->cost);
			}
			if(endOfLine) {
				printf("\n");
			}
		}
	}
}
#endif
