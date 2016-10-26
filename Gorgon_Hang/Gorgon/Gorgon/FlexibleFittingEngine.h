// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   An engine for flexible fitting of calpha atoms to volumes

// CVS Meta Information: 
//   $Source: /project/mm/cvs/graphics/ssa1/source/Gorgon/FlexibleFittingEngine.h,v $
//   $Revision: 1.9 $
//   $Date: 2012/06/18 13:27:00 $
//   $Author: burrowsd $
//   $State: Exp $
//
// History Log: 
//   $Log: FlexibleFittingEngine.h,v $
//   Revision 1.9  2010/06/18 13:27:00  coleman.r
//   Complete re-write of flexible fitting plugin
//
//   Revision 1.8  2010/07/29 20:21:59  coleman.r
//   gcc compile fix: gcc requires nested templates to end in "> >" not ">>"
//
//   Revision 1.7  2010/07/23 18:18:32  heiderp
//   Side chains now transform correctly.  PDB helices now color correctly and rigid initialization bug is fixed
//
//   Revision 1.6  2010/07/22 21:09:07  heiderp
//   Minor updates. Mostly commenting and removing extra material from CurveDeformer.h
//
//   Revision 1.5  2010/07/19 17:29:02  heiderp
//   LARGE update.  Added flexible fitting functionality, lots of logic in FlexibleFittingEngine.h
//
//   Revision 1.4  2010/06/17 19:31:47  ssa1
//   Visually displaying flexible fitting clusters.
//
//   Revision 1.3  2010/05/21 15:45:16  ssa1
//   Flexible fitting implemented in Gorgon
//
//   Revision 1.2  2010/05/20 21:55:53  ssa1
//   Rigid body alignment based on largest flexible cluster
//
//   Revision 1.1  2010/05/06 21:50:11  ssa1
//   Fixing performance bug when moving a volume
//

#ifndef FLEXIBLEFITTINGENGINE_H
#define FLEXIBLEFITTINGENGINE_H

#include <set>
#include <map>
#include <queue>
#include <vector>
#include <utility>
#include <math.h>
#include <tnt/tnt.h>
#include <tnt/jama_lu.h>
#include <tnt/jama_qr.h>
#include <tnt/jama_cholesky.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <MathTools/Matrix.h>
#include <MathTools/Vector3D.h>
#include <MathTools/LinearSolver.h>
#include <ProteinMorph/SSECorrespondenceNode.h>
#include <ProteinMorph/SSECorrespondenceFeature.h>
#include <ProteinMorph/SSECorrespondenceFinder.h>

#include <iomanip>
#include "Eigen/SVD"
#include "Eigen/Dense"
//#include "Eigen/SparseCore"
//#include "Eigen/SparseLU"
//#include "Eigen/SparseCholesky"
//#include "Eigen/SparseQR"
//#include "Eigen/IterativeLinearSolvers"

#if !defined(_WIN32) && !defined(_WIN64)
    #include <umfpack.h>
#endif

using namespace TNT;
using namespace JAMA;
using namespace wustl_mm::MathTools;

class FlexibleFittingEngine {
    public:
        // Constructors
        FlexibleFittingEngine();

        // Reset
        void resetCorrespondence();
        void resetDeformation();
        void resetDeformationVertices();
        void resetDeformationEdges();
        void resetDeformationHandles();
        void resetDeformationOptions();

        // Reset WSVD
		void resetWSVD();
        void resetIterationThreshold();

        // Setters
        void addCorrespondenceSourceFeature(const Vector3DFloat&, const Vector3DFloat&);
        void addCorrespondenceTargetFeature(const Vector3DFloat&, const Vector3DFloat&);
        void addCorrespondenceManual(const int&, const int&, const int&, const bool&);
        void addDeformationOriginalVertex(const Vector3DFloat&);
        void addDeformationVertex(const Vector3DFloat&);
        void addDeformationEdge(const int&, const int&);
        void addDeformationHandle(const int&, const Vector3DFloat&);
        void addDeformationRigidInitializationTransformation(const MatrixFloat&);
        void setCorrespondenceOptions(const int&, const float&, const float&, const float&, const float&);
        void setCorrespondenceDirection(const int&, const int&, const int&, const bool&);
        void setDeformationRigidInitialization(const bool&);
        void setDeformationNeighborhoodSizes(const int&, const int&, const int&);
        void setDeformationWeights(const double&, const double&);
        void setDeformationVerticesToOriginalVertices();
        void setDeformationVerticesToDeformedVertices();
        void mergeCorrespondenceClusters(const int&, const int&, const int&);

        void setIterationStopThreshold(const double & threshold);
        void setWSVDOptions(const double&);
        void addCorrespondencesByReadIn(const int&, const int&, const bool&);

        // Getters
        int getCorrespondenceAlignmentCount();
        int getCorrespondenceClusterCount(const int&);
        int getCorrespondenceFeatureCount(const int&, const int&);
        SSECorrespondenceNode getCorrespondence(const int&, const int&, const int&);
        MatrixFloat getCorrespondenceAlignmentTransformation(const int&);
        MatrixFloat getCorrespondenceClusterTransformation(const int&, const int&);
        MatrixFloat getCorrespondenceFeatureTransformation(const int&, const int&, const int&);
        vector<Vector3DFloat> getDeformedVertices();
        MatrixFloat getDeformationTransformation(const int &);

        double getShapeDistance() { return m_shapeDistance; }

        // Getters WSVD
        MatrixFloat getCorrespondenceFeatureTransformationWSVD(const int&, const int&, const int&);

        // Methods
        void calculateCorrespondences();
        void calculateNeighborVertexIndicesCaches(const int&);
        void calculateDeformationTransformations(const int&);
        void deformLaplacian();

        // Methods WSVD
        void computeCorrespondencePairWiseWeight(const int&);
        void deformWSVD();

    private:

		// For WSVD
		std::vector<std::vector<double>> m_correspondencePairWiseWeight;
		std::vector<int>  m_correspondenceSourceFeatureIndex;
		std::vector<int>  m_correspondenceTargetFeatureIndex;
        std::vector<bool> m_correspondenceFlipFlag;

        double m_bbRatio;
        double m_shapeDistance; // RMSD so far

        // Helpers
        Array2D<double> transformationForm(const Vector3DFloat&, const int&, const bool&);
        Array2D<double> transformationForm(const vector<Vector3DFloat>&, const int&, const bool&);
        
        // Matrix
        Array2D<double> inverse(const Array2D<double>&);
        Array2D<double> transpose(const Array2D<double>&);
        Array2D<double> linearSolve(const Array2D<double>&, const Array2D<double>&);
        double* sparseLinearSolve(const int&, const int&, const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&);

        // Correspondence Fields
        int maxCorrespondenceSolutions;

        float helixLengthDifference,
              helixCentroidDifference,
              jointAngleDifference,
              dihedralAngleDifference;

        SSECorrespondenceFinder sseCorrespondenceFinder;
        
        vector<SSECorrespondenceFeature> sourceFeatures,
                                         targetFeatures;

        vector< vector< vector<SSECorrespondenceNode> > > correspondences;

        // Deformation Fields
        bool rigidInitialization;

        int laplacianNeighborhoodSize,
            transformationNeighborhoodSize,
            deformationNeighborhoodSize;
       
        double fittingWeight,
              distortionWeight;

        vector<int> handleIndices;
        
        vector<Vector3DFloat> originalVertices,
                              vertices,
                              targetVertices,
                              deformedVertices;

        vector<MatrixFloat> rigidInitializationTransformations,
                            deformationTransformations;

        vector< vector<int> > laplacianNeighborVertexIndicesCache,
                              transformationNeighborVertexIndicesCache,
                              deformationNeighborVertexIndicesCache;

        vector< pair<int, int> > edges;
        
        LinearSolver linearSolver;
};

FlexibleFittingEngine::FlexibleFittingEngine() {
    resetCorrespondence();
    resetDeformation();
}

void FlexibleFittingEngine::resetWSVD()
{
    m_bbRatio = 0.3;

	m_correspondenceSourceFeatureIndex.clear();
	m_correspondenceTargetFeatureIndex.clear();    
    m_correspondenceFlipFlag.clear();
    m_correspondencePairWiseWeight.clear();
}

void FlexibleFittingEngine::resetCorrespondence() {
    maxCorrespondenceSolutions = 1;

    jointAngleDifference    = 0.0;
    dihedralAngleDifference = 0.0;
    helixLengthDifference   = 0.0;
    helixCentroidDifference = 0.0;
    
    sourceFeatures.clear();
    targetFeatures.clear();

    correspondences.clear();
}

void FlexibleFittingEngine::resetDeformation() {
    resetDeformationVertices();
    resetDeformationEdges();
    resetDeformationHandles();
    resetDeformationOptions();

	resetWSVD();
}

void FlexibleFittingEngine::resetDeformationVertices() {
    originalVertices.clear();
    vertices.clear();
    deformedVertices.clear();
    
    rigidInitializationTransformations.clear();
    deformationTransformations.clear();
}

void FlexibleFittingEngine::resetDeformationEdges() {
    edges.clear();
    
    laplacianNeighborVertexIndicesCache.clear();
    transformationNeighborVertexIndicesCache.clear();
    deformationNeighborVertexIndicesCache.clear();
}

void FlexibleFittingEngine::resetDeformationHandles() {
    handleIndices.clear();
    targetVertices.clear();
}

void FlexibleFittingEngine::resetDeformationOptions() {
    rigidInitialization = false;

    laplacianNeighborhoodSize      = 1;
    transformationNeighborhoodSize = 1;
    deformationNeighborhoodSize    = 1;

    fittingWeight    = 1.0;
    distortionWeight = 1.0;
}

void FlexibleFittingEngine::addCorrespondenceSourceFeature(const Vector3DFloat& startVertex, const Vector3DFloat& endVertex) {
    sourceFeatures.push_back(SSECorrespondenceFeature(startVertex, endVertex));
}

void FlexibleFittingEngine::addCorrespondenceTargetFeature(const Vector3DFloat& startVertex, const Vector3DFloat& endVertex) {
    targetFeatures.push_back(SSECorrespondenceFeature(startVertex, endVertex));
}

void FlexibleFittingEngine::addCorrespondenceManual(const int& alignmentIndex, const int& pIndex, const int& qIndex, const bool& isForward) {
    vector< vector<SSECorrespondenceNode> > newAlignment;
    for (int clusterIndex = 0; clusterIndex < getCorrespondenceClusterCount(alignmentIndex); clusterIndex++) {
        vector<SSECorrespondenceNode> newCluster;
        for (int featureIndex = 0; featureIndex < getCorrespondenceFeatureCount(alignmentIndex, clusterIndex); featureIndex++) {
            SSECorrespondenceNode correspondence = getCorrespondence(alignmentIndex, clusterIndex, featureIndex);
            if (correspondence.GetPIndex() != pIndex && correspondence.GetQIndex() != qIndex) {
                newCluster.push_back(correspondence);
            }
        }
        if (!newCluster.empty()) {
            newAlignment.push_back(newCluster);
        }
    }
    newAlignment.push_back(vector<SSECorrespondenceNode>(1, SSECorrespondenceNode(pIndex, qIndex, isForward))); 
    
    correspondences[alignmentIndex] = newAlignment;
}

void FlexibleFittingEngine::addDeformationOriginalVertex(const Vector3DFloat& vertex) {
    originalVertices.push_back(vertex);
}

void FlexibleFittingEngine::addDeformationVertex(const Vector3DFloat& vertex) {
    vertices.push_back(vertex);
}

void FlexibleFittingEngine::addDeformationEdge(const int& firstVertexIndex, const int& secondVertexIndex) {
    edges.push_back(make_pair(firstVertexIndex, secondVertexIndex));
}

void FlexibleFittingEngine::addDeformationHandle(const int& vertexIndex, const Vector3DFloat& targetVertex) {
    handleIndices.push_back(vertexIndex);
    targetVertices.push_back(targetVertex);
}

void FlexibleFittingEngine::addDeformationRigidInitializationTransformation(const MatrixFloat& rigidInitializationTransformation) {
    rigidInitializationTransformations.push_back(rigidInitializationTransformation);
}

void FlexibleFittingEngine::setCorrespondenceOptions(const int& newMaxCorrespondenceSolutions, const float& newHelixLengthDifference, const float& newHelixCentroidDifference, const float& newJointAngleDifference, const float& newDihedralAngleDifference) {
    maxCorrespondenceSolutions = newMaxCorrespondenceSolutions;
    
    helixLengthDifference      = newHelixLengthDifference;
    helixCentroidDifference    = newHelixCentroidDifference;

    jointAngleDifference       = (newJointAngleDifference    * PI) / 180.0;
    dihedralAngleDifference    = (newDihedralAngleDifference * PI) / 180.0;
}

void FlexibleFittingEngine::setWSVDOptions(const double & bbRatio) 
{
    m_bbRatio = bbRatio;
}

void FlexibleFittingEngine::setCorrespondenceDirection(const int& alignmentIndex, const int& clusterIndex, const int& featureIndex, const bool& isForward) {
    correspondences[alignmentIndex][clusterIndex][featureIndex].SetForward(isForward);
}

void FlexibleFittingEngine::setDeformationRigidInitialization(const bool& newRigidInitialization) {
    rigidInitialization = newRigidInitialization;
}

void FlexibleFittingEngine::setDeformationNeighborhoodSizes(const int& newLaplacianNeighborhoodSize, const int& newTransformationNeighborhoodSize, const int& newDeformationNeighborhoodSize) {
    laplacianNeighborhoodSize      = newLaplacianNeighborhoodSize;
    transformationNeighborhoodSize = newTransformationNeighborhoodSize;
    deformationNeighborhoodSize    = newDeformationNeighborhoodSize;
}

void FlexibleFittingEngine::setDeformationWeights(const double& newFittingWeight, const double& newDistortionWeight) {
    fittingWeight    = newFittingWeight;
    distortionWeight = newDistortionWeight;
}

void FlexibleFittingEngine::setDeformationVerticesToOriginalVertices() {
    vertices = originalVertices;
}

void FlexibleFittingEngine::setDeformationVerticesToDeformedVertices() {
    vertices = deformedVertices;
}

void FlexibleFittingEngine::mergeCorrespondenceClusters(const int& alignmentIndex, const int& firstClusterIndex, const int& secondClusterIndex) {
    // Add Second Cluster Features to First Cluster
    correspondences[alignmentIndex][firstClusterIndex].insert(correspondences[alignmentIndex][firstClusterIndex].end(), correspondences[alignmentIndex][secondClusterIndex].begin(), correspondences[alignmentIndex][secondClusterIndex].end());

    // Delete Second Cluster
    correspondences[alignmentIndex].erase(correspondences[alignmentIndex].begin() + secondClusterIndex);
}

int FlexibleFittingEngine::getCorrespondenceAlignmentCount() {
    return correspondences.size();
}

int FlexibleFittingEngine::getCorrespondenceClusterCount(const int& alignmentIndex) {
    return correspondences[alignmentIndex].size();
}

int FlexibleFittingEngine::getCorrespondenceFeatureCount(const int& alignmentIndex, const int& clusterIndex) {
    return correspondences[alignmentIndex][clusterIndex].size();
}

SSECorrespondenceNode FlexibleFittingEngine::getCorrespondence(const int& alignmentIndex, const int& clusterIndex, const int& featureIndex) {
    return correspondences[alignmentIndex][clusterIndex][featureIndex];
}

MatrixFloat FlexibleFittingEngine::getCorrespondenceAlignmentTransformation(const int& alignmentIndex) {
    vector<SSECorrespondenceNode> alignmentCorrespondences;
    for (int clusterIndex = 0; clusterIndex < getCorrespondenceClusterCount(alignmentIndex); clusterIndex++) {
        for (int featureIndex = 0; featureIndex < getCorrespondenceFeatureCount(alignmentIndex, clusterIndex); featureIndex++) {
            alignmentCorrespondences.push_back(getCorrespondence(alignmentIndex, clusterIndex, featureIndex));
        }
    }

    MatrixFloat alignmentTransformation = sseCorrespondenceFinder.GetTransform(alignmentCorrespondences, 10);

    return alignmentTransformation;
}

MatrixFloat FlexibleFittingEngine::getCorrespondenceClusterTransformation(const int& alignmentIndex, const int& clusterIndex) {
    vector<SSECorrespondenceNode> clusterCorrespondences;
    if (getCorrespondenceFeatureCount(alignmentIndex, clusterIndex) == 1) {
        // Get Correspondence
        SSECorrespondenceNode correspondence = getCorrespondence(alignmentIndex, clusterIndex, 0);

        int sourceFeatureIndex = correspondence.GetPIndex(),
            targetFeatureIndex = correspondence.GetQIndex();
 
        // Get Source Feature Centroid
        Vector3DFloat sourceFeatureCentroid = sourceFeatures[sourceFeatureIndex].GetCentroid();

        // Calculate Distances to Other Source Features
		vector< boost::tuple<double, int, int> > closestSourceFeatures;
        for (int otherClusterIndex = 0; otherClusterIndex < getCorrespondenceClusterCount(alignmentIndex); otherClusterIndex++) {
            if (otherClusterIndex == clusterIndex) {
                continue;
            }

            for (int otherFeatureIndex = 0; otherFeatureIndex < getCorrespondenceFeatureCount(alignmentIndex, otherClusterIndex); otherFeatureIndex++) {
                // Get Other Correspondence
                SSECorrespondenceNode otherCorrespondence = getCorrespondence(alignmentIndex, otherClusterIndex, otherFeatureIndex);

                // Get Other Source Feature Index
                int otherSourceFeatureIndex = otherCorrespondence.GetPIndex();
                
                // Get Other Source Feature Centroid
                Vector3DFloat otherSourceFeatureCentroid = sourceFeatures[otherSourceFeatureIndex].GetCentroid();

                // Calculate Distance
                double distance = (otherSourceFeatureCentroid - sourceFeatureCentroid).Length();

                // Store Distance and Correspondence
				closestSourceFeatures.push_back(boost::tuple<double, int, int>(distance, otherClusterIndex, otherFeatureIndex));
            }
        }

        // Sort Distances 
        sort(closestSourceFeatures.begin(), closestSourceFeatures.end());

        // Build Temporary Cluster
        int clusterCorrespondencesSize = 3;
        for (int clusterCorrespondenceIndex = 0; clusterCorrespondenceIndex < clusterCorrespondencesSize; clusterCorrespondenceIndex++) {
            double distance       = closestSourceFeatures[clusterCorrespondenceIndex].get<0>();

            int otherClusterIndex = closestSourceFeatures[clusterCorrespondenceIndex].get<1>(),
                otherFeatureIndex = closestSourceFeatures[clusterCorrespondenceIndex].get<2>();

            clusterCorrespondences.push_back(getCorrespondence(alignmentIndex, otherClusterIndex, otherFeatureIndex));
        }
    }
    else {
        clusterCorrespondences = correspondences[alignmentIndex][clusterIndex];
    }

    MatrixFloat clusterTransformation = sseCorrespondenceFinder.GetTransform(clusterCorrespondences, 10);

    return clusterTransformation;
}

MatrixFloat FlexibleFittingEngine::getCorrespondenceFeatureTransformation(const int& alignmentIndex, const int& clusterIndex, const int& featureIndex) {
    // Get Correspondence
    SSECorrespondenceNode correspondence = getCorrespondence(alignmentIndex, clusterIndex, featureIndex);

    int sourceFeatureIndex = correspondence.GetPIndex(),
        targetFeatureIndex = correspondence.GetQIndex();
    
    bool flip = !correspondence.IsForward();

    // Cluster Transformation
    MatrixFloat clusterTransformation = getCorrespondenceClusterTransformation(alignmentIndex, clusterIndex);

    // Translation
    Vector3DFloat sourceCentroid    = sourceFeatures[sourceFeatureIndex].GetCentroid().Transform(clusterTransformation),
                  targetCentroid    = targetFeatures[targetFeatureIndex].GetCentroid(),
                  translationVector = targetCentroid - sourceCentroid;

    MatrixFloat translationTransformation = MatrixFloat::Identity(4);
    translationTransformation.SetValue(translationVector[0], 0, 3);
    translationTransformation.SetValue(translationVector[1], 1, 3);
    translationTransformation.SetValue(translationVector[2], 2, 3);

    // Rotation
    Vector3DFloat sourceStartPosition = sourceFeatures[sourceFeatureIndex].GetEndPoint(0).Transform(translationTransformation * clusterTransformation),
                  sourceEndPosition   = sourceFeatures[sourceFeatureIndex].GetEndPoint(1).Transform(translationTransformation * clusterTransformation),
                  targetStartPosition = targetFeatures[targetFeatureIndex].GetEndPoint(0),
                  targetEndPosition   = targetFeatures[targetFeatureIndex].GetEndPoint(1),
                  sourceVector        = sourceEndPosition - sourceStartPosition,
                  targetVector        = targetEndPosition - targetStartPosition;

    sourceVector.Normalize();
    targetVector.Normalize();
    
    Vector3DFloat normal = sourceVector ^ targetVector;
    normal.Normalize();

    double angle = acos(sourceVector * targetVector);
    if (angle > (PI / 2)) {
        angle -= PI;
    }
    if (flip) {
        angle *= -1;
    }

    Vector3DFloat row0(1.0, 0.0, 0.0),
                  row1(0.0, 1.0, 0.0),
                  row2(0.0, 0.0, 1.0);

    row0.Rotate(normal, angle);
    row1.Rotate(normal, angle);
    row2.Rotate(normal, angle);

    MatrixFloat rotationTransformation = MatrixFloat::Identity(4);
    rotationTransformation.SetValue(row0[0], 0, 0); rotationTransformation.SetValue(row0[1], 0, 1); rotationTransformation.SetValue(row0[2], 0, 2);
    rotationTransformation.SetValue(row1[0], 1, 0); rotationTransformation.SetValue(row1[1], 1, 1); rotationTransformation.SetValue(row1[2], 1, 2);
    rotationTransformation.SetValue(row2[0], 2, 0); rotationTransformation.SetValue(row2[1], 2, 1); rotationTransformation.SetValue(row2[2], 2, 2);

    // Rotation Origin Offset
    sourceCentroid = sourceFeatures[sourceFeatureIndex].GetCentroid().Transform(translationTransformation * clusterTransformation);
    
    MatrixFloat translationToOriginTransformation = MatrixFloat::Identity(4);
    translationToOriginTransformation.SetValue(-sourceCentroid[0], 0, 3);
    translationToOriginTransformation.SetValue(-sourceCentroid[1], 1, 3);
    translationToOriginTransformation.SetValue(-sourceCentroid[2], 2, 3);

    MatrixFloat translationFromOriginTransformation = MatrixFloat::Identity(4);
    translationFromOriginTransformation.SetValue(sourceCentroid[0], 0, 3);
    translationFromOriginTransformation.SetValue(sourceCentroid[1], 1, 3);
    translationFromOriginTransformation.SetValue(sourceCentroid[2], 2, 3);

    return translationFromOriginTransformation * rotationTransformation * translationToOriginTransformation * translationTransformation * clusterTransformation;
}

MatrixFloat FlexibleFittingEngine::getCorrespondenceFeatureTransformationWSVD(const int& alignmentIndex, const int& clusterIndex, const int& featureIndex) 
{
    // Get Correspondence
    //SSECorrespondenceNode correspondence = getCorrespondence(alignmentIndex, clusterIndex, featureIndex);    

	// cluster1[f1,f2,f3], cluster2[f7,f8], cluster3[f4,f5], ...
	int corrIndex = 0;
	// Locate the correspondence index in weight matrix
    for (int i = 0; i < clusterIndex; ++i)
    {        
        corrIndex += getCorrespondenceFeatureCount(alignmentIndex, i);
    }
    corrIndex += featureIndex;

	//////// Construct the matrix for wsvd ////////

	// Compute the center of source points and target points
    Vector3DDouble sourceCenter(0, 0, 0);
    Vector3DDouble targetCenter(0, 0, 0);
    double weightSum = 0;
    std::vector<double> tempWeight;
    for (int i = 0; i < m_correspondenceSourceFeatureIndex.size(); ++i)
    {
        Vector3DDouble sourceEnd1(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).X(),
                                  sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).Y(),
                                  sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).Z());
        Vector3DDouble sourceEnd2(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).X(),
                                  sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).Y(),
                                  sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).Z());

        Vector3DDouble targetEnd1(targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).X(),
                                  targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Y(),
                                  targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Z());
        Vector3DDouble targetEnd2(targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).X(),
                                  targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Y(),
                                  targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Z());

        // Since we only have upper triangle of the pair wist weight matrix
        if (i > corrIndex)
        {
            sourceCenter += (sourceEnd1 + sourceEnd2) * m_correspondencePairWiseWeight[corrIndex][i];
            targetCenter += (targetEnd1 + targetEnd2) * m_correspondencePairWiseWeight[corrIndex][i];
            weightSum += m_correspondencePairWiseWeight[corrIndex][i];
            tempWeight.push_back(m_correspondencePairWiseWeight[corrIndex][i]);
        }
        else
        {
            sourceCenter += (sourceEnd1 + sourceEnd2) * m_correspondencePairWiseWeight[i][corrIndex];
            targetCenter += (targetEnd1 + targetEnd2) * m_correspondencePairWiseWeight[i][corrIndex];
            weightSum += m_correspondencePairWiseWeight[i][corrIndex];
            tempWeight.push_back(m_correspondencePairWiseWeight[i][corrIndex]);
        }
    }

    weightSum *= 2.0;
    sourceCenter = Vector3DDouble(sourceCenter.X()/weightSum, sourceCenter.Y()/weightSum, sourceCenter.Z()/weightSum);
    targetCenter = Vector3DDouble(targetCenter.X()/weightSum, targetCenter.Y()/weightSum, targetCenter.Z()/weightSum);

    //std::cout << "weight list:\n";
    //for (auto & i : tempWeight)
    //{
    //    std::cout << i << ", ";
    //}
    //std::cout << "\n";

    //std::cout << "sc:\n";
    //sourceCenter.Print();
    //std::cout << "tc:\n";
    //targetCenter.Print();
    //std::cout << "weight sum: " << weightSum << "\n";
    //std::cout << "correspondence.IsForward(): " << correspondence.IsForward() << "\n";

	// Construct the "weighted matrix" for source and target
	int pointDim = 3;

    // "m_correspondencePairWiseWeightMatrix" is the number of correspondences
    // each correspondence's line segment has two end points
    Eigen::MatrixXd sourcePointsMatrix(pointDim, 2*m_correspondencePairWiseWeight.size());
    Eigen::MatrixXd targetPointsMatrixTrans(2*m_correspondencePairWiseWeight.size(), pointDim);

    for (int i = 0; i < m_correspondencePairWiseWeight.size(); ++i)
	{
        double weight = i > corrIndex ? m_correspondencePairWiseWeight[corrIndex][i] : m_correspondencePairWiseWeight[i][corrIndex];

        // Source start point
        sourcePointsMatrix(0, 2*i) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).X() - sourceCenter.X()) * weight;
        sourcePointsMatrix(1, 2*i) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).Y() - sourceCenter.Y()) * weight;
        sourcePointsMatrix(2, 2*i) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(0).Z() - sourceCenter.Z()) * weight;
        // Source end point
        sourcePointsMatrix(0, 2*i + 1) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).X() - sourceCenter.X()) * weight;
        sourcePointsMatrix(1, 2*i + 1) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).Y() - sourceCenter.Y()) * weight;
        sourcePointsMatrix(2, 2*i + 1) = (double)(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetEndPoint(1).Z() - sourceCenter.Z()) * weight;

        // Make sure end points match correctly by checking matching direction
        if (m_correspondenceFlipFlag[i])
        {
            // Target start point
            targetPointsMatrixTrans(2*i, 0) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).X() - targetCenter.X();
            targetPointsMatrixTrans(2*i, 1) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Y() - targetCenter.Y();
            targetPointsMatrixTrans(2*i, 2) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Z() - targetCenter.Z();
            // Target end point
            targetPointsMatrixTrans(2*i + 1, 0) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).X() - targetCenter.X();
            targetPointsMatrixTrans(2*i + 1, 1) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Y() - targetCenter.Y();
            targetPointsMatrixTrans(2*i + 1, 2) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Z() - targetCenter.Z();
        }
        else
        {
            // Target start point
            targetPointsMatrixTrans(2*i, 0) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).X() - targetCenter.X();
            targetPointsMatrixTrans(2*i, 1) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Y() - targetCenter.Y();
            targetPointsMatrixTrans(2*i, 2) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(1).Z() - targetCenter.Z();
            // Target end point
            targetPointsMatrixTrans(2*i + 1, 0) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).X() - targetCenter.X();
            targetPointsMatrixTrans(2*i + 1, 1) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Y() - targetCenter.Y();
            targetPointsMatrixTrans(2*i + 1, 2) = (double)targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetEndPoint(0).Z() - targetCenter.Z();
        }
	}

    // Solve SVD
    //Eigen::MatrixXd upTRansM(3, 3);
    //upTRansM = sourcePointsMatrix*targetPointsMatrixTrans;    
    Eigen::JacobiSVD<Eigen::MatrixXd> svdM(sourcePointsMatrix*targetPointsMatrixTrans, Eigen::ComputeThinU | Eigen::ComputeThinV);

    //MatrixFloat svd = MatrixFloat(3, 3);
    //MatrixFloat u = MatrixFloat(3, 3);
    //MatrixFloat v = MatrixFloat(3, 3);
    //MatrixFloat w = MatrixFloat(3, 3);

    //svd.SetValue(upTRansM(0, 0), 0, 0); svd.SetValue(upTRansM(0, 1), 0, 1); svd.SetValue(upTRansM(0, 2), 0, 2);
    //svd.SetValue(upTRansM(1, 0), 1, 0); svd.SetValue(upTRansM(1, 1), 1, 1); svd.SetValue(upTRansM(1, 2), 1, 2);
    ////svd.SetValue(upTRansM(2, 0), 2, 0); svd.SetValue(upTRansM(2, 1), 2, 1); svd.SetValue(upTRansM(2, 2), 2, 2);
    //svd.SetValue(1, 0, 0); svd.SetValue(2, 0, 1); svd.SetValue(3, 0, 2);
    //svd.SetValue(4, 1, 0); svd.SetValue(5, 1, 1); svd.SetValue(6, 1, 2);
    //svd.SetValue(7, 2, 0); svd.SetValue(8, 2, 1); svd.SetValue(9, 2, 2);

    //svd.SingularValueDecomposition(u, w, v);

    //std::cout << "Mathtool U:";
    //u.Print(false);
    //std::cout << "Mathtool V:";
    //v.Print(false)

    //MatrixFloat rot = v * u.Transpose();
    //if (rot.Determinant() < 0)
    //{
    //    v.SetValue(-v.GetValue(0,2), 0, 2);
    //    v.SetValue(-v.GetValue(0, 2), 0, 2);
    //    v.SetValue(-v.GetValue(0, 2), 0, 2);
    //    rot = v * u.Transpose();
    //}

    //MatrixFloat rot4 = MatrixFloat(4, 4);
    //for (unsigned int i = 0; i < 3; i++) {
    //    for (unsigned int j = 0; j < 3; j++) {
    //        rot4.SetValue(rot.GetValue(i, j), i, j);
    //    }
    //}

    //std::cout << "rotation matrix second vesion:";
    //rot4.Print(false);

    // rot4.SetValue(1.0f, 3, 3);
    
    //std::cout << "upTRansM matrix one value:\n" << std::setprecision(30) << upTRansM(1, 1) << "\n";
    //std::cout << "svd matrix:\n" << std::setprecision(20) << sourcePointsMatrix*targetPointsMatrixTrans << "\n";
    //std::cout << "U matrix:\n" << svdM.matrixU() << "\n";
    //std::cout << "V matrix:\n" << svdM.matrixV() << "\n";
    //std::cout << "p:\n" << std::setprecision(20) << sourcePointsMatrix << "\n";
    //std::cout << "qTrans:\n" << std::setprecision(20) << targetPointsMatrixTrans << "\n";

    Eigen::Matrix3d matrixV = svdM.matrixV();
    Eigen::Matrix3d RotationM = matrixV * svdM.matrixU().transpose();
    // Avoid reflection from WSVD
    if (RotationM.determinant() < 0)
    {
        matrixV(0, 2) = -svdM.matrixV()(0, 2);
        matrixV(1, 2) = -svdM.matrixV()(1, 2);
        matrixV(2, 2) = -svdM.matrixV()(2, 2);
        // Re-compute the rotationMation
        RotationM = matrixV * svdM.matrixU().transpose();
    }

    //std::cout << "rotation matrix:\n" << std::setprecision(20) << RotationM << "\n";
    //std::cout << "source center:\n" << sourceCenter.X() << ", " << sourceCenter.Y() << ", " << sourceCenter.Z() << "\n";
    //std::cout << "target center:\n" << targetCenter.X() << ", " << targetCenter.Y() << ", " << targetCenter.Z() << "\n";

    MatrixFloat rotationTransformation = MatrixFloat::Identity(4);
    rotationTransformation.SetValue(RotationM(0, 0), 0, 0); rotationTransformation.SetValue(RotationM(0, 1), 0, 1); rotationTransformation.SetValue(RotationM(0, 2), 0, 2);
    rotationTransformation.SetValue(RotationM(1, 0), 1, 0); rotationTransformation.SetValue(RotationM(1, 1), 1, 1); rotationTransformation.SetValue(RotationM(1, 2), 1, 2);
    rotationTransformation.SetValue(RotationM(2, 0), 2, 0); rotationTransformation.SetValue(RotationM(2, 1), 2, 1); rotationTransformation.SetValue(RotationM(2, 2), 2, 2);

    // Rotation Origin Offset
    MatrixFloat translationToOriginTransformation = MatrixFloat::Identity(4);
    translationToOriginTransformation.SetValue(-sourceCenter.X(), 0, 3);
    translationToOriginTransformation.SetValue(-sourceCenter.Y(), 1, 3);
    translationToOriginTransformation.SetValue(-sourceCenter.Z(), 2, 3);

    MatrixFloat translationFromOriginTransformation = MatrixFloat::Identity(4);
    translationFromOriginTransformation.SetValue(targetCenter.X(), 0, 3);
    translationFromOriginTransformation.SetValue(targetCenter.Y(), 1, 3);
    translationFromOriginTransformation.SetValue(targetCenter.Z(), 2, 3);

    // Output deformed points
    //int sourceIndex = correspondence.GetPIndex();
    //int targetIndex = correspondence.GetQIndex();
    //std::cout << "s\n";
    //std::cout << std::setprecision(30) << sourceFeatures[sourceIndex].GetEndPoint(0).X() << ", " << sourceFeatures[sourceIndex].GetEndPoint(0).Y() << ", " << sourceFeatures[sourceIndex].GetEndPoint(0).Z() << "\n";
    //std::cout << std::setprecision(30) << sourceFeatures[sourceIndex].GetEndPoint(1).X() << ", " << sourceFeatures[sourceIndex].GetEndPoint(1).Y() << ", " << sourceFeatures[sourceIndex].GetEndPoint(1).Z() << "\n";

    //std::cout << "t\n";
    //std::cout << std::setprecision(30) << targetFeatures[targetIndex].GetEndPoint(0).X() << ", " << targetFeatures[targetIndex].GetEndPoint(0).Y() << ", " << targetFeatures[targetIndex].GetEndPoint(0).Z() << "\n";
    //std::cout << std::setprecision(30) << targetFeatures[targetIndex].GetEndPoint(1).X() << ", " << targetFeatures[targetIndex].GetEndPoint(1).Y() << ", " << targetFeatures[targetIndex].GetEndPoint(1).Z() << "\n";

    return translationFromOriginTransformation * rotationTransformation * translationToOriginTransformation;
    //return translationFromOriginTransformation * rot4 * translationToOriginTransformation;
}

vector<Vector3DFloat> FlexibleFittingEngine::getDeformedVertices() {
    return deformedVertices;
}

MatrixFloat FlexibleFittingEngine::getDeformationTransformation(const int& vertexIndex) {
    return deformationTransformations[vertexIndex];
}

// Test the helix-guided deformation
void FlexibleFittingEngine::deformWSVD() 
{
    deformedVertices.clear();
    deformedVertices.reserve(vertices.size());
    for (int vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex++)
    {
        MatrixFloat rigidInitializationTransformation(4, 4);
        if (rigidInitialization) 
            rigidInitializationTransformation = rigidInitializationTransformations[vertexIndex];

        // Lookup Vertex
        Vector3DFloat vertex = vertices[vertexIndex];
        if (rigidInitialization)
            vertex = vertex.Transform(rigidInitializationTransformation);

        deformedVertices.push_back(vertex);
    }
}

void FlexibleFittingEngine::deformLaplacian() {
    // Initialization
    int numVertices = vertices.size(),
        numHandles = handleIndices.size(),
        dimensions = 3,
        ataDim = numVertices * dimensions;

    // Initialization for Eigen linear solver
    //std::vector< Eigen::Triplet<double> > tripletList;
    Eigen::MatrixXd ata = Eigen::MatrixXd::Constant(numVertices * dimensions, numVertices * dimensions, 0.0);
    Eigen::VectorXd b;
    b.setConstant(numVertices * dimensions, 0.0);

    //map<int, double> ata;

    #if defined(_WIN32) || defined(_WIN64)
//        Array2D<double> AtA(numVertices * dimensions, numVertices * dimensions, 0.0),
//                        AtB(numVertices * dimensions, 1, 0.0);
    #else
        vector<int>    AtARowIndices,
                       AtAColumnIndices;
        vector<double> AtAValues,
                       AtB(numVertices * dimensions, 0.0);
    #endif

    //cout << "Done initializing matrices!" << endl;
   
    // Fitting Terms
    for (int handleIndex = 0; handleIndex < numHandles; handleIndex++) {
        int vertexIndex = handleIndices[handleIndex];
        Vector3DFloat targetVertex = targetVertices[handleIndex];

        for (int dimension = 0; dimension < dimensions; dimension++) {
            #if defined(_WIN32) || defined(_WIN64)
//                AtA[vertexIndex + (dimension * numVertices)][vertexIndex + (dimension * numVertices)] += fittingWeight * fittingWeight;
//                AtB[vertexIndex + (dimension * numVertices)][0] += fittingWeight * fittingWeight * targetVertex[dimension];
            #else
                AtARowIndices.push_back(vertexIndex + (dimension * numVertices));
                AtAColumnIndices.push_back(vertexIndex + (dimension * numVertices));
                AtAValues.push_back(fittingWeight * fittingWeight);
                AtB[vertexIndex + (dimension * numVertices)] += fittingWeight * fittingWeight * targetVertex[dimension];
            #endif

            // For Eigen linear solver
            //ata[(vertexIndex + (dimension * numVertices))*ataDim + vertexIndex + (dimension * numVertices)] += fittingWeight * fittingWeight;
            ata(vertexIndex + dimension * numVertices, vertexIndex + dimension * numVertices) += fittingWeight * fittingWeight;
            //tripletList.push_back(Eigen::Triplet<double>(vertexIndex + (dimension * numVertices), vertexIndex + (dimension * numVertices), fittingWeight * fittingWeight));
            b[vertexIndex + (dimension * numVertices)] += fittingWeight * fittingWeight * targetVertex[dimension];
        };
    }
    
    //cout << "Done adding fitting terms!" << endl;
    
    //std::cout << "LP neighbor vertices: " << laplacianNeighborVertexIndicesCache[handleIndices[numHandles/2]].size() << "\n";
    //std::cout << "Trans neighbor vertices: " << transformationNeighborVertexIndicesCache[handleIndices[numHandles/2]].size() << "\n";

    // Distortion Terms
    for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
        // Lookup Rigid Initialization Transformation
        MatrixFloat rigidInitializationTransformation(4, 4);
        if (rigidInitialization) {
            rigidInitializationTransformation = rigidInitializationTransformations[vertexIndex];
        }

        // Lookup Vertex
        Vector3DFloat vertex = vertices[vertexIndex]; 
        if (rigidInitialization) {
            vertex = vertex.Transform(rigidInitializationTransformation);
        }

        // Calculate Laplacian
        vector<int> laplacianNeighborVertexIndices = laplacianNeighborVertexIndicesCache[vertexIndex];
        int laplacianNumNeighbors = laplacianNeighborVertexIndices.size();
        if (laplacianNumNeighbors == 0) {
            cerr << "Unable to calculate laplacian, vertex with 0 neighbors!" << endl;
            exit(-1);
        }

        Vector3DFloat neighborSum(0.0, 0.0, 0.0);
        for (vector<int>::iterator neighborVertexIndexIterator = laplacianNeighborVertexIndices.begin(); neighborVertexIndexIterator != laplacianNeighborVertexIndices.end(); neighborVertexIndexIterator++) {
            Vector3DFloat neighborVertex = vertices[*neighborVertexIndexIterator];
            if (rigidInitialization) {
                neighborVertex = neighborVertex.Transform(rigidInitializationTransformation);
            }
            neighborSum += neighborVertex;
        }
        Vector3DFloat laplacian = vertex - (neighborSum * (1.0 / laplacianNumNeighbors));
  
        // Transformation Compensation
        vector<int> transformationNeighborVertexIndices = transformationNeighborVertexIndicesCache[vertexIndex];
        vector<Vector3DFloat> transformationNeighborhoodVertices(1, vertex);
        for (vector<int>::iterator neighborVertexIndexIterator = transformationNeighborVertexIndices.begin(); neighborVertexIndexIterator != transformationNeighborVertexIndices.end(); neighborVertexIndexIterator++) {
            Vector3DFloat neighborVertex = vertices[*neighborVertexIndexIterator];
            if (rigidInitialization) {
                neighborVertex = neighborVertex.Transform(rigidInitializationTransformation);
            }
            transformationNeighborhoodVertices.push_back(neighborVertex);
        }
        
        Array2D<double> C  = transformationForm(transformationNeighborhoodVertices, dimensions, true),
                        D  = transformationForm(laplacian, dimensions, false),
                        Ct = transpose(C),
                        T  = matmult(matmult(D, inverse(matmult(Ct, C))), Ct);

        // Construct Row
        for (int dimensionOuter = 0; dimensionOuter < dimensions; dimensionOuter++) {
            map<int, double> a;
            
            a[vertexIndex + (dimensionOuter * numVertices)] += distortionWeight;
            for (int dimensionInner = 0; dimensionInner < dimensions; dimensionInner++) {
                a[vertexIndex + (dimensionInner * numVertices)] -= distortionWeight * T[dimensionOuter][dimensionInner];
            }
           
            for (vector<int>::iterator neighborVertexIndexIterator = laplacianNeighborVertexIndices.begin(); neighborVertexIndexIterator != laplacianNeighborVertexIndices.end(); neighborVertexIndexIterator++) {
                int neighborVertexIndex = *neighborVertexIndexIterator;
                
                a[neighborVertexIndex + (dimensionOuter * numVertices)] -= distortionWeight / laplacianNumNeighbors;
            }

            int neighborIndex = 1;
            for (vector<int>::iterator neighborVertexIndexIterator = transformationNeighborVertexIndices.begin(); neighborVertexIndexIterator != transformationNeighborVertexIndices.end(); neighborVertexIndexIterator++) {
                int neighborVertexIndex = *neighborVertexIndexIterator;

                for (int dimensionInner = 0; dimensionInner < dimensions; dimensionInner++) {
                    a[neighborVertexIndex + (dimensionInner * numVertices)] -= distortionWeight * T[dimensionOuter][(neighborIndex * dimensions) + dimensionInner];
                }
                
                neighborIndex++;
            }
            
            for (map<int, double>::iterator atIterator = a.begin(); atIterator != a.end(); atIterator++) {
                //cout << "\n";
                for (map<int, double>::iterator aIterator = a.begin(); aIterator != a.end(); aIterator++) {
                    int atIndex = (*atIterator).first,
                        aIndex  = (*aIterator).first;
            
                    double atValue = (*atIterator).second,
                           aValue  = (*aIterator).second;
           
                    #if defined(_WIN32) || defined(_WIN64)
//                        AtA[atIndex][aIndex] += atValue * aValue;
                    #else
                        AtARowIndices.push_back(atIndex);
                        AtAColumnIndices.push_back(aIndex);
                        AtAValues.push_back(atValue * aValue);
                    #endif
                    //tripletList.push_back(Eigen::Triplet<double>(atIndex, aIndex, atValue * aValue));
                    //ata[atIndex*ataDim + aIndex] += atValue * aValue;
                    ata(atIndex, aIndex) += atValue * aValue;
                    //cout << (*atIterator).first << ", ";
                }
                //tripletList.push_back(Eigen::Triplet<double>(atIndex, aIndex, atValue * aValue));
            }
            //cout << "\n";
        }
    }

    //cout << "Done adding distortion terms!" << endl;
    //cout << "AtA fill rate: " << ata.size() << ", " << (numVertices * dimensions * numVertices * dimensions) << ", " << (float)ata.size() / (float)(numVertices * dimensions * numVertices * dimensions) << endl;

    // Solution
    /// Eigen sparse
    //Eigen::SparseMatrix<double, Eigen::ColMajor> A(numVertices * dimensions, numVertices * dimensions);
    //for (map<int, double>::iterator it = ata.begin(); it != ata.end(); ++it)
    //{
    //    tripletList.push_back(Eigen::Triplet<double>(it->first/ataDim, it->first % ataDim, it->second));
    //}    
    //A.setFromTriplets(tripletList.begin(), tripletList.end());

    //// Set up "x" for Ax=b
    //Eigen::VectorXd x(numVertices * dimensions);
    //// Setup solver
    //Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
    ////Eigen::BiCGSTAB <Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;

    //// Solve for Ax=b
    //A.makeCompressed();
    //solver.compute(A);
    //x = solver.solve(b);

    //if (solver.info() != Eigen::Success)
    //{
    //    // solving failed
    //    std::cout << "Eigen sparse solve AtAx=AtB fail...\n";
    //}

    //// Eigen dense
    Eigen::VectorXd x = ata.ldlt().solve(b); // ata.ldlt().solve(b) ata.llt().solve(b)


    #if defined(_WIN32) || defined(_WIN64)
//        Array2D<double> X = linearSolve(AtA, AtB);
    #else
        double* X = sparseLinearSolve(numVertices, dimensions, AtARowIndices, AtAColumnIndices, AtAValues, AtB);
    #endif

    //cout << "Done solving linear system!" << endl;

    // Convert to Vector
    deformedVertices.clear();
    deformedVertices.reserve(numVertices);
    for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
        #if defined(_WIN32) || defined(_WIN64)
//            deformedVertices.push_back(Vector3DFloat(
//                X[(vertexIndex + (0 * numVertices))][0], 
//                X[(vertexIndex + (1 * numVertices))][0], 
//                X[(vertexIndex + (2 * numVertices))][0]
//            ));
        #else
            deformedVertices.push_back(Vector3DFloat(
                X[(vertexIndex + (0 * numVertices))], 
                X[(vertexIndex + (1 * numVertices))], 
                X[(vertexIndex + (2 * numVertices))]
            ));
        #endif
        deformedVertices.push_back(Vector3DFloat(
            x[(vertexIndex + (0 * numVertices))], 
            x[(vertexIndex + (1 * numVertices))], 
            x[(vertexIndex + (2 * numVertices))]
        ));
    }

    // Compare the deformed vertices with vertices (by rmsd)
    m_shapeDistance = 0.0;
    for (int i = 0; i < numVertices; ++i)
    {
        m_shapeDistance += (deformedVertices[i][0] - vertices[i][0])*(deformedVertices[i][0] - vertices[i][0]) +
                         (deformedVertices[i][1]-vertices[i][1])*(deformedVertices[i][1]-vertices[i][1]) + 
                         (deformedVertices[i][2]-vertices[i][2])*(deformedVertices[i][2]-vertices[i][2]);
    }
    m_shapeDistance = sqrt(m_shapeDistance / (double)numVertices);

    //cout << "RMSD compared to last iteration: "<< shapeDistance << endl;
    
    // Cleanup
    #if !defined(_WIN32) && !defined(_WIN64)
        delete[] X;
    #endif
}

// Compute the weight matrix for pair wise correspondence
void FlexibleFittingEngine::computeCorrespondencePairWiseWeight(const int& selectedAlignmentIndex)
{
	//std::cout << "Feature sizes:\n" << sourceFeatures.size() << ", " << targetFeatures.size() << "\n";
    double sourceMinX =  MAX_DOUBLE, sourceMinY =  MAX_DOUBLE, sourceMinZ =  MAX_DOUBLE;
    double sourceMaxX = -MAX_DOUBLE, sourceMaxY = -MAX_DOUBLE, sourceMaxZ = -MAX_DOUBLE;
    double targetMinX =  MAX_DOUBLE, targetMinY =  MAX_DOUBLE, targetMinZ =  MAX_DOUBLE;
    double targetMaxX = -MAX_DOUBLE, targetMaxY = -MAX_DOUBLE, targetMaxZ = -MAX_DOUBLE;

    for (int clusterIndex = 0; clusterIndex < (int)correspondences[selectedAlignmentIndex].size(); ++clusterIndex)
        for (int featureIndex = 0; featureIndex < (int)correspondences[selectedAlignmentIndex][clusterIndex].size(); ++featureIndex)
		{
            // Cache the source feature, target feature and matching direction for each correspondence
            SSECorrespondenceNode correspondence = getCorrespondence(selectedAlignmentIndex, clusterIndex, featureIndex);
            int sourceIndex = correspondence.GetPIndex();
            int targetIndex = correspondence.GetQIndex();
            m_correspondenceSourceFeatureIndex.push_back(sourceIndex);
            m_correspondenceTargetFeatureIndex.push_back(targetIndex);
            m_correspondenceFlipFlag.push_back(correspondence.IsForward());
		}

    for (auto it = sourceFeatures.begin(); it != sourceFeatures.end(); ++it)
    {            
        // Upate the AABB bounding box
        sourceMinX = min(min(sourceMinX, it->GetEndPoint(0).X()), it->GetEndPoint(1).X());
        sourceMinY = min(min(sourceMinY, it->GetEndPoint(0).Y()), it->GetEndPoint(1).Y());
        sourceMinZ = min(min(sourceMinZ, it->GetEndPoint(0).Z()), it->GetEndPoint(1).Z());
        sourceMaxX = max(max(sourceMaxX, it->GetEndPoint(0).X()), it->GetEndPoint(1).X());
        sourceMaxY = max(max(sourceMaxY, it->GetEndPoint(0).Y()), it->GetEndPoint(1).Y());
        sourceMaxZ = max(max(sourceMaxZ, it->GetEndPoint(0).Z()), it->GetEndPoint(1).Z());

        //std::cout << "s\n";
        //std::cout << std::setprecision(30) << it->GetEndPoint(0).X() << ", " << it->GetEndPoint(0).Y() << ", " << it->GetEndPoint(0).Z() << "\n";
        //std::cout << std::setprecision(30) << it->GetEndPoint(1).X() << ", " << it->GetEndPoint(1).Y() << ", " << it->GetEndPoint(1).Z() << "\n";
    }

    for (auto it = targetFeatures.begin(); it != targetFeatures.end(); ++it)
    {
        targetMinX = min(min(targetMinX, it->GetEndPoint(0).X()), it->GetEndPoint(1).X());
        targetMinY = min(min(targetMinY, it->GetEndPoint(0).Y()), it->GetEndPoint(1).Y());
        targetMinZ = min(min(targetMinZ, it->GetEndPoint(0).Z()), it->GetEndPoint(1).Z());
        targetMaxX = max(max(targetMaxX, it->GetEndPoint(0).X()), it->GetEndPoint(1).X());
        targetMaxY = max(max(targetMaxY, it->GetEndPoint(0).Y()), it->GetEndPoint(1).Y());
        targetMaxZ = max(max(targetMaxZ, it->GetEndPoint(0).Z()), it->GetEndPoint(1).Z());

        //std::cout << "t\n";
        //std::cout << std::setprecision(30) << it->GetEndPoint(0).X() << ", " << it->GetEndPoint(0).Y() << ", " << it->GetEndPoint(0).Z() << "\n";
        //std::cout << std::setprecision(30) << it->GetEndPoint(1).X() << ", " << it->GetEndPoint(1).Y() << ", " << it->GetEndPoint(1).Z() << "\n";
    }
    
    //double ratioBB = 0.3;
    double sourceBBDiagonal = (Vector3DDouble(sourceMaxX, sourceMaxY, sourceMaxZ) - Vector3DDouble(sourceMinX, sourceMinY, sourceMinZ)).Length();
    double targetBBDiagonal = (Vector3DDouble(targetMaxX, targetMaxY, targetMaxZ) - Vector3DDouble(targetMinX, targetMinY, targetMinZ)).Length();

    // Gaussian weight --- 3*sigma == ratio*boundingBoxDiagonal
    double sigma = m_bbRatio * min(sourceBBDiagonal, targetBBDiagonal) / 3.0;

	// For all the paired features, compute their weight based on midpoint distance
    // Since we will fill in all the upper triangle, it is safe to use "std::vector<>().resize" here
    m_correspondencePairWiseWeight.resize(m_correspondenceSourceFeatureIndex.size(), std::vector<double>(m_correspondenceSourceFeatureIndex.size()));
    for (int i = 0; i < (int)m_correspondencePairWiseWeight.size(); ++i)
        for (int j = i; j < (int)m_correspondencePairWiseWeight.size(); ++j)
        {
            if (i == j)
            {
                m_correspondencePairWiseWeight[i][j] = 1.0;
                continue;
            }

            Vector3DDouble sourceCentroidI = Vector3DDouble(sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetCentroid().X(), sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetCentroid().Y(), sourceFeatures[m_correspondenceSourceFeatureIndex[i]].GetCentroid().Z());
            Vector3DDouble sourceCentroidJ = Vector3DDouble(sourceFeatures[m_correspondenceSourceFeatureIndex[j]].GetCentroid().X(), sourceFeatures[m_correspondenceSourceFeatureIndex[j]].GetCentroid().Y(), sourceFeatures[m_correspondenceSourceFeatureIndex[j]].GetCentroid().Z());

            Vector3DDouble targetCentroidI = Vector3DDouble(targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetCentroid().X(), targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetCentroid().Y(), targetFeatures[m_correspondenceTargetFeatureIndex[i]].GetCentroid().Z());
            Vector3DDouble targetCentroidJ = Vector3DDouble(targetFeatures[m_correspondenceTargetFeatureIndex[j]].GetCentroid().X(), targetFeatures[m_correspondenceTargetFeatureIndex[j]].GetCentroid().Y(), targetFeatures[m_correspondenceTargetFeatureIndex[j]].GetCentroid().Z());

            // sourceCenterDistance and targetCenterDistance
            double scDist = (sourceCentroidI - sourceCentroidJ).Length();
            double tcDist = (targetCentroidI - targetCentroidJ).Length();

            m_correspondencePairWiseWeight[i][j] = std::exp(-0.5 * (scDist*scDist + tcDist*tcDist) / (sigma*sigma));
        }

        //std::cout << "sigma: " << sigma << "\n";
        //std::cout << "sourceBBDiagonal: " << sourceBBDiagonal << "\n";
        //std::cout << "targetBBDiagonal: " << targetBBDiagonal << "\n";
        //std::cout << "m_correspondencePairWiseWeight:\n";
        //for (const auto &i : m_correspondencePairWiseWeight)
        //{
        //    for (const auto &j : i)
        //    {
        //        std::cout << j << ",";
        //    }
        //    std::cout << "\n";
        //}        
}

void FlexibleFittingEngine::calculateCorrespondences() {
    sseCorrespondenceFinder.InitializeFeatures(sourceFeatures, targetFeatures); 
    sseCorrespondenceFinder.InitializeConstants(0, helixLengthDifference, 0, 0, 0, 0, 0, jointAngleDifference, dihedralAngleDifference, helixCentroidDifference, maxCorrespondenceSolutions, 10);
    
    correspondences = sseCorrespondenceFinder.GetAStarTriangleBasedCliqueDistanceFeatureCorrespondence(false, true, 4);
}

void FlexibleFittingEngine::addCorrespondencesByReadIn(const int& pIndex, const int& qIndex, const bool& isForward) {

    // When reading in the correspondence, we assume there is only one alignment and one correspondence in one cluster
    if (correspondences.empty())
        correspondences.push_back(vector< vector<SSECorrespondenceNode> >(1, vector<SSECorrespondenceNode>(1, SSECorrespondenceNode(pIndex, qIndex, isForward))));
    else
        correspondences[0].push_back(vector<SSECorrespondenceNode>(1, SSECorrespondenceNode(pIndex, qIndex, isForward)));
}

void FlexibleFittingEngine::calculateNeighborVertexIndicesCaches(const int& numVertices) {
    // Initialize Immediate Neighbor Cache
    vector< vector<int> > neighborVertexIndicesCache(numVertices, vector<int>());

    // Calculate Immediate Neighbor Cache
    for (vector< pair<int, int> >::iterator edgeIterator = edges.begin(); edgeIterator != edges.end(); edgeIterator++) {
        pair<int, int> edge = *edgeIterator;
    
        neighborVertexIndicesCache[edge.first].push_back(edge.second);
        neighborVertexIndicesCache[edge.second].push_back(edge.first);
    }

    // Calculate Extended Neighbor Caches
    laplacianNeighborVertexIndicesCache      = vector< vector<int> >(numVertices, vector<int>());
    transformationNeighborVertexIndicesCache = vector< vector<int> >(numVertices, vector<int>());
    deformationNeighborVertexIndicesCache    = vector< vector<int> >(numVertices, vector<int>());
    for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
        queue< pair<int, int> > frontier;
        set<int> explored;

        frontier.push(make_pair(vertexIndex, 0));
        explored.insert(vertexIndex);
        while (!frontier.empty()) {
            pair<int, int> exploration = frontier.front();

            int explorationVertexIndex = exploration.first,
                explorationDistance    = exploration.second;

            if (explorationVertexIndex != vertexIndex) {
                if (explorationDistance <= laplacianNeighborhoodSize) {
                    laplacianNeighborVertexIndicesCache[vertexIndex].push_back(explorationVertexIndex);
                }

                if (explorationDistance <= transformationNeighborhoodSize) {
                    transformationNeighborVertexIndicesCache[vertexIndex].push_back(explorationVertexIndex);
                }

                if (explorationDistance <= deformationNeighborhoodSize) {
                    deformationNeighborVertexIndicesCache[vertexIndex].push_back(explorationVertexIndex);
                }

                if (explorationDistance > laplacianNeighborhoodSize && explorationDistance > transformationNeighborhoodSize && explorationDistance > deformationNeighborhoodSize) {
                    break;
                }
            }
            frontier.pop();

            vector<int> neighborVertexIndices = neighborVertexIndicesCache[explorationVertexIndex];
            for (vector<int>::iterator neighborVertexIndexIterator = neighborVertexIndices.begin(); neighborVertexIndexIterator != neighborVertexIndices.end(); neighborVertexIndexIterator++) {
                int neighborVertexIndex = *neighborVertexIndexIterator;

                if (explored.count(neighborVertexIndex) == 0) {
                    frontier.push(make_pair(neighborVertexIndex, explorationDistance + 1));
                    explored.insert(neighborVertexIndex);
                }
            }
        }
        
        sort(laplacianNeighborVertexIndicesCache[vertexIndex].begin(), laplacianNeighborVertexIndicesCache[vertexIndex].end());
        sort(transformationNeighborVertexIndicesCache[vertexIndex].begin(), transformationNeighborVertexIndicesCache[vertexIndex].end());
        sort(deformationNeighborVertexIndicesCache[vertexIndex].begin(), deformationNeighborVertexIndicesCache[vertexIndex].end());
    }

    //cout << "Done calculating neighbor vertex indices caches!" << endl;
}

void FlexibleFittingEngine::calculateDeformationTransformations(const int& numVertices) {
    deformationTransformations.reserve(numVertices);
    for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
        vector<Vector3DFloat> sourceNeighborhood,
                              targetNeighborhood;

        sourceNeighborhood.push_back(originalVertices[vertexIndex]); 
        targetNeighborhood.push_back(deformedVertices[vertexIndex]); 

        vector<int> neighborVertexIndices = deformationNeighborVertexIndicesCache[vertexIndex];
        for (vector<int>::iterator neighborVertexIndexIterator = neighborVertexIndices.begin(); neighborVertexIndexIterator != neighborVertexIndices.end(); neighborVertexIndexIterator++) {
            sourceNeighborhood.push_back(originalVertices[*neighborVertexIndexIterator]); 
            targetNeighborhood.push_back(deformedVertices[*neighborVertexIndexIterator]); 
        }

        deformationTransformations.push_back(linearSolver.FindRotationTranslation(sourceNeighborhood, targetNeighborhood));
    }

    //cout << "Done calculating deformation transformations!" << endl;
}

Array2D<double> FlexibleFittingEngine::transformationForm(const Vector3DFloat& vertex, const int& dimensions, const bool& oneFlag) {
    vector<Vector3DFloat> vertices;
    vertices.push_back(vertex);

    return transformationForm(vertices, dimensions, oneFlag);
}

Array2D<double> FlexibleFittingEngine::transformationForm(const vector<Vector3DFloat>& workVertices, const int& dimensions, const bool& oneFlag) {
    int numVertices = workVertices.size();
    
    Array2D<double> form(numVertices * dimensions, (2 * dimensions) + 1, 0.0);
    for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
        // X  0  Z -Y  1  0  0
        // Y -Z  0  X  0  1  0
        // Z  Y -X  0  0  0  1

        Vector3DFloat vertex = workVertices[vertexIndex];

        form[vertexIndex * dimensions + 0][0] =  vertex.X();
        form[vertexIndex * dimensions + 0][2] =  vertex.Z();
        form[vertexIndex * dimensions + 0][3] = -vertex.Y();
       
        form[vertexIndex * dimensions + 1][0] =  vertex.Y();
        form[vertexIndex * dimensions + 1][1] = -vertex.Z();
        form[vertexIndex * dimensions + 1][3] =  vertex.X();
        
        form[vertexIndex * dimensions + 2][0] =  vertex.Z();
        form[vertexIndex * dimensions + 2][1] =  vertex.Y();
        form[vertexIndex * dimensions + 2][2] = -vertex.X();

        if (oneFlag) {
            form[vertexIndex * dimensions + 0][4] = 1.0;
            form[vertexIndex * dimensions + 1][5] = 1.0;
            form[vertexIndex * dimensions + 2][6] = 1.0;
        }
    }

    return form;
}

Array2D<double> FlexibleFittingEngine::inverse(const Array2D<double>& array) {
    if (array.dim1() != array.dim2()) {
        cerr << "Unable to invert matrix, not square!" << endl;
        exit(-1);
    }

    Array2D<double> identity(array.dim1(), array.dim1(), 0.0);
    for (int i = 0; i < array.dim1(); i++) {
        identity[i][i] = 1.0;
    }

    return linearSolve(array, identity);
}

Array2D<double> FlexibleFittingEngine::transpose(const Array2D<double>& array) {
    Array2D<double> arrayTranspose(array.dim2(), array.dim1());

    for (int row = 0; row < array.dim1(); row++) {
        for (int col = 0; col < array.dim2(); col++) {
            arrayTranspose[col][row] = array[row][col];
        }
    }

    return arrayTranspose;
}

Array2D<double> FlexibleFittingEngine::linearSolve(const Array2D<double>& A, const Array2D<double>& B) {
    Array2D<double> solution;
    
    LU<double> lu(A);
    if (lu.isNonsingular()) {
        solution = lu.solve(B);
    }
    else {
        cerr << "Unable to solve linear system of equations using LU!" << endl;
        exit(-1);
    }

    //QR<double> qr(A);
    //if (qr.isFullRank()) {
    //    solution = qr.solve(B);
    //}
    //else {
    //    cerr << "Unable to solve linear system of equations using QR!" << endl;
    //    exit(-1);
    //}
            
    //Cholesky<double> ch(A);
    //if (ch.is_spd()) {
    //    solution = ch.solve(B);
    //}
    //else {
    //    cerr << "Unable to solve linear system of equations using Cholesky!" << endl;
    //    exit(-1);
    //}

    return solution;
}

#if !defined(_WIN32) && !defined(_WIN64)
    double* FlexibleFittingEngine::sparseLinearSolve(const int& numVertices, const int& dimensions, const vector<int>& ATripletRowIndices, const vector<int>& ATripletColumnIndices, const vector<double>& ATripletValues, const vector<double>& B) {
        // Convert to UMFPACK Format
        int *AColumnCount      = new int [(numVertices * dimensions) + 1],
            *AColumnRowIndices = new int [ATripletValues.size()];

        double *AColumnValues = new double [ATripletValues.size()];

        umfpack_di_triplet_to_col(numVertices * dimensions, numVertices * dimensions, ATripletValues.size(), &ATripletRowIndices[0], &ATripletColumnIndices[0], &ATripletValues[0], AColumnCount, AColumnRowIndices, AColumnValues, NULL);

        // Solve
        void *symbolic,
             *numeric;
        
        double *X = new double [numVertices * dimensions];

        umfpack_di_symbolic(numVertices * dimensions, numVertices * dimensions, AColumnCount, AColumnRowIndices, AColumnValues, &symbolic, NULL, NULL);
        umfpack_di_numeric(AColumnCount, AColumnRowIndices, AColumnValues, symbolic, &numeric, NULL, NULL);
        umfpack_di_solve(UMFPACK_A, AColumnCount, AColumnRowIndices, AColumnValues, X, &B[0], numeric, NULL, NULL);

        // Cleanup
        delete [] AColumnCount;
        delete [] AColumnRowIndices;
        delete [] AColumnValues;
        umfpack_di_free_symbolic(&symbolic);
        umfpack_di_free_numeric(&numeric);

        return X;
    }
#endif

#endif
