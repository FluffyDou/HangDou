// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Rendering engine responsible for rendering C-Alpha atoms.

// CVS Meta Information: 
//   $Source: /project/mm/cvs/graphics/ssa1/source/Gorgon/CAlphaRenderer.h,v $
//   $Revision: 1.69 $
//   $Date: 2011/10/12 21:08:47 $
//   $Author: coleman.r $
//   $State: Exp $
//
// History Log: 
//   $Log: CAlphaRenderer.h,v $
//   Revision 1.69  2011/10/12 21:08:47  coleman.r
//   commented out SSEHunter C++ code that is not in use. Some of it has been replaced with Python code in Gorgon/src_py/sse_hunter_engine.py.
//
//   Revision 1.68  2011/08/20 20:03:21  coleman.r
//   Reverting the algorithm for GetSelectedAtom, fixing a logical error I made in my last commit.
//
//   Revision 1.67  2011/06/07 16:03:23  coleman.r
//   We had been using memory addresses for "names" in glLoadName() in GL_SELECT mode. Now, we are storing atom hash keys in a vector (hash keys are 64 bit) and using the indices of these hash keys as "names". (The indices should be 32 bit and can be cast to GLuint safely.) This avoids a bug on 64 bit systems where memory addresses were too large to fit in a GLuint type. Thus, a 64 bit OS segmentation fault (we previously cast int "names" back to pointers with SelectionToggle()) is avoided.
//
//   Revision 1.66  2010/10/11 23:24:36  coleman.r
//   Fixing last commit (I had made similar changes to the code on two different computers)
//
//   Revision 1.65  2010/10/11 23:18:49  coleman.r
//   added GetAtomHashes()
//
//   Revision 1.64  2010/08/20 13:52:43  coleman.r
//   gcc compile fix: gcc requires nested templates to end in "> >" not ">>"
//
//   Revision 1.63  2010/08/19 23:05:07  chenb
//   Cleaned and commented ribbon diagram code
//
//   Revision 1.62  2010/08/13 21:20:16  coleman.r
//   AutoHelixBuilder changes
//
//   Revision 1.61  2010/07/29 20:21:58  coleman.r
//   gcc compile fix: gcc requires nested templates to end in "> >" not ">>"
//
//   Revision 1.60  2010/07/27 23:18:58  chenb
//   Ribbon diagram code now merged with flexible fitting code
//
//   Revision 1.59  2010/07/23 18:18:32  heiderp
//   Side chains now transform correctly.  PDB helices now color correctly and rigid initialization bug is fixed
//
//   Revision 1.58  2010/07/22 21:09:07  heiderp
//   Minor updates. Mostly commenting and removing extra material from CurveDeformer.h
//
//   Revision 1.57  2010/07/19 17:29:02  heiderp
//   LARGE update.  Added flexible fitting functionality, lots of logic in FlexibleFittingEngine.h
//
//   Revision 1.51  2010/06/23 19:11:51  ssa1
//   Adding simple ribbon rendering and associated events for flexible fitting
//
//   Revision 1.50  2010/06/23 13:02:56  ssa1
//   Allowing users to reset a flexible fitting if need be.
//
//   Revision 1.49  2010/05/27 18:28:46  ssa1
//   Better color control for all atom visualization
//
//   Revision 1.48  2010/05/27 17:10:18  ssa1
//   Better color control for all atom visualization
//
//   Revision 1.47  2010/05/27 05:08:49  ssa1
//   Side chain visualization on Gorgon
//
//   Revision 1.46  2010/05/27 04:41:54  ssa1
//   Side chain visualization on Gorgon
//
//   Revision 1.45  2010/05/26 20:17:35  ssa1
//   Adding in display styles for atom rendering.
//
//   Revision 1.44  2010/05/21 15:45:16  ssa1
//   Flexible fitting implemented in Gorgon
//
//   Revision 1.43  2010/05/20 21:55:53  ssa1
//   Rigid body alignment based on largest flexible cluster
//
//   Revision 1.42  2010/02/11 23:19:11  ssa1
//   Allowing the ability to save pseudoatoms generated from SSEHunter
//
//   Revision 1.41  2010/01/10 05:31:43  colemanr
//   PDBAtoms now store their correlation, skeleton, and geometry scores. Changing the weighting for these three scores in the GUI now changes the total score for each pseudoatom.
//
//   Revision 1.40  2009/10/13 18:09:34  ssa1
//   Refactoring Volume.h
//
//   Revision 1.39  2009/08/10 20:03:40  ssa1
//   SSEHunter interfaced into Gorgon
//
//   Revision 1.38  2009/08/10 13:54:38  ssa1
//   Adding initial ssehunter program
//
//   Revision 1.37  2009/07/01 21:25:13  ssa1
//   Centering the volume cropped using a radius around the point selected by the atom selection tool.
//
//   Revision 1.36  2009/06/30 21:23:24  ssa1
//   SSEHunter results have range between -3 and 3, not -1 and 1
//
//   Revision 1.35  2009/06/22 20:17:27  ssa1
//   Adding in SSEBuilder Functionality: Selection to Helix functionality
//
//   Revision 1.34  2009/06/19 18:51:05  ssa1
//   Adding in SSEBuilder Functionality
//
//   Revision 1.33  2009/03/30 21:36:12  ssa1
//   Interactive loop building
//
//   Revision 1.32  2008/12/07 07:11:36  ssa1
//   Coloring bonds with red and blue if they exceed maximum or minimum length restrictions
//
//   Revision 1.31  2008/12/03 21:58:25  ssa1
//   Selection rotations for atoms and helices.
//
//   Revision 1.30  2008/12/02 23:55:43  colemanr
//   Fixed logic for GetBondIndex().
//
//   Revision 1.29  2008/12/02 21:25:44  ssa1
//   adding getBondIndex method to give access to bonds
//
//   Revision 1.28  2008/11/13 20:54:40  ssa1
//   Using the correct scale when loading volumes
//
//   Revision 1.27  2008/11/10 16:15:43  ssa1
//   Making python and C++ use the same PDBAtom objects
//
//   Revision 1.26  2008/11/07 21:32:21  ssa1
//   Fixing returning of the actual c++ pdbatom object instead of a copy
//
//   Revision 1.25  2008/10/10 14:25:55  ssa1
//   Setting the cost functions to scale with the edge length
//
//   Revision 1.24  2008/10/07 23:48:14  colemanr
//   added a function which returns the PDBAtom for a given hitStack
//
//   Revision 1.23  2008/09/29 20:36:35  ssa1
//   Drawing skeletal curves as cylinders and spheres
//
//   Revision 1.22  2008/09/29 16:01:17  ssa1
//   Adding in CVS meta information
//

#ifndef GORGON_CALPHA_RENDERER_H
#define GORGON_CALPHA_RENDERER_H


#include <glut.h>
#include <cstdlib>
#include <cstdio>
#include <ProteinMorph/NonManifoldMesh.h>
#include <ProteinMorph/SSEHunter.h>
#include <SkeletonMaker/volume.h>
#include <GraphMatch/PDBReader.h>
#include <GraphMatch/PDBAtom.h>
#include <GraphMatch/PDBBond.h>
#include "Renderer.h"
#include <map>
#include <list>
#include <GraphMatch/VectorMath.h>
#include <boost/tuple/tuple.hpp>
#include <set>
#include <iostream>
#include <fstream>
#include <boost/python.hpp>
#include <math.h>
#include <tuple>
#include <Eigen/Dense>

using namespace std;
using namespace wustl_mm::Protein_Morph;
using namespace wustl_mm::GraphMatch;
using namespace wustl_mm::SkeletonMaker;
using namespace boost::tuples;

namespace wustl_mm {
	namespace Visualization {	
        typedef map<unsigned long long, PDBAtom> AtomMapType;
        struct SerialAndHashType {
            unsigned int serial;
            unsigned long long hashKey;
        };

        class SerialAndHashTypePredicate {
        public:
            bool operator() (const SerialAndHashType& lhs, const SerialAndHashType& rhs) {
                return lhs.serial < rhs.serial;
            }
        };

        const int CALPHA_DISPLAY_STYLE_BACKBONE = 3;
        const int CALPHA_DISPLAY_STYLE_RIBBON = 4;
        const int CALPHA_DISPLAY_STYLE_SIDE_CHAIN = 5;
        const int CALPHA_DISPLAY_STYLE_BACKBONE_PATHWALKER = 6;

        /**
        Begin Hermite Curve code, to be moved into another file after testing
        -this code based on molscript's hermite_curve.c file, and produced with the help
        of wikipedia's article on the cubic hermite spline
        */
        class HermiteCurve{
        public:
            Vector3DFloat p0, p1, m0, m1;

            void setCurve(Vector3DFloat pstart, Vector3DFloat pend, Vector3DFloat tstart, Vector3DFloat tend);
            Vector3DFloat getPos(double t);
            Vector3DFloat getTangent(double t);
        };

        void HermiteCurve::setCurve(Vector3DFloat pstart, Vector3DFloat pend, Vector3DFloat tstart, Vector3DFloat tend){
            p0 = pstart;
            p1 = pend;
            m0 = tstart;
            m1 = tend;
        }

        Vector3DFloat HermiteCurve::getPos(double t){
            double tsquared = t*t;
            double tcubed = tsquared * t;

            double cp0 = 2 * tcubed - 3 * tsquared + 1;
            double cm0 = tcubed - 2 * tsquared + t;
            double cp1 = (cp0 - 1)*(-1);
            double cm1 = tcubed - tsquared;

            double xt = cp0*p0.X() + cm0*m0.X() + cp1*p1.X() + cm1*m1.X();
            double yt = cp0*p0.Y() + cm0*m0.Y() + cp1*p1.Y() + cm1*m1.Y();
            double zt = cp0*p0.Z() + cm0*m0.Z() + cp1*p1.Z() + cm1*m1.Z();

            return Vector3DFloat(xt, yt, zt);
        }

        // I don't know how this method works, but it is a part of the entirely functional
        // molscript code - BC
        Vector3DFloat HermiteCurve::getTangent(double t){
            double t2 = t * t;
            double cp0 = 6.0 * (t2 - t);
            double cp1 = 6.0 * (-t2 + t);
            double cm0 = 3.0 * t2 - 4.0 * t + 1.0;
            double cm1 = 3.0 * t2 - 2.0 * t;
            double vxt = p0.X()*cp0 + p1.X() * cp1 + m0.X() * cm0 + m1.X() * cm1;
            double vyt = p0.Y()*cp0 + p1.Y() * cp1 + m0.Y() * cm0 + m1.Y() * cm1;
            double vzt = p0.Z()*cp0 + p1.Z() * cp1 + m0.Z() * cm0 + m1.Z() * cm1;

            return Vector3DFloat(vxt, vyt, vzt);
        }
        /**
        End Hermite Curve code
        */

        struct Quad {
            float p1[3];
            float p2[3];
            float p3[3];
            float p4[3];
            float normal[3];
            float saliency[3];
            float localIntensity;
            float eigenValue;
            float eigenValue0;
            float eigenValue1;
            int type;
            int index1;
            int index2;
            int index3;
            int index4;
            unsigned long long hash1;
            unsigned long long hash2;
            unsigned long long hash3;
            unsigned long long hash4;
            float eigenVector0[3];
            float eigenVector1[3];
            float eigenVector2[3];
        };

        struct fittedHelix {
            std::vector<int> fittedHelixBondIndices;
            std::vector<Vector3DFloat> fittedHelixBondPositions;
            float coilHeight;
            float coils;
            float orthoDist;
            float helixError;
        };

        class CAlphaRenderer : public Renderer {
        public:
            struct Secel{
                std::vector<unsigned long long> atomHashes;
                bool selected;
            };

            CAlphaRenderer();
            ~CAlphaRenderer();

            void Draw(int subSceneIndex, bool selectEnabled);
            void LoadFile(string fileName);
            void LoadSSEHunterFile(string fileName);
            bool SaveSSEHunterFile(string fileName);
            //			void GetSSEHunterAtoms(Volume * vol, NonManifoldMesh_Annotated * skeleton, float resolution, float threshold, float correlationCoeff, float skeletonCoeff, float geometryCoeff);
            void UpdateTotalScoreSSEHunterAtoms(float correlationCoeff, float skeletonCoeff, float geometryCoeff);
            void ColorSSEHunterAtoms();
            int SelectionObjectCount();
            int SelectionAtomCount();
            Vector3DFloat SelectionCenterOfMass();
            bool SelectionRotate(Vector3DFloat centerOfMass, Vector3DFloat rotationAxis, float angle);
            bool SelectionMove(Vector3DFloat moveDirection);
            bool SelectionClear();
            void SelectionToggle(int subsceneIndex, bool forceTrue, int ix0, int ix1 = -1, int ix2 = -1, int ix3 = -1, int ix4 = -1);
            void Unload();
            string GetSupportedLoadFileFormats();
            string GetSupportedSaveFileFormats();
            Vector3DFloat Get3DCoordinates(int subsceneIndex, int ix0, int ix1 = -1, int ix2 = -1, int ix3 = -1, int ix4 = -1);
            void TransformAllAtomLocations(MatrixFloat transform);
            void UpdateBoundingBox();
            //void DrawSurface(int quad1, int quad2, int quad3, int quad4, float nx, float ny, float nz);
            void DrawSurface(Quad quad);
            void setMaxCurveEigenMaxesAndMins(float maxCurveEigen0Max, float maxCurveEigen1Max, float maxCurveEigen2Max, float maxCurveEigen0Min, float maxCurveEigen1Min, float maxCurveEigen2Min);
            Matrix4 alignWithAxis(float v1[3], float axis[3], Vector3DFloat minPt, Vector3DFloat maxPt);

            // Controlling the atom std::vector
            PDBAtom * AddAtom(PDBAtom atom);
            PDBAtom * GetAtom(unsigned long long index);
            PDBAtom * GetAtomFromHitStack(int subsceneIndex, bool forceTrue, int ix0, int ix1, int ix2, int ix3, int ix4);
            PDBAtom * GetSelectedAtom(unsigned int selectionId);
            void DeleteAtom(unsigned long long index);
            int GetAtomCount();
            std::vector<unsigned long long> GetAtomHashes();

            //Controlling the bond std::vector
            void AddBond(PDBBond bond);
            PDBBond * GetBond(int index);
            int GetBondIndex(unsigned long long atom0, unsigned long long atom1);
            void DeleteBond(int index);
            int GetBondCount();

            //Controlling the bond std::vector
            void AddSideChainBond(PDBBond bond);
            PDBBond * GetSideChainBond(int index);
            int GetSideChainBondIndex(unsigned long long atom0, unsigned long long atom1);
            void DeleteSideChainBond(int index);
            int GetSideChainBondCount();

            void SetNumSegments(int segments);
            void SetNumSlices(int slices);
            void DrawDeletedBond(int atom1, int atom2);
            void UndrawDeletedBond(int atom1, int atom2);
            void addMaxEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2);
            void addMinEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2);
            void addSaddleEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2);

            void DrawAddedBond(int atom1, int atom2);
            void UndrawAddedBond(int atom1, int atom2);

            void AddMaxBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1);
            void AddMinBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1);
            void AddSaddleBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1);
            void addOrd(float v1, float v2, float v3);
            void setEllipsoidDisplay(bool displayEllipsoid);
            void setOrthoDistMin(float min);
            void setOrthoDistMax(float max);



            string FindDistance(int atom1, int atom2);

            int StartHelix(); //StartHelix creates a new helix element in aHelices and returns its index
            void AddHelixElement(int, unsigned long long); //adds a helix element to the helix indexed at param 1

            int StartStrand(); //StartStrand creates a new strand element in bStrands and returns its index
            void AddStrandElement(int, unsigned long long); //adds a strand element to the strand indexed at param 1

            int StartLoop(); //StartLoop creates a new strand element in loops and returns its index
            void AddLoopElement(int, unsigned long long); //adds a loop element to the loop indexed at param 1

            //void DrawHelices();
            //void DrawStrands();
            //void DrawLoops();

            bool CleanSecondaryStructures(); //empties the aHelices, bStrands and loops variables
            //what should really happen is that the code should check if it is
            //trying to reload the same one, and then if it did return false
            std::vector<Vector3DFloat> CreatePointVector(PDBAtom first, PDBAtom last); // functionality mirrored in previously implemented method,
            // will try to refactor
            std::vector<Vector3DFloat> LaplacianSmoothing(std::vector<Vector3DFloat> points, int steps); // applies Laplacian smoothing to a std::vector of
            // Vector3DFloats
            std::vector<Vector3DFloat> CreateStrandNormals(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next); // create line segment normals to be used in drawing Beta
            // strands
            std::vector<PDBBond> getDeletedBonds();
            std::vector< std::vector<Vector3DFloat> > originalHelices;
            std::vector<float> strandCycles;
            std::vector<float> orthoDistances;

            void CreateHelixAxesTangentsAndPoints(std::vector<Vector3DFloat>& axes, std::vector<Vector3DFloat>& tangents, std::vector<Vector3DFloat>& interpPoints, std::vector<Vector3DFloat> points,
                Vector3DFloat previous, Vector3DFloat next, double HELIX_ALPHA, double HELIX_BETA, double HELIX_HERMITE_FACTOR);
            void DrawOpenBox(std::vector<Vector3DFloat> points, std::vector<Vector3DFloat> normals); // takes a std::vector of 8 points and draws a rectangular prism with two of its six sides not
            // filled in; the first 4 points are from the beggining edge of the box, with the second four
            // forming the end
            void DrawTube(std::vector<Vector3DFloat> points, std::vector<Vector3DFloat> normals, int stacks, int slices);
            std::vector<Vector3DFloat> InterpolateLoopPoints(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS); // creates interpolated points for loops
            //std::vector<Vector3DFloat> InterpolateStrandPoints(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS);
            //std::vector<Vector3DFloat> InterpolateHelixPoints(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS);

            // for testing purposes only; allow changing of highlight color values
            void SetHltRValue(int col);
            void SetHltGValue(int col);
            void SetHltBValue(int col);
            void SetHltAValue(int col);


            std::vector<int> GetSelectedHelixIndices();
            void SetHelixCorrs(std::vector < int > flatCorrespondences);
            void SetSelectedSSEHelices(std::vector<int>);
            void ClearOtherHighlights();
            void SetFeatureVecs(std::vector<Vector3DFloat> flatFeatureVecs);
            void SetHelixColor(int helixNum, float r, float g, float b);
            void ClearHelixColors();

            string getDeletedBondAtoms();
            std::vector<unsigned long long> getDeletedBonds1Ix();
            void RemoveSelectedBonds(string nobonds);
            void addSelectedBonds(string newBonds);
            void DrawBackboneModelPathwalker(int subSceneIndex, bool selectEnabled);
            void DeleteAtomFromVisualization(unsigned long long deletedAtom);
            void setMaxOn(bool display);
            void setMinOn(bool display);
            void setSaddleOn(bool display);
            void SetExtremalMode(bool extremal);
            void sortMaxCurveDimensions();
            bool maxOn = false, minOn = false, saddleOn = false;

            void toggleSaddleOn(int atom1, int atom2, bool display);
            void toggleMaxOn(int atom1, int atom2, bool display);
            void toggleMinOn(int atom1, int atom2, bool display);
            void addQuad(int q1, int q2, int q3, int q4, float nx, float ny, float nz);
            //void AddQuadSurface(string quad1, string quad2, string quad3, string quad4, float nx, float ny, float nz, float saliency, float intensity, float eigenvalue, int type);
            //void AddQuadSurface(string quad1, string quad2, string quad3, string quad4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type);

            //void AddQuadSurface(int index1, int index2, int index3, int index4, string quad1, string quad2, string quad3, string quad4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type);
            void AddQuadSurface(int index1, int index2, int index3, int index4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type, float eigenValue0, float eigenValue1);
            void addQuadEigenvectors(int index1, int index2, int index3, int index4, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8);
            void interpolateHelixPoints();

            void clearQuads();
            py::list getQuadPts();
            py::list getQuadNormals();
            py::list findMaxCurveStrand(int i);

            int findHelices();
            void makeMax(int index, bool maxB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1);
            void makeMin(int index, bool minB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1);
            void makeSaddle(int index, bool minB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1);

            void addMaxPtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8);
            void addMinPtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8);
            void addSaddlePtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8);


            void setExtremalParams(float curve, float point, float surface, float minG, float maxG, float eigen);
            void setExtremalChecks(bool saliency, bool intensity, bool eigenvalue);
            void setExtremalHide(bool isHide);
            void setExtremalSurfaceHide(bool isHide);
            void setMaxSurfaceOn(bool maxOn);
            void setMinSurfaceOn(bool minOn);
            void setMinMaxIntensitiesEigenvalues(float minI, float maxI, float minE, float maxE);
            void setMinDisplay(bool minPt);
            void setMaxDisplay(bool maxPt);
            void setSaddleDisplay(bool saddlePt);
            void addVertexPos(float px, float py, float pz);
            void writeMaxCurveToPDB();
            bool CheckMinPointHide(float saliency[3], float intensity, float eigenvalue);
            bool CheckMaxPointHide(float saliency[3], float intensity, float eigenvalue);
            bool CheckSaddlePointHide(float saliency[3], float intensity, float eigenvalue);
            void findMaxCurves();
            int findMaxHelices();
            //void addCurveCoeff(float coefA, float coefB, float coefC, float coefD);
            void addCurveCoeff(float coefA, float coefB, float coefC);
            void findCurveMinsMaxes();
            void inRangeOfSkeleton();
            void scoreHelix(int strandIndex);
            //void scoreHelix(std::vector<Vector3DFloat> currentStrand, );
            void normalizedHelicesMinsAndMaxes();
            void findEigenMinMaxes();
            void setHelixDisplay(bool display);
            void setHelixDisplayPoints(bool display);

            void setHelixDisplayThreshold(float helixThresh);
            void setScaleEllipsoid(float eScale);
            void setHelixCoilHeight(float coilHeight);
            void setSegmentThreshold(int threshold);
            void averageEllipsoidScale();

            //float helixThreshold = 30.0;

            bool ellipsoidVisible;
            float ellipsoidScale;
            int helixCoilHeight;
            int helixSegmentThreshold;
            float orthoDistMin;
            float orthoDistMax;
            int ellipsoidScaleFactor;
            std::vector<Vector3DFloat> ordLines;
            std::vector<float> helixScoreThresholds;
            std::vector<int> quadPts;
            std::vector<float> quadNormals;
            std::vector<Quad> quads;
            std::vector<float> vertexPositions;
            std::vector<Vector3DFloat> maxCoord;
            std::vector<Vector3DFloat> minCoord;
            std::vector<Vector3DFloat> saddleCoord;

            std::vector<float> helixTurnHeights;

            std::vector<float> maxPointSaliencies;
            std::vector< fittedHelix > fittedHelices;

            std::vector<float> maxPointIntensities;
            std::vector<float> maxPointEigenvalues;
            std::vector<float> maxPointEigenvalues0;
            std::vector<float> maxPointEigenvalues1;

            std::vector<float> minPointSaliencies;

            std::vector<float> minPointIntensities;
            std::vector<float> minPointEigenvalues;
            std::vector<float> minPointEigenvalues0;
            std::vector<float> minPointEigenvalues1;
            std::vector< std::vector<Vector3DFloat> > normalizedHelices;

            std::vector<float> saddlePointSaliencies;

            std::vector<float> saddlePointIntensities;
            std::vector<float> saddlePointEigenvalues;
            std::vector<float> saddlePointEigenvalues0;
            std::vector<float> saddlePointEigenvalues1;

            std::vector< std::vector<unsigned long long> > helices;
            std::vector<unsigned long long> bondAtm1s;
            std::vector<unsigned long long> bondAtm2s;

            std::vector<int> maxCurveIndices1;
            std::vector<int> maxCurveIndices2;

            std::vector<float> maxCurveSaliencies;
            std::vector<float> maxCurveIntensities;
            std::vector<float> maxCurveEigenvalues;
            std::vector<Vector3DFloat> axisPts;
            std::vector< std::vector<float> > maxCurveSepPts;

            std::vector<float> coefsA;
            std::vector<float> coefsB;
            std::vector<float> coefsC;
            //std::vector<float> coefsD;
            std::vector<float> curveMinsX;
            std::vector<float> curveMinsY;
            std::vector<float> curveMinsZ;
            std::vector<float> curveMaxsX;
            std::vector<float> curveMaxsY;
            std::vector<float> curveMaxsZ;

            std::vector<float>curveMeanXs;
            std::vector<float>curveMeanYs;
            std::vector<float>curveMeanZs;
            //std::vector<Vector3DFloat> pseudoatoms;
            //void toggleHideOn(bool display);
            std::vector<float> eigenValues0;
            std::vector<float> eigenValues1;
            std::vector<float> eigenValues2;

            std::vector<Vector3DFloat> helixMaxs;
            std::vector<Vector3DFloat> helixMins;
            std::vector< std::vector<Vector3DFloat> > foundHelices;
            std::vector<int> foundHelixIndices;

            std::vector<float> normalizedMinXs;
            std::vector<float>	normalizedMinYs;
            std::vector<float>	normalizedMinZs;
            std::vector<float>	normalizedMaxXs;
            std::vector<float>	normalizedMaxYs;
            std::vector<float>	normalizedMaxZs;

            std::vector<Matrix4> rotationMatrices;
            std::vector<float> rotAngles;

            bool extremalMode;

            float curveRatio;
            float pointRatio;
            float surfaceRatio;
            float minGeo;
            float maxGeo;
            float eigenValue;

            float logEigenMin0;
            float logEigenMin1;
            float logEigenMin2;

            float logEigenMax0;
            float logEigenMax1;
            float logEigenMax2;

            bool extremalHide;
            bool extremalSurfaceHide;

            bool saliencyCheck;
            bool intensityCheck;
            bool eigenvalueCheck;

            bool maxSurfaceOn = false;
            bool minSurfaceOn = false;

            bool minPointOn = false;
            bool maxPointOn = false;
            bool saddlePointOn = false;
            bool displayHelix = false;
            bool displayHelixPoints = false;

            float minDataIntensity, maxDataIntensity, minDataEigenvalue, maxDataEigenvalue;
            float maxCurveEigenVMax0, maxCurveEigenVMax1, maxCurveEigenVMax2, maxCurveEigenVMin0, maxCurveEigenVMin1, maxCurveEigenVMin2;

            float minPointX = -1000.0;
            float minPointY = -1000.0;
            float minPointZ = -1000.0;
            float maxPointX = -1000.0;
            float maxPointY = -1000.0;
            float maxPointZ = -1000.0;

            float helixThreshold;
        private:
            void DrawBackboneModel(int subSceneIndex, bool selectEnabled);

            void DrawRibbonModel(int subSceneIndex, bool selectEnabled);
            void DrawSideChainModel(int subSceneIndex, bool selectEnabled);
        private:
            AtomMapType atoms;

            //TODO: possibly implement mouse picking using ray intersection
            std::vector<unsigned long long> atomHashKeys; //glLoadName(index of this std::vector)... used for selection
            std::vector<unsigned long long> ix0s;
            std::vector<unsigned long long> ix1s;

            std::vector<PDBBond> bonds;
            std::vector<PDBBond> sidechainBonds;
            std::vector<PDBBond> bondsToDelete;

            std::vector<Secel> aHelices;
            std::vector<Secel> bStrands;
            std::vector<Secel> loops;

            std::vector<int> selectedHelixIndices;
            //std::vector<int> selectedSecelIndices; //unsure if I can just keep track of secels as one structure or not
            std::vector<int> selectedStrandIndices;
            std::vector<int> selectedLoopIndices;
            std::vector < boost::tuple<int, int> > corrs;
            std::vector<int> selectedSSEHelices;
            std::vector< boost::tuple<Vector3DFloat, Vector3DFloat> > featureVecs;

            map<int, boost::tuple<float, float, float> > helixColors;

            int renderingType;
            float thinRibbThickness;

            int NUM_SEGMENTS;
            int NUM_SLICES;

            /* These three constants used in rendering alpha helices */
            float HELIX_ALPHA;
            float HELIX_BETA;
            float HELIX_HERMITE_FACTOR;
            float HELIX_WIDTH;

            float STRAND_HERMITE_FACTOR;

            float LOOP_RADIUS;

            float hlt_r, hlt_g, hlt_b, hlt_a;
        };


        CAlphaRenderer::CAlphaRenderer() {
            atoms.clear();
            bonds.clear();
            sidechainBonds.clear();
            selectedHelixIndices.clear();
            selectedStrandIndices.clear();
            featureVecs.clear();

            NUM_SEGMENTS = 10;
            NUM_SLICES = 10;
            HELIX_HERMITE_FACTOR = 4.7;
            HELIX_ALPHA = 32.0 * PI / 180.0;
            HELIX_BETA = -11.0 * PI / 180.0; // these three values taken from molscript code
            HELIX_WIDTH = 4.0;

            STRAND_HERMITE_FACTOR = .5;
            LOOP_RADIUS = .25;

            renderingType = 1;
            thinRibbThickness = .05;

            hlt_r = 1.0;
            hlt_g = 1.0;
            hlt_b = 1.0;
            hlt_a = 1.0;
        }


        CAlphaRenderer::~CAlphaRenderer() {
            atoms.clear();
            bonds.clear();
            sidechainBonds.clear();
            selectedHelixIndices.clear();
            selectedStrandIndices.clear();
            selectedLoopIndices.clear();
            featureVecs.clear();
        }

        void CAlphaRenderer::addQuad(int q1, int q2, int q3, int q4, float nx, float ny, float nz) {
            quadPts.push_back(q1);
            quadPts.push_back(q2);
            quadPts.push_back(q3);
            quadPts.push_back(q4);
            quadNormals.push_back(nx);
            quadNormals.push_back(ny);
            quadNormals.push_back(nz);
        }

        float findDistance(Vector3DFloat linePt1, Vector3DFloat linePt2, Vector3DFloat pt){
            Vector3DFloat sub1 = pt - linePt1;
            Vector3DFloat sub2 = linePt2 - linePt1;
            float mag = (sub1.X() * sub2.X()) + (sub1.Y() * sub2.Y()) + (sub1.Z() * sub2.Z());
            mag = mag / (sub2.Length() * sub2.Length());
            sub2 = Vector3DFloat(mag*sub2.X(), mag*sub2.Y(), mag*sub2.Z());
            mag = (sub1 - sub2).Length();
            return mag;
        }

        void CAlphaRenderer::setOrthoDistMin(float min) {
            orthoDistMin = min;
        }

        void CAlphaRenderer::setOrthoDistMax(float max) {
            orthoDistMax = max;
        }

        void CAlphaRenderer::SetExtremalMode(bool extremal) {
            extremalMode = extremal;
        }

        void CAlphaRenderer::setSegmentThreshold(int threshold) {
            helixSegmentThreshold = threshold;
        }

        std::vector<int> CAlphaRenderer::GetSelectedHelixIndices(){
            return selectedHelixIndices;
        }

        PDBAtom * CAlphaRenderer::AddAtom(PDBAtom atom) {
            atoms[atom.GetHashKey()] = atom;
            UpdateBoundingBox();
            return &atoms[atom.GetHashKey()];
        }

        void CAlphaRenderer::AddBond(PDBBond bond) {
            bonds.push_back(bond);
        }

        void CAlphaRenderer::clearQuads() {
            quadPts.clear();
            quadNormals.clear();
        }

        int CAlphaRenderer::findMaxHelices() {
            return (int)helices.size();
        }

        void CAlphaRenderer::setScaleEllipsoid(float eScale) {
            ellipsoidScale = eScale;
        }

        void CAlphaRenderer::setEllipsoidDisplay(bool displayEllipsoid) {
            ellipsoidVisible = displayEllipsoid;
        }

        void CAlphaRenderer::setMinDisplay(bool minPt) {
            minPointOn = minPt;
        }

        void CAlphaRenderer::setMaxDisplay(bool maxPt){
            maxPointOn = maxPt;
        }

        void CAlphaRenderer::addMaxPtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8) {
            unsigned long long maxHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            atoms[maxHash].eigenVector0[0] = ev0;
            atoms[maxHash].eigenVector0[1] = ev1;
            atoms[maxHash].eigenVector0[2] = ev2;

            atoms[maxHash].eigenVector1[0] = ev3;
            atoms[maxHash].eigenVector1[1] = ev4;
            atoms[maxHash].eigenVector1[2] = ev5;

            atoms[maxHash].eigenVector2[0] = ev6;
            atoms[maxHash].eigenVector2[1] = ev7;
            atoms[maxHash].eigenVector2[2] = ev8;
        }

        void CAlphaRenderer::addMinPtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8) {
            unsigned long long minHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            atoms[minHash].eigenVector0[0] = ev0;
            atoms[minHash].eigenVector0[1] = ev1;
            atoms[minHash].eigenVector0[2] = ev2;

            atoms[minHash].eigenVector1[0] = ev3;
            atoms[minHash].eigenVector1[1] = ev4;
            atoms[minHash].eigenVector1[2] = ev5;

            atoms[minHash].eigenVector2[0] = ev6;
            atoms[minHash].eigenVector2[1] = ev7;
            atoms[minHash].eigenVector2[2] = ev8;
        }

        void CAlphaRenderer::setHelixCoilHeight(float coilHeight) {
            helixCoilHeight = 20.0 * (float)coilHeight / 100.0;
        }

        void CAlphaRenderer::addSaddlePtEigenVector(int index, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8) {
            unsigned long long saddleHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            atoms[saddleHash].eigenVector0[0] = ev0;
            atoms[saddleHash].eigenVector0[1] = ev1;
            atoms[saddleHash].eigenVector0[2] = ev2;

            atoms[saddleHash].eigenVector1[0] = ev3;
            atoms[saddleHash].eigenVector1[1] = ev4;
            atoms[saddleHash].eigenVector1[2] = ev5;

            atoms[saddleHash].eigenVector2[0] = ev6;
            atoms[saddleHash].eigenVector2[1] = ev7;
            atoms[saddleHash].eigenVector2[2] = ev8;
        }


        void CAlphaRenderer::makeMax(int index, bool maxB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1) {
            unsigned long long maxHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            if (!atoms.count(maxHash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index * 3];
                float py = vertexPositions[index * 3 + 1];
                float pz = vertexPositions[index * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[maxHash] = newAtom;

            }
            float px = atoms[maxHash].GetPosition().X();
            float py = atoms[maxHash].GetPosition().Y();
            float pz = atoms[maxHash].GetPosition().Z();
            atoms[maxHash].isMax = maxB;
            maxPointSaliencies.push_back(saliency1);
            maxPointSaliencies.push_back(saliency2);
            maxPointSaliencies.push_back(saliency3);
            maxPointIntensities.push_back(intensity);
            maxPointEigenvalues.push_back(eigenvalue);
            maxPointEigenvalues0.push_back(eigenValue0);
            maxPointEigenvalues1.push_back(eigenValue1);
            maxCoord.push_back(Vector3DFloat(px, py, pz));
        }

        void CAlphaRenderer::setSaddleDisplay(bool saddlePt){
            saddlePointOn = saddlePt;
        }

        void CAlphaRenderer::addMaxEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2) {
            unsigned long long atom1Hash = std::stoll(atm1);
            unsigned long long atom2Hash = std::stoll(atm2);
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            Vector3DFloat v1 = Vector3DFloat(ev0, ev1, ev2);
            v1.Normalize();
            Vector3DFloat v2 = Vector3DFloat(ev3, ev4, ev5);
            v2.Normalize();
            Vector3DFloat v3 = Vector3DFloat(ev6, ev7, ev8);
            v3.Normalize();
            bonds[bondIndex].eigenVector0[0] = v1.X();
            bonds[bondIndex].eigenVector0[1] = v1.Y();
            bonds[bondIndex].eigenVector0[2] = v1.Z();
            bonds[bondIndex].eigenVector1[0] = v2.X();
            bonds[bondIndex].eigenVector1[1] = v2.Y();
            bonds[bondIndex].eigenVector1[2] = v2.Z();
            bonds[bondIndex].eigenVector2[0] = v3.X();
            bonds[bondIndex].eigenVector2[1] = v3.Y();
            bonds[bondIndex].eigenVector2[2] = v3.Z();
        }

        void CAlphaRenderer::addMinEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2) {
            unsigned long long atom1Hash = std::stoll(atm1);
            unsigned long long atom2Hash = std::stoll(atm2);
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            Vector3DFloat v1 = Vector3DFloat(ev0, ev1, ev2);
            v1.Normalize();
            Vector3DFloat v2 = Vector3DFloat(ev3, ev4, ev5);
            v2.Normalize();
            Vector3DFloat v3 = Vector3DFloat(ev6, ev7, ev8);
            v3.Normalize();
            bonds[bondIndex].eigenVector0[0] = v1.X();
            bonds[bondIndex].eigenVector0[1] = v1.Y();
            bonds[bondIndex].eigenVector0[2] = v1.Z();
            bonds[bondIndex].eigenVector1[0] = v2.X();
            bonds[bondIndex].eigenVector1[1] = v2.Y();
            bonds[bondIndex].eigenVector1[2] = v2.Z();
            bonds[bondIndex].eigenVector2[0] = v3.X();
            bonds[bondIndex].eigenVector2[1] = v3.Y();
            bonds[bondIndex].eigenVector2[2] = v3.Z();
        }

        void CAlphaRenderer::setHelixDisplayThreshold(float helixThresh) {
            //cout << "helixThresh set " << helixThresh << endl;
            helixThreshold = helixThresh;
        }

        void CAlphaRenderer::addSaddleEigenVector(float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8, string atm1, string atm2) {
            unsigned long long atom1Hash = std::stoll(atm1);
            unsigned long long atom2Hash = std::stoll(atm2);
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            Vector3DFloat v1 = Vector3DFloat(ev0, ev1, ev2);
            v1.Normalize();
            Vector3DFloat v2 = Vector3DFloat(ev3, ev4, ev5);
            v2.Normalize();
            Vector3DFloat v3 = Vector3DFloat(ev6, ev7, ev8);
            v3.Normalize();
            bonds[bondIndex].eigenVector0[0] = v1.X();
            bonds[bondIndex].eigenVector0[1] = v1.Y();
            bonds[bondIndex].eigenVector0[2] = v1.Z();
            bonds[bondIndex].eigenVector1[0] = v2.X();
            bonds[bondIndex].eigenVector1[1] = v2.Y();
            bonds[bondIndex].eigenVector1[2] = v2.Z();
            bonds[bondIndex].eigenVector2[0] = v3.X();
            bonds[bondIndex].eigenVector2[1] = v3.Y();
            bonds[bondIndex].eigenVector2[2] = v3.Z();
        }

        void CAlphaRenderer::makeSaddle(int index, bool saddleB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1) {
            unsigned long long saddleHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            if (!atoms.count(saddleHash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index * 3];
                float py = vertexPositions[index * 3 + 1];
                float pz = vertexPositions[index * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[saddleHash] = newAtom;

            }
            float px = atoms[saddleHash].GetPosition().X();
            float py = atoms[saddleHash].GetPosition().Y();
            float pz = atoms[saddleHash].GetPosition().Z();
            atoms[saddleHash].isSaddle = saddleB;
            saddlePointSaliencies.push_back(saliency1);
            saddlePointSaliencies.push_back(saliency2);
            saddlePointSaliencies.push_back(saliency3);
            saddlePointIntensities.push_back(intensity);
            saddlePointEigenvalues.push_back(eigenvalue);
            saddlePointEigenvalues0.push_back(eigenValue0);
            saddlePointEigenvalues1.push_back(eigenValue1);
            saddleCoord.push_back(Vector3DFloat(px, py, pz));
        }


        void CAlphaRenderer::setMinMaxIntensitiesEigenvalues(float minI, float maxI, float minE, float maxE) {
            minDataIntensity = minI;
            maxDataIntensity = maxI;
            minDataEigenvalue = minE;
            maxDataEigenvalue = maxE;
        }

        void CAlphaRenderer::setMaxCurveEigenMaxesAndMins(float maxCurveEigen0Max, float maxCurveEigen1Max, float maxCurveEigen2Max, float maxCurveEigen0Min, float maxCurveEigen1Min, float maxCurveEigen2Min) {
            maxCurveEigenVMax0 = maxCurveEigen0Max;
            maxCurveEigenVMax1 = maxCurveEigen1Max;
            maxCurveEigenVMax2 = maxCurveEigen2Max;

            maxCurveEigenVMin0 = maxCurveEigen0Min;
            maxCurveEigenVMin1 = maxCurveEigen1Min;
            maxCurveEigenVMin2 = maxCurveEigen2Min;
        }




        void CAlphaRenderer::findCurveMinsMaxes() {
            for (int i = 0; i < maxCurveSepPts.size(); i++) {
                std::vector<float> currentSeg = maxCurveSepPts[i];
                float xMin = 10000.0;
                float yMin = 10000.0;
                float zMin = 10000.0;
                float xMax = -10000.0;
                float yMax = -10000.0;
                float zMax = -10000.0;

                float xMean = 0.0;
                float yMean = 0.0;
                float zMean = 0.0;

                for (int j = 0; j < currentSeg.size() / 3; j++) {
                    if (currentSeg[3 * j] < xMin) {
                        xMin = currentSeg[3 * j];
                    }
                    if (currentSeg[3 * j + 1] < yMin) {
                        yMin = currentSeg[3 * j + 1];
                    }
                    if (currentSeg[3 * j + 2] < zMin) {
                        zMin = currentSeg[3 * j + 1];
                    }
                    if (currentSeg[3 * j] > xMax) {
                        xMax = currentSeg[3 * j];
                    }
                    if (currentSeg[3 * j + 1] > yMax) {
                        yMax = currentSeg[3 * j + 1];
                    }
                    if (currentSeg[3 * j + 2] > zMax) {
                        zMax = currentSeg[3 * j + 2];
                    }
                    xMean += currentSeg[3 * j];
                    yMean += currentSeg[3 * j + 1];
                    zMean += currentSeg[3 * j + 2];

                }
                xMean /= (currentSeg.size() / 3);
                yMean /= (currentSeg.size() / 3);
                zMean /= (currentSeg.size() / 3);
                curveMinsX.push_back(xMin);
                curveMinsY.push_back(yMin);
                curveMinsZ.push_back(zMin);
                curveMaxsX.push_back(xMax);
                curveMaxsY.push_back(yMax);
                curveMaxsZ.push_back(zMax);

                curveMeanXs.push_back(xMean);
                curveMeanYs.push_back(yMean);
                curveMeanZs.push_back(zMean);
            }
        }

        void CAlphaRenderer::setHelixDisplay(bool display) {
            displayHelix = display;
        }

        void CAlphaRenderer::setHelixDisplayPoints(bool display) {
            displayHelixPoints = display;
        }

        void CAlphaRenderer::inRangeOfSkeleton() {
            for (int i = 0; i < maxCurveSepPts.size(); i++) {
                std::vector<float> currentStrand = maxCurveSepPts[i];
                for (int j = 0; j < currentStrand.size() / 6; j++) {
                    Vector3DFloat atm1 = Vector3DFloat(currentStrand[6 * j], currentStrand[6 * j + 1], currentStrand[6 * j + 2]);
                    Vector3DFloat atm2 = Vector3DFloat(currentStrand[6 * j + 3], currentStrand[6 * j + 4], currentStrand[6 * j + 5]);

                    for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                        Vector3DFloat currentPos = it->second.GetPosition();
                        if (findDistance(atm1, atm2, currentPos) < 1.0) {
                            it = atoms.erase(it);
                            //atomHashes.push_back(it->second.GetHashKey());
                        }


                    }
                }



                //Vector3DFloat atm1 = atoms[bonds[i].GetAtom0Ix()].GetPosition();
                //Vector3DFloat atm2 = atoms[bonds[i].GetAtom1Ix()].GetPosition();
                //for(int j = 0; j < pseudoatoms.size(); j++) {
                //if(findDistance(atm1, atm2, pseudoatoms[j]) < 1.0) {
                //pseudoatoms.erase(pseudoatoms.begin() + j);
                //}
                //}
            }
            for (int i = 0; i < maxCurveSepPts.size(); i++) {

                std::vector<float> currentStrand = maxCurveSepPts[i];
                for (int j = 0; j < currentStrand.size() / 6; j++) {
                    Vector3DFloat atm1 = Vector3DFloat(currentStrand[6 * j], currentStrand[6 * j + 1], currentStrand[6 * j + 2]);
                    Vector3DFloat atm2 = Vector3DFloat(currentStrand[6 * j + 3], currentStrand[6 * j + 4], currentStrand[6 * j + 5]);

                    PDBAtom newAtm = PDBAtom("___3", 'A', (unsigned int)atoms.size(), "CA");

                    newAtm.SetPosition(atm1);
                    AddAtom(newAtm);

                    PDBAtom newAtm1 = PDBAtom("___3", 'A', (unsigned int)atoms.size(), "CA");
                    newAtm1.SetPosition(atm2);
                    AddAtom(newAtm1);

                }
            }
        }

        //void CAlphaRenderer::addCurveCoeff(float coefA, float coefB, float coefC, float coefD) {
        void CAlphaRenderer::addCurveCoeff(float coefA, float coefB, float coefC) {

            coefsA.push_back(coefA);
            coefsB.push_back(coefB);
            coefsC.push_back(coefC);
            //coefsD.push_back(coefD);
        }

        void CAlphaRenderer::makeMin(int index, bool minB, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenValue0, float eigenValue1) {
            unsigned long long minHash = ConstructHashKey("___3", 'A', (unsigned int)index, "CA");
            if (!atoms.count(minHash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index * 3];
                float py = vertexPositions[index * 3 + 1];
                float pz = vertexPositions[index * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[minHash] = newAtom;

            }
            float px = atoms[minHash].GetPosition().X();
            float py = atoms[minHash].GetPosition().Y();
            float pz = atoms[minHash].GetPosition().Z();
            atoms[minHash].isMin = minB;
            minPointSaliencies.push_back(saliency1);
            minPointSaliencies.push_back(saliency2);
            minPointSaliencies.push_back(saliency3);
            minPointIntensities.push_back(intensity);
            minPointEigenvalues.push_back(eigenvalue);
            minPointEigenvalues0.push_back(eigenValue0);
            minPointEigenvalues1.push_back(eigenValue1);
            //atoms[minHash].saliency[0] = saliency1;
            //atoms[minHash].saliency[1] = saliency2;
            //atoms[minHash].saliency[2] = saliency3;
            //atoms[minHash].saliency[2] = saliency3;
            //atoms[minHash].intensity = intensity;
            //atoms[minHash].eigenvalue = eigenvalue;
            minCoord.push_back(Vector3DFloat(px, py, pz));


            //unsigned long long atom2Hash = 1;
            /**
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == index) {
            //atom2Hash = it->second.GetHashKey();
            float px = it->second.GetPosition().X();
            float py = it->second.GetPosition().Y();
            float pz = it->second.GetPosition().Z();
            it->second.isMin = minB;
            minPointX = px;
            minPointY = py;
            minPointZ = pz;
            //OpenGLUtils::SetColor(0.1, 0.2, 0.8, 1.0);
            //DrawSphere(it->second.GetPosition(), 0.3);
            //atomHashes.push_back(it->second.GetHashKey());
            }
            }
            **/
        }

        py::list CAlphaRenderer::getQuadPts() {
            return std_vector_to_py_list(quadPts);
        }

        py::list CAlphaRenderer::getQuadNormals() {
            return std_vector_to_py_list(quadNormals);
        }

        template<class T>
        py::list std_vector_to_py_list(const std::vector<T>& v)
        {
            py::object get_iter = py::iterator<std::vector<T> >();
            py::object iter = get_iter(v);
            py::list l(iter);
            return l;
        }

        py::list CAlphaRenderer::findMaxCurveStrand(int i) {
            return std_vector_to_py_list(maxCurveSepPts[i]);
        }

        void CAlphaRenderer::AddSideChainBond(PDBBond bond) {
            sidechainBonds.push_back(bond);
        }


        //void CAlphaRenderer::AddQuadSurface(string quad1, string quad2, string quad3, string quad4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type) {

        //void CAlphaRenderer::AddQuadSurface(int index1, int index2, int index3, int index4, string quad1, string quad2, string quad3, string quad4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type) {
        static unsigned long long GetCharIndex(char c) {
            unsigned long long value;
            if ((c >= 48) && (c <= 57)) {		// 0..9
                value = (unsigned long long)c - 48;
            }
            else if ((c >= 65) && (c <= 90))  {  // A..Z
                value = (unsigned long long)c - 55;
            }
            else if ((c >= 97) && (c <= 122))  {  // a..z  (same as A..Z)
                value = (unsigned long long)c - 87;
            }
            else {
                value = 36;
            }
            return value;
        }

        static unsigned long long GetPDBIdIndex(string pdbId) {
            return GetCharIndex(pdbId.c_str()[0]) * 37 * 37 * 37 +
                GetCharIndex(pdbId.c_str()[1]) * 37 * 37 +
                GetCharIndex(pdbId.c_str()[2]) * 37 +
                GetCharIndex(pdbId.c_str()[3]);
        }

        unsigned long long GetChainIdIndex(char chainId) {
            return GetCharIndex(chainId);
        }


        void CAlphaRenderer::writeMaxCurveToPDB() {
            bool displayMax = false;
            bool display = false;
            ofstream myfile;

            myfile.open("1.pdb");
            int atomIndex = 0;
            for (int i = 0; i < bonds.size(); i++) {
                if (bonds[i].maxOn) {

                    displayMax = true;

                    display = true;
                    float pointR = pointRatio;
                    float curveR = curveRatio;
                    float surfaceR = surfaceRatio;
                    float minG = minGeo;
                    float maxG = maxGeo;
                    float eigenV = eigenValue;
                    if (saliencyCheck) {
                        float *saliencies = bonds[i].saliencies;
                        if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])){
                            display = false;
                            displayMax = false;
                        }
                    }
                    if (intensityCheck) {
                        float maxI = maxDataIntensity;
                        float minI = minDataIntensity;
                        float localI = bonds[i].intensity;
                        if (((localI - minI) / (maxI - minI) < maxG)) {
                            display = false;
                            displayMax = false;
                        }
                    }
                    if (eigenvalueCheck) {
                        float maxE = maxDataEigenvalue;
                        float minE = minDataEigenvalue;
                        float localE = bonds[i].eigenvalue;
                        if ((localE - minE) / (maxE - minE) < eigenV){
                            display = false;
                            displayMax = false;
                        }
                    }

                }
                else {
                    display = false;
                }
                if (display) {
                    string vIndex = std::to_string(atomIndex);
                    padTo(vIndex, 5);

                    float xPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().X();
                    float yPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().Y();
                    float zPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().Z();

                    string posx0 = std::to_string(xPos0);
                    string posy0 = std::to_string(yPos0);
                    string posz0 = std::to_string(zPos0);

                    myfile << "ATOM  " << vIndex << "  CA  ALA " << vIndex << "     " << posx0.substr(0, 7) << " " << posy0.substr(0, 7) << " " << posz0.substr(0, 7) << "  1.00  1.00      S_00  0 " << endl;

                    atomIndex++;
                    vIndex = std::to_string(atomIndex);
                    padTo(vIndex, 5);

                    float xPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().X();
                    float yPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().Y();
                    float zPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().Z();

                    string posx1 = std::to_string(xPos1);
                    string posy1 = std::to_string(yPos1);
                    string posz1 = std::to_string(zPos1);

                    myfile << "ATOM  " << vIndex << "  CA  ALA " << vIndex << "     " << posx1.substr(0, 7) << " " << posy1.substr(0, 7) << " " << posz1.substr(0, 7) << "  1.00  1.00      S_00  0 " << endl;

                    atomIndex++;
                }
            }
            myfile.close();


        }



        static unsigned long long ConstructHashKey(string pdbId, char chainId, unsigned int resSeq, string name) {
            unsigned long long chainIDCount = 37;
            unsigned long long residueNumCount = 10000;
            unsigned long long atomTypeCount = 45;

            //GetAtomTypeIndex(name)

            return GetPDBIdIndex(pdbId) * chainIDCount * residueNumCount * atomTypeCount +
                GetChainIdIndex(chainId) * residueNumCount * atomTypeCount +
                (unsigned long long)resSeq * atomTypeCount +
                1;
        }

        void CAlphaRenderer::addQuadEigenvectors(int index1, int index2, int index3, int index4, float ev0, float ev1, float ev2, float ev3, float ev4, float ev5, float ev6, float ev7, float ev8) {
            unsigned long long quad1Hash = ConstructHashKey("___3", 'A', (unsigned int)index1, "CA");
            unsigned long long quad2Hash = ConstructHashKey("___3", 'A', (unsigned int)index2, "CA");
            unsigned long long quad3Hash = ConstructHashKey("___3", 'A', (unsigned int)index3, "CA");
            unsigned long long quad4Hash = ConstructHashKey("___3", 'A', (unsigned int)index4, "CA");
            for (int i = 0; i < quads.size(); i++) {
                if (quads[i].index1 == index1 && quads[i].index2 == index2 && quads[i].index3 == index3 && quads[i].index4 == index4) {
                    //if (quads[i].hash1 == quad1Hash && quads[i].hash2 == quad2Hash && quads[i].hash3 == quad3Hash && quads[i].hash4 == quad4Hash) {
                    quads[i].eigenVector0[0] = ev0;
                    quads[i].eigenVector0[1] = ev1;
                    quads[i].eigenVector0[2] = ev2;

                    quads[i].eigenVector1[0] = ev3;
                    quads[i].eigenVector1[1] = ev4;
                    quads[i].eigenVector1[2] = ev5;

                    quads[i].eigenVector2[0] = ev6;
                    quads[i].eigenVector2[1] = ev7;
                    quads[i].eigenVector2[2] = ev8;
                }
                //}
            }
        }

        void CAlphaRenderer::findEigenMinMaxes() {
            float logEigen0Min = 1000.0;
            float logEigen1Min = 1000.0;
            float logEigen2Min = 1000.0;

            float logEigen0Max = -1000.0;
            float logEigen1Max = -1000.0;
            float logEigen2Max = -1000.0;

            for (int i = 0; i < bonds.size(); i++) {
                if (bonds[i].logEigen0 < logEigen0Min) {
                    logEigen0Min = bonds[i].logEigen0;
                }
                if (bonds[i].logEigen1 < logEigen1Min) {
                    logEigen1Min = bonds[i].logEigen1;
                }
                if (bonds[i].logEigen2 < logEigen2Min) {
                    logEigen2Min = bonds[i].logEigen2;
                }

                if (bonds[i].logEigen0 > logEigen0Max) {
                    logEigen0Max = bonds[i].logEigen0;
                }
                if (bonds[i].logEigen1 > logEigen1Max) {
                    logEigen1Max = bonds[i].logEigen1;
                }
                if (bonds[i].logEigen2 > logEigen2Max) {
                    logEigen2Max = bonds[i].logEigen2;
                }
            }
            logEigenMin0 = logEigen0Min;
            logEigenMin1 = logEigen1Min;
            logEigenMin2 = logEigen2Min;

            logEigenMax0 = logEigen0Max;
            logEigenMax1 = logEigen1Max;
            logEigenMax2 = logEigen2Max;
        }

        void CAlphaRenderer::AddQuadSurface(int index1, int index2, int index3, int index4, float nx, float ny, float nz, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, int type, float eigenValue0, float eigenValue1) {
            unsigned long long quad1Hash = ConstructHashKey("___3", 'A', (unsigned int)index1, "CA");
            unsigned long long quad2Hash = ConstructHashKey("___3", 'A', (unsigned int)index2, "CA");
            unsigned long long quad3Hash = ConstructHashKey("___3", 'A', (unsigned int)index3, "CA");
            unsigned long long quad4Hash = ConstructHashKey("___3", 'A', (unsigned int)index4, "CA");

            /**
            unsigned long long quad1Str = std::stoll(quad1);
            unsigned long long quad2Str = std::stoll(quad2);
            unsigned long long quad3Str= std::stoll(quad3);
            unsigned long long quad4Str = std::stoll(quad4);
            unsigned long long atom1Hash = 1;
            float p1x = -1000.0, p1y = -1000.0, p1z = -1000.0, p2x = -1000.0, p2y = -1000.0,
            p2z = -1000.0, p3x = -1000.0, p3y = -1000.0, p3z = -1000.0, p4x = -1000.0, p4y = -1000.0, p4z = -1000.0;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

            if(it->second.GetSerial() == quad1Str) {
            atom1Hash = it->second.GetHashKey();
            p1x = it->second.GetPosition().X();
            p1y = it->second.GetPosition().Y();
            p1z = it->second.GetPosition().Z();
            }


            }
            unsigned long long atom2Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == quad2Str) {
            atom2Hash = it->second.GetHashKey();
            p2x = it->second.GetPosition().X();
            p2y = it->second.GetPosition().Y();
            p2z = it->second.GetPosition().Z();
            //OpenGLUtils::SetColor(1.0, 0.0, 0.1, 1.0);
            //DrawSphere(it->second.GetPosition(), 0.2);
            //atomHashes.push_back(it->second.GetHashKey());
            }
            }
            unsigned long long atom3Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == quad3Str) {
            atom3Hash = it->second.GetHashKey();
            p3x = it->second.GetPosition().X();
            p3y = it->second.GetPosition().Y();
            p3z = it->second.GetPosition().Z();
            //OpenGLUtils::SetColor(1.0, 0.85, 0.1, 1.0);
            //DrawSphere(it->second.GetPosition(), 0.15);
            //atomHashes.push_back(it->second.GetHashKey());
            }
            }
            unsigned long long atom4Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == quad4Str) {
            atom4Hash = it->second.GetHashKey();
            p4x = it->second.GetPosition().X();
            p4y = it->second.GetPosition().Y();
            p4z = it->second.GetPosition().Z();

            }
            }
            if (atom1Hash == 1 || atom2Hash == 1 || atom3Hash == 1 || atom4Hash == 1) {
            return;
            }
            else {
            **/

            //unsigned long long quad1Hash = std::stoll(quad1);
            //unsigned long long quad2Hash = std::stoll(quad2);
            //unsigned long long quad3Hash = std::stoll(quad3);
            //unsigned long long quad4Hash = std::stoll(quad4);
            if (!atoms.count(quad1Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index1 * 3];
                float py = vertexPositions[index1 * 3 + 1];
                float pz = vertexPositions[index1 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[quad1Hash] = newAtom;
            }

            if (!atoms.count(quad2Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index2 * 3];
                float py = vertexPositions[index2 * 3 + 1];
                float pz = vertexPositions[index2 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[quad2Hash] = newAtom;
            }

            if (!atoms.count(quad3Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index3 * 3];
                float py = vertexPositions[index3 * 3 + 1];
                float pz = vertexPositions[index3 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[quad3Hash] = newAtom;
            }

            if (!atoms.count(quad4Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index4 * 3];
                float py = vertexPositions[index4 * 3 + 1];
                float pz = vertexPositions[index4 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[quad4Hash] = newAtom;
            }


            Vector3DFloat pos1 = atoms[quad1Hash].GetPosition();
            Vector3DFloat pos2 = atoms[quad2Hash].GetPosition();
            Vector3DFloat pos3 = atoms[quad3Hash].GetPosition();
            Vector3DFloat pos4 = atoms[quad4Hash].GetPosition();
            float qNorm[3] = { nx, ny, nz };
            float p1[3] = { pos1.X(), pos1.Y(), pos1.Z() };
            float p2[3] = { pos2.X(), pos2.Y(), pos2.Z() };
            float p3[3] = { pos3.X(), pos3.Y(), pos3.Z() };
            float p4[3] = { pos4.X(), pos4.Y(), pos4.Z() };
            Quad quad;
            quad.hash1 = quad1Hash;
            quad.hash2 = quad2Hash;
            quad.hash3 = quad3Hash;
            quad.hash4 = quad4Hash;
            quad.p1[0] = p1[0];
            quad.p1[1] = p1[1];
            quad.p1[2] = p1[2];
            quad.p2[0] = p2[0];
            quad.p2[1] = p2[1];
            quad.p2[2] = p2[2];
            quad.p3[0] = p3[0];
            quad.p3[1] = p3[1];
            quad.p3[2] = p3[2];
            quad.p4[0] = p4[0];
            quad.p4[1] = p4[1];
            quad.p4[2] = p4[2];
            quad.normal[0] = qNorm[0];
            quad.normal[1] = qNorm[1];
            quad.normal[2] = qNorm[2];
            quad.localIntensity = intensity;
            quad.saliency[0] = saliency1;
            quad.saliency[1] = saliency2;
            quad.saliency[2] = saliency3;
            quad.eigenValue = eigenvalue;
            quad.type = type;
            quad.eigenValue0 = eigenValue0;
            quad.eigenValue1 = eigenValue1;
            quad.index1 = index1;
            quad.index2 = index2;
            quad.index3 = index3;
            quad.index4 = index4;
            quads.push_back(quad);
            //}
            //}
        }

        float findAngle(float x1, float y1, float z1, float x2, float y2, float z2) {
            float dot = x1*x2 + y1*y2 + z1*z2;
            float lenSq1 = x1*x1 + y1*y1 + z1*z1;
            float lenSq2 = x2*x2 + y2*y2 + z2*z2;
            float angle = acos(dot / sqrt(lenSq1 * lenSq2));
            return angle;
        }

        static float getDotP(float v1[3], float v2[3]) {
            return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        }

        Vector3DFloat findMean(std::vector<Vector3DFloat> points) {
            float sumX = 0.0;
            float sumY = 0.0;
            float sumZ = 0.0;
            for (int i = 0; i < points.size(); i++) {
                sumX += points[i].X();
                sumY += points[i].Y();
                sumZ += points[i].Z();
            }
            float avgX = sumX / points.size();
            float avgY = sumY / points.size();
            float avgZ = sumZ / points.size();
            return Vector3DFloat(avgX, avgY, avgZ);
        }


        float * crossProduct(float v1[3], float v2[3]) {
            float cross[3];
            cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
            cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
            cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return cross;
        }

        Matrix4 CAlphaRenderer::alignWithAxis(float v1[3], float axis[3], Vector3DFloat minPt, Vector3DFloat maxPt) {

            float rot_angle = acos(getDotP(v1, axis));

            Vector3 v1Vector = Vector3(v1[0], v1[1], v1[2]);


            //cos(rot_angle) * v1Vector.Length();

            float * cross_p = crossProduct(v1, axis);
            Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
            crossPNormalized.Normalize();
            //float avgX = (maxPt.X() + minPt.X()) / 2.0;
            //float avgY = (maxPt.Y() + minPt.Y()) / 2.0;
            //float avgZ = (maxPt.Z() + minPt.Z()) / 2.0;
            Vector3 crossNormal = Vector3(crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
            float conversion = 180.0 / M_PI;
            //rotAngles.push_back(rot_angle);
            Matrix4 rotMatrix = Matrix4::rotation(crossNormal, rot_angle);
            Matrix4 inverseRotMatrix = Matrix4::rotation(crossNormal, -rot_angle);
            rotationMatrices.push_back(inverseRotMatrix);
            //Vector3DFloat zAxisAdj = Vector3DFloat(0.0, 1.0, 1.0);
            //zAxisAdj.Normalize();
            //Vector3 xy = Vector3(zAxisAdj.X(), zAxisAdj.Y(), zAxisAdj.Z());
            //Matrix4 rotMatrix2 = Matrix4::rotation(xy, M_PI);

            //Vec3f Translation = Vec3f(0.0, -1.0, 0.0);
            //Matrix4 transMatrix  = Matrix4::Translation(Translation);

            return rotMatrix;
            //return rotMatrix2 * rotMatrix;
            //Vector3DFloat segmentDirection = maxPt - minPt;
            //segmentDirection.Normalize();



            //Vector3DFloat segDirMinPt = Vector3DFloat(minPt.X(), minPt.Y(), minPt.Z());
            //segDirMinPt.Normalize();
            //Vector3 segDir = Vector3(segDirMinPt.X(), segDirMinPt.Y(), segDirMinPt.Z());

            //Vector3 segDir = rotMatrix * v1Vector;

            //cout << "segdir " << segDir[0] << " " << segDir[1] << " " << segDir[2] << endl;

            //return segDir;
            //Vector3DFloat segDirMaxPt = Vector3DFloat(maxPt.X(), maxPt.Y(), maxPt.Z());
            //segDirMaxPt.Normalize();
            //Vector3 segDir1 = Vector3(segDirMaxPt.X(), segDirMaxPt.Y(), segDirMaxPt.Z());
            //segDir1 = rotMatrix * segDir1;

            //cout << "segdir1 " << segDir1[0] << " " << segDir1[1] << " " << segDir1[2] << endl;

            //segmentDirection.Normalize();
            //glPushMatrix();
            //glRotatef( rot_angle * 57.2957795, crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

            //glTranslatef(avgX, avgY, avgZ);
            //CAlphaRenderer::DrawCylinder(minPt, maxPt, 0.1);
            //CAlphaRenderer::DrawCylinder(Vector3DFloat(0.0, 0.0, 0.0), segmentDirection, 0.1);
            //glPopMatrix();

        }
        //need translation
        void DrawEllipsoid(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC, float x, float y, float z, float target_dirx, float target_diry, float target_dirz, Vector3DFloat dir1, Vector3DFloat dir2)
        {
            GLfloat ambient_light[] = { 1.0f, 1.0f, 2.0f, 1.0f };
            glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

            GLfloat	diffuse_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);

            glEnable(GL_LIGHT0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            float tStep = (M_PI) / (float)uiSlices;
            float sStep = (M_PI) / (float)uiStacks;

            float target_dir[3];
            target_dir[0] = target_dirx;
            target_dir[1] = target_diry;
            target_dir[2] = target_dirz;

            float xaxis[3];
            xaxis[0] = 1.0;
            xaxis[1] = 0.0;
            xaxis[2] = 0.0;

            float yaxis[3];
            yaxis[0] = 0.0;
            yaxis[1] = 1.0;
            yaxis[2] = 0.0;

            float zaxis[3];
            zaxis[0] = 0.0;
            zaxis[1] = 0.0;
            zaxis[2] = 1.0;

            float target_dir_max[3];
            target_dir_max[0] = dir2.X();
            target_dir_max[1] = dir2.Y();
            target_dir_max[2] = dir2.Z();

            float target_dir_mid[3];
            target_dir_mid[0] = dir1.X();
            target_dir_mid[1] = dir1.Y();
            target_dir_mid[2] = dir1.Z();
            glPushMatrix();
            glTranslatef(x, y, z);
            float rot_angle = acos(getDotP(target_dir, xaxis));
            cout << "rot_angle " << rot_angle << endl;

            if (fA < fB && fA < fC) {
                float rot_angle = acos(getDotP(target_dir, xaxis));
                cout << "rot_angle " << rot_angle << endl;
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, xaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                    //z greatest dir
                    /**
                    if(fB < fC) {
                    rot_angle = acos(getDotP(target_dir_max , zaxis));
                    cross_p = crossProduct(target_dir_max, zaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , yaxis));
                    cross_p = crossProduct(target_dir_mid, yaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    else {
                    rot_angle = acos(getDotP(target_dir_max , yaxis));
                    cross_p = crossProduct(target_dir_max, yaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , zaxis));
                    cross_p = crossProduct(target_dir_mid, zaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    **/
                }

            }


            //y smallest
            else if (fB < fC && fB < fA) {
                float rot_angle = acos(getDotP(target_dir, yaxis));
                cout << "rot_angle " << rot_angle << endl;
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, yaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                    /**
                    if(fA < fC) {
                    rot_angle = acos(getDotP(target_dir_max , zaxis));
                    cross_p = crossProduct(target_dir_max, zaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , xaxis));
                    cross_p = crossProduct(target_dir_mid, xaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    else {
                    rot_angle = acos(getDotP(target_dir_max , xaxis));
                    cross_p = crossProduct(target_dir_max, xaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , zaxis));
                    cross_p = crossProduct(target_dir_mid, zaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    **/
                }

            }



            //z smallest

            else if (fC < fA && fC < fB) {
                float rot_angle = acos(getDotP(target_dir, zaxis));
                cout << "rot_angle " << rot_angle << endl;
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, zaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                    /**
                    if(fA < fB) {
                    rot_angle = acos(getDotP(target_dir_max , yaxis));
                    cross_p = crossProduct(target_dir_max, yaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , xaxis));
                    cross_p = crossProduct(target_dir_mid, xaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    else {
                    rot_angle = acos(getDotP(target_dir_max , xaxis));
                    cross_p = crossProduct(target_dir_max, xaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    rot_angle = acos(getDotP(target_dir_mid , yaxis));
                    cross_p = crossProduct(target_dir_mid, yaxis);
                    crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1],  cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef( rot_angle * (180.0 / M_PI), crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z() );

                    }
                    **/
                }
            }




            for (float t = -M_PI / 2; t <= (M_PI / 2) + .0001; t += tStep)
            {


                //float angle = findAngle(fA, fB, fC, xangle, yangle, zangle);

                //glRotatef(yangle * 57.2958, 0.0f, 1.0f, 0.0f);

                //glRotatef(zangle * 57.2958, 0.0f, 0.0f, 1.0f);

                //glPushMatrix();
                //glRotatef(angle * 57.2958, fA, fB, fC);
                glBegin(GL_TRIANGLE_STRIP);


                for (float s = -M_PI; s <= M_PI + .0001; s += sStep)
                {
                    //glNormal3f(x + fA * cos(t) * cos(s), y + fB * cos(t) * sin(s), z + fC * sin(t));

                    float xPos1 = fA * cos(t) * cos(s);
                    float yPos1 = fB * cos(t) * sin(s);
                    float zPos1 = fC * sin(t);
                    /**
                    // x smallest
                    if(fA < fB && fA < fC) {

                    //y axis
                    xPos1 = cos(yangle) * xPos1 + sin(yangle) * zPos1;
                    zPos1 = (-sin(yangle) * xPos1) + cos(yangle) * zPos1;

                    //z axis
                    xPos1 = cos(zangle) * xPos1 + (-sin(zangle) * yPos1);
                    yPos1 = sin(zangle) * xPos1 + (cos(zangle) * yPos1);


                    }
                    //y smallest
                    else if(fB < fC && fB < fA) {
                    //x axis
                    yPos1 = cos(xangle) * yPos1 + (-(sin(xangle)) * zPos1);
                    zPos1 = sin(xangle) * yPos1 + cos(xangle) * zPos1;

                    //z axis
                    xPos1 = cos(zangle) * xPos1 + (-sin(zangle) * yPos1);
                    yPos1 = sin(zangle) * xPos1 + (cos(zangle) * yPos1);

                    }
                    //z smallest
                    else if(fC < fA && fC < fB) {
                    //x axis
                    yPos1 = cos(xangle) * yPos1 + (-(sin(xangle)) * zPos1);
                    zPos1 = sin(xangle) * yPos1 + cos(xangle) * zPos1;

                    //y axis
                    xPos1 = cos(yangle) * xPos1 + sin(yangle) * zPos1;
                    zPos1 = (-sin(yangle) * xPos1) + cos(yangle) * zPos1;
                    }
                    **/

                    glNormal3f(xPos1, yPos1, zPos1);
                    glVertex3f(xPos1, yPos1, zPos1);




                    float xPos2 = fA * cos(t + tStep) * cos(s);
                    float yPos2 = fB * cos(t + tStep) * sin(s);
                    float zPos2 = fC * sin(t + tStep);
                    /**
                    // x smallest
                    if(fA < fB && fA < fC) {
                    //y axis
                    xPos2 = cos(yangle) * xPos2 + sin(yangle) * zPos2;
                    zPos2 = (-sin(yangle) * xPos2) + cos(yangle) * zPos2;

                    //z axis
                    xPos2 = cos(zangle) * xPos2 + (-sin(zangle) * yPos2);
                    yPos2 = sin(zangle) * xPos2 + (cos(zangle) * yPos2);


                    }
                    //y smallest
                    else if(fB < fC && fB < fA) {
                    //x axis
                    yPos2 = cos(xangle) * yPos2 + (-(sin(xangle)) * zPos2);
                    zPos2 = sin(xangle) * yPos2 + cos(xangle) * zPos2;

                    //z axis
                    xPos2 = cos(zangle) * xPos2 + (-sin(zangle) * yPos2);
                    yPos2 = sin(zangle) * xPos2 + (cos(zangle) * yPos2);

                    }
                    //z smallest
                    else if(fC < fA && fC < fB) {
                    //x axis
                    yPos2 = cos(xangle) * yPos2 + (-(sin(xangle)) * zPos2);
                    zPos2 = sin(xangle) * yPos2 + cos(xangle) * zPos2;

                    //y axis
                    xPos2 = cos(yangle) * xPos2 + sin(yangle) * zPos2;
                    zPos2 = (-sin(yangle) * xPos2) + cos(yangle) * zPos2;
                    }
                    **/


                    glVertex3f(xPos2, yPos2, zPos2);


                }

                glEnd();


            }

            glPopMatrix();
            //glClear();
        }

        void DrawEllipsoidSurface(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC, float x, float y, float z, float target_dirx, float target_diry, float target_dirz)
        {
            GLfloat ambient_light[] = { 1.0f, 1.0f, 2.0f, 1.0f };
            glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

            GLfloat	diffuse_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);

            glEnable(GL_LIGHT0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            float tStep = (M_PI) / (float)uiSlices;
            float sStep = (M_PI) / (float)uiStacks;

            float target_dir[3];
            target_dir[0] = target_dirx;
            target_dir[1] = target_diry;
            target_dir[2] = target_dirz;

            float xaxis[3];
            xaxis[0] = 1.0;
            xaxis[1] = 0.0;
            xaxis[2] = 0.0;

            float yaxis[3];
            yaxis[0] = 0.0;
            yaxis[1] = 1.0;
            yaxis[2] = 0.0;

            float zaxis[3];
            zaxis[0] = 0.0;
            zaxis[1] = 0.0;
            zaxis[2] = 1.0;
            glPushMatrix();
            glTranslatef(x, y, z);


            if (fA > fB && fA > fC) {
                float rot_angle = acos(getDotP(target_dir, xaxis));
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, xaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * 57.2957795, crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                }

            }


            //y smallest
            else if (fB > fC && fB > fA) {
                float rot_angle = acos(getDotP(target_dir, yaxis));
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, yaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * 57.2957795, crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                }

            }



            //z smallest

            else if (fC > fA && fC > fB) {
                float rot_angle = acos(getDotP(target_dir, zaxis));
                if (fabs(rot_angle) > 0.0001)
                {
                    float * cross_p = crossProduct(target_dir, zaxis);
                    Vector3DFloat crossPNormalized = Vector3DFloat(cross_p[0], cross_p[1], cross_p[2]);
                    crossPNormalized.Normalize();
                    glRotatef(rot_angle * 57.2957795, crossPNormalized.X(), crossPNormalized.Y(), crossPNormalized.Z());
                }
            }




            for (float t = -M_PI / 2; t <= (M_PI / 2) + .0001; t += tStep)
            {

                glBegin(GL_TRIANGLE_STRIP);


                for (float s = -M_PI; s <= M_PI + .0001; s += sStep)
                {

                    float xPos1 = fA * cos(t) * cos(s);
                    float yPos1 = fB * cos(t) * sin(s);
                    float zPos1 = fC * sin(t);


                    glNormal3f(xPos1, yPos1, zPos1);
                    glVertex3f(xPos1, yPos1, zPos1);




                    float xPos2 = fA * cos(t + tStep) * cos(s);
                    float yPos2 = fB * cos(t + tStep) * sin(s);
                    float zPos2 = fC * sin(t + tStep);


                    glVertex3f(xPos2, yPos2, zPos2);


                }

                glEnd();


            }

            glPopMatrix();
            //glClear();
        }

        void DrawEllipsoidPoint(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC, float x, float y, float z)
        {
            GLfloat ambient_light[] = { 1.0f, 1.0f, 2.0f, 1.0f };
            glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

            GLfloat	diffuse_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);

            glEnable(GL_LIGHT0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            float tStep = (M_PI) / (float)uiSlices;
            float sStep = (M_PI) / (float)uiStacks;

            glPushMatrix();
            glTranslatef(x, y, z);
            for (float t = -M_PI / 2; t <= (M_PI / 2) + .0001; t += tStep)
            {

                glBegin(GL_TRIANGLE_STRIP);


                for (float s = -M_PI; s <= M_PI + .0001; s += sStep)
                {

                    float xPos1 = fA * cos(t) * cos(s);
                    float yPos1 = fB * cos(t) * sin(s);
                    float zPos1 = fC * sin(t);


                    glNormal3f(xPos1, yPos1, zPos1);
                    glVertex3f(xPos1, yPos1, zPos1);




                    float xPos2 = fA * cos(t + tStep) * cos(s);
                    float yPos2 = fB * cos(t + tStep) * sin(s);
                    float zPos2 = fC * sin(t + tStep);


                    glVertex3f(xPos2, yPos2, zPos2);


                }

                glEnd();


            }

            glPopMatrix();
            //glClear();
        }

        void CAlphaRenderer::DrawSurface(Quad quad) {
            float pointR = pointRatio;
            float curveR = curveRatio;
            float surfaceR = surfaceRatio;
            float minG = minGeo;
            float maxG = maxGeo;
            float eigenV = eigenValue;

            int type = quad.type;
            bool hide = false;
            //if(extremalSurfaceHide) {
            if (saliencyCheck) {
                if ((surfaceR * quad.saliency[0] < curveR * quad.saliency[1]) ||
                    (surfaceR * quad.saliency[0] < pointR * quad.saliency[2])) hide = true;
                //if ( ( curveR * quad.saliency[1] < surfaceR * quad.saliency[0] ) || ( curveR * quad.saliency[1] < pointR * quad.saliency[2] ) ) hide = true;
            }
            if (intensityCheck) {
                float maxI = maxDataIntensity;
                float minI = minDataIntensity;
                float localI = quad.localIntensity;
                if ((type == 1) && ((localI - minI) / (maxI - minI) < maxG)) hide = true;
                if ((type == 2) && ((localI - minI) / (maxI - minI) > minG)) hide = true;
            }
            if (eigenvalueCheck) {
                float maxE = maxDataEigenvalue;
                float minE = minDataEigenvalue;
                float localE = quad.eigenValue;
                if ((localE - minE) / (maxE - minE) < eigenV) hide = true;
            }
            //}

            if (!hide) {
                //GLuint textureName;
                //glEnable(GL_TEXTURE_3D);

                //glBindTexture(GL_TEXTURE_3D, textureName);
                //glEnable(GL_CULL_FACE);
                GLfloat ambient_light[] = { 1.0f, 1.0f, 2.0f, 1.0f };
                glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
                glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

                GLfloat	diffuse_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
                glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);

                glEnable(GL_LIGHT0);


                glBegin(GL_QUADS);
                //glBegin(GL_POLYGON);
                //glTexCoord3d(quad.p1[0], quad.p1[1], quad.p1[2]);
                cout << "quadp1 " << quad.p1[0] << " " << quad.p2[0] << " " << quad.p3[0] << " " << quad.p4[0] << endl;

                glNormal3f(quad.normal[0], quad.normal[1], quad.normal[2]);
                glVertex3f(quad.p1[0], quad.p1[1], quad.p1[2]);

                //glTexCoord3d(quad.p2[0], quad.p2[1], quad.p2[2]);
                glVertex3f(quad.p2[0], quad.p2[1], quad.p2[2]);
                //glTexCoord3d(quad.p3[0], quad.p3[1], quad.p3[2]);
                glVertex3f(quad.p3[0], quad.p3[1], quad.p3[2]);

                //glTexCoord3d(quad.p4[0], quad.p4[1], quad.p4[2]);
                glVertex3f(quad.p4[0], quad.p4[1], quad.p4[2]);

                glEnd();

                float v1[3];
                v1[0] = quad.p2[0] - quad.p1[0];
                v1[1] = quad.p2[1] - quad.p1[1];
                v1[2] = quad.p2[2] - quad.p1[2];
                float v2[3];
                v2[0] = quad.p3[0] - quad.p1[0];
                v2[1] = quad.p3[1] - quad.p1[1];
                v2[2] = quad.p3[2] - quad.p1[2];
                float* quadNormal = crossProduct(v1, v2);
                Vector3DFloat normalizedQuadNormal = Vector3DFloat(quadNormal[0], quadNormal[1], quadNormal[2]);
                normalizedQuadNormal.Normalize();
                quadNormal[0] = normalizedQuadNormal.X();
                quadNormal[1] = normalizedQuadNormal.Y();
                quadNormal[2] = normalizedQuadNormal.Z();

                float e1 = 0.0;
                float e2 = 0.0;
                float e3 = 0.0;
                float bondDirection1[3];
                float bondDirection2[3];
                if (quad.eigenValue0 < quad.eigenValue1 && quad.eigenValue0 < quad.eigenValue) {

                    e1 = quad.eigenVector0[0];
                    e2 = quad.eigenVector0[1];
                    e3 = quad.eigenVector0[2];
                    if (quad.eigenValue1 < quad.eigenValue) {
                        bondDirection1[0] = quad.eigenVector1[0];
                        bondDirection1[1] = quad.eigenVector1[1];
                        bondDirection1[2] = quad.eigenVector1[2];

                        bondDirection2[0] = quad.eigenVector2[0];
                        bondDirection2[1] = quad.eigenVector2[1];
                        bondDirection2[2] = quad.eigenVector2[2];
                    }
                    else {
                        bondDirection1[0] = quad.eigenVector2[0];
                        bondDirection1[1] = quad.eigenVector2[1];
                        bondDirection1[2] = quad.eigenVector2[2];

                        bondDirection2[0] = quad.eigenVector1[0];
                        bondDirection2[1] = quad.eigenVector1[1];
                        bondDirection2[2] = quad.eigenVector1[2];
                    }
                }
                if (quad.eigenValue1 < quad.eigenValue0 && quad.eigenValue1 < quad.eigenValue) {
                    e1 = quad.eigenVector1[0];
                    e2 = quad.eigenVector1[1];
                    e3 = quad.eigenVector1[2];
                    if (quad.eigenValue0 < quad.eigenValue) {
                        bondDirection1[0] = quad.eigenVector0[0];
                        bondDirection1[1] = quad.eigenVector0[1];
                        bondDirection1[2] = quad.eigenVector0[2];

                        bondDirection2[0] = quad.eigenVector2[0];
                        bondDirection2[1] = quad.eigenVector2[1];
                        bondDirection2[2] = quad.eigenVector2[2];
                    }
                    else {
                        bondDirection1[0] = quad.eigenVector2[0];
                        bondDirection1[1] = quad.eigenVector2[1];
                        bondDirection1[2] = quad.eigenVector2[2];

                        bondDirection2[0] = quad.eigenVector0[0];
                        bondDirection2[1] = quad.eigenVector0[1];
                        bondDirection2[2] = quad.eigenVector0[2];
                    }
                }
                if (quad.eigenValue < quad.eigenValue0 && quad.eigenValue < quad.eigenValue1) {
                    e1 = quad.eigenVector2[0];
                    e2 = quad.eigenVector2[1];
                    e3 = quad.eigenVector2[2];
                    if (quad.eigenValue0 < quad.eigenValue1) {
                        bondDirection1[0] = quad.eigenVector0[0];
                        bondDirection1[1] = quad.eigenVector0[1];
                        bondDirection1[2] = quad.eigenVector0[2];

                        bondDirection2[0] = quad.eigenVector1[0];
                        bondDirection2[1] = quad.eigenVector1[1];
                        bondDirection2[2] = quad.eigenVector1[2];
                    }
                    else {
                        bondDirection1[0] = quad.eigenVector1[0];
                        bondDirection1[1] = quad.eigenVector1[1];
                        bondDirection1[2] = quad.eigenVector1[2];

                        bondDirection2[0] = quad.eigenVector0[0];
                        bondDirection2[1] = quad.eigenVector0[1];
                        bondDirection2[2] = quad.eigenVector0[2];
                    }
                }

                //float e0Size = 1.5 / (-log(quad.eigenValue0));
                //float e1Size = 1.5 / (-log(quad.eigenValue0));
                //float e2Size = 1.5 / (-log(quad.eigenValue0));
                if (ellipsoidVisible) {
                    float e0Size = 20.0*quad.eigenValue0;
                    float e1Size = 20.0*quad.eigenValue1;
                    float e2Size = 20.0*quad.eigenValue;
                    Vector3DFloat bondDir1 = Vector3DFloat(bondDirection1[0], bondDirection1[1], bondDirection1[2]);
                    bondDir1.Normalize();
                    Vector3DFloat bondDir2 = Vector3DFloat(bondDirection2[0], bondDirection2[1], bondDirection2[2]);
                    bondDir2.Normalize();
                    OpenGLUtils::SetColor(1.0, 0.0, 0.1, 1.0);
                    //cout << "e1 e2 e3 " << quad.eigenValue0 << " " << quad.eigenValue1 << " " << quad.eigenValue << " " << quad.eigenVector0[0] << " " << quad.eigenVector0[1] << " " << quad.eigenVector0[2] << " " << quad.eigenVector1[0] << " " << quad.eigenVector1[1] << " " << quad.eigenVector1[2] << " " << quad.eigenVector2[0] << " " << quad.eigenVector2[1] << " " << quad.eigenVector2[2] << " " << e1 << " " << e2 << " " << e3 << endl;
                    DrawEllipsoid(10, 10, ellipsoidScale * e0Size, ellipsoidScale * e1Size, ellipsoidScale * e2Size, (quad.p1[0] + quad.p2[0] + quad.p3[0] + quad.p4[0]) / 4.0, (quad.p1[1] + quad.p2[1] + quad.p3[1] + quad.p4[1]) / 4.0, (quad.p1[2] + quad.p2[2] + quad.p3[2] + quad.p4[2]) / 4.0, e1, e2, e3, bondDir1, bondDir2);
                }
            }
        }
        /**
        void CAlphaRenderer::DrawSurface(int quad1, int quad2, int quad3, int quad4, float nx, float ny, float nz) {
        unsigned long long atom1Hash = 1;
        float p1x = -1000.0, p1y = -1000.0, p1z = -1000.0, p2x = -1000.0, p2y = -1000.0,
        p2z = -1000.0, p3x = -1000.0, p3y = -1000.0, p3z = -1000.0, p4x = -1000.0, p4y = -1000.0, p4z = -1000.0;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

        if(it->second.GetSerial() == quad1) {
        atom1Hash = it->second.GetHashKey();
        p1x = it->second.GetPosition().X();
        p1y = it->second.GetPosition().Y();
        p1z = it->second.GetPosition().Z();
        //OpenGLUtils::SetColor(0.0, 1.0, 0.5, 1.0);
        //DrawSphere(it->second.GetPosition(), 0.25);
        //atomHashes.push_back(it->second.GetHashKey());
        }


        }
        unsigned long long atom2Hash = 1;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
        if(it->second.GetSerial() == quad2) {
        atom2Hash = it->second.GetHashKey();
        p2x = it->second.GetPosition().X();
        p2y = it->second.GetPosition().Y();
        p2z = it->second.GetPosition().Z();
        //OpenGLUtils::SetColor(1.0, 0.0, 0.1, 1.0);
        //DrawSphere(it->second.GetPosition(), 0.2);
        //atomHashes.push_back(it->second.GetHashKey());
        }
        }
        unsigned long long atom3Hash = 1;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
        if(it->second.GetSerial() == quad3) {
        atom3Hash = it->second.GetHashKey();
        p3x = it->second.GetPosition().X();
        p3y = it->second.GetPosition().Y();
        p3z = it->second.GetPosition().Z();
        //OpenGLUtils::SetColor(1.0, 0.85, 0.1, 1.0);
        //DrawSphere(it->second.GetPosition(), 0.15);
        //atomHashes.push_back(it->second.GetHashKey());
        }
        }
        unsigned long long atom4Hash = 1;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
        if(it->second.GetSerial() == quad4) {
        atom4Hash = it->second.GetHashKey();
        p4x = it->second.GetPosition().X();
        p4y = it->second.GetPosition().Y();
        p4z = it->second.GetPosition().Z();
        //OpenGLUtils::SetColor(0.6, 0.1, 0.85, 1.0);
        //DrawSphere(it->second.GetPosition(), 0.1);
        //atomHashes.push_back(it->second.GetHashKey());
        }
        }
        if (atom1Hash == 1 || atom2Hash == 1 || atom3Hash == 1 || atom4Hash == 1) {
        return;
        }
        else {
        glBegin(GL_QUADS);
        // top
        //glColor3f(1.0f, 0.0f, 0.0f);
        glNormal3f(nx, ny, nz);
        glVertex3f(p1x, p1y, p1z);
        glVertex3f(p2x, p2y, p2z);
        glVertex3f(p3x, p3y, p3z);
        glVertex3f(p4x, p4y, p4z);

        glEnd();


        //glDisable ( GL_COLOR_MATERIAL );


        }
        }
        **/
        static bool msg_sort_criteria0(boost::tuple<float, float, float, int, int> lhs,
            boost::tuple<float, float, float, int, int> rhs)
        {
            return boost::get<0>(lhs) < boost::get<0>(rhs);
        }

        static bool msg_sort_criteria1(boost::tuple<float, float, float, int, int> lhs,
            boost::tuple<float, float, float, int, int> rhs)
        {
            return boost::get<1>(lhs) < boost::get<1>(rhs);
        }

        static bool msg_sort_criteria2(boost::tuple<float, float, float, int, int> lhs,
            boost::tuple<float, float, float, int, int> rhs)
        {
            return boost::get<2>(lhs) < boost::get<2>(rhs);
        }

        static bool sortZValue(Vector3DFloat lhs, Vector3DFloat rhs)
        {
            return lhs.Z() < rhs.Z();
        }

        void CAlphaRenderer::sortMaxCurveDimensions() {
            //cout << "called " << endl;
            std::vector< boost::tuple<float, float, float, int, int> > eigenvalues;
            //std::vector<float> eigenV0;
            //std::vector<float> eigenV1;
            //std::vector<float> eigenV2;
            for (int i = 0; i < bonds.size(); i++) {
                //if(bonds[i].maxOn) {
                boost::tuple<float, float, float, int, int> eigenValueTuple(bonds[i].eigenValue0, bonds[i].eigenValue1, bonds[i].eigenvalue, i, 0);
                eigenvalues.push_back(eigenValueTuple);
                //}
                //eigenV0.push_back(bonds[i].eigenValue0);
                //eigenV1.push_back(bonds[i].eigenValue1);
                //eigenV2.push_back(bonds[i].eigenvalue);
            }
            for (int i = 0; i < quads.size(); i++) {
                boost::tuple<float, float, float, int, int> eigenValueTuple(quads[i].eigenValue0, quads[i].eigenValue1, quads[i].eigenValue, i, 1);
                eigenvalues.push_back(eigenValueTuple);
            }

            for (int i = 0; i < saddlePointEigenvalues.size() / 3; i++) {
                boost::tuple<float, float, float, int, int> eigenValueTuple(saddlePointEigenvalues[3 * i], saddlePointEigenvalues[3 * i + 1], saddlePointEigenvalues[3 * i + 2], i, 2);
                eigenvalues.push_back(eigenValueTuple);
            }

            for (int i = 0; i < maxPointEigenvalues.size() / 3; i++) {
                boost::tuple<float, float, float, int, int> eigenValueTuple(maxPointEigenvalues[3 * i], maxPointEigenvalues[3 * i + 1], maxPointEigenvalues[3 * i + 2], i, 3);
                eigenvalues.push_back(eigenValueTuple);
            }

            for (int i = 0; i < minPointEigenvalues.size() / 3; i++) {
                boost::tuple<float, float, float, int, int> eigenValueTuple(minPointEigenvalues[3 * i], minPointEigenvalues[3 * i + 1], minPointEigenvalues[3 * i + 2], i, 4);
                eigenvalues.push_back(eigenValueTuple);
            }

            std::sort(eigenvalues.begin(), eigenvalues.end(), msg_sort_criteria0);
            //for(int i = 0; i < bondEigenvalues.size(); i++) {
            //cout << "bond eigenvalues " << get<3>(bondEigenvalues[i]) << " " << i << endl;
            //}
            for (int i = 0; i < bonds.size(); i++) {
                for (int j = 0; j < eigenvalues.size(); j++) {
                    //cout << get<3>(bondEigenvalues[j]) << " " << i << " " << j << " " << bondEigenvalues.size() << endl;
                    if (get<3>(eigenvalues[j]) == i && get<4>(eigenvalues[j]) == 0) {
                        bonds[i].eigenValue0Rank = (float)(j) / (float)eigenvalues.size();
                    }
                }
            }



            std::sort(eigenvalues.begin(), eigenvalues.end(), msg_sort_criteria1);
            for (int i = 0; i < bonds.size(); i++) {
                for (int j = 0; j < eigenvalues.size(); j++) {
                    if (get<3>(eigenvalues[j]) == i && get<4>(eigenvalues[j]) == 0) {
                        bonds[i].eigenValue1Rank = (float)(j) / (float)eigenvalues.size();
                    }
                }
            }

            std::sort(eigenvalues.begin(), eigenvalues.end(), msg_sort_criteria2);
            for (int i = 0; i < bonds.size(); i++) {
                for (int j = 0; j < eigenvalues.size(); j++) {
                    if (get<3>(eigenvalues[j]) == i && get<4>(eigenvalues[j]) == 0) {
                        bonds[i].eigenValue2Rank = (float)(j) / (float)eigenvalues.size();
                    }
                }
            }

            //
            //std::sort(eigenV0.begin(), eigenV0.end());
            //eigenValues0 = eigenV0;

            //std::sort(eigenV1.begin(), eigenV1.end());
            //eigenValues1 = eigenV1;

            //std::sort(eigenV2.begin(), eigenV2.end());
            //eigenValues2 = eigenV2;

            /**
            for(int j = 0; j < bonds.size(); j++) {
            for(int i = 0; i < eigenV0.size(); i++) {
            if(bonds[j].eigenValue0 == eigenV0[i]) {
            bonds[j].eigenValue0Rank = (float)i / (float)bonds.size();
            }
            if(bonds[j].eigenValue1 == eigenV1[i]) {
            bonds[j].eigenValue1Rank = (float)i / (float)bonds.size();
            }
            if(bonds[j].eigenValue2 == eigenV1[i]) {
            bonds[j].eigenValue1Rank = (float)i / (float)bonds.size();
            }
            }
            }
            **/


        }


        int CAlphaRenderer::findHelices() {
            helices.clear();

            //AddMaxBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1) {

            std::vector< boost::tuple<unsigned long long, unsigned long long> > maxCurveSegments;
            //for(int i = 0; i < bonds.size(); i++) {
            //if(bonds[i].maxOn) {
            //boost::tuple<unsigned long long, unsigned long long> segmentTuple(bonds[i].GetAtom0Ix(), bonds[i].GetAtom1Ix());
            //maxCurveSegments.push_back(segmentTuple);
            //}
            //}
            for (int i = 0; i < bondAtm1s.size(); i++) {

                //AddMaxBond(maxCurveIndices1[i], maxCurveIndices2[i], std::to_string(bondAtm1s[i]), std::to_string(bondAtm2s[i]), maxCurveSaliencies[3*i], maxCurveSaliencies[3*i+1], maxCurveSaliencies[3*i+2], maxCurveIntensities[i], maxCurveEigenvalues[3*i], maxCurveEigenvalues[3*i+1], maxCurveEigenvalues[3*i+2]);
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[maxCurveIndices1[i] * 3];
                float py = vertexPositions[maxCurveIndices1[i] * 3 + 1];
                float pz = vertexPositions[maxCurveIndices1[i] * 3 + 2];
                //cout << "px py pz " << px << " " << py << " " << pz << endl;
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[bondAtm1s[i]] = newAtom;

                PDBAtom newAtom2 = PDBAtom();
                float px2 = vertexPositions[maxCurveIndices2[i] * 3];
                float py2 = vertexPositions[maxCurveIndices2[i] * 3 + 1];
                float pz2 = vertexPositions[maxCurveIndices2[i] * 3 + 2];
                newAtom2.SetPosition(Vector3DFloat(px2, py2, pz2));
                atoms[bondAtm2s[i]] = newAtom2;

                boost::tuple<unsigned long long, unsigned long long> segmentTuple(bondAtm1s[i], bondAtm2s[i]);
                maxCurveSegments.push_back(segmentTuple);
            }
            int i = 0;
            while (maxCurveSegments.size() > 0) {
                std::vector<unsigned long long> currentHelix;
                int j = 0;
                boost::tuple<unsigned long long, unsigned long long> currentSegment = maxCurveSegments[j];
                //int bondIndex = max(GetBondIndex(get<0>(currentSegment), get<1>(currentSegment) ), GetBondIndex( get<1>(currentSegment), get<0>(currentSegment) ));
                //bonds[bondIndex].helixIndex = i;
                currentHelix.push_back(get<0>(currentSegment));
                currentHelix.push_back(get<1>(currentSegment));
                maxCurveSegments.erase(maxCurveSegments.begin() + j);
                unsigned long long segment1 = get<0>(currentSegment);
                unsigned long long segment2 = get<1>(currentSegment);
                //cout << "segment " << std::to_string(segment1) << " " << std::to_string(segment2) << endl;
                while (segment2 != -1) {
                    bool found = false;
                    for (int k = 0; k < maxCurveSegments.size(); k++) {
                        boost::tuple<unsigned long long, unsigned long long> nextSegment = maxCurveSegments[k];
                        if (get<0>(nextSegment) == segment2) {
                            segment2 = get<1>(nextSegment);
                            currentHelix.push_back(segment2);
                            segment1 = segment2;
                            found = true;
                            maxCurveSegments.erase(maxCurveSegments.begin() + k);
                            break;
                        }

                        if (get<1>(nextSegment) == segment2) {
                            segment2 = get<0>(nextSegment);
                            currentHelix.push_back(segment2);
                            segment1 = segment2;
                            found = true;
                            maxCurveSegments.erase(maxCurveSegments.begin() + k);
                            break;
                        }
                    }
                    if (found == false) {
                        segment2 = -1;
                    }
                }
                //cout << "max curve segments size " << maxCurveSegments.size() << endl;
                if (currentHelix.size() > 4) {
                    fittedHelix newHelix;

                    for (int i = 0; i < currentHelix.size(); i++) {
                        int bondIndex = max(GetBondIndex(currentHelix[i], currentHelix[i + 1]), GetBondIndex(currentHelix[i + 1], currentHelix[i]));
                        //bonds[bondIndex].helixIndex = 1;
                        //bonds[bondIndex].helixIndices = currentHelix;
                        newHelix.fittedHelixBondIndices.push_back(bondIndex);

                    }
                    fittedHelices.push_back(newHelix);
                    //boost::tuple<unsigned long long, unsigned long long> currentSegment = maxCurveSegments[j];
                    //int bondIndex = max(GetBondIndex(get<0>(currentSegment), get<1>(currentSegment) ), GetBondIndex( get<1>(currentSegment), get<0>(currentSegment) ));


                    /**
                    std::vector<Vector3DFloat> firstOrderDifferences;
                    for(int i = 0; i < currentHelix.size()-1; i++) {
                    Vector3DFloat firstOrderDiff = atoms[currentHelix[i+1]].GetPosition() - atoms[currentHelix[i]].GetPosition();
                    firstOrderDifferences.push_back(firstOrderDiff);
                    }
                    for(int i = 0; i < firstOrderDifferences.size()-1; i++) {
                    Vector3DFloat v1 = firstOrderDifferences[i];
                    Vector3DFloat v2 = firstOrderDifferences[i+1];
                    float bondAngle = findAngle(v1.X(), v1.Y(), v1.Z(), v2.X(), v2.Y(), v2.Z());
                    float multipliedBondAngle = bondAngle * 57.2957795;
                    cout << "bond angle " << bondAngle << " " << multipliedBondAngle << endl;
                    }
                    **/
                    /**
                    std::vector<Vector3DFloat> secondOrderDifferences;
                    for(int i = 0; i < firstOrderDifferences.size(); i++) {
                    Vector3DFloat secondOrderDiff = firstOrderDifferences[i+1] - firstOrderDifferences[i];

                    secondOrderDifferences.push_back(secondOrderDiff);
                    }
                    **/
                    /**
                    std::vector<Vector3DFloat> atomPositions;
                    for(int i = 0; i < currentHelix.size(); i++) {
                    atomPositions.push_back(atoms[currentHelix[i]].GetPosition());
                    }

                    //fitting central axis
                    Vector3DFloat mean = findMean(atomPositions);
                    Vector3DFloat qtr1 = atoms[currentHelix[0]].GetPosition() + mean;
                    qtr1 = Vector3DFloat(qtr1.X()/2.0, qtr1.Y()/2.0, qtr1.Z()/2.0);

                    Vector3DFloat qtr2 = atoms[currentHelix[currentHelix.size()-1]].GetPosition() + mean;
                    qtr2 = Vector3DFloat(qtr2.X()/2.0, qtr2.Y()/2.0, qtr2.Z()/2.0);


                    Vector3DFloat axisDirection = qtr2 - qtr1;
                    axisDirection.Normalize();

                    float axisLength = (atoms[currentHelix[currentHelix.size()-1]].GetPosition() - atoms[currentHelix[0]].GetPosition()).Length();

                    Vector3DFloat axisPt1 = mean + (axisDirection * (axisLength / 2.0) );
                    Vector3DFloat axisPt2 = mean - (axisDirection * (axisLength / 2.0) );

                    std::vector<Vector3DFloat> meanDifferences;

                    for(int i = 0; i < currentHelix.size(); i++) {
                    Vector3DFloat currentDifference = atoms[currentHelix[i]].GetPosition() - mean;
                    meanDifferences.push_back(currentDifference);
                    }
                    for(int i = 0; i < currentHelix.size(); i++) {
                    float dist = findDistance(axisPt1, axisPt2, atoms[currentHelix[i]].GetPosition());
                    cout << "dist " << dist << endl;
                    }
                    **/

                    helices.push_back(currentHelix);
                }
                currentHelix.clear();
            }
            /**
            for(int i = 0; i < helices.size(); i++) {
            std::vector<unsigned long long> currentH = helices[i];
            for(int j = 0; j < currentH.size(); j++) {
            cout << "helix number " << i << " helixAtm " << currentH[j] << endl;
            }
            }
            **/
            for (int i = 0; i < helices.size(); i++) {
                std::vector<unsigned long long> currentH = helices[i];
                std::vector<float> helixPositions;
                for (int j = 0; j < currentH.size(); j++) {
                    Vector3DFloat currentA = atoms[currentH[j]].GetPosition();
                    helixPositions.push_back(currentA.X());
                    helixPositions.push_back(currentA.Y());
                    helixPositions.push_back(currentA.Z());
                }
                maxCurveSepPts.push_back(helixPositions);
            }
            return (int)helices.size();
        }

        //Divide into buckets


        void CAlphaRenderer::findMaxCurves() {
            helices.clear();

            //AddMaxBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1) {

            std::vector< boost::tuple<unsigned long long, unsigned long long> > maxCurveSegments;
            //for(int i = 0; i < bonds.size(); i++) {
            //if(bonds[i].maxOn) {
            //boost::tuple<unsigned long long, unsigned long long> segmentTuple(bonds[i].GetAtom0Ix(), bonds[i].GetAtom1Ix());
            //maxCurveSegments.push_back(segmentTuple);
            //}
            //}
            for (int i = 0; i < bonds.size(); i++) {

                //AddMaxBond(maxCurveIndices1[i], maxCurveIndices2[i], std::to_string(bondAtm1s[i]), std::to_string(bondAtm2s[i]), maxCurveSaliencies[3*i], maxCurveSaliencies[3*i+1], maxCurveSaliencies[3*i+2], maxCurveIntensities[i], maxCurveEigenvalues[3*i], maxCurveEigenvalues[3*i+1], maxCurveEigenvalues[3*i+2]);

                bool display = false;
                bool displayMax = false;

                if (bonds[i].maxOn) {
                    /**
                    displayMax = true;

                    display = true;
                    float pointR = pointRatio;
                    float curveR = curveRatio;
                    float surfaceR = surfaceRatio;
                    float minG = minGeo;
                    float maxG = maxGeo;
                    float eigenV = eigenValue;
                    if (saliencyCheck) {
                    float *saliencies = bonds[i].saliencies;
                    if ( ( curveR * saliencies[1] < surfaceR * saliencies[0] ) || ( curveR * saliencies[1] < pointR * saliencies[2] ) ){
                    display = false;
                    displayMax = false;
                    }
                    }
                    if(intensityCheck) {
                    float maxI = maxDataIntensity;
                    float minI = minDataIntensity;
                    float localI = bonds[i].intensity;
                    if (( ( localI - minI ) / ( maxI - minI ) < maxG ) ) {
                    display = false;
                    displayMax = false;
                    }
                    }
                    if(eigenvalueCheck) {
                    float maxE = maxDataEigenvalue;
                    float minE = minDataEigenvalue;
                    float localE = bonds[i].eigenvalue;
                    if ( ( localE - minE ) / ( maxE - minE ) < eigenV ){
                    display = false;
                    displayMax = false;
                    }
                    }

                    }
                    else {
                    display = false;
                    }
                    **/
                    //if(display) {
                    boost::tuple<unsigned long long, unsigned long long> segmentTuple(bondAtm1s[i], bondAtm2s[i]);
                    maxCurveSegments.push_back(segmentTuple);
                    //}


                }
            }
            int i = 0;
            while (maxCurveSegments.size() > 0) {
                std::vector<unsigned long long> currentHelix;
                int j = 0;
                boost::tuple<unsigned long long, unsigned long long> currentSegment = maxCurveSegments[j];
                currentHelix.push_back(get<0>(currentSegment));
                currentHelix.push_back(get<1>(currentSegment));
                maxCurveSegments.erase(maxCurveSegments.begin() + j);
                unsigned long long segment1 = get<0>(currentSegment);
                unsigned long long segment2 = get<1>(currentSegment);
                //cout << "segment " << std::to_string(segment1) << " " << std::to_string(segment2) << endl;
                while (segment2 != -1) {
                    bool found = false;
                    for (int k = 0; k < maxCurveSegments.size(); k++) {
                        boost::tuple<unsigned long long, unsigned long long> nextSegment = maxCurveSegments[k];
                        if (get<0>(nextSegment) == segment2) {
                            segment2 = get<1>(nextSegment);
                            currentHelix.push_back(segment2);
                            segment1 = segment2;
                            found = true;
                            maxCurveSegments.erase(maxCurveSegments.begin() + k);
                            break;
                        }
                    }
                    if (found == false) {
                        segment2 = -1;
                    }
                }
                //cout << "max curve segments size " << maxCurveSegments.size() << endl;
                if (currentHelix.size() > 3) {

                    helices.push_back(currentHelix);
                }
                currentHelix.clear();
            }
            /**
            for(int i = 0; i < helices.size(); i++) {
            std::vector<unsigned long long> currentH = helices[i];
            for(int j = 0; j < currentH.size(); j++) {
            cout << "helix number " << i << " helixAtm " << currentH[j] << endl;
            }
            }
            **/
            ofstream myfile;
            myfile.open("maxCurves.pdb");
            int l = 0;
            for (int j = 0; j < helices.size(); j++) {

                string currentIndex = std::to_string(j);

                std::vector<unsigned long long> currentMaxCurve = helices[j];

                for (int k = 0; k < currentMaxCurve.size(); k++) {


                    string vIndex = std::to_string(l);
                    padTo(vIndex, 5);

                    float xPos0 = atoms[currentMaxCurve[k]].GetPosition().X();
                    float yPos0 = atoms[currentMaxCurve[k]].GetPosition().Y();
                    float zPos0 = atoms[currentMaxCurve[k]].GetPosition().Z();

                    string posx0 = std::to_string(xPos0);
                    string posy0 = std::to_string(yPos0);
                    string posz0 = std::to_string(zPos0);

                    myfile << "ATOM  " << vIndex << "  CA  ALA " << vIndex << "     " << posx0.substr(0, 7) << " " << posy0.substr(0, 7) << " " << posz0.substr(0, 7) << "  1.00  1.00      S_00  0 " << endl;
                    l++;
                }

            }
            myfile.close();
        }

        bool CAlphaRenderer::CheckMaxPointHide(float saliencies[3], float intensity, float eigenvalue) {
            bool hide = false;


            if (saliencyCheck) {
                //float *saliencies = point.relativeSaliencies;
                //cout << "saliencies " << saliencies[2] << " " << saliencies[1] << " " << saliencies[0] << endl;
                //cout << pointRatio * saliencies[0] << " " << surfaceRatio * saliencies[2] << " " << curveRatio * saliencies[1] << endl;
                if ((pointRatio * saliencies[2] < surfaceRatio * saliencies[0]) ||
                    (pointRatio * saliencies[2] < curveRatio * saliencies[1])) hide = true;
            }
            if (intensityCheck){
                float maxI = maxDataIntensity;
                float minI = minDataIntensity;
                float localI = intensity;
                //cout << "localI " << localI << "minI " << minI << "maxI " << maxI << endl;
                //cout << ( localI - minI ) / ( maxI - minI ) << " " << maxGeo << endl;
                if ((localI - minI) / (maxI - minI) < maxGeo) hide = true;

            }
            if (eigenvalueCheck) {
                float maxE = maxDataEigenvalue;
                float minE = minDataEigenvalue;
                float localE = eigenvalue;
                //cout << "localE " << localE << " maxE " << maxE << " minE " << minE << " eigenValue " << endl;
                //cout << ( localE - minE ) / ( maxE - minE ) << endl;
                if ((localE - minE) / (maxE - minE) < eigenValue) hide = true;

            }
            //cout << "hide " << hide << endl;
            return hide;
        }

        bool CAlphaRenderer::CheckMinPointHide(float saliencies[3], float intensity, float eigenvalue) {
            bool hide = false;


            if (saliencyCheck) {
                if ((pointRatio * saliencies[2] < surfaceRatio * saliencies[0]) ||
                    (pointRatio * saliencies[2] < curveRatio * saliencies[1])) hide = true;
            }
            if (intensityCheck){
                float maxI = maxDataIntensity;
                float minI = minDataIntensity;
                float localI = intensity;
                //cout << "localI " << localI << "minI " << minI << "maxI " << maxI << endl;
                //cout << ( localI - minI ) / ( maxI - minI ) << " " << maxGeo << endl;
                if (((localI - minI) / (maxI - minI) > minGeo)) hide = true;

            }
            if (eigenvalueCheck) {
                float maxE = maxDataEigenvalue;
                float minE = minDataEigenvalue;
                float localE = eigenvalue;
                //cout << "localE " << localE << " maxE " << maxE << " minE " << minE << " eigenValue " << endl;
                //cout << ( localE - minE ) / ( maxE - minE ) << endl;
                if ((localE - minE) / (maxE - minE) < eigenValue) hide = true;

            }
            //cout << "hide " << hide << endl;
            return hide;
        }

        bool CAlphaRenderer::CheckSaddlePointHide(float saliencies[3], float intensity, float eigenvalue) {
            bool hide = false;


            if (saliencyCheck) {
                if ((pointRatio * saliencies[2] < surfaceRatio * saliencies[0]) ||
                    (pointRatio * saliencies[2] < curveRatio * saliencies[1])) hide = true;
            }
            if (eigenvalueCheck) {
                float maxE = maxDataEigenvalue;
                float minE = minDataEigenvalue;
                float localE = eigenvalue;
                //cout << "localE " << localE << " maxE " << maxE << " minE " << minE << " eigenValue " << endl;
                //cout << ( localE - minE ) / ( maxE - minE ) << endl;
                if ((localE - minE) / (maxE - minE) < eigenValue) hide = true;

            }
            //cout << "hide " << hide << endl;
            return hide;
        }

        //float max(float a, float b) 
        //{
        //	if(a > b) 
        //    {
        //		return a;
        //	}
        //
        //	return b;
        //}

        /**
        void CAlphaRenderer::fitHelix(std::vector<unsigned long long> currentStrand) {
        for(int i = 0; i < currentStrand.size()-1; i++) {
        //currentStrand[i]
        }
        }


        void CAlphaRenderer::fitHelices() {
        for(int i = 0; i < helices.size(); i++) {
        std::vector<unsigned long long> currentHelix = helices[i];
        fitHelix(currentHelix);
        }
        }
        **/

        void CAlphaRenderer::interpolateHelixPoints() {
            //ofstream myfile;
            //myfile.open("origPoints.txt");
            for (int i = 0; i < helices.size(); i++) {
                //	std::vector<float> currentSeg = maxCurveSepPts[i];
                //cout << "curveminsx " << curveMinsX.size() << " segments " << maxCurveSepPts.size() << endl;
                std::vector<unsigned long long> currentHelix = helices[i];
                //DrawSphere(atoms[currentHelix[0]].GetPosition(), 0.3);
                //DrawSphere(atoms[currentHelix[currentHelix.size()-1]].GetPosition(), 0.3);
                float currentMinX = curveMinsX[i];
                float currentMinY = curveMinsY[i];
                float currentMaxX = curveMaxsX[i];
                float currentMaxY = curveMaxsY[i];
                float currentMinZ = curveMinsZ[i];
                float currentMaxZ = curveMaxsZ[i];
                float currentMeanX = curveMeanXs[i];
                float currentMeanY = curveMeanYs[i];
                float currentMeanZ = curveMeanZs[i];


                Vector3DFloat firstAtomPos = atoms[currentHelix[0]].GetPosition();
                Vector3DFloat lastAtomPos = atoms[currentHelix[currentHelix.size() - 1]].GetPosition();

                float maxDistance = Vector3DFloat(firstAtomPos.X() - lastAtomPos.X(), firstAtomPos.Y() - lastAtomPos.Y(), firstAtomPos.Z() - lastAtomPos.Z()).Length();
                maxDistance = maxDistance / 2.0;
                Vector3DFloat minPt = Vector3DFloat(currentMeanX, currentMeanY, currentMeanZ) - Vector3DFloat(maxDistance * ordLines[i].X(), maxDistance * ordLines[i].Y(), maxDistance * ordLines[i].Z());
                Vector3DFloat maxPt = Vector3DFloat(currentMeanX, currentMeanY, currentMeanZ) + Vector3DFloat(maxDistance * ordLines[i].X(), maxDistance * ordLines[i].Y(), maxDistance * ordLines[i].Z());
                //if(maxPt.X() < minPt.X() && maxPt.Y() < minPt.Y() && maxPt.Z() < minPt.Z() ) {
                //Vector3DFloat tempPt = minPt;
                //minPt = maxPt;
                //maxPt = tempPt;
                //}

                float currentA = coefsA[i];
                float currentB = coefsB[i];
                float currentC = coefsC[i];
                //float currentD = coefsD[i];
                //float predictedZMin = currentA * currentMinX + currentMinY * currentB + currentC + currentD;
                //float predictedZMax = currentA * currentMaxX + currentMaxY * currentB + currentC + currentD;
                float predictedZMin = currentA + currentMinX  * currentB + currentC * currentMinY;
                float predictedZMax = currentA + currentMaxX  * currentB + currentC * currentMaxY;

                //float predictedZMin = currentA + currentSeg[0]*currentB + currentC*currentSeg[1];
                //float predictedZMax = currentA + currentSeg[currentSeg.size()-3]*currentB + currentSeg[currentSeg.size()-2]*currentC;

                //Vector3DFloat minPt = Vector3DFloat(currentMinX, currentMinY, predictedZMin);
                //Vector3DFloat maxPt = Vector3DFloat(currentMaxX, currentMaxY, predictedZMax);

                //Vector3DFloat minPt = Vector3DFloat(currentSeg[0], currentSeg[1], predictedZMin);
                //Vector3DFloat maxPt = Vector3DFloat(currentSeg[currentSeg.size()-3], currentSeg[currentSeg.size()-2], predictedZMin);
                //cout << "minpt " << minPt.X() << " " << minPt.Y() << " " << minPt.Z() << endl;
                //cout << "maxpt " << maxPt.X() << " " << maxPt.Y() << " " << maxPt.Z() << endl;

                //DrawCylinder(minPt, maxPt, 0.2);
                Vector3DFloat segmentDirection = maxPt - minPt;
                //orthogonal distances

                segmentDirection.Normalize();
                float currentDirection[3];
                currentDirection[0] = segmentDirection.X();
                currentDirection[1] = segmentDirection.Y();
                currentDirection[2] = segmentDirection.Z();

                float meanOrthoDistance = 0.0;
                for (int i = 0; i < currentHelix.size(); i++) {
                    Vector3DFloat PtoP = atoms[currentHelix[i]].GetPosition() - minPt;
                    float pTop[3];
                    pTop[0] = PtoP.X();
                    pTop[1] = PtoP.Y();
                    pTop[2] = PtoP.Z();

                    float * crossP = crossProduct(pTop, currentDirection);
                    Vector3DFloat proj = Vector3DFloat(crossP[0], crossP[1], crossP[2]);
                    float dist = proj.Length() / segmentDirection.Length();
                    meanOrthoDistance += dist;

                }
                meanOrthoDistance /= currentHelix.size();
                //cout << "mean ortho dist " << meanOrthoDistance << endl;
                orthoDistances.push_back(fabs(meanOrthoDistance));




                float zaxis[3];
                zaxis[0] = 0.0;
                zaxis[1] = 0.0;
                zaxis[2] = 1.0;


                helixMaxs.push_back(maxPt);
                helixMins.push_back(minPt);
                //DrawCylinder(minPt, maxPt, 0.03);
                Matrix4 rotation = alignWithAxis(currentDirection, zaxis, minPt, maxPt);
                //cout << "origDir " << segmentDirection.X() << " " << segmentDirection.Y() << " " << segmentDirection.Z() << endl;
                Vector3 newDir = rotation * Vector3(segmentDirection.X(), segmentDirection.Y(), segmentDirection.Z());
                //cout << "newDir " << newDir[0] << " "  << newDir[1] << " " << newDir[2] << endl;
                std::vector<Vector3DFloat> helixNormalizedPos;

                for (int j = 0; j < currentHelix.size(); j++) {
                    Vector3DFloat currentPos = atoms[currentHelix[j]].GetPosition();
                    //myfile << "helix Index " << i << " atom num " << j << " " << currentPos.X() << " " << currentPos.Y() << " " << currentPos.Z() << endl;
                    //float interpolatedX = (currentPos.X() - minPt.X() ) / (maxPt.X() - minPt.X() );
                    //float interpolatedY = (currentPos.Y() - minPt.Y() ) / (maxPt.Y() - minPt.Y() );
                    //float interpolatedZ = (currentPos.Z() - minPt.Z() ) / (maxPt.Z() - minPt.Z() );
                    float interpolatedX = (currentPos.X() - minPt.X());
                    float interpolatedY = (currentPos.Y() - minPt.Y());
                    float interpolatedZ = (currentPos.Z() - minPt.Z());

                    Vector3DFloat interpPos = Vector3DFloat(interpolatedX, interpolatedY, interpolatedZ);
                    //Vector3DFloat interpPos = currentPos - Vector3DFloat(currentMinX, currentMinY, currentMinZ);
                    //Vector3DFloat maxMinusMin = Vector3DFloat(currentMaxX, currentMaxY, currentMaxZ) - Vector3DFloat(currentMinX, currentMinY, currentMinZ);
                    //interpPos = Vector3DFloat(interpPos.X()  / maxMinusMin.X(), interpPos.Y()  / maxMinusMin.Y(), interpPos.Z()  / maxMinusMin.Z());
                    //interpPos.Normalize();
                    //interpPos = Vector3DFloat(interpPos.X() / segmentDirection.X(), interpPos.Y() / segmentDirection.Y(), interpPos.Z() / segmentDirection.Z());
                    //interpPos.Normalize();
                    Vector3 inPos = Vector3(interpPos.X(), interpPos.Y(), interpPos.Z());
                    //cout << "originalPos " << currentPos.X() << " " << currentPos.Y() << " " << currentPos.Z() << endl;
                    Vector3 adjustedCurrentPos = Vector3(currentPos.X() - minPt.X(), currentPos.Y() - minPt.Y(), currentPos.Z() - minPt.Z());
                    //Vector3 adjustedCurrentPos = Vector3(currentPos.X(),  currentPos.Y(), currentPos.Z());

                    //cout << "currentPos " << adjustedCurrentPos[0] << " " << adjustedCurrentPos[1] << " " << adjustedCurrentPos[2] << endl;
                    //Vector3 rotatedPos = rotation * inPos;
                    Vector3 rotatedPos = rotation * adjustedCurrentPos;
                    //cout << "currentPos after" << rotatedPos[0] << " " << rotatedPos[1] << " " << rotatedPos[2] << endl;


                    //rotatedPos.normalize();

                    helixNormalizedPos.push_back(Vector3DFloat(rotatedPos[0], rotatedPos[1], rotatedPos[2]));
                    //cout << "rotated pos " << rotatedPos[0] << " " << rotatedPos[1] << " " << rotatedPos[2] << endl;
                    //float inPos[3];

                    //inPos[0] = interpPos.X();
                    //inPos[1] = interpPos.Y();
                    //inPos[2] = interpPos.Z();
                    //alignWithAxis(inPos, zaxis, minPt, maxPt);
                }
                normalizedHelices.push_back(helixNormalizedPos);

                //rotationMatrices.push_back(rotation.inverse());

            }
            normalizedHelicesMinsAndMaxes();
            for (int i = 0; i < normalizedHelices.size(); i++) {
                //std::vector< Vector3DFloat > currentStrand = normalizedHelices[i];
                //foundHelixIndices.push_back(strandIndex);
                //helixScoreThresholds.push_back(currentZDist / thisStrand.size());
                //strandCycles.push_back(helixCycles);

                scoreHelix(i);
                fittedHelices[i].fittedHelixBondPositions = foundHelices[i];
                fittedHelices[i].orthoDist = orthoDistances[i];
                fittedHelices[i].coils = strandCycles[i];
                fittedHelices[i].coilHeight = helixTurnHeights[i];
                fittedHelices[i].helixError = helixScoreThresholds[i];
            }
            //myfile.close();	
        }

        void CAlphaRenderer::normalizedHelicesMinsAndMaxes() {
            for (int i = 0; i < normalizedHelices.size(); i++) {
                std::vector< Vector3DFloat > currentStrand = normalizedHelices[i];
                float minX = 1000.0;
                float minY = 1000.0;
                float minZ = 1000.0;
                float maxX = -1000.0;
                float maxY = -1000.0;
                float maxZ = -1000.0;
                for (int j = 0; j < currentStrand.size(); j++) {
                    if (currentStrand[j].X() > maxX) {
                        maxX = currentStrand[j].X();
                    }
                    if (currentStrand[j].X() < minX) {
                        minX = currentStrand[j].X();
                    }

                    if (currentStrand[j].Y() > maxY) {
                        maxY = currentStrand[j].Y();
                    }
                    if (currentStrand[j].Y() < minY) {
                        minY = currentStrand[j].Y();
                    }

                    if (currentStrand[j].Z() > maxZ) {
                        maxZ = currentStrand[j].Z();
                    }
                    if (currentStrand[j].Z() < minZ) {
                        minZ = currentStrand[j].Z();
                    }
                }
                normalizedMinXs.push_back(minX);
                normalizedMinYs.push_back(minY);
                normalizedMinZs.push_back(minZ);
                normalizedMaxXs.push_back(maxX);
                normalizedMaxYs.push_back(maxY);
                normalizedMaxZs.push_back(maxZ);
            }
        }

        void CAlphaRenderer::addOrd(float v1, float v2, float v3) {
            ordLines.push_back(Vector3DFloat(v1, v2, v3));
        }

        void CAlphaRenderer::scoreHelix(int strandIndex) {
            //ofstream ofs;
            //ofs.open("helixPos.txt", ofstream::app);
            std::vector<Vector3DFloat> thisStrand = normalizedHelices[strandIndex];
            float avgA = 0.0;
            float avgB = 0.0;
            float avgW = 0.0;
            float avgSigma = 0.0;
            //std::vector<Vector3DFloat> thisStrand = currentStrand;
            std::sort(thisStrand.begin(), thisStrand.end(), sortZValue);

            float strandSize = thisStrand.size();
            float helixCycles = 0;
            float avgRadiusX = 0.0;
            float avgRadiusY = 0.0;
            float semiCircleTracker = 0.0;
            int semiCircleCounts = 0;
            Vector3DFloat strand0 = thisStrand[0];
            for (int i = 0; i < thisStrand.size() - 1; i++) {
                //DrawSphere(Vector3DFloat(thisStrand[i].X() * 100.0, thisStrand[i].Y() * 100.0, thisStrand[i].Z() * 100.0), 0.05);
                //cout << "strand coord " << thisStrand[i].X() << " " << thisStrand[i].Y() << " " << thisStrand[i].Z() << endl;
                /**
                float x0 = thisStrand[i].X() / 10.0;
                float x1 = thisStrand[i+1].X() / 10.0;
                float y0 = thisStrand[i].Y() / 10.0;
                float y1 = thisStrand[i+1].Y() / 10.0;
                float delta = ((x0 * x0) * (y1 * y1)) - ((x1 * x1) * (y0 * y0));
                if (delta == 0) {
                strandSize -= 1;
                continue;
                }
                float aSquared = abs( (((y1 * y1) - (y0 * y0))) / delta );
                float bSquared = abs( (((x0 * x0) - (x1 * x1))) / delta );
                if (aSquared < 0 || bSquared < 0) {
                strandSize -= 1;
                continue;
                }
                float a = sqrt( aSquared );
                float b = sqrt( bSquared );
                float angle0 = atan( (a * y0) / (b * x0) );
                float angle1 = atan( (a * y1) / (b * x1) );
                float z0 = thisStrand[i].Z() / 10.0;
                float z1 = thisStrand[i+1].Z() / 10.0;
                **/
                Vector3DFloat strand1 = Vector3DFloat(thisStrand[i].X(), thisStrand[i].Y(), 0.0);
                strand1.Normalize();
                float str1[3];
                str1[0] = strand1.X();
                str1[1] = strand1.Y();
                str1[2] = 0.0;
                Vector3DFloat strand2 = Vector3DFloat(thisStrand[i + 1].X(), thisStrand[i + 1].Y(), 0.0);
                strand2.Normalize();
                float str2[3];
                str2[0] = strand2.X();
                str2[1] = strand2.Y();
                str2[2] = 0.0;
                float rot_angle = acos(getDotP(str1, str2));
                //cout << "rot_angle " << rot_angle << endl;
                helixCycles += (rot_angle * (180.0 / M_PI));
                semiCircleTracker += (rot_angle * (180.0 / M_PI));
                if (semiCircleTracker > 180.0) {
                    semiCircleCounts++;
                    //float remaining = semiCircleTracker - 180.0;
                    float currentXRadius = fabs(strand2.X() - strand0.X());
                    float proportion = sin((360.0 - semiCircleTracker) / 2.0 / (180.0 / M_PI));
                    currentXRadius = (currentXRadius / proportion) / 2.0;

                    float currentYRadius = fabs(strand2.Y() - strand0.Y());
                    //avgRadiusX += currentXRadius;
                    //avgRadiusX += (remaining / 180.0) * avgRadiusX;
                    //avgRadiusY += currentYRadius;
                    //avgRadiusY += (remaining / 180.0) * avgRadiusY;
                    currentYRadius = (currentYRadius / proportion) / 2.0;
                    avgRadiusX += currentXRadius;
                    avgRadiusY += currentYRadius;
                    strand0 = strand2;
                    semiCircleTracker -= 180.0;
                    //cout << "semi circle tracker inc " << avgRadiusX << " " << avgRadiusY << " " << remaining << " " << strand0.X() << " " << strand2.X() << " " << semiCircleCounts << endl;
                }

                //float w = (angle1 - angle0) / (z1 - z0);
                //float sigma = ( (angle0 * z1) - (angle1 * z0) ) / (z1 - z0);
                //cout << "debug " << x0 << " " << x1 << " " << y0 << " " << y1 << " " << z0 << " " << z1 << " " << a << " " << b << " " << angle0 << " " << angle1 << " " << z0 << " " << z1 << " " << w << " " << sigma << endl; 
                //avgA += a;
                //avgB += b;
                //avgW += w;
                //avgSigma += sigma;
            }
            if (semiCircleCounts == 0) {
                //float proportionIncrease = 180.0 / semiCircleTracker;
                avgRadiusX = abs(thisStrand[thisStrand.size() - 1].X() - strand0.X());
                avgRadiusY = abs(thisStrand[thisStrand.size() - 1].Y() - strand0.Y());
                float proportion = sin((360.0 - semiCircleTracker) / 2.0 / (180.0 / M_PI));
                avgRadiusX = (avgRadiusX / proportion) / 2.0;
                avgRadiusY = (avgRadiusY / proportion) / 2.0;


                semiCircleCounts = 1;
            }
            avgRadiusX /= semiCircleCounts;
            avgRadiusY /= semiCircleCounts;
            //cout << "radii " << avgRadiusX << " " << avgRadiusY << " " << semiCircleCounts << endl;
            helixCycles /= 360.0;

            cout << "helix cycles " << helixCycles << endl;
            //if(helixCycles < 1.5) {
            //cout << "cycles less than 1.5 " << endl;
            //return;
            //}
            //avgA = avgA / (strandSize - 1);
            //avgB = avgB / (strandSize - 1);
            //avgW = avgW / (strandSize - 1) ;
            //avgSigma = avgSigma / (strandSize - 1);
            double xRange = normalizedMaxXs[strandIndex] - normalizedMinXs[strandIndex];
            double yRange = normalizedMaxYs[strandIndex] - normalizedMinYs[strandIndex];
            double zRange = normalizedMaxZs[strandIndex] - normalizedMinZs[strandIndex];
            //cout << "zrange " << zRange << endl;
            //OpenGLUtils::SetColor(0.5, 0.5, 1.0, 1.0);
            //DrawCylinder(Vector3DFloat(0.0, 0.0, 0.0), Vector3DFloat(0.0, 0.0, 10.0), 0.05);

            avgA = avgRadiusX;
            avgB = avgRadiusY;
            //avgA = xRange * 0.48;
            //avgB = yRange * 0.48;
            float turnHeight = zRange / helixCycles;
            cout << "ratios " << avgA << " " << avgB / zRange << endl;
            helixTurnHeights.push_back(turnHeight);
            cout << "turnHeight " << turnHeight << endl;
            float c = turnHeight / (2.0 * M_PI);
            avgW = 1.0;
            avgSigma = 0.0;
            //for(float i = 0.0; i < 2.0 * M_PI; i = i + 0.01) {
            float minErrorRotation = 0.0;
            float minError = 10000000.0;

            float minErrorRotationInverse = 0.0;
            bool inverseOrientation = false;

            //float zTrans = 0.0;
            //float zSearchRange = turnHeight / 2.0;
            //float currentZDist = 1000.0;

            for (float d = 0.0; d < 2.0 * M_PI; d = d + 0.01) {
                float currentRotationSum = 0.0;
                float currentRotationZInverseSum = 0.0;
                //float thisZDist = 10000.0;
                //for(float i = -zSearchRange; i < zSearchRange; i = i + 0.01) {

                for (int j = 0; j < thisStrand.size(); j++) {
                    float zPos = (thisStrand[j].Z());
                    //float finalDist = 100000.0;
                    //for(float h = -1.0 * M_PI; h < 1.0 * M_PI; h = h + 0.01) {

                    float currentT = zPos / c;

                    float predictedX = avgA * cos(avgW * currentT + avgSigma);
                    float predictedY = avgB * sin(avgW * currentT + avgSigma);
                    Vector3 origPrediction = Vector3(predictedX, predictedY, zPos);

                    Vector3 rotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), d) * origPrediction;

                    predictedX = rotatedPrediction[0];
                    predictedY = rotatedPrediction[1];
                    //zPos = rotatedPrediction[2];

                    float xPos = (thisStrand[j].X());
                    float yPos = (thisStrand[j].Y());

                    float dist = (xPos - predictedX) * (xPos - predictedX) + (yPos - predictedY) * (yPos - predictedY) + (zPos - rotatedPrediction[2]) * (zPos - rotatedPrediction[2]);

                    Vector3 inverseRotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), d + M_PI) * origPrediction;
                    float inversePredictedX = inverseRotatedPrediction[0];
                    float inversePredictedY = inverseRotatedPrediction[1];
                    float inverseZPos = inverseRotatedPrediction[2];

                    float inverseDist = (xPos - inversePredictedX) * (xPos - inversePredictedX) + (yPos - inversePredictedY) * (yPos - inversePredictedY) + (zPos - inverseZPos) * (zPos - inverseZPos);

                    //if (dist < finalDist) {
                    //finalDist = dist;
                    //}
                    //}
                    currentRotationZInverseSum += inverseDist;
                    currentRotationSum += dist;
                }
                if (currentRotationZInverseSum < currentRotationSum) {
                    currentRotationSum = currentRotationZInverseSum;
                    inverseOrientation = true;
                }
                else {
                    inverseOrientation = false;
                }

                if (currentRotationSum < minError) {
                    minError = currentRotationSum;
                    minErrorRotation = d;
                    //zTrans = i;
                }
                //}


            }

            float zTrans = 0.0;
            float zSearchRange = turnHeight / 2.0;
            float currentZDist = 1000.0;
            //minError = 1000.0;
            for (float i = -zSearchRange; i < zSearchRange; i = i + 0.01) {
                float currentPredT = (thisStrand[0].Z() + i) / c;

                float predictedX = avgA * cos(avgW * currentPredT + avgSigma);
                float predictedY = avgB * sin(avgW * currentPredT + avgSigma);
                float xPos = (thisStrand[0].X());
                float yPos = (thisStrand[0].Y());
                Vector3 origPrediction = Vector3(predictedX, predictedY, thisStrand[0].Z() + i);
                float inverseAdd = 0.0;
                if (inverseOrientation) {
                    inverseAdd = M_PI;
                }
                Vector3 rotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), minErrorRotation + inverseAdd) * origPrediction;
                predictedX = rotatedPrediction[0];
                predictedY = rotatedPrediction[1];


                float dist = (xPos - predictedX) * (xPos - predictedX) + (yPos - predictedY) * (yPos - predictedY);


                for (int j = 1; j < thisStrand.size(); j++) {
                    currentPredT = (thisStrand[j].Z() + i) / c;

                    predictedX = avgA * cos(avgW * currentPredT + avgSigma);
                    predictedY = avgB * sin(avgW * currentPredT + avgSigma);
                    xPos = (thisStrand[j].X());
                    yPos = (thisStrand[j].Y());
                    origPrediction = Vector3(predictedX, predictedY, thisStrand[j].Z() + i);
                    inverseAdd = 0.0;
                    if (inverseOrientation) {
                        inverseAdd = M_PI;
                    }
                    rotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), minErrorRotation + inverseAdd) * origPrediction;
                    predictedX = rotatedPrediction[0];
                    predictedY = rotatedPrediction[1];


                    dist += (xPos - predictedX) * (xPos - predictedX) + (yPos - predictedY) * (yPos - predictedY);

                }
                if (dist < currentZDist) {
                    currentZDist = dist;
                    zTrans = i;
                }

            }



            cout << "dist thresh " << currentZDist / thisStrand.size() << endl;
            //cout << "minErrors " << minErrorRotation << " " << minError << endl;
            //cout << "inverse " << inverseOrientation << endl;
            //if(currentZDist / thisStrand.size() > 0.3) {
            //return;
            //}
            //OpenGLUtils::SetColor(0.9, 0.6, 0.3, 1.0);
            //std::vector<Vector3DFloat> originalPositions;
            for (int i = 0; i < thisStrand.size(); i++) {
                float xPos = (thisStrand[i].X());
                float yPos = (thisStrand[i].Y());
                float zPos = (thisStrand[i].Z());
                Vector3 unrotatedPos = rotationMatrices[strandIndex] * Vector3(xPos, yPos, zPos);
                //Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0]+helixMins[strandIndex].X(), unrotatedPos[1]+helixMins[strandIndex].Y(), unrotatedPos[2]+helixMins[strandIndex].Z()); 
                Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0], unrotatedPos[1], unrotatedPos[2]);
                //DrawSphere(Vector3DFloat(helixMins[strandIndex].X() + originalPos.X(), helixMins[strandIndex].Y() + originalPos.Y(), helixMins[strandIndex].Z() + originalPos.Z()), 0.1);
                //ofs << "rotatedPosOrig " << strandIndex << " " << i << " " << helixMins[strandIndex].X() + originalPos.X() << " " << helixMins[strandIndex].Y() + originalPos.Y() << " " << helixMins[strandIndex].Z() + originalPos.Z() << endl;
                //DrawSphere(Vector3DFloat(originalPos.X(), originalPos.Y(), originalPos.Z()), 0.1);

                //float xPos = (thisStrand[i].X() - normalizedMinXs[strandIndex]) / xRange;
                //float yPos = (thisStrand[i].Y() - normalizedMinYs[strandIndex]) / yRange;
                //float zPos = (thisStrand[i].Z() - normalizedMinZs[strandIndex]) / zRange;
                //cout << "sphere pos converted back " << helixMins[strandIndex].X() + originalPos.X() << " " << helixMins[strandIndex].Y() + originalPos.Y() << " " << helixMins[strandIndex].Z() + originalPos.Z() << endl;
                //DrawSphere(Vector3DFloat(10.0*strandIndex + xPos, 10.0*strandIndex + yPos, 10.0*strandIndex + zPos), 0.1);
                //originalPositions.push_back(originalPos);

            }
            //originalHelices.push_back(originalPositions);
            std::vector<Vector3DFloat> validHelixPts;
            for (float i = thisStrand[0].Z(); i < thisStrand[thisStrand.size() - 1].Z(); i = i + 0.01) {

                float currentT = (i + zTrans) / c;
                float xPos = avgA * cos(avgW * currentT + avgSigma);
                float yPos = avgB * sin(avgW * currentT + avgSigma);
                //float zPos = i;

                Vector3 rotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), minErrorRotation) * Vector3(xPos, yPos, i);
                if (inverseOrientation) {
                    rotatedPrediction = Matrix4::rotation(Vector3(0, 0, 1.0), M_PI) * rotatedPrediction;
                }
                //cout << "xPos " << xPos << " " << yPos << " " << zPos << endl;
                //OpenGLUtils::SetColor(0.6, 0.9, 0.1, 1.0);
                Vector3 unrotatedPos = rotationMatrices[strandIndex] * Vector3(rotatedPrediction[0], rotatedPrediction[1], rotatedPrediction[2]);
                //Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0]+helixMins[strandIndex].X(), unrotatedPos[1]+helixMins[strandIndex].Y(), unrotatedPos[2]+helixMins[strandIndex].Z()); 
                Vector3DFloat originalPos = Vector3DFloat(helixMins[strandIndex].X() + unrotatedPos[0], helixMins[strandIndex].Y() + unrotatedPos[1], helixMins[strandIndex].Z() + unrotatedPos[2]);
                //Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0], unrotatedPos[1], unrotatedPos[2]); 

                validHelixPts.push_back(originalPos);

                //DrawSphere(originalPos, 0.03);
                //cout << "pushing back " << endl;
                //DrawSphere(Vector3DFloat(10.0*strandIndex + rotatedPrediction[0], 10.0*strandIndex +rotatedPrediction[1], 10.0*strandIndex + rotatedPrediction[2]), 0.02);
            }
            foundHelices.push_back(validHelixPts);
            foundHelixIndices.push_back(strandIndex);
            helixScoreThresholds.push_back(currentZDist / thisStrand.size());
            strandCycles.push_back(helixCycles);

            /**
            for(float i = thisStrand[0].Z(); i < thisStrand[0].Z() + (2.0 * M_PI * helixCycles); i = i + 0.01) {

            float xPos = avgA * cos(avgW * i + avgSigma);
            float yPos = avgB * sin(avgW * i + avgSigma);
            //float zPos = i;
            float zPos = c * i;
            Vector3 rotatedPrediction = Matrix4::rotation(Vector3(0,0,1.0), minErrorRotation) * Vector3(xPos, yPos, zPos);
            if(inverseOrientation) {
            rotatedPrediction = Matrix4::rotation(Vector3(0,0,1.0), M_PI) * rotatedPrediction;
            }
            //cout << "xPos " << xPos << " " << yPos << " " << zPos << endl;
            OpenGLUtils::SetColor(0.6, 0.9, 0.1, 1.0);
            Vector3 unrotatedPos = rotationMatrices[strandIndex] * Vector3(rotatedPrediction[0], rotatedPrediction[1], rotatedPrediction[2]);
            //Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0]+helixMins[strandIndex].X(), unrotatedPos[1]+helixMins[strandIndex].Y(), unrotatedPos[2]+helixMins[strandIndex].Z());
            Vector3DFloat originalPos = Vector3DFloat(helixMins[strandIndex].X() + unrotatedPos[0], helixMins[strandIndex].Y() + unrotatedPos[1], helixMins[strandIndex].Z() + unrotatedPos[2]);
            //Vector3DFloat originalPos = Vector3DFloat(unrotatedPos[0], unrotatedPos[1], unrotatedPos[2]);



            DrawSphere(originalPos, 0.03);

            //DrawSphere(Vector3DFloat(10.0*strandIndex + rotatedPrediction[0], 10.0*strandIndex +rotatedPrediction[1], 10.0*strandIndex + rotatedPrediction[2]), 0.02);
            }
            **/
            //cout << "helix params " << avgA << " " << avgB << " " << avgW << " " << avgSigma << endl;
        }

        void CAlphaRenderer::averageEllipsoidScale() {
            float avgEScale = 0.0;
            for (int i = 0; i < bonds.size(); i++) {
                avgEScale += ((bonds[i].eigenValue0 + bonds[i].eigenValue1 + bonds[i].eigenvalue) / 3.0);
            }
            avgEScale /= bonds.size();
            int scale = 0;
            if (avgEScale < 0.1) {
                while (avgEScale < 0.1) {
                    avgEScale *= 10.0;
                    scale += 1;
                }
            }
            else if (avgEScale > 10.0) {
                while (avgEScale > 10.0) {
                    scale -= 1;
                }
            }
            ellipsoidScaleFactor = scale;
        }

        void CAlphaRenderer::DrawBackboneModel(int subSceneIndex, bool selectEnabled) {
            GLfloat emissionColor[4] = { 1.0, 1.0, 1.0, 1.0 };

            if (subSceneIndex == 0) { // Drawing Atoms
                if (selectEnabled) {
                    atomHashKeys.clear();
                    glPushName(0);
                    glPushName(0);
                }

                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if (it->second.GetName() == "CA") {
                        glPushAttrib(GL_LIGHTING_BIT);
                        if (it->second.GetSelected()) {
                            glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                            glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                        }
                        else {
                            OpenGLUtils::SetColor(it->second.GetColorR(), it->second.GetColorG(), it->second.GetColorB(), it->second.GetColorA());
                        }

                        if (selectEnabled){
                            //TODO: possibly implement mouse picking using ray intersection
                            atomHashKeys.push_back(it->first); // adding the atom hash key as an element
                            glLoadName(static_cast<GLuint>(atomHashKeys.size() - 1)); // the index of the element just added
                        }
                        if (it->second.GetVisible()) {
                            if (!extremalMode) {
                                DrawSphere(it->second.GetPosition(), it->second.GetAtomRadius() * 0.3);
                            }
                        }

                        glPopAttrib();
                    }

                }
                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 1) { // Drawing Bonds
                if (selectEnabled) {
                    glPushName(1);
                    glPushName(0);
                }
                //OpenGLUtils::SetColor(0.6, 0.9, 0.1, 1.0);
                /**
                if(displayHelix) {
                for(int i = 0; i < foundHelices.size(); i++) {
                //cout << helixScoreThresholds[i] << " thr " <<  helixThreshold << endl;
                //cout << helixTurnHeights[i] << " coil height " << helixCoilHeight << endl;
                //cout << "seg size " <<  foundHelices[i].size() << " " << helixSegmentThreshold << endl;
                int helixIndex = foundHelixIndices[i];
                std::vector<float> currentStrand = maxCurveSepPts[helixIndex];
                if(orthoDistMin < orthoDistances[i] && orthoDistMax > orthoDistances[i]) {

                if(helixScoreThresholds[i] < helixThreshold && helixTurnHeights[i] < helixCoilHeight) {
                std::vector<Vector3DFloat> foundHelix = foundHelices[i];
                if((currentStrand.size()/3)-1  >= helixSegmentThreshold) {
                for(int j = 0; j <  foundHelix.size(); j++) {
                DrawSphere(foundHelix[j], 0.03);
                }
                }
                }
                }
                }
                }
                **/
                //OpenGLUtils::SetColor(0.1, 0.6, 0.9, 1.0);
                /**
                if(displayHelixPoints) {

                for(int i = 0; i < originalHelices.size(); i++) {
                std::vector<Vector3DFloat> thisStrand = originalHelices[i];
                if(helixScoreThresholds[i] < helixThreshold && helixTurnHeights[i] > helixCoilHeight) {
                for(int j = 0; j < thisStrand.size(); j++) {
                DrawSphere(thisStrand[j], 0.03);
                }
                }
                }


                for(int i = 0; i < foundHelixIndices.size(); i++) {
                int helixIndex = foundHelixIndices[i];
                //cout << helixScoreThresholds[helixIndex] << " thr1 " <<  helixTurnHeights[helixIndex] << endl;
                if(orthoDistMin < orthoDistances[i] && orthoDistMax > orthoDistances[i]) {

                if(helixScoreThresholds[helixIndex] < helixThreshold && helixTurnHeights[helixIndex] < helixCoilHeight) {
                std::vector<float> currentStrand = maxCurveSepPts[helixIndex];
                //cout << "seg size 1 " <<  currentStrand.size()/3 << " " << helixSegmentThreshold << endl;
                if( (currentStrand.size()/3)-1 >= helixSegmentThreshold) {
                for(int j = 0; j < (currentStrand.size() / 3)-1; j ++) {
                float xPos = currentStrand[3*j];
                float yPos = currentStrand[3*j+1];
                float zPos = currentStrand[3*j+2];

                float xPos1 = currentStrand[3*(j+1)];
                float yPos1 = currentStrand[3*(j+1)+1];
                float zPos1 = currentStrand[3*(j+1)+2];
                DrawCylinder(Vector3DFloat(xPos, yPos, zPos), Vector3DFloat(xPos1, yPos1, zPos1), 0.04);
                }
                }
                }
                }


                }

                }
                **/
                for (int i = 0; i < (int)bonds.size(); i++) {

                    glPushAttrib(GL_LIGHTING_BIT);
                    if (bonds[i].GetSelected()) {
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }

                    if (selectEnabled){
                        glLoadName(i);
                    }
                    float length = (atoms[bonds[i].GetAtom0Ix()].GetPosition() - atoms[bonds[i].GetAtom1Ix()].GetPosition()).Length();
                    if (length > 4.2) {
                        OpenGLUtils::SetColor(1.0, 0, 0, 1.0);
                    }

                    if (length < 3.3) {
                        OpenGLUtils::SetColor(0, 0, 1.0, 1.0);
                    }
                    if (bonds[i].tempDeleted == true) {
                        OpenGLUtils::SetColor(0.3, 0.86, 0.72, 1.0);
                    }
                    if (bonds[i].tempNew == true) {
                        OpenGLUtils::SetColor(0.7, 0.2, 0.4, 1.0);
                    }




                    float xPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().X();
                    float yPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().Y();
                    float zPos0 = atoms[bonds[i].GetAtom0Ix()].GetPosition().Z();

                    float xPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().X();
                    float yPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().Y();
                    float zPos1 = atoms[bonds[i].GetAtom1Ix()].GetPosition().Z();


                    bool display = true;
                    Vector3DFloat segmentDirection = atoms[bonds[i + 1].GetAtom0Ix()].GetPosition() - atoms[bonds[i].GetAtom1Ix()].GetPosition();
                    if (extremalMode) {
                        segmentDirection.Normalize();
                        display = true;

                        /**
                        float e0Size = bonds[i].eigenValue0;
                        float e1Size = bonds[i].eigenValue1;
                        float e2Size = bonds[i].eigenvalue;
                        cout << "ellipsoid sizes " << e0Size << " " << e1Size << " " << e2Size << endl;

                        if(ellipsoidVisible) {
                        DrawEllipsoid(10, 10, 20.0*ellipsoidScale * e0Size, 20.0*ellipsoidScale * e1Size, 20.0*ellipsoidScale * e2Size, (xPos0 + xPos1)/2.0, (yPos0 + yPos1)/2.0, (zPos0 + zPos1)/2.0, bondDirection[0], bondDirection[1], bondDirection[2], bondDNormalized1, bondDNormalized2);
                        }
                        **/
                        float pointR = pointRatio;
                        float curveR = curveRatio;
                        float surfaceR = surfaceRatio;
                        float minG = minGeo;
                        float maxG = maxGeo;
                        float eigenV = eigenValue;
                        if (saliencyCheck) {
                            float *saliencies = bonds[i].saliencies;
                            if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])){
                                display = false;
                            }
                        }
                        if (intensityCheck) {
                            float maxI = maxDataIntensity;
                            float minI = minDataIntensity;
                            float localI = bonds[i].intensity;
                            if (((localI - minI) / (maxI - minI) < maxG)) {
                                display = false;
                            }
                        }
                        if (eigenvalueCheck) {
                            float maxE = maxDataEigenvalue;
                            float minE = minDataEigenvalue;
                            float localE = bonds[i].eigenvalue;
                            if ((localE - minE) / (maxE - minE) < eigenV){
                                display = false;
                            }
                        }


                        if (bonds[i].maxOn) {
                            bonds[i].isDisplay = display;
                        }


                        display = false;
                        bool displayMax = false;
                        if (maxOn) {

                            //cout << "max curve eigen min " << maxCurveEigenVMin0 << " " << maxCurveEigenVMin1 << " " << maxCurveEigenVMin2 << " " << maxCurveEigenVMax0 << " " << maxCurveEigenVMax1 << " " << maxCurveEigenVMax2 << endl;
                            /**
                            float sqrtMin0 = log(sqrt(maxCurveEigenVMin0));
                            float sqrtMin1 = log(sqrt(maxCurveEigenVMin1));
                            float sqrtMin2 = log(sqrt(maxCurveEigenVMin2));
                            float sqrtMax0 = log(sqrt(maxCurveEigenVMax0));
                            float sqrtMax1 = log(sqrt(maxCurveEigenVMax1));
                            float sqrtMax2 = log(sqrt(maxCurveEigenVMax2));
                            float sqrtMinEigen = sqrtMin0;
                            if (sqrtMin1 < sqrtMinEigen) {
                            sqrtMinEigen = sqrtMin1;
                            }
                            if (sqrtMin2 < sqrtMinEigen) {
                            sqrtMinEigen = sqrtMin2;
                            }
                            float sqrtMaxEigen = sqrtMax0;
                            if (sqrtMax1 > sqrtMaxEigen) {
                            sqrtMaxEigen = sqrtMax1;
                            }
                            if (sqrtMax2 > sqrtMaxEigen) {
                            sqrtMaxEigen = sqrtMax2;
                            }
                            float maxLogEigen = sqrtMaxEigen;
                            float minLogEigen = sqrtMinEigen;
                            **/
                            //bonds[i].isDisplay = display;
                            if (bonds[i].maxOn) {
                                //OpenGLUtils::SetColor(0.5, 0.5, 1.0, 0.6);
                                displayMax = true;
                                /**
                                float bondDirection[3];
                                float bondDirection1[3];
                                float bondDirection2[3];
                                if(bonds[i].eigenValue0 < bonds[i].eigenValue1 && bonds[i].eigenValue0 < bonds[i].eigenvalue) {
                                bondDirection[0] = bonds[i].eigenVector0[0];
                                bondDirection[1] = bonds[i].eigenVector0[1];
                                bondDirection[2] = bonds[i].eigenVector0[2];
                                if(bonds[i].eigenValue1 < bonds[i].eigenvalue) {
                                bondDirection1[0] = bonds[i].eigenVector1[0];
                                bondDirection1[1] = bonds[i].eigenVector1[1];
                                bondDirection1[2] = bonds[i].eigenVector1[2];

                                bondDirection2[0] = bonds[i].eigenVector2[0];
                                bondDirection2[1] = bonds[i].eigenVector2[1];
                                bondDirection2[2] = bonds[i].eigenVector2[2];
                                }
                                else {
                                bondDirection1[0] = bonds[i].eigenVector2[0];
                                bondDirection1[1] = bonds[i].eigenVector2[1];
                                bondDirection1[2] = bonds[i].eigenVector2[2];

                                bondDirection2[0] = bonds[i].eigenVector1[0];
                                bondDirection2[1] = bonds[i].eigenVector1[1];
                                bondDirection2[2] = bonds[i].eigenVector1[2];
                                }
                                }
                                if(bonds[i].eigenValue1 < bonds[i].eigenValue0 && bonds[i].eigenValue1 < bonds[i].eigenvalue) {
                                bondDirection[0] = bonds[i].eigenVector1[0];
                                bondDirection[1] = bonds[i].eigenVector1[1];
                                bondDirection[2] = bonds[i].eigenVector1[2];
                                if(bonds[i].eigenValue0 < bonds[i].eigenvalue) {
                                bondDirection1[0] = bonds[i].eigenVector0[0];
                                bondDirection1[1] = bonds[i].eigenVector0[1];
                                bondDirection1[2] = bonds[i].eigenVector0[2];

                                bondDirection2[0] = bonds[i].eigenVector2[0];
                                bondDirection2[1] = bonds[i].eigenVector2[1];
                                bondDirection2[2] = bonds[i].eigenVector2[2];
                                }
                                else {
                                bondDirection1[0] = bonds[i].eigenVector2[0];
                                bondDirection1[1] = bonds[i].eigenVector2[1];
                                bondDirection1[2] = bonds[i].eigenVector2[2];

                                bondDirection2[0] = bonds[i].eigenVector0[0];
                                bondDirection2[1] = bonds[i].eigenVector0[1];
                                bondDirection2[2] = bonds[i].eigenVector0[2];
                                }
                                }
                                if(bonds[i].eigenvalue < bonds[i].eigenValue0 && bonds[i].eigenvalue < bonds[i].eigenValue1) {
                                bondDirection[0] = bonds[i].eigenVector2[0];
                                bondDirection[1] = bonds[i].eigenVector2[1];
                                bondDirection[2] = bonds[i].eigenVector2[2];
                                if(bonds[i].eigenValue0 < bonds[i].eigenValue1) {
                                bondDirection1[0] = bonds[i].eigenVector0[0];
                                bondDirection1[1] = bonds[i].eigenVector0[1];
                                bondDirection1[2] = bonds[i].eigenVector0[2];

                                bondDirection2[0] = bonds[i].eigenVector1[0];
                                bondDirection2[1] = bonds[i].eigenVector1[1];
                                bondDirection2[2] = bonds[i].eigenVector1[2];
                                }
                                else {
                                bondDirection1[0] = bonds[i].eigenVector1[0];
                                bondDirection1[1] = bonds[i].eigenVector1[1];
                                bondDirection1[2] = bonds[i].eigenVector1[2];

                                bondDirection2[0] = bonds[i].eigenVector0[0];
                                bondDirection2[1] = bonds[i].eigenVector0[1];
                                bondDirection2[2] = bonds[i].eigenVector0[2];
                                }
                                }


                                Vector3DFloat bondDNormalized = Vector3DFloat(bondDirection[0], bondDirection[1], bondDirection[2]);
                                bondDNormalized.Normalize();
                                bondDirection[0] = bondDNormalized.X();
                                bondDirection[1] = bondDNormalized.Y();
                                bondDirection[2] = bondDNormalized.Z();

                                Vector3DFloat bondDNormalized1 = Vector3DFloat(bondDirection1[0], bondDirection1[1], bondDirection1[2]);
                                bondDNormalized1.Normalize();

                                Vector3DFloat bondDNormalized2 = Vector3DFloat(bondDirection2[0], bondDirection2[1], bondDirection2[2]);
                                bondDNormalized2.Normalize();
                                **/
                                display = true;

                                /**
                                float e0Size = bonds[i].eigenValue0;
                                float e1Size = bonds[i].eigenValue1;
                                float e2Size = bonds[i].eigenvalue;
                                cout << "ellipsoid sizes " << e0Size << " " << e1Size << " " << e2Size << endl;

                                if(ellipsoidVisible) {
                                DrawEllipsoid(10, 10, 20.0*ellipsoidScale * e0Size, 20.0*ellipsoidScale * e1Size, 20.0*ellipsoidScale * e2Size, (xPos0 + xPos1)/2.0, (yPos0 + yPos1)/2.0, (zPos0 + zPos1)/2.0, bondDirection[0], bondDirection[1], bondDirection[2], bondDNormalized1, bondDNormalized2);
                                }
                                **/
                                float pointR = pointRatio;
                                float curveR = curveRatio;
                                float surfaceR = surfaceRatio;
                                float minG = minGeo;
                                float maxG = maxGeo;
                                float eigenV = eigenValue;
                                if (saliencyCheck) {
                                    float *saliencies = bonds[i].saliencies;
                                    if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])){
                                        display = false;
                                        displayMax = false;
                                    }
                                }
                                if (intensityCheck) {
                                    float maxI = maxDataIntensity;
                                    float minI = minDataIntensity;
                                    float localI = bonds[i].intensity;
                                    if (((localI - minI) / (maxI - minI) < maxG)) {
                                        display = false;
                                        displayMax = false;
                                    }
                                }
                                if (eigenvalueCheck) {
                                    float maxE = maxDataEigenvalue;
                                    float minE = minDataEigenvalue;
                                    float localE = bonds[i].eigenvalue;
                                    if ((localE - minE) / (maxE - minE) < eigenV){
                                        display = false;
                                        displayMax = false;
                                    }
                                }

                            }
                            else {
                                display = false;
                            }

                        }
                        bool displayMin = false;
                        if (minOn) {
                            if (bonds[i].minOn) {
                                OpenGLUtils::SetColor(0.6, 0.1, 0.85, 1.0);
                                display = true;
                                displayMin = true;

                                float xaxis[3];
                                xaxis[0] = 1.0;
                                xaxis[1] = 0.0;
                                xaxis[2] = 0.0;

                                float yaxis[3];
                                yaxis[0] = 0.0;
                                yaxis[1] = 1.0;
                                yaxis[2] = 0.0;

                                float zaxis[3];
                                zaxis[0] = 0.0;
                                zaxis[1] = 0.0;
                                zaxis[2] = 1.0;

                                float bondDirection[3];

                                bondDirection[0] = fabs(xPos1 - xPos0);
                                bondDirection[1] = fabs(yPos1 - yPos0);
                                bondDirection[2] = fabs(zPos1 - zPos0);
                                Vector3DFloat bondDNormalized = Vector3DFloat(bondDirection[0], bondDirection[1], bondDirection[2]);
                                bondDNormalized.Normalize();
                                bondDirection[0] = bondDNormalized.X();
                                bondDirection[1] = bondDNormalized.Y();
                                bondDirection[2] = bondDNormalized.Z();

                                if (ellipsoidVisible) {
                                    DrawEllipsoid(10, 10, 20.0 * ellipsoidScale * bonds[i].eigenValue0, 20.0* ellipsoidScale * bonds[i].eigenValue1, 20.0 * ellipsoidScale * bonds[i].eigenvalue, (xPos0 + xPos1) / 2.0, (yPos0 + yPos1) / 2.0, (zPos0 + zPos1) / 2.0, segmentDirection[0], segmentDirection[1], segmentDirection[2], Vector3DFloat(1.0, 1.0, 1.0), Vector3DFloat(1.0, 1.0, 1.0));
                                }

                                float pointR = pointRatio;
                                float curveR = curveRatio;
                                float surfaceR = surfaceRatio;
                                float minG = minGeo;
                                float maxG = maxGeo;
                                float eigenV = eigenValue;
                                if (saliencyCheck) {
                                    float *saliencies = bonds[i].saliencies;
                                    if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])){
                                        display = false;
                                        displayMin = false;
                                    }
                                }
                                if (intensityCheck) {
                                    float maxI = maxDataIntensity;
                                    float minI = minDataIntensity;
                                    float localI = bonds[i].intensity;
                                    //cout << "localI " << localI << " " << minI << " " << maxI << " " << minG << endl;
                                    if (((localI - minI) / (maxI - minI) > minG)){
                                        display = false;
                                        displayMin = false;
                                    }
                                }
                                if (eigenvalueCheck) {
                                    float maxE = maxDataEigenvalue;
                                    float minE = minDataEigenvalue;
                                    float localE = bonds[i].eigenvalue;
                                    if ((localE - minE) / (maxE - minE) < eigenV) {
                                        display = false;
                                        displayMin = false;
                                    }
                                }

                            }
                            else {
                                if (!displayMax) {
                                    display = false;
                                }
                            }
                        }
                        if (saddleOn) {
                            if (bonds[i].saddleOn) {
                                OpenGLUtils::SetColor(0.3, 0.3, 0.3, 1.0);
                                display = true;

                                float xaxis[3];
                                xaxis[0] = 1.0;
                                xaxis[1] = 0.0;
                                xaxis[2] = 0.0;

                                float yaxis[3];
                                yaxis[0] = 0.0;
                                yaxis[1] = 1.0;
                                yaxis[2] = 0.0;

                                float zaxis[3];
                                zaxis[0] = 0.0;
                                zaxis[1] = 0.0;
                                zaxis[2] = 1.0;

                                float bondDirection[3];
                                bondDirection[0] = fabs(xPos1 - xPos0);
                                bondDirection[1] = fabs(yPos1 - yPos0);
                                bondDirection[2] = fabs(zPos1 - zPos0);
                                Vector3DFloat bondDNormalized = Vector3DFloat(bondDirection[0], bondDirection[1], bondDirection[2]);
                                bondDNormalized.Normalize();
                                bondDirection[0] = bondDNormalized.X();
                                bondDirection[1] = bondDNormalized.Y();
                                bondDirection[2] = bondDNormalized.Z();

                                if (ellipsoidVisible) {
                                    DrawEllipsoid(10, 10, 20.0 * ellipsoidScale * bonds[i].eigenValue0, 20.0*ellipsoidScale * bonds[i].eigenValue1, 20.0 * ellipsoidScale * bonds[i].eigenvalue, (xPos0 + xPos1) / 2.0, (yPos0 + yPos1) / 2.0, (zPos0 + zPos1) / 2.0, bondDirection[0], bondDirection[1], bondDirection[2], Vector3DFloat(1.0, 1.0, 1.0), Vector3DFloat(1.0, 1.0, 1.0));
                                }
                                float pointR = pointRatio;
                                float curveR = curveRatio;
                                float surfaceR = surfaceRatio;
                                float minG = minGeo;
                                float maxG = maxGeo;
                                float eigenV = eigenValue;
                                if (saliencyCheck) {
                                    float *saliencies = bonds[i].saliencies;
                                    if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])){
                                        display = false;
                                    }
                                }
                                if (intensityCheck) {
                                    float maxI = maxDataIntensity;
                                    float minI = minDataIntensity;
                                    float localI = bonds[i].intensity;
                                }
                                if (eigenvalueCheck) {
                                    float maxE = maxDataEigenvalue;
                                    float minE = minDataEigenvalue;
                                    float localE = bonds[i].eigenvalue;
                                    if ((localE - minE) / (maxE - minE) < eigenV) display = false;
                                }

                            }
                            else {
                                if (!displayMax && !displayMin){
                                    display = false;
                                }
                            }
                        }
                    }


                    if (display && atoms[bonds[i].GetAtom0Ix()].GetVisible() && atoms[bonds[i].GetAtom1Ix()].GetVisible()) {
                        if (bonds[i].saddleOn == true) {
                            OpenGLUtils::SetColor(0.3, 0.3, 0.3, 1.0);
                        }
                        if (bonds[i].maxOn == true) {
                            OpenGLUtils::SetColor(0.1, 0.85, 0.6, 1.0);
                        }
                        if (bonds[i].minOn == true) {
                            OpenGLUtils::SetColor(0.6, 0.1, 0.85, 1.0);
                        }

                        DrawCylinder(atoms[bonds[i].GetAtom0Ix()].GetPosition(), atoms[bonds[i].GetAtom1Ix()].GetPosition(), 0.1, 10, 2);




                        OpenGLUtils::SetColor(0.5, 0.5, 1.0, 0.6);
                        float sqrtMin0 = log(sqrt(maxCurveEigenVMin0));
                        float sqrtMin1 = log(sqrt(maxCurveEigenVMin1));
                        float sqrtMin2 = log(sqrt(maxCurveEigenVMin2));
                        float sqrtMax0 = log(sqrt(maxCurveEigenVMax0));
                        float sqrtMax1 = log(sqrt(maxCurveEigenVMax1));
                        float sqrtMax2 = log(sqrt(maxCurveEigenVMax2));
                        float sqrtMinEigen = sqrtMin0;
                        if (sqrtMin1 < sqrtMinEigen) {
                            sqrtMinEigen = sqrtMin1;
                        }
                        if (sqrtMin2 < sqrtMinEigen) {
                            sqrtMinEigen = sqrtMin2;
                        }
                        float sqrtMaxEigen = sqrtMax0;
                        if (sqrtMax1 > sqrtMaxEigen) {
                            sqrtMaxEigen = sqrtMax1;
                        }
                        if (sqrtMax2 > sqrtMaxEigen) {
                            sqrtMaxEigen = sqrtMax2;
                        }
                        float maxLogEigen = sqrtMaxEigen;
                        float minLogEigen = sqrtMinEigen;

                        float bondDirection[3];
                        float bondDirection1[3];
                        float bondDirection2[3];
                        if (bonds[i].eigenValue0 < bonds[i].eigenValue1 && bonds[i].eigenValue0 < bonds[i].eigenvalue) {
                            bondDirection[0] = bonds[i].eigenVector0[0];
                            bondDirection[1] = bonds[i].eigenVector0[1];
                            bondDirection[2] = bonds[i].eigenVector0[2];

                            if (bonds[i].eigenValue1 < bonds[i].eigenvalue) {
                                bondDirection1[0] = bonds[i].eigenVector1[0];
                                bondDirection1[1] = bonds[i].eigenVector1[1];
                                bondDirection1[2] = bonds[i].eigenVector1[2];

                                bondDirection2[0] = bonds[i].eigenVector2[0];
                                bondDirection2[1] = bonds[i].eigenVector2[1];
                                bondDirection2[2] = bonds[i].eigenVector2[2];
                            }
                            else {
                                bondDirection1[0] = bonds[i].eigenVector2[0];
                                bondDirection1[1] = bonds[i].eigenVector2[1];
                                bondDirection1[2] = bonds[i].eigenVector2[2];

                                bondDirection2[0] = bonds[i].eigenVector1[0];
                                bondDirection2[1] = bonds[i].eigenVector1[1];
                                bondDirection2[2] = bonds[i].eigenVector1[2];
                            }

                        }
                        if (bonds[i].eigenValue1 < bonds[i].eigenValue0 && bonds[i].eigenValue1 < bonds[i].eigenvalue) {
                            bondDirection[0] = bonds[i].eigenVector1[0];
                            bondDirection[1] = bonds[i].eigenVector1[1];
                            bondDirection[2] = bonds[i].eigenVector1[2];
                            if (bonds[i].eigenValue0 < bonds[i].eigenvalue) {
                                bondDirection1[0] = bonds[i].eigenVector0[0];
                                bondDirection1[1] = bonds[i].eigenVector0[1];
                                bondDirection1[2] = bonds[i].eigenVector0[2];

                                bondDirection2[0] = bonds[i].eigenVector2[0];
                                bondDirection2[1] = bonds[i].eigenVector2[1];
                                bondDirection2[2] = bonds[i].eigenVector2[2];
                            }
                            else {
                                bondDirection1[0] = bonds[i].eigenVector2[0];
                                bondDirection1[1] = bonds[i].eigenVector2[1];
                                bondDirection1[2] = bonds[i].eigenVector2[2];

                                bondDirection2[0] = bonds[i].eigenVector0[0];
                                bondDirection2[1] = bonds[i].eigenVector0[1];
                                bondDirection2[2] = bonds[i].eigenVector0[2];
                            }
                        }
                        if (bonds[i].eigenvalue < bonds[i].eigenValue0 && bonds[i].eigenvalue < bonds[i].eigenValue1) {
                            bondDirection[0] = bonds[i].eigenVector2[0];
                            bondDirection[1] = bonds[i].eigenVector2[1];
                            bondDirection[2] = bonds[i].eigenVector2[2];
                            if (bonds[i].eigenValue0 < bonds[i].eigenValue1) {
                                bondDirection1[0] = bonds[i].eigenVector0[0];
                                bondDirection1[1] = bonds[i].eigenVector0[1];
                                bondDirection1[2] = bonds[i].eigenVector0[2];

                                bondDirection2[0] = bonds[i].eigenVector1[0];
                                bondDirection2[1] = bonds[i].eigenVector1[1];
                                bondDirection2[2] = bonds[i].eigenVector1[2];
                            }
                            else {
                                bondDirection1[0] = bonds[i].eigenVector1[0];
                                bondDirection1[1] = bonds[i].eigenVector1[1];
                                bondDirection1[2] = bonds[i].eigenVector1[2];

                                bondDirection2[0] = bonds[i].eigenVector0[0];
                                bondDirection2[1] = bonds[i].eigenVector0[1];
                                bondDirection2[2] = bonds[i].eigenVector0[2];
                            }
                        }


                        Vector3DFloat bondDNormalized = Vector3DFloat(bondDirection[0], bondDirection[1], bondDirection[2]);
                        //bondDNormalized.Normalize();
                        //bondDirection[0] = bondDNormalized.X();
                        //bondDirection[1] = bondDNormalized.Y();
                        //bondDirection[2] = bondDNormalized.Z();

                        bondDirection[0] = fabs(xPos1 - xPos0);
                        bondDirection[1] = fabs(yPos1 - yPos0);
                        bondDirection[2] = fabs(zPos1 - zPos0);
                        //cout << "bond direction " << bondDirection[0] << " " << bondDirection[1] << " " << bondDirection[2] << endl;

                        Vector3DFloat bondDNormalized1 = Vector3DFloat(bondDirection1[0], bondDirection1[1], bondDirection1[2]);
                        bondDNormalized1.Normalize();

                        Vector3DFloat bondDNormalized2 = Vector3DFloat(bondDirection2[0], bondDirection2[1], bondDirection2[2]);
                        bondDNormalized2.Normalize();

                        //display = true;


                        float e0Size = bonds[i].eigenValue0;
                        float e1Size = bonds[i].eigenValue1;
                        float e2Size = bonds[i].eigenvalue;
                        /**
                        float eAvg = (e0Size + e1Size + e2Size) / 2.0;
                        float scale = 0;
                        if (eAvg < 0.1) {
                        while(eAvg < 0.1) {
                        eAvg = 10.0 * eAvg;
                        scale++;
                        }
                        }
                        else if(eAvg > 10.0) {
                        while(eAvg > 10.0) {
                        eAvg /= 10.0;
                        scale--;
                        }
                        }
                        if(scale != 0 ){
                        if(scale > 0) {
                        for(int i = 0; i < scale; i++) {
                        e0Size = 10.0 * e0Size;
                        e1Size = 10.0 * e1Size;
                        e2Size = 10.0 * e2Size;
                        }
                        }
                        if(scale < 0) {
                        for(int i = 0; i > scale; i--) {
                        e0Size /= 10.0;
                        e1Size /= 10.0;
                        e2Size /= 10.0;
                        }
                        }
                        }
                        **/

                        float drawScale = 1.0;
                        if (ellipsoidScaleFactor != 0) {
                            if (ellipsoidScaleFactor > 0) {
                                for (int j = 0; j < ellipsoidScaleFactor; j++) {
                                    drawScale *= 10.0;
                                }
                            }
                            else if (ellipsoidScaleFactor < 0) {
                                for (int j = 0; j < ellipsoidScaleFactor; j++) {
                                    drawScale /= 10.0;
                                }
                            }
                        }

                        //cout << "ellipsoid sizes " << e0Size << " " << e1Size << " " << e2Size << endl;

                        if (ellipsoidVisible) {
                            DrawEllipsoid(10, 10, drawScale * ellipsoidScale * e0Size, drawScale * ellipsoidScale * e1Size, drawScale * ellipsoidScale * e2Size, (xPos0 + xPos1) / 2.0, (yPos0 + yPos1) / 2.0, (zPos0 + zPos1) / 2.0, bondDirection[0], bondDirection[1], bondDirection[2], bondDNormalized1, bondDNormalized2);
                        }


                    }
                    //else {
                    //display = false;
                    //}



                    //}

                    //Here
                    //OpenGLUtils::SetColor(1.0, 0.1, 0.85, 1.0);
                    //for(int i = 0; i < helices.size(); i++) {
                    //std::vector<unsigned long long> currentHelix = helices[i];
                    //for(int j = 0; j < currentHelix.size()-1; j++) {
                    //DrawCylinder(atoms[currentHelix[j]].GetPosition(), atoms[currentHelix[j+1]].GetPosition(), 0.1, 10, 2);
                    //}
                    //}

                    //OpenGLUtils::SetColor(0.1, 1.0, 0.45, 1.0);
                    //for(int i = 0; i < axisPts.size()/2; i++) {
                    //DrawCylinder(axisPts[2*i], axisPts[2*i+1], 0.05, 10, 1);
                    //}

                    /**
                    OpenGLUtils::SetColor(1.0, 0.0, 0.1, 0.5);
                    cout << "quadpts size " << std::to_string(quadPts.size()) << endl;
                    cout << std::to_string(quadNormals.size()) << endl;
                    if (quadPts.size() > 0) {
                    for(int i = 0; i < quadPts.size()-4; i+=4) {
                    int j = i/4;
                    DrawSurface(quadPts[i], quadPts[i+1], quadPts[i+2], quadPts[i+3], quadNormals[3*j], quadNormals[(3*j)+1],quadNormals[(3*j)+2] );
                    }
                    }

                    for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if(it->second.isMax) {
                    OpenGLUtils::SetColor(0.1, 0.2, 0.85, 1.0);
                    DrawSphere(it->second.GetPosition(), 0.3);
                    }
                    }
                    for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if(it->second.isMin) {
                    OpenGLUtils::SetColor(0.1, 0.5, 0.1, 1.0);
                    DrawSphere(it->second.GetPosition(), 0.3);
                    }
                    }
                    **/
                    glPopAttrib();
                }
                //glEnable( GL_BLEND ) ;
                //glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA ) ; // works fine, except cracks

                //glBlendFunc( GL_SRC_ALPHA_SATURATE, GL_ONE ) ; // works only on black background...

                //glPolygonMode( GL_FRONT, GL_FILL ) ;
                //glEnable( GL_POLYGON_SMOOTH ) ;
                //glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST ) ;
                //glDisable(GL_LIGHTING);
                //glDisable(GL_TEXTURE_2D);

                //GLfloat mcolor[] = {0.1, 0.9, 0.4, 1.0};
                //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mcolor);
                //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mcolor);
                //glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);
                //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


                //glColor3fv(mcolor); 

                if (maxSurfaceOn && minSurfaceOn){
                    for (int i = 0; i < quads.size(); i++) {
                        if (quads[i].type == 1) {
                            OpenGLUtils::SetColor(1.0, 0.0, 0.1, 0.6);
                            //cout << "drawing surface " << quads.size() << endl;
                            DrawSurface(quads[i]);

                        }
                        if (quads[i].type == 2) {
                            OpenGLUtils::SetColor(0.0, 1.0, 0.1, 1.0);
                            DrawSurface(quads[i]);
                        }
                    }
                }
                else if (maxSurfaceOn && !minSurfaceOn) {
                    for (int i = 0; i < quads.size(); i++) {
                        if (quads[i].type == 1) {
                            OpenGLUtils::SetColor(1.0, 0.0, 0.1, 1.0);
                            DrawSurface(quads[i]);
                        }
                    }
                }
                else if (!maxSurfaceOn && minSurfaceOn) {
                    for (int i = 0; i < quads.size(); i++) {
                        if (quads[i].type == 2) {
                            OpenGLUtils::SetColor(0.0, 1.0, 0.1, 1.0);
                            DrawSurface(quads[i]);
                        }
                    }
                }
                //glEnable(GL_LIGHTING);
                //glEnable(GL_TEXTURE_2D);

                //glDisable( GL_BLEND ) ;

                if (maxPointOn) {
                    //if(maxPointX != -1000.0 && maxPointY != -1000.0 && maxPointZ != -1000.0) {
                    OpenGLUtils::SetColor(0.1, 0.2, 0.85, 1.0);
                    for (int i = 0; i < maxCoord.size(); i++) {
                        float saliencies[3];
                        saliencies[0] = maxPointSaliencies[3 * i];
                        saliencies[1] = maxPointSaliencies[3 * i + 1];
                        saliencies[2] = maxPointSaliencies[3 * i + 2];
                        if (!CheckMaxPointHide(saliencies, maxPointIntensities[i], maxPointEigenvalues[i])) {
                            DrawSphere(Vector3DFloat(maxCoord[i].X(), maxCoord[i].Y(), maxCoord[i].Z()), 0.3);

                            DrawEllipsoidPoint(10, 10, 10.0*maxPointEigenvalues0[i], 10.0*maxPointEigenvalues1[i], 10.0*maxPointEigenvalues[i], maxCoord[i].X(), maxCoord[i].Y(), maxCoord[i].Z());
                        }
                    }
                    //}
                    //for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    //if(it->second.isMax) {
                    //OpenGLUtils::SetColor(0.1, 0.2, 0.85, 1.0);

                    //DrawSphere(it->second.GetPosition(), 0.3);
                    //}
                    //}
                }
                if (minPointOn) {
                    //if(minPointX != -1000.0 && minPointY != -1000.0 && minPointZ != -1000.0) {

                    OpenGLUtils::SetColor(0.1, 0.5, 0.1, 1.0);
                    for (int i = 0; i < minCoord.size(); i++) {
                        float saliencies[3];
                        saliencies[0] = minPointSaliencies[3 * i];
                        saliencies[1] = minPointSaliencies[3 * i + 1];
                        saliencies[2] = minPointSaliencies[3 * i + 2];
                        if (!CheckMinPointHide(saliencies, minPointIntensities[i], minPointEigenvalues[i])) {
                            DrawSphere(Vector3DFloat(minCoord[i].X(), minCoord[i].Y(), minCoord[i].Z()), 0.3);
                            DrawEllipsoidPoint(10, 10, 10.0*minPointEigenvalues0[i], 10.0*minPointEigenvalues1[i], 10.0*minPointEigenvalues[i], minCoord[i].X(), minCoord[i].Y(), minCoord[i].Z());
                        }

                    }
                }

                if (saddlePointOn) {
                    //if(minPointX != -1000.0 && minPointY != -1000.0 && minPointZ != -1000.0) {

                    OpenGLUtils::SetColor(0.5, 0.1, 0.5, 1.0);
                    for (int i = 0; i < saddleCoord.size(); i++) {
                        float saliencies[3];
                        saliencies[0] = saddlePointSaliencies[3 * i];
                        saliencies[1] = saddlePointSaliencies[3 * i + 1];
                        saliencies[2] = saddlePointSaliencies[3 * i + 2];
                        if (!CheckSaddlePointHide(saliencies, saddlePointIntensities[i], saddlePointEigenvalues[i])) {
                            DrawSphere(Vector3DFloat(saddleCoord[i].X(), saddleCoord[i].Y(), saddleCoord[i].Z()), 0.3);
                            DrawEllipsoidPoint(10, 10, 10.0*saddlePointEigenvalues0[i], 10.0*saddlePointEigenvalues1[i], 10.0*saddlePointEigenvalues[i], saddleCoord[i].X(), saddleCoord[i].Y(), saddleCoord[i].Z());
                        }

                    }
                }
                //if( coefsA.size() > 0 && helixMins.size() > 0) {
                //OpenGLUtils::SetColor(1.0, 0.5, 0.5, 1.0);
                //for(int i = 0; i < maxCurveSepPts.size(); i++) {
                //for(int i = 0; i < helixMins.size(); i++ ) {
                //DrawCylinder(helixMins[i], helixMaxs[i], 0.2);
                //}
                //}

                //for(int i = 0; i < normalizedHelices.size(); i++) {
                //std::vector< Vector3DFloat > currentStrand = normalizedHelices[i];
                //scoreHelix(i);
                //}


                OpenGLUtils::SetColor(0.4, 0.5, 0.5, 1.0);
                //for(int i = 0; i < ordLines.size(); i++) {
                //DrawCylinder(Vector3DFloat(0.0,0.0,0.0), Vector3DFloat(10.0*ordLines[i].X(), 10.0*ordLines[i].Y(), 10.0*ordLines[i].Z()), 0.2);
                //}
                //}
                //for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                //if(it->second.isMin) {
                //OpenGLUtils::SetColor(0.1, 0.5, 0.1, 1.0);
                //DrawSphere(it->second.GetPosition(), 0.3);
                //}
                //}
                //}


                //if(displayHelix) {
                for (int i = 0; i < fittedHelices.size(); i++) {
                    bool displayCurrentHelix = false;
                    for (int j = 0; j < fittedHelices[i].fittedHelixBondIndices.size(); j++) {
                        if (bonds[j].isDisplay) {
                            displayCurrentHelix = true;
                        }
                    }

                    if (displayCurrentHelix) {
                        if (displayHelix) {
                            OpenGLUtils::SetColor(0.6, 0.9, 0.1, 1.0);
                            // && strandCycles[i] > helixSegmentThreshold

                            if (orthoDistMin < orthoDistances[i] && orthoDistMax > orthoDistances[i]) {
                                //cout << "helix const " << " " << strandCycles[i] << " " << helixScoreThresholds[i] << " " << helixTurnHeights[i] << " " << displayCurrentHelix << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].helixError << " " << fittedHelices[i].coilHeight << endl;
                                if (helixTurnHeights[i] < helixCoilHeight && helixScoreThresholds[i] < helixThreshold) {
                                    if (fittedHelices[i].fittedHelixBondIndices.size() >= helixSegmentThreshold) {
                                        //if(helixScoreThresholds[i] < helixThreshold && helixTurnHeights[i] < helixCoilHeight) {
                                        std::vector<Vector3DFloat> helPositions = fittedHelices[i].fittedHelixBondPositions;
                                        //cout << "helBondSize " << fittedHelices[i].fittedHelixBondIndices.size() << endl;
                                        for (int k = 0; k < helPositions.size(); k++) {
                                            //cout << "helPositions " << helPositions[k].X() << " " << helPositions[k].Y() << " " << helPositions[k].Z() << endl;
                                            DrawSphere(helPositions[k], 0.03);
                                        }
                                    }
                                }
                                //}
                            }
                        }
                        if (displayHelixPoints) {
                            OpenGLUtils::SetColor(0.1, 0.6, 0.9, 1.0);
                            if (orthoDistMin < orthoDistances[i] && orthoDistMax > orthoDistances[i]) {
                                //cout << "helix const " << helixScoreThresholds[i] << " " << helixTurnHeights[i] << " " << displayCurrentHelix << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].helixError << " " << fittedHelices[i].coilHeight << endl;
                                //cout << "helix const " << " " << strandCycles[i] << " " << helixScoreThresholds[i] << " " << helixTurnHeights[i] << " " << displayCurrentHelix << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].orthoDist << " " << fittedHelices[i].helixError << " " << fittedHelices[i].coilHeight << endl;
                                if (helixTurnHeights[i] < helixCoilHeight) {
                                    if (fittedHelices[i].fittedHelixBondIndices.size() >= helixSegmentThreshold) {
                                        for (int k = 0; k < fittedHelices[i].fittedHelixBondIndices.size(); k++) {
                                            DrawCylinder(atoms[bonds[fittedHelices[i].fittedHelixBondIndices[k]].GetAtom0Ix()].GetPosition(), atoms[bonds[fittedHelices[i].fittedHelixBondIndices[k]].GetAtom1Ix()].GetPosition(), 0.03);
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
                /**
                if(bonds[i].helixIndex != -1) {
                int helixIndex = bonds[i].helixIndex;
                std::vector<float> currentStrand = maxCurveSepPts[helixIndex];
                if(orthoDistMin < orthoDistances[bonds[i].helixIndex] && orthoDistMax > orthoDistances[bonds[i].helixIndex]) {

                if(helixScoreThresholds[bonds[i].helixIndex] < helixThreshold && helixTurnHeights[bonds[i].helixIndex] < helixCoilHeight) {
                std::vector<Vector3DFloat> foundHelix = foundHelices[bonds[i].helixIndex];
                if((currentStrand.size()/3)-1  >= helixSegmentThreshold) {
                for(int j = 0; j <  foundHelix.size(); j++) {
                DrawSphere(foundHelix[j], 0.03);
                }
                }
                }
                }
                }
                //}
                OpenGLUtils::SetColor(0.1, 0.6, 0.9, 1.0);
                if(displayHelixPoints) {


                //for(int i = 0; i < foundHelixIndices.size(); i++) {
                if(bonds[i].helixIndex != -1) {
                std::vector<unsigned long long> bonds[bondIndex].helixIndices;
                int helixIndex = foundHelixIndices[bonds[i].helixIndex];
                //cout << helixScoreThresholds[helixIndex] << " thr1 " <<  helixTurnHeights[helixIndex] << endl;
                if(orthoDistMin < orthoDistances[bonds[i].helixIndex] && orthoDistMax > orthoDistances[bonds[i].helixIndex]) {

                if(helixScoreThresholds[helixIndex] < helixThreshold && helixTurnHeights[helixIndex] < helixCoilHeight) {
                std::vector<float> currentStrand = maxCurveSepPts[helixIndex];
                //cout << "seg size 1 " <<  currentStrand.size()/3 << " " << helixSegmentThreshold << endl;
                if( (currentStrand.size()/3)-1 >= helixSegmentThreshold) {
                for(int j = 0; j < (currentStrand.size() / 3)-1; j ++) {
                float xPos = currentStrand[3*j];
                float yPos = currentStrand[3*j+1];
                float zPos = currentStrand[3*j+2];

                float xPos1 = currentStrand[3*(j+1)];
                float yPos1 = currentStrand[3*(j+1)+1];
                float zPos1 = currentStrand[3*(j+1)+2];
                DrawCylinder(Vector3DFloat(xPos, yPos, zPos), Vector3DFloat(xPos1, yPos1, zPos1), 0.04);
                }
                }
                }
                }
                }


                //}

                }
                **/



                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 2) { // Drawing spheres to cover up the cylinder edges				
                //for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                //if(i->second.GetName() == "CA") {
                //DrawSphere(i->second.GetPosition(), 0.1);
                //}
                //}
            }
        }

        void CAlphaRenderer::toggleSaddleOn(int atom1, int atom2, bool display) {
            if (display == true) {
                unsigned long long atom1Hash = 1;

                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                    if (it->second.GetSerial() == atom1) {
                        atom1Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }


                }
                unsigned long long atom2Hash = 1;
                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if (it->second.GetSerial() == atom2) {
                        atom2Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }
                }
                if (atom1Hash == 1 || atom2Hash == 1) {
                    return;
                }
                else {
                    int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                    if (bondIndex == -1) {
                        PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                        tempBond.saddleOn = true;
                        AddBond(tempBond);
                    }
                    else {
                        bonds[bondIndex].saddleOn = true;
                    }
                }
            }
        }

        void CAlphaRenderer::toggleMinOn(int atom1, int atom2, bool display) {
            if (display == true) {
                unsigned long long atom1Hash = 1;

                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                    if (it->second.GetSerial() == atom1) {
                        atom1Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }


                }
                unsigned long long atom2Hash = 1;
                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if (it->second.GetSerial() == atom2) {
                        atom2Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }
                }
                if (atom1Hash == 1 || atom2Hash == 1) {
                    return;
                }
                else {
                    int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                    if (bondIndex == -1) {
                        PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                        tempBond.tempNew = true;
                        tempBond.original = false;
                        tempBond.minOn = true;
                        AddBond(tempBond);
                    }
                    else {
                        bonds[bondIndex].minOn = true;
                        bonds[bondIndex].tempNew = true;
                    }
                }
            }
        }

        void CAlphaRenderer::toggleMaxOn(int atom1, int atom2, bool display) {
            if (display == true) {
                unsigned long long atom1Hash = 1;

                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                    if (it->second.GetSerial() == atom1) {
                        atom1Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }


                }
                unsigned long long atom2Hash = 1;
                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if (it->second.GetSerial() == atom2) {
                        atom2Hash = it->second.GetHashKey();
                        //atomHashes.push_back(it->second.GetHashKey());
                    }
                }
                if (atom1Hash == 1 || atom2Hash == 1) {
                    return;
                }
                else {
                    int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                    if (bondIndex == -1) {
                        PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                        tempBond.tempNew = true;
                        tempBond.original = false;
                        tempBond.maxOn = true;
                        AddBond(tempBond);
                    }
                    else {
                        bonds[bondIndex].maxOn = true;
                        bonds[bondIndex].tempNew = true;
                    }
                }
            }
        }




        void CAlphaRenderer::setMaxOn(bool display) {
            maxOn = display;
        }

        void CAlphaRenderer::setMinOn(bool display) {
            minOn = display;
        }

        void CAlphaRenderer::setSaddleOn(bool display) {
            saddleOn = display;
        }

        void CAlphaRenderer::DeleteAtomFromVisualization(unsigned long long deletedAtom) {
            atoms[deletedAtom].SetVisible(false);
        }

        string CAlphaRenderer::FindDistance(int atom1, int atom2) {
            unsigned long long atom1Hash = 1;
            PDBAtom atm1;
            PDBAtom atm2;

            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                if (it->second.GetSerial() == atom1) {
                    atom1Hash = it->second.GetHashKey();
                    atm1 = it->second;
                    //atomHashes.push_back(it->second.GetHashKey());
                }


            }
            unsigned long long atom2Hash = 1;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSerial() == atom2) {
                    atom2Hash = it->second.GetHashKey();
                    atm2 = it->second;
                    //atomHashes.push_back(it->second.GetHashKey());
                }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
                return "unknown";
            }
            else {
                return std::to_string((atm1.GetPosition() - atm2.GetPosition()).Length());

            }
        }

        void CAlphaRenderer::AddSaddleBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1) {
            /**
            unsigned long long atom1Hash = 1;

            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

            if(it->second.GetSerial() == atom1) {
            atom1Hash = it->second.GetHashKey();
            }


            }
            unsigned long long atom2Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == atom2) {
            atom2Hash = it->second.GetHashKey();
            }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
            return;
            }
            else {
            **/
            unsigned long long atom1Hash = std::stoll(atom1);
            unsigned long long atom2Hash = std::stoll(atom2);
            if (!atoms.count(atom1Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index1 * 3];
                float py = vertexPositions[index1 * 3 + 1];
                float pz = vertexPositions[index1 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[atom1Hash] = newAtom;
            }
            if (!atoms.count(atom2Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index2 * 3];
                float py = vertexPositions[index2 * 3 + 1];
                float pz = vertexPositions[index2 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[atom2Hash] = newAtom;
            }
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            if (bondIndex == -1) {
                PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                tempBond.saddleOn = true;
                tempBond.saliencies[0] = saliency1;
                tempBond.saliencies[1] = saliency2;
                tempBond.saliencies[2] = saliency3;
                tempBond.intensity = intensity;
                tempBond.eigenvalue = eigenvalue;
                tempBond.eigenValue0 = eigenvalue0;
                tempBond.eigenValue1 = eigenvalue1;
                //tempBond.original = false;
                AddBond(tempBond);
            }
            else {
                bonds[bondIndex].saddleOn = true;
                bonds[bondIndex].saliencies[0] = saliency1;
                bonds[bondIndex].saliencies[1] = saliency2;
                bonds[bondIndex].saliencies[2] = saliency3;
                bonds[bondIndex].intensity = intensity;
                bonds[bondIndex].eigenvalue = eigenvalue;
                bonds[bondIndex].eigenValue0 = eigenvalue0;
                bonds[bondIndex].eigenValue1 = eigenvalue1;
            }
            //}
        }

        template <class T>
        inline std::string to_string(const T& t)
        {
            std::stringstream ss;
            ss << t;
            return ss.str();
        }

        void CAlphaRenderer::addVertexPos(float px, float py, float pz) {
            vertexPositions.push_back(px);
            vertexPositions.push_back(py);
            vertexPositions.push_back(pz);
        }

        void CAlphaRenderer::AddMaxBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1) {

            /**unsigned long long atom1Hash = 1;

            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

            if(it->second.GetSerial() == atom1) {
            atom1Hash = it->second.GetHashKey();
            }


            }
            unsigned long long atom2Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == atom2) {
            atom2Hash = it->second.GetHashKey();
            }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
            return;
            }
            **/
            //else {

            //cout << "casted unsigned long long atom1" << " " << std::stoll(atom1) << endl;

            unsigned long long atom1Hash = std::stoll(atom1);
            unsigned long long atom2Hash = std::stoll(atom2);
            bondAtm1s.push_back(atom1Hash);
            bondAtm2s.push_back(atom2Hash);
            maxCurveIndices1.push_back(index1);
            maxCurveIndices2.push_back(index2);
            //maxCurveSaliencies.push_back(saliency1);
            //maxCurveSaliencies.push_back(saliency2);
            //maxCurveSaliencies.push_back(saliency3);
            //maxCurveIntensities.push_back(intensity);
            //maxCurveEigenvalues.push_back(eigenvalue);
            //maxCurveEigenvalues.push_back(eigenvalue0);
            //maxCurveEigenvalues.push_back(eigenvalue1);

            //if(!atoms.count(atom1Hash)) {
            PDBAtom newAtom = PDBAtom();
            float px = vertexPositions[index1 * 3];
            float py = vertexPositions[index1 * 3 + 1];
            float pz = vertexPositions[index1 * 3 + 2];
            newAtom.SetPosition(Vector3DFloat(px, py, pz));
            atoms[atom1Hash] = newAtom;
            //}
            //if(!atoms.count(atom2Hash)) {
            PDBAtom newAtom1 = PDBAtom();
            float px1 = vertexPositions[index2 * 3];
            float py1 = vertexPositions[index2 * 3 + 1];
            float pz1 = vertexPositions[index2 * 3 + 2];
            newAtom1.SetPosition(Vector3DFloat(px1, py1, pz1));
            atoms[atom2Hash] = newAtom1;
            //}
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            //cout << atom1Hash << " " << atom2Hash << endl;
            if (bondIndex == -1) {
                PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                tempBond.maxOn = true;
                tempBond.saliencies[0] = saliency1;
                tempBond.saliencies[1] = saliency2;
                tempBond.saliencies[2] = saliency3;
                tempBond.intensity = intensity;
                tempBond.eigenvalue = eigenvalue;
                tempBond.eigenValue0 = eigenvalue0;
                tempBond.eigenValue1 = eigenvalue1;

                tempBond.logEigen0 = log(eigenvalue0);
                tempBond.logEigen1 = log(eigenvalue1);
                tempBond.logEigen2 = log(eigenvalue);
                AddBond(tempBond);
            }
            else {
                bonds[bondIndex].maxOn = true;
                bonds[bondIndex].SetAtom0Ix(atom1Hash);
                bonds[bondIndex].SetAtom1Ix(atom2Hash);
                bonds[bondIndex].saliencies[0] = saliency1;
                bonds[bondIndex].saliencies[1] = saliency2;
                bonds[bondIndex].saliencies[2] = saliency3;
                bonds[bondIndex].intensity = intensity;
                bonds[bondIndex].eigenvalue = eigenvalue;
                bonds[bondIndex].eigenValue0 = eigenvalue0;
                bonds[bondIndex].eigenValue1 = eigenvalue1;

                bonds[bondIndex].logEigen0 = log(eigenvalue0);
                bonds[bondIndex].logEigen1 = log(eigenvalue1);
                bonds[bondIndex].logEigen2 = log(eigenvalue);
            }
            //}
        }

        void CAlphaRenderer::setExtremalParams(float curve, float point, float surface, float minG, float maxG, float eigen) {
            curveRatio = curve;
            pointRatio = point;
            surfaceRatio = surface;
            minGeo = minG;
            maxGeo = maxG;
            eigenValue = eigen;
        }

        void CAlphaRenderer::setExtremalChecks(bool saliency, bool intensity, bool eigenvalue) {
            saliencyCheck = saliency;
            intensityCheck = intensity;
            eigenvalueCheck = eigenvalue;
        }

        void CAlphaRenderer::setExtremalHide(bool isHide) {
            extremalHide = isHide;
        }

        void CAlphaRenderer::setExtremalSurfaceHide(bool isHide) {
            extremalSurfaceHide = isHide;
        }

        void CAlphaRenderer::setMaxSurfaceOn(bool maxOn) {
            maxSurfaceOn = maxOn;
        }

        void CAlphaRenderer::setMinSurfaceOn(bool minOn) {
            minSurfaceOn = minOn;
        }

        void CAlphaRenderer::AddMinBond(int index1, int index2, string atom1, string atom2, float saliency1, float saliency2, float saliency3, float intensity, float eigenvalue, float eigenvalue0, float eigenvalue1) {
            /**unsigned long long atom1Hash = 1;

            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

            if(it->second.GetSerial() == atom1) {
            atom1Hash = it->second.GetHashKey();
            //atomHashes.push_back(it->second.GetHashKey());
            }


            }
            unsigned long long atom2Hash = 1;
            for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
            if(it->second.GetSerial() == atom2) {
            atom2Hash = it->second.GetHashKey();
            //atomHashes.push_back(it->second.GetHashKey());
            }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
            return;
            }
            else {
            **/
            unsigned long long atom1Hash = std::stoll(atom1);
            unsigned long long atom2Hash = std::stoll(atom2);
            if (!atoms.count(atom1Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index1 * 3];
                float py = vertexPositions[index1 * 3 + 1];
                float pz = vertexPositions[index1 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[atom1Hash] = newAtom;
            }
            if (!atoms.count(atom2Hash)) {
                PDBAtom newAtom = PDBAtom();
                float px = vertexPositions[index2 * 3];
                float py = vertexPositions[index2 * 3 + 1];
                float pz = vertexPositions[index2 * 3 + 2];
                newAtom.SetPosition(Vector3DFloat(px, py, pz));
                atoms[atom2Hash] = newAtom;
            }
            int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
            if (bondIndex == -1) {
                PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                tempBond.minOn = true;
                tempBond.saliencies[0] = saliency1;
                tempBond.saliencies[1] = saliency2;
                tempBond.saliencies[2] = saliency3;
                tempBond.intensity = intensity;
                tempBond.eigenvalue = eigenvalue;
                tempBond.eigenValue0 = eigenvalue0;
                tempBond.eigenValue1 = eigenvalue1;
                AddBond(tempBond);
            }
            else {
                bonds[bondIndex].minOn = true;
                bonds[bondIndex].saliencies[0] = saliency1;
                bonds[bondIndex].saliencies[1] = saliency2;
                bonds[bondIndex].saliencies[2] = saliency3;
                bonds[bondIndex].intensity = intensity;
                bonds[bondIndex].eigenvalue = eigenvalue;
                bonds[bondIndex].eigenValue0 = eigenvalue0;
                bonds[bondIndex].eigenValue1 = eigenvalue1;
            }
            //}
        }

        void CAlphaRenderer::DrawAddedBond(int atom1, int atom2) {
            unsigned long long atom1Hash = 1;

            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                if (it->second.GetSerial() == atom1) {
                    atom1Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }


            }
            unsigned long long atom2Hash = 1;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSerial() == atom2) {
                    atom2Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
                return;
            }
            else {
                int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                if (bondIndex == -1) {
                    PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                    tempBond.tempNew = true;
                    tempBond.original = false;
                    AddBond(tempBond);
                }
                else {
                    bonds[bondIndex].tempNew = true;
                }
            }
        }

        void CAlphaRenderer::UndrawAddedBond(int atom1, int atom2) {
            unsigned long long atom1Hash = 1;

            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                if (it->second.GetSerial() == atom1) {
                    atom1Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }


            }
            unsigned long long atom2Hash = 1;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSerial() == atom2) {
                    atom2Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
                return;
            }
            else {
                int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                if (bondIndex == -1) {
                    return;
                }
                if (bonds[bondIndex].original == false) {
                    DeleteBond(bondIndex);
                }
                else {
                    bonds[bondIndex].tempNew = false;
                }
            }
        }

        void CAlphaRenderer::DrawDeletedBond(int atom1, int atom2) {
            //std::vector<unsigned long long> atomHashes;
            unsigned long long atom1Hash = 1;

            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                if (it->second.GetSerial() == atom1) {
                    atom1Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }


            }
            unsigned long long atom2Hash = 1;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSerial() == atom2) {
                    atom2Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
                return;
            }
            else {
                int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                if (bondIndex == -1) {
                    cout << "no bond found " << endl;
                    PDBBond tempBond = PDBBond(atom1Hash, atom2Hash, false);
                    tempBond.tempDeleted = true;
                    tempBond.original = false;
                    AddBond(tempBond);
                }
                else {
                    bonds[bondIndex].tempDeleted = true;
                }
                //OpenGLUtils::SetColor(0.3, 1.0, 0.5, 1.0);
                //atoms[atom1Hash].Print();
                //atoms[atom2Hash].Print();
                //DrawCylinder(atoms[atom1Hash].GetPosition(), atoms[atom2Hash].GetPosition(), 0.1, 10, 2);
            }

        }

        void CAlphaRenderer::UndrawDeletedBond(int atom1, int atom2) {
            //std::vector<unsigned long long> atomHashes;
            unsigned long long atom1Hash = 1;

            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {

                if (it->second.GetSerial() == atom1) {
                    atom1Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }


            }
            unsigned long long atom2Hash = 1;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSerial() == atom2) {
                    atom2Hash = it->second.GetHashKey();
                    //atomHashes.push_back(it->second.GetHashKey());
                }
            }
            if (atom1Hash == 1 || atom2Hash == 1) {
                return;
            }
            else {
                int bondIndex = max(GetBondIndex(atom1Hash, atom2Hash), GetBondIndex(atom2Hash, atom1Hash));
                if (bondIndex == -1) {
                    return;
                }
                if (bonds[bondIndex].original == false) {
                    DeleteBond(bondIndex);
                }
                else {
                    bonds[bondIndex].tempDeleted = false;
                }
            }

        }


        void CAlphaRenderer::DrawBackboneModelPathwalker(int subSceneIndex, bool selectEnabled) {
            GLfloat emissionColor[4] = { 1.0, 1.0, 1.0, 1.0 };

            if (subSceneIndex == 0) { // Drawing Atoms
                if (selectEnabled) {
                    atomHashKeys.clear();
                    glPushName(0);
                    glPushName(0);
                }
                for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                    if (it->second.GetName() == "CA") {
                        glPushAttrib(GL_LIGHTING_BIT);
                        if (it->second.GetSelected()) {
                            glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                            glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                        }
                        else {
                            OpenGLUtils::SetColor(it->second.GetColorR(), it->second.GetColorG(), it->second.GetColorB(), it->second.GetColorA());
                        }

                        if (selectEnabled){
                            //TODO: possibly implement mouse picking using ray intersection
                            atomHashKeys.push_back(it->first); // adding the atom hash key as an element
                            glLoadName(static_cast<GLuint>(atomHashKeys.size() - 1)); // the index of the element just added
                        }
                        if (it->second.GetVisible()) {
                            Vector3DFloat currentPosition = it->second.GetPosition();
                            DrawSphere(currentPosition, it->second.GetAtomRadius() * 0.3);
                        }

                        glPopAttrib();
                    }

                }
                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 1) { // Drawing Bonds
                if (selectEnabled) {
                    glPushName(1);
                    glPushName(0);
                }
                for (int i = 0; i < (int)bonds.size(); i++) {
                    glPushAttrib(GL_LIGHTING_BIT);
                    if (bonds[i].GetSelected()) {
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }

                    if (selectEnabled){
                        glLoadName(i);
                    }
                    float length = (atoms[bonds[i].GetAtom0Ix()].GetPosition() - atoms[bonds[i].GetAtom1Ix()].GetPosition()).Length();
                    if (length > 4.2) {
                        OpenGLUtils::SetColor(1.0, 0, 0, 1.0);
                    }

                    if (length < 3.3) {
                        OpenGLUtils::SetColor(0, 0, 1.0, 1.0);
                    }

                    if (atoms[bonds[i].GetAtom0Ix()].GetVisible() && atoms[bonds[i].GetAtom1Ix()].GetVisible()) {
                        Vector3DFloat currentPosition1 = atoms[bonds[i].GetAtom0Ix()].GetPosition();
                        Vector3DFloat currentPosition2 = atoms[bonds[i].GetAtom1Ix()].GetPosition();
                        DrawCylinder(currentPosition1, currentPosition2, 0.1, 10, 2);
                    }
                    glPopAttrib();
                }
                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 2) { // Drawing spheres to cover up the cylinder edges				
                for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                    if (i->second.GetName() == "CA") {
                        Vector3DFloat currentPosition = i->second.GetPosition();
                        DrawSphere(currentPosition, 0.1);
                    }
                }
            }
        }

        //void CAlphaRenderer::DrawBackboneModelPathwalker(int subSceneIndex, bool selectEnabled)
        //{
        //    return;
        //}

        void CAlphaRenderer::DrawRibbonModel(int subSceneIndex, bool selectEnabled) {
            if (selectEnabled) {
                glPushName(subSceneIndex);
                glPushName(0);
            }
            //GLfloat emissionColor[4] = {1.0, 1.0, 1.0, 1.0};
            //GLfloat frontColor[4] = {1.0, 0.0, 0.0, 1.0};
            //GLfloat backColor[4] = {0.0, 0.0, 1.0, 1.0};
            GLfloat emissionColor[4] = { hlt_r, hlt_g, hlt_b, hlt_a };
            GLfloat frontColor[4] = { 1.0, 0.0, 0.0, 1.0 };
            GLfloat backColor[4] = { 0.0, 0.0, 1.0, 1.0 };
            std::vector<int> PDBIndices;

            switch (subSceneIndex) {
            case 0: // Helices

                //if(subSceneIndex == 0){
                for (unsigned int i = 0; i < corrs.size(); ++i){
                    int SSEIndex = get<1>(corrs[i]);
                    for (unsigned int k = 0; k < selectedSSEHelices.size(); ++k){
                        if (selectedSSEHelices[k] == SSEIndex){
                            PDBIndices.push_back(get<0>(corrs[i]));
                        }
                    }
                }

                for (int i = 0; i < aHelices.size(); i++) {
                    if (selectEnabled){
                        glLoadName(i);
                    }

                    glPushAttrib(GL_LIGHTING_BIT);

                    Secel currentSecel = aHelices[i];

                    if (currentSecel.selected == true){
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }
                    map<int, boost::tuple<float, float, float> >::iterator iter = helixColors.begin();
                    iter = helixColors.find(i);
                    if (iter != helixColors.end()){

                        OpenGLUtils::SetColor(get<0>(helixColors[i]), get<1>(helixColors[i]), get<2>(helixColors[i]), 1.0);
                    }
                    else{
                        //OpenGLUtils::SetColor(0.8,0.8,0.8,1.0);
                    }

                    if (currentSecel.atomHashes.size() > 0){

                        PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
                        PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size() - 1])->second;
                        Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
                        Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

                        std::vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
                        std::vector<Vector3DFloat> tangents = std::vector<Vector3DFloat>(points);
                        std::vector<Vector3DFloat> axes = std::vector<Vector3DFloat>(points);
                        std::vector<Vector3DFloat> interpPoints = std::vector<Vector3DFloat>((points.size() - 1)*NUM_SEGMENTS + 1);
                        int flatSlices = 2;
                        int rptsize = interpPoints.size() * 4;
                        switch (renderingType){
                        case 0:
                            rptsize = interpPoints.size()*flatSlices;
                            break;
                        case 1:
                            rptsize = interpPoints.size() * 4;
                            break;
                        default:
                            rptsize = interpPoints.size() * 4;
                            break;
                        }
                        std::vector<Vector3DFloat> renderingPoints(rptsize);
                        std::vector<Vector3DFloat> renderingNormals(renderingPoints.size());
                        /*std::vector<Vector3DFloat> renderingPoints(interpPoints.size()*NUM_SLICES);
                        std::vector<Vector3DFloat> renderingNormals(renderingPoints.size());*/

                        HermiteCurve curve;
                        Vector3DFloat m0, m1;

                        CreateHelixAxesTangentsAndPoints(axes, tangents, interpPoints, points, preSecelAtomPos, postSecelAtomPos, HELIX_ALPHA, HELIX_BETA, HELIX_HERMITE_FACTOR);

                        for (unsigned int x = 0; x < points.size() - 1; ++x){

                            m0 = tangents[x];
                            m1 = tangents[x + 1];

                            curve.setCurve(points[x], points[x + 1], m0, m1);

                            // used in rendering a helix as glowing if selected
                            std::vector<Vector3DFloat> selectedBoxPositions(8);

                            float halfwidth = HELIX_WIDTH / 2.0;
                            float halfthickness = LOOP_RADIUS;
                            Vector3DFloat lastPos = points[x];
                            //int NUM_SECTIONS = 10;
                            for (int sect = 0; sect <= NUM_SEGMENTS; ++sect){
                                if (sect == 0 && x != 0){
                                    continue;
                                }
                                double tsect = ((double)sect) / ((double)NUM_SEGMENTS);
                                Vector3DFloat nextPos = curve.getPos(tsect);

                                Vector3DFloat currentAxis = axes[x] * (1.0 - tsect) + axes[x + 1] * tsect;
                                currentAxis.Normalize();

                                Vector3DFloat curnormal = curve.getTangent(tsect);
                                curnormal = curnormal^currentAxis;
                                curnormal.Normalize();

                                if (x == 0){
                                    halfwidth = LOOP_RADIUS + (0.5 * HELIX_WIDTH - LOOP_RADIUS) * 0.5 * (-1 * cos(PI*tsect) + 1.0);
                                }
                                else if (x == points.size() - 2){
                                    halfwidth = LOOP_RADIUS + (0.5 * HELIX_WIDTH - LOOP_RADIUS) * 0.5 * (cos(PI*tsect) + 1.0);
                                }

                                //for(int y = 0; y < NUM_SLICES; ++y){
                                //	renderingPoints[(x*NUM_SEGMENTS + sect) + y*interpPoints.size()] = nextPos + currentAxis*halfwidth*cos(y*2*PI/NUM_SLICES) + curnormal*halfthickness*sin(y*2*PI/NUM_SLICES);
                                //	renderingNormals[(x*NUM_SEGMENTS + sect) + y*interpPoints.size()] = currentAxis*halfwidth*cos(y*2*PI/NUM_SLICES)*sin(y*2*PI/NUM_SLICES)/(cos(y*2*PI/NUM_SLICES)+sin(y*2*PI/NUM_SLICES))
                                //		+ (curnormal*halfthickness*sin(y*2*PI/NUM_SLICES)*cos(y*2*PI/NUM_SLICES))/(cos(y*2*PI/NUM_SLICES)+sin(y*2*PI/NUM_SLICES));
                                //	renderingNormals[(x*NUM_SEGMENTS + sect) + y*interpPoints.size()].Normalize();
                                //}


                                switch (renderingType) {
                                case 0:
                                    for (int q = 0; q < flatSlices; ++q){
                                        renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*(q % 3 == 0 ? -1 : 1); //+ curnormal*hlt_r*pow(-1.0, q/2);
                                        renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal;//*(q%3 == 0 ? -1 : 1);
                                        renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()].Normalize();
                                    }
                                    break;
                                case 1:
                                    for (int q = 0; q < 4; ++q){
                                        renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*pow(-1.0, q / 2)
                                            + curnormal*halfthickness*(q % 3 == 0 ? -1 : 1);
                                        renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal*(q % 3 == 0 ? -1 : 1);
                                    }
                                    break;
                                default:
                                    cout << "should not have reached default case in ahelices rendering method" << endl;


                                }

                                //for(int q = 0; q < 4; ++q){
                                //	switch (renderingType){
                                //		case 0:
                                //			renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*pow(-1.0, q/2); //+ curnormal*thinRibbThickness*(q%3 == 0 ? -1 : 1);
                                //			renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal*(q%3 == 0 ? -1 : 1);
                                //			break;
                                //		case 1:
                                //			renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*pow(-1.0, q/2) + curnormal*halfthickness*(q%3 == 0 ? -1 : 1);
                                //			renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal*(q%3 == 0 ? -1 : 1);
                                //			break;
                                //		default:
                                //			break;
                                //	}

                                //}

                                //double HLT_FACTOR = 1.05;
                            }
                        }
                        switch (renderingType){
                        case 0:
                            DrawTube(renderingPoints, renderingNormals, interpPoints.size() - 1, flatSlices);
                            break;
                        case 1:
                            DrawTube(renderingPoints, renderingNormals, interpPoints.size() - 1, 4);
                            break;
                        default:
                            break;
                        }
                        //DrawTube(renderingPoints, renderingNormals, interpPoints.size() - 1, 4);
                        //DrawTube(renderingPoints, renderingNormals, interpPoints.size() - 1, NUM_SLICES);
                    }
                    glPopAttrib();

                    if (aHelices[i].selected == true){
                        glPushAttrib(GL_LIGHTING_BIT);

                        if (featureVecs.size() > 0){
                            OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                            DrawSphere(featureVecs[i].get<0>(), 1.0);
                            OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                            DrawSphere(featureVecs[i].get<1>(), 1.0);
                        }
                        else{
                            OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                            DrawSphere(atoms[aHelices[i].atomHashes[0]].GetPosition(), 1.0);
                            OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                            DrawSphere(atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size() - 1]].GetPosition(), 1.0);
                        }

                        glPopAttrib();
                        Vector3DFloat pos1 = atoms[aHelices[i].atomHashes[0]].GetPosition();
                        Vector3DFloat pos2 = atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size() - 1]].GetPosition();
                        printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i + 1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

                        fflush(stdout);
                    }

                    for (unsigned int j = 0; j < PDBIndices.size(); ++j){
                        if (PDBIndices[j] == i){
                            glPushAttrib(GL_LIGHTING_BIT);
                            OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                            DrawSphere(atoms[aHelices[i].atomHashes[0]].GetPosition(), 1.0);
                            OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                            DrawSphere(atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size() - 1]].GetPosition(), 1.0);
                            glPopAttrib();
                        }
                    }
                }
                //}
                break;
            case 1: // Strands
                for (unsigned int i = 0; i < bStrands.size(); ++i){
                    if (selectEnabled){
                        glLoadName(i);
                    }

                    int atom_counter = 0;
                    PDBAtom lastEnd;

                    HermiteCurve curve;
                    Vector3DFloat m0, m1, dir1, dir2;

                    glPushAttrib(GL_LIGHTING_BIT);

                    Secel currentSecel = bStrands[i];

                    if (currentSecel.selected == true){
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }

                    if (currentSecel.atomHashes.size() > 0){
                        PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
                        PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size() - 1])->second;
                        Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
                        Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

                        std::vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
                        int num_interp_points = (points.size() - 1)*NUM_SEGMENTS + 1;
                        int num_rendering_points = num_interp_points * 4;
                        switch (renderingType){
                        case 0:
                            num_rendering_points = num_interp_points * 2;
                            break;
                        case 1:
                            num_rendering_points = num_interp_points * 4;
                            break;
                        default:
                            cout << "bstrands" << endl;
                        }
                        std::vector<Vector3DFloat> renderingPoints(num_rendering_points);
                        std::vector<Vector3DFloat> normals = CreateStrandNormals(points, preSecelAtomPos, postSecelAtomPos);
                        std::vector<Vector3DFloat> renderingNormals(renderingPoints);
                        double arrowhead_factor = 1.0;
                        //std::vector<Vector3DFloat> boxpositions(8);
                        //std::vector<Vector3DFloat> boxnormals(8);

                        bool LAPLACIAN_SMOOTHING = true;
                        int SMOOTHING_STEPS = 1;
                        if (LAPLACIAN_SMOOTHING){
                            points = LaplacianSmoothing(points, SMOOTHING_STEPS);
                        }

                        for (unsigned int i = 0; i < points.size() - 1; ++i){
                            if (i == 0){
                                m0 = points[i + 1] - points[i];
                            }
                            else {
                                if (i + 2 < points.size()){
                                    m0 = points[i + 2] - points[i];
                                }
                                else {
                                    m0 = postSecelAtomPos - points[i];
                                }
                            }

                            if (i + 3 >= points.size()){
                                m1 = postSecelAtomPos - points[i];
                            }
                            else {
                                m1 = points[i + 3] - points[i + 1];
                            }

                            m0 = m0*STRAND_HERMITE_FACTOR;
                            m1 = m1*STRAND_HERMITE_FACTOR;

                            if (i == 0){
                                dir2 = points[1] - points[0];
                                dir2.Normalize();
                            }

                            dir1 = dir2;
                            if (i + 2 < points.size()){
                                dir2 = points[i + 2] - points[i];
                            }
                            else {
                                dir2 = postSecelAtomPos - points[i];
                            }

                            dir2.Normalize();

                            curve.setCurve(points[i], points[i + 1], m0, m1);

                            float WIDTH = 1.3;
                            float THICKNESS = LOOP_RADIUS;
                            Vector3DFloat lastPos = points[i];
                            Vector3DFloat direction = dir1;
                            Vector3DFloat currentNormal = normals[i];
                            Vector3DFloat side = currentNormal^direction;
                            side.Normalize();

                            for (int sect = 0; sect <= NUM_SEGMENTS; ++sect){
                                if (sect == 0 && i != 0){
                                    continue;
                                }

                                double tsect = ((double)sect) / ((double)NUM_SEGMENTS);
                                if (i > points.size() - 3){
                                    arrowhead_factor = 2.0*(1 - tsect);
                                }
                                direction = dir1*(1.0 - tsect) + dir2*tsect;
                                currentNormal = normals[i] * (1.0 - tsect) + normals[i + 1] * (tsect);
                                side = currentNormal^direction;
                                side.Normalize();
                                Vector3DFloat nextPos = curve.getPos(tsect);

                                switch (renderingType){
                                case 0:
                                    for (int q = 0; q < 2; ++q){
                                        renderingPoints[i*NUM_SEGMENTS + sect + q*num_interp_points] = nextPos + side*(0.5*WIDTH*arrowhead_factor + LOOP_RADIUS / 2.0)*(q % 3 == 0 ? -1 : 1);
                                        renderingNormals[i*NUM_SEGMENTS + sect + q*num_interp_points] = currentNormal;//*pow(-1.0, q/2);
                                    }
                                    break;
                                case 1:
                                    for (int q = 0; q < 4; ++q){
                                        renderingPoints[i*NUM_SEGMENTS + sect + q*num_interp_points] = nextPos + side*(0.5*WIDTH*arrowhead_factor + LOOP_RADIUS / 2.0)*(q % 3 == 0 ? -1 : 1)
                                            + currentNormal*0.5*THICKNESS*pow(-1.0, q / 2);
                                        renderingNormals[i*NUM_SEGMENTS + sect + q*num_interp_points] = currentNormal*pow(-1.0, q / 2);
                                    }
                                    break;
                                default:
                                    std::cout << "You've reached the default.  You should never reach the default." << std::endl;
                                    break;
                                }


                            }

                            if (i == points.size() - 3){
                                arrowhead_factor = 2.0;
                            }
                        }
                        switch (renderingType){
                        case 0:
                            DrawTube(renderingPoints, renderingNormals, num_interp_points - 1, 2);
                            break;
                        case 1:
                            DrawTube(renderingPoints, renderingNormals, num_interp_points - 1, 4);
                            break;
                        default:
                            cout << "should not have reached bstrand drawtube default" << endl;
                        }

                    }
                    glPopAttrib();

                    if (currentSecel.selected == true){
                        glPushAttrib(GL_LIGHTING_BIT);

                        Vector3DFloat pos1 = atoms[currentSecel.atomHashes[0]].GetPosition();
                        Vector3DFloat pos2 = atoms[currentSecel.atomHashes[currentSecel.atomHashes.size() - 1]].GetPosition();
                        printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i + 1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

                        /*if(featureVecs.size() > 0){
                        OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                        DrawSphere(featureVecs[i].get<0>(), 1.0);
                        OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                        DrawSphere(featureVecs[i].get<1>(), 1.0);
                        }else{*/
                        OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                        DrawSphere(pos1, 1.0);
                        OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                        DrawSphere(pos2, 1.0);
                        /*}*/

                        glPopAttrib();

                        fflush(stdout);
                    }
                }
                break;
                //}
                //else if (subSceneIndex == 2){
            case 2: // Loops
                for (unsigned int i = 0; i < loops.size(); ++i){
                    if (selectEnabled){
                        glLoadName(i);
                    }
                    //GLfloat emissionColor[4] = {1.0, 1.0, 1.0, 1.0};

                    int atom_counter = 0;
                    PDBAtom lastEnd;

                    HermiteCurve curve;
                    Vector3DFloat m0, m1;

                    double HERMITE_FACTOR = 0.5;

                    glPushAttrib(GL_LIGHTING_BIT);
                    Secel currentSecel = loops[i];

                    if (currentSecel.selected == true){
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }

                    if (currentSecel.atomHashes.size() > 1){
                        PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
                        PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size() - 1])->second;
                        Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
                        Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

                        std::vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
                        std::vector<Vector3DFloat> normals(points.size());
                        normals = CreateStrandNormals(points, preSecelAtomPos, postSecelAtomPos);

                        // generate smoothed interpolated points
                        std::vector<Vector3DFloat> interpolatedPoints = InterpolateLoopPoints(points, preSecelAtomPos, postSecelAtomPos, NUM_SEGMENTS);
                        int ptsize = interpolatedPoints.size();

                        // create vectors to hold the vertices and normals of our loop polygon
                        int renderingsize = (ptsize)*NUM_SLICES;
                        switch (renderingType){
                        case 0:
                            renderingsize = ptsize * 2;
                            break;
                        case 1:
                            renderingsize = (ptsize)*NUM_SLICES;
                            break;
                        default:
                            break;
                        }
                        std::vector<Vector3DFloat> renderingPoints(renderingsize);
                        std::vector<Vector3DFloat> renderingNormals(renderingPoints.size());

                        Vector3DFloat nextVector = interpolatedPoints[1] - interpolatedPoints[0];
                        Vector3DFloat previousVector = preSecelAtomPos - interpolatedPoints[0];
                        if (previousVector.Length() < .001){
                            if (2 < interpolatedPoints.size()){
                                previousVector = interpolatedPoints[2] - interpolatedPoints[0];
                            }
                            else {
                                previousVector = postSecelAtomPos - interpolatedPoints[0];
                            }
                        }
                        Vector3DFloat curAxis = nextVector^previousVector;
                        Vector3DFloat curNormal = nextVector^curAxis;
                        curNormal.Normalize();
                        Vector3DFloat curPos = interpolatedPoints[0];
                        nextVector.Normalize();

                        // generate first stack of points
                        Vector3DFloat outward = normals[0] ^ (interpolatedPoints[1] - interpolatedPoints[0]);
                        outward.Normalize();
                        switch (renderingType){
                        case 0:
                            for (unsigned int k = 0; k < 2; ++k){
                                renderingPoints[k*ptsize] = curPos + outward*LOOP_RADIUS*(k % 2 == 0 ? 1 : -1);
                                renderingNormals[k*ptsize] = normals[0];
                            }
                            break;
                        case 1:
                            for (unsigned int k = 0; k < NUM_SLICES; ++k){
                                Vector3DFloat outwardNormal = (curNormal.Rotate(nextVector, ((double)k) * 2 * PI / ((double)NUM_SLICES)));
                                outwardNormal.Normalize();
                                renderingPoints[k*(ptsize)] = curPos + outwardNormal*LOOP_RADIUS;
                                renderingNormals[k*(ptsize)] = outwardNormal; //renderingPoints[j*(ptsize)]-curPos;
                            }
                            break;
                        default:
                            break;
                        }

                        for (unsigned int j = 1; j < ptsize; ++j){
                            Vector3DFloat nextVector, previousVector;
                            curPos = interpolatedPoints[j];
                            if (j == ptsize - 1){
                                nextVector = postSecelAtomPos - interpolatedPoints[j];
                            }
                            else {
                                nextVector = interpolatedPoints[j + 1] - interpolatedPoints[j];
                            }

                            if (j == 0){
                                previousVector = interpolatedPoints[0] - preSecelAtomPos;
                            }
                            else {
                                previousVector = interpolatedPoints[j - 1] - interpolatedPoints[j];
                            }
                            Vector3DFloat newAxis = nextVector^previousVector;
                            if (newAxis.Length() < .0001){
                                int pix = j - 2;
                                int nix = j + 2;
                                if (pix < 0){
                                    pix = 0;
                                }
                                if (nix >= ptsize){
                                    nix = ptsize - 1;
                                }

                                newAxis = (interpolatedPoints[nix] - interpolatedPoints[j]) ^ (interpolatedPoints[pix] - interpolatedPoints[j]);
                            }
                            curAxis = newAxis;
                            double alpha = asin(curAxis.Length() / (nextVector.Length()*previousVector.Length()));
                            curAxis.Normalize();
                            curNormal.Normalize();
                            if (nextVector.Length() > .001 && previousVector.Length() > .001){
                                curNormal = curNormal.Rotate(curAxis, -1 * alpha);
                            }

                            nextVector.Normalize();

                            Vector3DFloat dirtemp;
                            if (j + 1 < ptsize){
                                dirtemp = (interpolatedPoints[j + 1] - interpolatedPoints[j]);
                            }
                            else {
                                dirtemp = interpolatedPoints[j] - interpolatedPoints[j - 1];
                            }
                            float tsect = ((float)(j%NUM_SEGMENTS)) / ((float)NUM_SEGMENTS);
                            Vector3DFloat outward = (normals[j / NUM_SEGMENTS] * (1.0 - tsect) + normals[j / NUM_SEGMENTS + 1] * (tsect)) ^ dirtemp;
                            outward.Normalize();
                            switch (renderingType){
                            case 0:
                                for (unsigned int k = 0; k < 2; ++k){
                                    renderingPoints[j + k*(ptsize)] = curPos + outward*LOOP_RADIUS*(k % 2 == 0 ? 1 : -1);
                                    renderingNormals[j + k*(ptsize)] = normals[j / NUM_SEGMENTS];
                                }
                                break;
                            case 1:
                                for (unsigned int k = 0; k < NUM_SLICES; ++k){
                                    Vector3DFloat outwardNormal = (curNormal.Rotate(nextVector, ((double)k * 2 * PI) / NUM_SLICES));
                                    outwardNormal.Normalize();
                                    renderingPoints[j + k*(ptsize)] = curPos + outwardNormal*LOOP_RADIUS;
                                    renderingNormals[j + k*(ptsize)] = outwardNormal;
                                }
                                break;
                            default:
                                cout << "shouldn't've gotten to the loop rendering default case" << endl;
                                break;
                            }


                        }
                        float r, g, b, a;
                        if (currentSecel.selected == true){
                            OpenGLUtils::GetColor(r, g, b, a);
                            OpenGLUtils::SetColor(hlt_r, hlt_g, hlt_b, hlt_a);
                        }
                        switch (renderingType){
                        case 0:
                            DrawTube(renderingPoints, renderingNormals, interpolatedPoints.size() - 1, 2);
                            break;
                        case 1:
                            DrawTube(renderingPoints, renderingNormals, interpolatedPoints.size() - 1, NUM_SLICES);
                            break;
                        default:
                            cout << "shouldn't have reached loop drawtube default case" << endl;
                        }
                        if (currentSecel.selected == true){
                            OpenGLUtils::SetColor(r, g, b, a);
                        }
                    }

                    glPopAttrib();

                    if (currentSecel.selected == true){
                        glPushAttrib(GL_LIGHTING_BIT);

                        Vector3DFloat pos1 = atoms[currentSecel.atomHashes[0]].GetPosition();
                        Vector3DFloat pos2 = atoms[currentSecel.atomHashes[currentSecel.atomHashes.size() - 1]].GetPosition();
                        printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i + 1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

                        /*if(featureVecs.size() > 0){
                        OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                        DrawSphere(featureVecs[i].get<0>(), 1.0);
                        OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                        DrawSphere(featureVecs[i].get<1>(), 1.0);
                        }else{*/
                        OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
                        DrawSphere(pos1, 1.0);
                        OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
                        DrawSphere(pos2, 1.0);
                        /*}*/

                        glPopAttrib();

                        fflush(stdout);
                    }
                }
                break;
            }

            if (selectEnabled) {
                glPopName();
                glPopName();
            }
        }

        void CAlphaRenderer::DrawSideChainModel(int subSceneIndex, bool selectEnabled) {
            GLfloat emissionColor[4] = { 1.0, 1.0, 1.0, 1.0 };
            float r, g, b, a;

            if (subSceneIndex == 0) { // Drawing Atoms				
                if (selectEnabled) {
                    atomHashKeys.clear();
                    glPushName(0);
                    glPushName(0);
                }
                for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                    glPushAttrib(GL_LIGHTING_BIT);
                    if (i->second.GetSelected()) {
                        glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                        glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                    }
                    else {
                        i->second.GetColor(r, g, b, a);
                        OpenGLUtils::SetColor(r, g, b, a);
                    }

                    if (selectEnabled){
                        //TODO: possibly implement mouse picking using ray intersection
                        atomHashKeys.push_back(i->first); // adding the hash key
                        glLoadName(static_cast<GLuint>(atomHashKeys.size() - 1)); // using the index of the element just added
                    }
                    if (i->second.GetVisible()) {
                        DrawSphere(i->second.GetPosition(), i->second.GetAtomRadius() * 0.3);
                    }

                    glPopAttrib();

                }
                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 1) { // Drawing Bonds
                if (selectEnabled) {
                    glPushName(1);
                    glPushName(0);
                }
                Vector3DFloat v1, vc, v2;


                for (int i = 0; i < (int)sidechainBonds.size(); i++) {
                    glPushAttrib(GL_LIGHTING_BIT);
                    if (selectEnabled){
                        glLoadName(i);
                    }

                    if (atoms[sidechainBonds[i].GetAtom0Ix()].GetVisible() && atoms[sidechainBonds[i].GetAtom1Ix()].GetVisible()) {
                        v1 = atoms[sidechainBonds[i].GetAtom0Ix()].GetPosition();
                        v2 = atoms[sidechainBonds[i].GetAtom1Ix()].GetPosition();
                        if (sidechainBonds[i].GetSelected()) {
                            glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
                            glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
                            DrawCylinder(v1, v2, 0.1, 6, 2);
                        }
                        else {
                            vc = (v1 + v2) * 0.5;
                            atoms[sidechainBonds[i].GetAtom0Ix()].GetColor(r, g, b, a);
                            OpenGLUtils::SetColor(r, g, b, a);
                            DrawCylinder(v1, vc, 0.1, 6, 2);
                            atoms[sidechainBonds[i].GetAtom1Ix()].GetColor(r, g, b, a);
                            OpenGLUtils::SetColor(r, g, b, a);
                            DrawCylinder(vc, v2, 0.1, 6, 2);
                        }
                    }
                    glPopAttrib();
                }
                if (selectEnabled) {
                    glPopName();
                    glPopName();
                }
            }
            else if (subSceneIndex == 2) { // Drawing spheres to cover up the cylinder edges					
                for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                    glPushAttrib(GL_LIGHTING_BIT);
                    i->second.GetColor(r, g, b, a);
                    OpenGLUtils::SetColor(r, g, b, a);
                    DrawSphere(i->second.GetPosition(), 0.1);
                    glPopAttrib();
                }
            }
        }

        void CAlphaRenderer::Draw(int subSceneIndex, bool selectEnabled) {
            switch (displayStyle) {
            case CALPHA_DISPLAY_STYLE_BACKBONE: // Backbone only
                DrawBackboneModel(subSceneIndex, selectEnabled);
                break;
            case CALPHA_DISPLAY_STYLE_RIBBON: // Ribbon mode
                DrawRibbonModel(subSceneIndex, selectEnabled);
                break;
            case CALPHA_DISPLAY_STYLE_SIDE_CHAIN: // Side chains
                DrawSideChainModel(subSceneIndex, selectEnabled);
                break;
            }
        }

        PDBAtom * CAlphaRenderer::GetAtomFromHitStack(int subsceneIndex, bool forceTrue, int ix0, int ix1, int ix2, int ix3, int ix4) {
            if (subsceneIndex == 0) {
                //TODO: possibly implement mouse picking using ray intersection
                AtomMapType::iterator it = atoms.find(atomHashKeys.at(ix0));
                if (it == atoms.end())
                    return NULL;
                else
                    return &(it->second);
            }
            return NULL;
        }



        void CAlphaRenderer::LoadFile(string fileName) {
            Renderer::LoadFile(fileName);
            atoms.clear();
            bonds.clear();
            atoms = PDBReader::ReadAtomPositions(fileName);

            // Keeping only C-Alpha atoms
            std::vector<unsigned long long> eraseKeys;
            eraseKeys.clear();

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                if (i->second.GetName().compare("CA") != 0) {
                    eraseKeys.push_back(i->first);
                }
            }

            for (unsigned int i = 0; i < eraseKeys.size(); i++) {
                atoms.erase(atoms.find(eraseKeys[i]));
            }

            eraseKeys.clear();

            std::list<SerialAndHashType> sortedSerials;
            SerialAndHashType elem;
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                elem.hashKey = i->first;
                elem.serial = i->second.GetSerial();

                sortedSerials.push_back(elem);
            }
            //sortedSerials.sort(SerialAndHashTypePredicate());


            std::list<SerialAndHashType>::iterator oldAtom = sortedSerials.begin();
            std::list<SerialAndHashType>::iterator startAtom = sortedSerials.begin();

            startAtom++;
            for (std::list<SerialAndHashType>::iterator i = startAtom; i != sortedSerials.end(); i++) {
                bonds.push_back(PDBBond(oldAtom->hashKey, i->hashKey, false));
                oldAtom = i;
            }
            sortedSerials.clear();
            UpdateBoundingBox();

        }

        void CAlphaRenderer::LoadSSEHunterFile(string fileName) {
            Renderer::LoadFile(fileName);
            atoms.clear();
            bonds.clear();
            atoms = PDBReader::ReadAtomPositions(fileName);

            float maxTempFactor = -10000.0f, minTempFactor = 10000.0f;
            float tempFactor;

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                tempFactor = i->second.GetTempFactor();
                if (tempFactor > maxTempFactor) {
                    maxTempFactor = tempFactor;
                }
                if (tempFactor < minTempFactor) {
                    minTempFactor = tempFactor;
                }
            }
            float r, g, b;

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                i->second.SetAtomRadius(3.0);
                tempFactor = i->second.GetTempFactor();
                if (tempFactor < 0) {
                    tempFactor = (tempFactor / minTempFactor);
                    r = 1.0f - tempFactor;
                    g = 1.0f - tempFactor;
                    b = 1.0f;
                }
                else {
                    tempFactor = (tempFactor / maxTempFactor);
                    r = 1.0f;
                    g = 1.0f - tempFactor;
                    b = 1.0f - tempFactor;
                }

                i->second.SetColor(r, g, b, 1.0f);
            }
            UpdateBoundingBox();

        }
        bool CAlphaRenderer::SaveSSEHunterFile(string fileName) {
            return PDBReader::WriteAtomPositions(atoms, fileName);
        }

        //		void CAlphaRenderer::GetSSEHunterAtoms(Volume * vol, NonManifoldMesh_Annotated * skeleton, float resolution, float threshold, float correlationCoeff, float skeletonCoeff, float geometryCoeff) {
        //			Renderer::LoadFile("");
        //			atoms.clear();
        //			bonds.clear();
        //
        //			SSEHunter * hunter = new SSEHunter();
        //			atoms = hunter->GetScoredAtoms(vol, skeleton, resolution, threshold, correlationCoeff, skeletonCoeff, geometryCoeff);
        //			delete hunter;
        //
        //			ColorSSEHunterAtoms();
        //		}

        void CAlphaRenderer::UpdateTotalScoreSSEHunterAtoms(float correlationCoeff, float skeletonCoeff, float geometryCoeff) {
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                i->second.SetTempFactor(i->second.GetTotalScore(correlationCoeff, skeletonCoeff, geometryCoeff));
            }
            ColorSSEHunterAtoms();
        }

        void CAlphaRenderer::ColorSSEHunterAtoms() {
            float maxTempFactor = -10000.0f, minTempFactor = 10000.0f;
            float tempFactor;

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                tempFactor = i->second.GetTempFactor();
                if (tempFactor > maxTempFactor) {
                    maxTempFactor = tempFactor;
                }
                if (tempFactor < minTempFactor) {
                    minTempFactor = tempFactor;
                }
            }
            float r, g, b;

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                i->second.SetAtomRadius(3.0);
                tempFactor = i->second.GetTempFactor();
                if (tempFactor < 0) {
                    tempFactor = (tempFactor / minTempFactor);
                    r = 1.0f - tempFactor;
                    g = 1.0f - tempFactor;
                    b = 1.0f;
                }
                else {
                    tempFactor = (tempFactor / maxTempFactor);
                    r = 1.0f;
                    g = 1.0f - tempFactor;
                    b = 1.0f - tempFactor;
                }

                i->second.SetColor(r, g, b, 1.0f);
            }

            UpdateBoundingBox();
        }

        int CAlphaRenderer::SelectionObjectCount(){
            int count = SelectionAtomCount();
            for (unsigned int i = 0; i < bonds.size(); i++) {
                if (bonds[i].GetSelected()) {
                    count++;
                }
            }
            return count;
        }

        int CAlphaRenderer::SelectionAtomCount(){
            int count = 0;
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                if (i->second.GetSelected()) {
                    count++;
                }
            }
            return count;
        }


        Vector3DFloat CAlphaRenderer::SelectionCenterOfMass() {
            int count = 0;
            Vector3DFloat centerOfMass = Vector3DFloat(0, 0, 0);
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                if (i->second.GetSelected()) {
                    count++;
                    centerOfMass = centerOfMass + i->second.GetPosition();
                }
            }

            for (unsigned int i = 0; i < bonds.size(); i++) {
                if (bonds[i].GetSelected()) {
                    count++;
                    centerOfMass = centerOfMass + (atoms[bonds[i].GetAtom0Ix()].GetPosition() + atoms[bonds[i].GetAtom1Ix()].GetPosition()) * 0.5;
                }
            }
            if (count == 0) {
                centerOfMass = Renderer::SelectionCenterOfMass();
            }
            else {
                centerOfMass = centerOfMass * (1.0f / (float)count);
            }
            return centerOfMass;
        }

        bool CAlphaRenderer::SelectionRotate(Vector3DFloat centerOfMass, Vector3DFloat rotationAxis, float angle) {
            bool rotated = false;
            Vector3 centerOfMassP3 = Vector3(centerOfMass.X(), centerOfMass.Y(), centerOfMass.Z());
            Vector3 rotationV3 = Vector3(rotationAxis.X(), rotationAxis.Y(), rotationAxis.Z());

            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                if (i->second.GetSelected()) {
                    rotated = true;
                    Vector3DFloat move = centerOfMass - i->second.GetPosition();
                    Vector3 moveV3 = Vector3(move.X(), move.Y(), move.Z());
                    Matrix4 rotMatrix = Matrix4::rotation(rotationV3, angle);
                    Vector3 newMove = rotMatrix * moveV3;
                    newMove = centerOfMassP3 - newMove;
                    i->second.SetPosition(Vector3DFloat(newMove[0], newMove[1], newMove[2]));
                }
            }
            return rotated;
        }

        bool CAlphaRenderer::SelectionMove(Vector3DFloat moveDirection) {
            bool moved = false;
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                if (i->second.GetSelected()) {
                    i->second.SetPosition(i->second.GetPosition() + moveDirection);
                    i->second.SetFlag(1);
                    moved = true;
                }
                else {
                    i->second.SetFlag(0);
                }
            }

            for (unsigned int i = 0; i < bonds.size(); i++) {
                if (bonds[i].GetSelected()) {
                    PDBAtom a = atoms[bonds[i].GetAtom0Ix()];
                    if (a.GetFlag() == 0) {
                        a.SetPosition(a.GetPosition() + moveDirection);
                        a.SetFlag(1);
                        moved = true;
                    }

                    a = atoms[bonds[i].GetAtom1Ix()];
                    if (a.GetFlag() == 0) {
                        a.SetPosition(a.GetPosition() + moveDirection);
                        a.SetFlag(1);
                        moved = true;
                    }
                }
            }
            UpdateBoundingBox();
            return moved;

        }


        bool CAlphaRenderer::SelectionClear() {
            if (Renderer::SelectionClear()) {
                for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                    i->second.SetSelected(false);
                }

                for (unsigned int i = 0; i < bonds.size(); i++) {
                    bonds[i].SetSelected(false);
                }
                for (unsigned int i = 0; i < aHelices.size(); i++) {
                    aHelices[i].selected = false;
                }
                for (unsigned int i = 0; i < bStrands.size(); ++i) {
                    bStrands[i].selected = false;
                }
                for (unsigned int i = 0; i < loops.size(); ++i) {
                    loops[i].selected = false;
                }
                selectedHelixIndices.clear();
                selectedSSEHelices.clear();
                selectedStrandIndices.clear();
                selectedLoopIndices.clear();
                return true;
            }
            return false;
        }

        void CAlphaRenderer::ClearOtherHighlights(){
            selectedSSEHelices.clear();
        }

        void CAlphaRenderer::SelectionToggle(int subsceneIndex, bool forceTrue, int ix0, int ix1, int ix2, int ix3, int ix4) {
            Renderer::SelectionToggle(subsceneIndex, forceTrue, ix0, ix1, ix2, ix3, ix4);
            AtomMapType::iterator it;
            PDBAtom * a;
            if (subsceneIndex == 0) {
                //TODO: possibly implement mouse picking using ray intersection
                switch (displayStyle) {
                case CALPHA_DISPLAY_STYLE_BACKBONE:
                case CALPHA_DISPLAY_STYLE_SIDE_CHAIN:
                    it = atoms.find(atomHashKeys.at(ix0));
                    if (it != atoms.end()) {
                        a = &(it->second);
                        a->SetSelected(forceTrue || !a->GetSelected());
                    }
                    break;
                case CALPHA_DISPLAY_STYLE_RIBBON:

                    if (aHelices[ix0].selected == true && !forceTrue){
                        aHelices[ix0].selected = false;
                    }
                    else{
                        cout << "Updating selectedHelix" << " ix0=" << ix0 << " forceTrue=" << forceTrue << endl;
                        aHelices[ix0].selected = true;
                        selectedHelixIndices.push_back(ix0);
                    }
                    break;
                }
            }
            else if ((subsceneIndex == 1) && (ix0 >= 0) && (ix0 <= (int)bonds.size())) {
                switch (displayStyle) {
                case CALPHA_DISPLAY_STYLE_BACKBONE:
                    bonds[ix0].SetSelected(forceTrue || !bonds[ix0].GetSelected());
                    break;
                case CALPHA_DISPLAY_STYLE_RIBBON:
                    if (bStrands[ix0].selected == true && !forceTrue){
                        bStrands[ix0].selected = false;
                    }
                    else {
                        cout << "Updating selectedStrand ix0=" << ix0 << " forceTrue=" << forceTrue << endl;
                        bStrands[ix0].selected = true;
                        selectedStrandIndices.push_back(ix0);
                    }
                    //cout << "A Ribbon was selected and subscene is 1" << endl;
                    break;
                case CALPHA_DISPLAY_STYLE_SIDE_CHAIN:
                    sidechainBonds[ix0].SetSelected(forceTrue || !sidechainBonds[ix0].GetSelected());
                    break;
                }
            }
            else if ((subsceneIndex == 2) && (ix0 != NULL)) {
                switch (displayStyle) {
                case CALPHA_DISPLAY_STYLE_BACKBONE:
                    break;
                case CALPHA_DISPLAY_STYLE_RIBBON:
                    if (loops[ix0].selected == true && !forceTrue){
                        loops[ix0].selected = false;
                    }
                    else {
                        cout << "Updating selectedStrand ix0=" << ix0 << " forceTrue=" << forceTrue << endl;
                        loops[ix0].selected = true;
                        selectedLoopIndices.push_back(ix0);
                    }
                    //cout << "A Ribbon was selected and subscene is 2" << endl;
                    break;
                case CALPHA_DISPLAY_STYLE_SIDE_CHAIN:
                    break;
                }
            }
            //cout << "Finished updating selected calpha helix" << endl;
            cout << "Finished updating selected secel" << endl;
        }

        void CAlphaRenderer::Unload() {
            atoms.clear();
            bonds.clear();
            sidechainBonds.clear();
            UpdateBoundingBox();
        }

        void CAlphaRenderer::UpdateBoundingBox() {
            if (atoms.size() > 0) {
                for (int i = 0; i < 3; i++) {
                    minPts[i] = atoms.begin()->second.GetPosition().values[i];
                    maxPts[i] = atoms.begin()->second.GetPosition().values[i];
                }

                for (AtomMapType::iterator j = atoms.begin(); j != atoms.end(); j++) {
                    for (int i = 0; i < 3; i++) {
                        minPts[i] = min(minPts[i], j->second.GetPosition().values[i]);
                        maxPts[i] = max(maxPts[i], j->second.GetPosition().values[i]);
                    }
                }
            }
            else {
                for (int i = 0; i < 3; i++) {
                    minPts[i] = -0.5;
                    maxPts[i] = 0.5;
                }
            }
        }


        string CAlphaRenderer::GetSupportedLoadFileFormats() {
            return "Atom Positions (*.pdb)";
        }

        string CAlphaRenderer::GetSupportedSaveFileFormats() {
            return "Atom Positions (*.atom)";
        }
        PDBAtom * CAlphaRenderer::GetAtom(unsigned long long index) {
            return &atoms[index];
        }

        PDBBond * CAlphaRenderer::GetBond(int index) {
            return &bonds[index];
        }

        PDBBond * CAlphaRenderer::GetSideChainBond(int index) {
            return &sidechainBonds[index];
        }

        PDBAtom * CAlphaRenderer::GetSelectedAtom(unsigned int selectionId) {
            //Python uses this with SelectionAtomCount() to get all the selected atoms
            int count = 0;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSelected()) {
                    if (count == selectionId) {
                        return &it->second;
                    }
                    count++;
                }
            }
            return NULL;
        }

        std::vector<unsigned long long> CAlphaRenderer::GetAtomHashes() {
            std::vector<unsigned long long> atomHashes;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                atomHashes.push_back(it->first);
            }
            return atomHashes;
        }

        int CAlphaRenderer::GetBondIndex(unsigned long long atom0, unsigned long long atom1) {
            for (unsigned int i = 0; i < bonds.size(); i++) {
                if (((bonds[i].GetAtom0Ix() == atom0) && (bonds[i].GetAtom1Ix() == atom1)) ||
                    ((bonds[i].GetAtom0Ix() == atom1) && (bonds[i].GetAtom1Ix() == atom0))) {
                    return i;
                }
            }
            return -1;
        }

        int CAlphaRenderer::GetSideChainBondIndex(unsigned long long atom0, unsigned long long atom1) {
            for (unsigned int i = 0; i < sidechainBonds.size(); i++) {
                if (((sidechainBonds[i].GetAtom0Ix() == atom0) && (sidechainBonds[i].GetAtom1Ix() == atom1)) ||
                    ((sidechainBonds[i].GetAtom0Ix() == atom1) && (sidechainBonds[i].GetAtom1Ix() == atom0))) {
                    return i;
                }
            }
            return -1;
        }

        int CAlphaRenderer::GetAtomCount() {
            return atoms.size();
        }

        int CAlphaRenderer::GetBondCount() {
            return bonds.size();
        }

        int CAlphaRenderer::GetSideChainBondCount() {
            return sidechainBonds.size();
        }

        void CAlphaRenderer::DeleteAtom(unsigned long long index) {
            /**
            for(auto it = begin(atoms); it != end(atoms);)
            {
            if (it->first == index)
            {
            it = atoms.erase(it);
            }
            else
            {
            ++it;
            }
            }
            **/
            atoms.erase(atoms.find(index));
        }

        void CAlphaRenderer::DeleteBond(int index) {
            bonds.erase(bonds.begin() + index);
        }

        void CAlphaRenderer::DeleteSideChainBond(int index) {
            sidechainBonds.erase(sidechainBonds.begin() + index);
        }


        Vector3DFloat CAlphaRenderer::Get3DCoordinates(int subsceneIndex, int ix0, int ix1, int ix2, int ix3, int ix4) {
            Vector3DFloat position;
            switch (subsceneIndex) {
            case(0) :
                if ((ix0 >= 0) && (ix0 <= (int)atoms.size())) {
                    PDBAtom * a = &(atoms[ix0]);
                    position = a->GetPosition();
                }
                break;
            case(1) :
                if ((ix0 >= 0) && (ix0 <= (int)bonds.size())) {
                    position = (atoms[bonds[ix0].GetAtom0Ix()].GetPosition() + atoms[bonds[ix0].GetAtom1Ix()].GetPosition()) * 0.5;
                }
                break;
            default:
                position = Vector3DFloat(0, 0, 0);
                break;
            }
            return position;
        }

        void CAlphaRenderer::TransformAllAtomLocations(MatrixFloat transform) {
            for (AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
                i->second.Transform(transform);
            }
        }

        int CAlphaRenderer::StartHelix() {
            aHelices.push_back(Secel());
            return aHelices.size() - 1;
        }

        void CAlphaRenderer::AddHelixElement(int index, unsigned long long hashKey){
            aHelices[index].atomHashes.push_back(hashKey);
        }

        int CAlphaRenderer::StartStrand() {
            bStrands.push_back(Secel());
            return bStrands.size() - 1;
        }

        void CAlphaRenderer::AddStrandElement(int index, unsigned long long hashKey){
            bStrands[index].atomHashes.push_back(hashKey);
        }

        int CAlphaRenderer::StartLoop() {
            loops.push_back(Secel());
            return loops.size() - 1;
        }

        void CAlphaRenderer::AddLoopElement(int index, unsigned long long hashKey){
            loops[index].atomHashes.push_back(hashKey);
        }

        bool CAlphaRenderer::CleanSecondaryStructures(){
            aHelices.clear();
            bStrands.clear();
            loops.clear();
            return true;
        }

        void CAlphaRenderer::SetHelixCorrs(std::vector < int > flatCorrespondences){
            if (flatCorrespondences.size() % 2 != 0)
                return;
            else
                corrs.clear();
            for (int i = 0; i < flatCorrespondences.size(); i = i + 2){
                corrs.push_back(boost::tuple<int, int>(flatCorrespondences[i], flatCorrespondences[i + 1]));
            }
        }

        void CAlphaRenderer::SetFeatureVecs(std::vector<Vector3DFloat> flatFeatureVecs){
            if (flatFeatureVecs.size() % 2 != 0)
                return;
            else
                featureVecs.clear();
            for (int i = 0; i < flatFeatureVecs.size(); i = i + 2){
                featureVecs.push_back(boost::tuple<Vector3DFloat, Vector3DFloat>(flatFeatureVecs[i], flatFeatureVecs[i + 1]));
            }

        }
        void CAlphaRenderer::SetSelectedSSEHelices(std::vector<int> indices){
            selectedSSEHelices.clear();
            selectedSSEHelices = indices;
        }

        void CAlphaRenderer::SetHelixColor(int helixNum, float r, float g, float b){
            cout << "setting helix color " << helixNum << " to (" << r << ", " << g << ", " << b << ")" << endl;
            helixColors.erase(helixNum);
            helixColors.insert(pair<int, boost::tuple<float, float, float> >(helixNum, boost::tuple<float, float, float>(r, g, b)));
        }

        void CAlphaRenderer::ClearHelixColors() {
            helixColors.clear();
        }

        // creates a std::vector of Vector3DFloats that represents the locations of all the PDBAtoms
        // starting with start and ending with end; it does not error check, so incorrectly
        // ordered points will break this method.  there are more efficient ways to handle this
        // functionality, but this seems simple and flexible enough
        std::vector<Vector3DFloat> CAlphaRenderer::CreatePointVector(PDBAtom start, PDBAtom end){
            std::vector<Vector3DFloat> points;

            PDBAtom current = start;
            while (current.GetHashKey() != end.GetHashKey()){
                points.push_back(current.GetPosition());
                if (current.GetHashKey() == current.GetNextCAHash()){
                    break;
                }
                current = atoms.find(current.GetNextCAHash())->second;
            }

            points.push_back(end.GetPosition());
            return points;
        }

        // implementation of Laplacian smoothing for a std::vector of Vector3DFloats (treats them like points)
        // creating copies of "points" twice seems unnecessary, but I am unsure about the performance cost,
        // so I am leaving it for simplicity of implementation
        std::vector<Vector3DFloat> CAlphaRenderer::LaplacianSmoothing(std::vector<Vector3DFloat> points, int steps){
            std::vector<Vector3DFloat> pointsTemp(points);
            std::vector<Vector3DFloat> smoothedPoints(points);

            for (int i = 0; i < steps; ++i){
                for (int j = 1; j < points.size() - 1; ++j){
                    smoothedPoints[j] = (pointsTemp[j - 1] + pointsTemp[j + 1])*.5;
                    smoothedPoints[j] = (smoothedPoints[j] + pointsTemp[j])*.5;
                }
                pointsTemp = smoothedPoints;
            }
            return pointsTemp;
        }

        // unsure of what behavior should be if points.size() < 3; in molscript the strand is skipped in this case
        std::vector<Vector3DFloat> CAlphaRenderer::CreateStrandNormals(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next){
            std::vector<Vector3DFloat> normals(points);
            int ptsSize = points.size();

            for (int i = 1, length = ptsSize - 1; i < length; ++i){
                Vector3DFloat newPos = (points[i - 1] + points[i + 1])*.5;
                normals[i] = points[i] - newPos;
                normals[i].Normalize();
            }

            normals[0] = (points[1] + previous)*.5 - points[0];
            if ((points[0] - previous).Length() < .0001){
                normals[0] = normals[1];
            }

            normals[ptsSize - 1] = (points[ptsSize - 2] + next)*.5 - points[ptsSize - 1];
            if ((points[ptsSize - 2] - next).Length() < .0001){
                normals[ptsSize - 1] = normals[ptsSize - 2];
            }

            // "normals must point the same way" - molscript/graphics.c
            for (int j = 0, size = ptsSize - 1; j < size; ++j){
                if (normals[j] * normals[j + 1] < 0){
                    normals[j + 1] = normals[j + 1] * -1;
                }
            }

            // "smooth normals, one iteration" - molscript/graphics.c
            std::vector<Vector3DFloat> smoothedNormals(normals);

            for (int k = 1, size = ptsSize - 1; k < size; ++k){
                smoothedNormals[k] = normals[k - 1] + normals[k] + normals[k + 1];
                smoothedNormals[k].Normalize();
            }

            // "normals exactly perpendicular to strand" - molscript/graphics.c
            Vector3DFloat direction = points[1] - points[0];
            Vector3DFloat side = direction^smoothedNormals[0];
            smoothedNormals[0] = side ^ direction;
            smoothedNormals[0].Normalize();

            for (int i = 1, size = ptsSize - 1; i < size; ++i){
                direction = points[i + 1] - points[i - 1];
                side = direction^smoothedNormals[i];
                smoothedNormals[i] = side^direction;
                smoothedNormals[i].Normalize();
            }

            direction = points[ptsSize - 1] - points[ptsSize - 2];
            side = direction^smoothedNormals[ptsSize - 1];
            smoothedNormals[ptsSize - 1] = side^direction;
            smoothedNormals[ptsSize - 1].Normalize();
            return smoothedNormals;
        }

        void CAlphaRenderer::CreateHelixAxesTangentsAndPoints(std::vector<Vector3DFloat>& axes, std::vector<Vector3DFloat>& tangents, std::vector<Vector3DFloat>& interpPoints, std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, double HELIX_ALPHA, double HELIX_BETA, double HELIX_HERMITE_FACTOR){
            if (points.size() > 2){

                for (int i = 0; i < points.size() - 1; ++i){

                    if (i > 0){
                        Vector3DFloat cvec = points[i + 1] - points[i - 1];
                        cvec.Normalize();

                        Vector3DFloat rvec = (points[i] - points[i - 1]) ^ (points[i + 1] - points[i]);
                        rvec.Normalize();

                        axes[i] = rvec*sin(HELIX_ALPHA) + cvec*cos(HELIX_ALPHA);
                        tangents[i] = rvec*sin(HELIX_BETA) + cvec*cos(HELIX_BETA);
                        tangents[i] = tangents[i] * HELIX_HERMITE_FACTOR;
                    }
                }
                axes[0] = axes[1];
                axes[axes.size() - 1] = axes[axes.size() - 2];

                tangents[0] = previous - points[1];
                tangents[0].Normalize();
                tangents[0] = tangents[0] * HELIX_HERMITE_FACTOR;
                tangents[tangents.size() - 1] = next - points[points.size() - 2];
                tangents[tangents.size() - 1].Normalize();
                tangents[tangents.size() - 1] = tangents[tangents.size() - 1] * HELIX_HERMITE_FACTOR;
            }
        }

        // method works like drawing the side of a cylinder with only one stack and 4 slices
        void CAlphaRenderer::DrawOpenBox(std::vector<Vector3DFloat> points, std::vector<Vector3DFloat> normals){
            glBegin(GL_TRIANGLE_STRIP);

            for (int j = 0, runlength = points.size() + 2; j < runlength; ++j){
                glNormal3f(normals[j%normals.size()].X(), normals[j%normals.size()].Y(), normals[j%normals.size()].Z());
                glVertex3f(points[j%points.size()].X(), points[j%points.size()].Y(), points[j%points.size()].Z());
            }

            glEnd();
        }



        // renders a set of points and normals assuming that they are laid out like the side of a cylinder's points and normals
        void CAlphaRenderer::DrawTube(std::vector<Vector3DFloat> points, std::vector<Vector3DFloat> normals, int stacks, int slices){
            //glLightModeli ( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
            //glDisable(GL_CULL_FACE);
            switch (renderingType){
            case 0:
                glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR | GL_AMBIENT_AND_DIFFUSE);
                glEnable(GL_COLOR_MATERIAL);
                break;
            default:
                break;
            }
            for (int i = 0, runlength = points.size(); i < runlength; ++i){
                if (i % (stacks + 1) == 0){
                    glBegin(GL_TRIANGLE_STRIP);
                }

                int nextSliceIx = (i + stacks + 1) % runlength;


                glNormal3f(normals[i].X(), normals[i].Y(), normals[i].Z());
                glVertex3f(points[i].X(), points[i].Y(), points[i].Z());

                glNormal3f(normals[nextSliceIx].X(), normals[nextSliceIx].Y(), normals[nextSliceIx].Z());
                glVertex3f(points[nextSliceIx].X(), points[nextSliceIx].Y(), points[nextSliceIx].Z());

                if ((i + 1) % (stacks + 1) == 0){
                    glEnd();
                }
            }
        }

        std::vector<Vector3DFloat> CAlphaRenderer::InterpolateLoopPoints(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS){
            HermiteCurve curve;
            Vector3DFloat m0, m1;
            std::vector<Vector3DFloat> pointstemp(points);
            bool LAPLACIAN_SMOOTHING = true;
            int SMOOTHING_STEPS = 1;
            double HERMITE_FACTOR = 0.5;
            int LOOP_SLICES = 10;
            if (LAPLACIAN_SMOOTHING){
                pointstemp = LaplacianSmoothing(points, SMOOTHING_STEPS);
            }

            std::vector<Vector3DFloat> interpolatedPoints((pointstemp.size() - 1)*(NUM_SEGMENTS));

            for (unsigned int i = 0; i < points.size() - 1; ++i){
                if (i == 0){
                    m0 = pointstemp[i + 1] - previous;
                }
                else {
                    m0 = pointstemp[i + 1] - pointstemp[i - 1];
                    m0 = m0*HERMITE_FACTOR;
                }

                if (i + 2 > pointstemp.size() - 1){
                    m1 = next - pointstemp[i];
                }
                else {
                    m1 = pointstemp[i + 2] - pointstemp[i];
                    m1 = m1*HERMITE_FACTOR;
                }

                curve.setCurve(pointstemp[i], pointstemp[i + 1], m0, m1);
                interpolatedPoints[i*(NUM_SEGMENTS)] = pointstemp[i];
                for (int sect = 1; sect < NUM_SEGMENTS; ++sect){
                    double tsect = ((double)sect) / ((double)NUM_SEGMENTS);
                    interpolatedPoints[i*(NUM_SEGMENTS)+sect] = curve.getPos(tsect);
                }
            }
            interpolatedPoints[interpolatedPoints.size() - 1] = points[points.size() - 1];
            return interpolatedPoints;
        }

        void CAlphaRenderer::SetNumSegments(int segments){
            NUM_SEGMENTS = segments;
        }

        void CAlphaRenderer::SetNumSlices(int slices){
            NUM_SLICES = slices;
        }

        void CAlphaRenderer::SetHltRValue(int col){
            //hlt_r = ((double)col)/100.0;
            hlt_r = ((double)col) / 1000.0;
            //thinRibbThickness = hlt_r;
            cout << "hlt_r: " << hlt_r << endl;
        }

        void CAlphaRenderer::SetHltGValue(int col){
            hlt_g = ((double)col) / 100.0;
            cout << "hlt_g: " << hlt_g << endl;
        }

        void CAlphaRenderer::SetHltBValue(int col){
            hlt_b = ((double)col) / 100.0;
            cout << "hlt_b: " << hlt_b << endl;
        }

        void CAlphaRenderer::SetHltAValue(int col){
            hlt_a = ((double)col) / 100.0;
            cout << "hlt_a: " << hlt_a << endl;
        }


        std::vector<PDBBond> CAlphaRenderer::getDeletedBonds() {
            return bondsToDelete;
        }

        string CAlphaRenderer::getDeletedBondAtoms() {

            ofstream myfile;
            myfile.open("noBondConstraints");

            string atomsToNotBond;
            int atomIndex = 0;
            unsigned int atom1;
            unsigned int atom2;
            for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
                if (it->second.GetSelected()) {
                    atomsToNotBond += std::to_string(it->second.GetResSeq());
                    if (atomIndex % 2 == 0) {
                        atom1 = it->second.GetResSeq();
                        atomsToNotBond += " ";
                    }
                    else {
                        atom2 = it->second.GetResSeq();
                        atomsToNotBond += " ";

                    }
                    atomIndex += 1;
                }

            }
            myfile << atomsToNotBond;
            myfile.close();
            return atomsToNotBond;
        }

        std::vector<unsigned long long> CAlphaRenderer::getDeletedBonds1Ix() {
            return ix1s;
        }

        /**
        void CAlphaRenderer::RemoveSelectedBonds(string nobonds) {
        string sLine = nobonds;
        std::stringstream ss(sLine);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> strs(begin, end);
        std::vector<unsigned int> atomNums;
        for(int i = 0; i < strs.size(); i++) {
        unsigned int currentAtomNum = (unsigned int)std::stoi(strs[i]);

        atomNums.push_back(currentAtomNum);
        }
        std::vector<unsigned long long> atomHashes;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
        for(int i = 0; i < atomNums.size(); i++) {
        unsigned int currentAtomNum = atomNums[i];
        if(it->second.GetResSeq() == currentAtomNum) {
        atomHashes.push_back(it->second.GetHashKey());
        }
        }
        }
        for(int i = 0; i < atomHashes.size()-1; i = i+2) {
        int currentBondIndex = GetBondIndex(atomHashes[i], atomHashes[i+1]);
        if (currentBondIndex != -1) {
        DeleteBond(currentBondIndex);
        }
        }

        }




        void CAlphaRenderer::addSelectedBonds(string newBonds) {

        string sLine = newBonds;
        std::stringstream ss(sLine);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> strs(begin, end);
        std::vector<unsigned int> atomNums;
        for(int i = 0; i < strs.size(); i++) {
        unsigned int currentAtomNum = (unsigned int)std::stoi(strs[i]);

        atomNums.push_back(currentAtomNum);
        }
        std::vector<unsigned long long> atomHashes;
        for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
        for(int i = 0; i < atomNums.size(); i++) {
        unsigned int currentAtomNum = atomNums[i];
        if(it->second.GetResSeq() == currentAtomNum) {
        atomHashes.push_back(it->second.GetHashKey());
        }
        }
        }
        for(int i = 0; i < atomHashes.size()-1; i++) {
        AddBond(PDBBond(atomHashes[i], atomHashes[i+1], true));
        }

        }
        **/




    }

}

#endif
