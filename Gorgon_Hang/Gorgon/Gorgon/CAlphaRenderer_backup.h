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

			double cp0 = 2*tcubed - 3*tsquared + 1;
			double cm0 = tcubed - 2*tsquared + t;
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

		class CAlphaRenderer : public Renderer{
		public:
			struct Secel{
				vector<unsigned long long> atomHashes;
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

			// Controlling the atom vector
			PDBAtom * AddAtom(PDBAtom atom);
			PDBAtom * GetAtom(unsigned long long index);
			PDBAtom * GetAtomFromHitStack(int subsceneIndex, bool forceTrue, int ix0, int ix1, int ix2, int ix3, int ix4);
			PDBAtom * GetSelectedAtom(unsigned int selectionId);
			void DeleteAtom(unsigned long long index);
			int GetAtomCount();
			vector<unsigned long long> GetAtomHashes();

			//Controlling the bond vector
			void AddBond(PDBBond bond);
			PDBBond * GetBond(int index);
			int GetBondIndex(unsigned long long atom0, unsigned long long atom1);
			void DeleteBond(int index);
			int GetBondCount();

			//Controlling the bond vector
			void AddSideChainBond(PDBBond bond);
			PDBBond * GetSideChainBond(int index);
			int GetSideChainBondIndex(unsigned long long atom0, unsigned long long atom1);
			void DeleteSideChainBond(int index);
			int GetSideChainBondCount();

			void SetNumSegments(int segments);
			void SetNumSlices(int slices);

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
			vector<Vector3DFloat> CreatePointVector(PDBAtom first, PDBAtom last); // functionality mirrored in previously implemented method,
			// will try to refactor
			vector<Vector3DFloat> LaplacianSmoothing(vector<Vector3DFloat> points, int steps); // applies Laplacian smoothing to a vector of
			// Vector3DFloats
			vector<Vector3DFloat> CreateStrandNormals(vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next); // create line segment normals to be used in drawing Beta
			// strands
			void CreateHelixAxesTangentsAndPoints(vector<Vector3DFloat>& axes, vector<Vector3DFloat>& tangents, vector<Vector3DFloat>& interpPoints, vector<Vector3DFloat> points, 
				Vector3DFloat previous, Vector3DFloat next, double HELIX_ALPHA, double HELIX_BETA, double HELIX_HERMITE_FACTOR);
			void DrawOpenBox(vector<Vector3DFloat> points, vector<Vector3DFloat> normals); // takes a vector of 8 points and draws a rectangular prism with two of its six sides not
			// filled in; the first 4 points are from the beggining edge of the box, with the second four
			// forming the end
			void DrawTube(vector<Vector3DFloat> points, vector<Vector3DFloat> normals, int stacks, int slices);
			vector<Vector3DFloat> InterpolateLoopPoints(vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS); // creates interpolated points for loops
			//vector<Vector3DFloat> InterpolateStrandPoints(vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS);
			//vector<Vector3DFloat> InterpolateHelixPoints(vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS);

			// for testing purposes only; allow changing of highlight color values
			void SetHltRValue(int col);
			void SetHltGValue(int col);
			void SetHltBValue(int col);
			void SetHltAValue(int col);


			vector<int> GetSelectedHelixIndices();
			void SetHelixCorrs( vector < int > flatCorrespondences);
			void SetSelectedSSEHelices(vector<int>);
			void ClearOtherHighlights();
			void SetFeatureVecs(vector<Vector3DFloat> flatFeatureVecs);
			void SetHelixColor(int helixNum, float r, float g, float b);
            void ClearHelixColors();
		private:
			void DrawBackboneModel(int subSceneIndex, bool selectEnabled);
			void DrawRibbonModel(int subSceneIndex, bool selectEnabled);
			void DrawSideChainModel(int subSceneIndex, bool selectEnabled);
		private:
			AtomMapType atoms;

			//TODO: possibly implement mouse picking using ray intersection
			vector<unsigned long long> atomHashKeys; //glLoadName(index of this vector)... used for selection

			vector<PDBBond> bonds;
			vector<PDBBond> sidechainBonds;

			vector<Secel> aHelices;
			vector<Secel> bStrands;
			vector<Secel> loops;

			vector<int> selectedHelixIndices;
			//vector<int> selectedSecelIndices; //unsure if I can just keep track of secels as one structure or not
			vector<int> selectedStrandIndices;
			vector<int> selectedLoopIndices;
			vector < boost::tuple<int, int> > corrs;
			vector<int> selectedSSEHelices;
			vector< boost::tuple<Vector3DFloat, Vector3DFloat> > featureVecs;

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
			HELIX_ALPHA = 32.0 * PI/180.0; 
			HELIX_BETA = -11.0 * PI/180.0; // these three values taken from molscript code
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

		vector<int> CAlphaRenderer::GetSelectedHelixIndices(){
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

		void CAlphaRenderer::AddSideChainBond(PDBBond bond) {
			sidechainBonds.push_back(bond);
		}

		void CAlphaRenderer::DrawBackboneModel(int subSceneIndex, bool selectEnabled) {
			GLfloat emissionColor[4] = {1.0, 1.0, 1.0, 1.0};

			if(subSceneIndex == 0) { // Drawing Atoms
				if(selectEnabled) {
					atomHashKeys.clear();
					glPushName(0);
					glPushName(0);
				}
				for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
					if(it->second.GetName() == "CA") {
						glPushAttrib(GL_LIGHTING_BIT);
						if(it->second.GetSelected()) {
							glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
							glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
						} else {
							OpenGLUtils::SetColor(it->second.GetColorR(), it->second.GetColorG(), it->second.GetColorB(), it->second.GetColorA());
						}

						if(selectEnabled){
							//TODO: possibly implement mouse picking using ray intersection
							atomHashKeys.push_back(it->first); // adding the atom hash key as an element
							glLoadName(static_cast<GLuint>( atomHashKeys.size() - 1)); // the index of the element just added
						}
						if(it->second.GetVisible()) {
							DrawSphere(it->second.GetPosition(), it->second.GetAtomRadius() * 0.3);
						}

						glPopAttrib();
					}

				}
				if(selectEnabled) {
					glPopName();
					glPopName();
				}
			} else if(subSceneIndex == 1) { // Drawing Bonds
				if(selectEnabled) {
					glPushName(1);
					glPushName(0);
				}
				for(int i=0; i < (int)bonds.size(); i++) {
					glPushAttrib(GL_LIGHTING_BIT);
					if(bonds[i].GetSelected()) {
						glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
						glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
					}

					if(selectEnabled){
						glLoadName(i);
					}
					float length = (atoms[bonds[i].GetAtom0Ix()].GetPosition() - atoms[bonds[i].GetAtom1Ix()].GetPosition()).Length();
					if(length > 4.2) {
						OpenGLUtils::SetColor(1.0, 0, 0, 1.0);
					}

					if(length < 3.3) {
						OpenGLUtils::SetColor(0, 0, 1.0, 1.0);
					}

					if(atoms[bonds[i].GetAtom0Ix()].GetVisible() && atoms[bonds[i].GetAtom1Ix()].GetVisible()) {
						DrawCylinder(atoms[bonds[i].GetAtom0Ix()].GetPosition(), atoms[bonds[i].GetAtom1Ix()].GetPosition(), 0.1, 10, 2);
					}
					glPopAttrib();
				}
				if(selectEnabled) {
					glPopName();
					glPopName();
				}
			} else if(subSceneIndex == 2) { // Drawing spheres to cover up the cylinder edges				
				for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
					if(i->second.GetName() == "CA") {
						DrawSphere(i->second.GetPosition(), 0.1);
					}
				}
			}
		}

		void CAlphaRenderer::DrawRibbonModel(int subSceneIndex, bool selectEnabled) {
			if(selectEnabled) {
				glPushName(subSceneIndex);
				glPushName(0);
			}
			//GLfloat emissionColor[4] = {1.0, 1.0, 1.0, 1.0};
			//GLfloat frontColor[4] = {1.0, 0.0, 0.0, 1.0};
			//GLfloat backColor[4] = {0.0, 0.0, 1.0, 1.0};
			GLfloat emissionColor[4] = {hlt_r, hlt_g, hlt_b, hlt_a};
			GLfloat frontColor[4] = {1.0, 0.0, 0.0, 1.0};
			GLfloat backColor[4] = {0.0, 0.0, 1.0, 1.0};
			vector<int> PDBIndices;

			switch(subSceneIndex) {
			case 0: // Helices

				//if(subSceneIndex == 0){
				for(unsigned int i = 0; i < corrs.size(); ++i){
					int SSEIndex = get<1> (corrs[i]);
					for(unsigned int k = 0; k < selectedSSEHelices.size(); ++k){
						if(selectedSSEHelices[k] == SSEIndex){
							PDBIndices.push_back( get<0>( corrs[i]) );
						}
					}
				}

				for(int i = 0; i < aHelices.size(); i++) {
					if(selectEnabled){
						glLoadName(i);
					}

					glPushAttrib(GL_LIGHTING_BIT);

					Secel currentSecel = aHelices[i];

					if(currentSecel.selected == true){
						glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
						glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
					}
					map<int, boost::tuple<float, float, float> >::iterator iter = helixColors.begin();
					iter = helixColors.find(i);
					if(iter != helixColors.end()){

						OpenGLUtils::SetColor(get<0>(helixColors[i]), get<1>(helixColors[i]), get<2>(helixColors[i]), 1.0);
					}else{
						//OpenGLUtils::SetColor(0.8,0.8,0.8,1.0);
					}

					if(currentSecel.atomHashes.size() > 0){

						PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
						PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size()-1])->second;
						Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
						Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

						vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
						vector<Vector3DFloat> tangents = vector<Vector3DFloat>(points);
						vector<Vector3DFloat> axes = vector<Vector3DFloat>(points);
						vector<Vector3DFloat> interpPoints = vector<Vector3DFloat>((points.size()-1)*NUM_SEGMENTS + 1);
						int flatSlices = 2;
						int rptsize = interpPoints.size()*4;
						switch (renderingType){
							case 0:
								rptsize = interpPoints.size()*flatSlices;
								break;
							case 1:
								rptsize = interpPoints.size()*4;
								break;
							default:
								rptsize = interpPoints.size()*4;
								break;
						}
						vector<Vector3DFloat> renderingPoints(rptsize);
						vector<Vector3DFloat> renderingNormals(renderingPoints.size());
						/*vector<Vector3DFloat> renderingPoints(interpPoints.size()*NUM_SLICES);
						vector<Vector3DFloat> renderingNormals(renderingPoints.size());*/

						HermiteCurve curve;
						Vector3DFloat m0, m1;

						CreateHelixAxesTangentsAndPoints(axes, tangents, interpPoints, points, preSecelAtomPos, postSecelAtomPos, HELIX_ALPHA, HELIX_BETA, HELIX_HERMITE_FACTOR);

						for(unsigned int x = 0; x < points.size()-1; ++x){

							m0 = tangents[x];
							m1 = tangents[x+1];

							curve.setCurve(points[x], points[x+1], m0, m1);

							// used in rendering a helix as glowing if selected
							vector<Vector3DFloat> selectedBoxPositions(8);

							float halfwidth = HELIX_WIDTH/2.0;
							float halfthickness = LOOP_RADIUS;
							Vector3DFloat lastPos = points[x];
							//int NUM_SECTIONS = 10;
							for (int sect = 0; sect <= NUM_SEGMENTS; ++sect){
								if(sect == 0 && x != 0){
									continue;
								}
								double tsect = ((double)sect)/((double)NUM_SEGMENTS);
								Vector3DFloat nextPos = curve.getPos(tsect);

								Vector3DFloat currentAxis = axes[x]*(1.0-tsect) + axes[x+1]*tsect;
								currentAxis.Normalize();

								Vector3DFloat curnormal = curve.getTangent(tsect);
								curnormal = curnormal^currentAxis;
								curnormal.Normalize();

								if(x == 0){
									halfwidth = LOOP_RADIUS + (0.5 * HELIX_WIDTH - LOOP_RADIUS) * 0.5 * (-1*cos(PI*tsect) + 1.0);
								} else if (x == points.size() - 2){
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
											renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*(q%3 == 0 ? -1 : 1); //+ curnormal*hlt_r*pow(-1.0, q/2);
											renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal;//*(q%3 == 0 ? -1 : 1);
											renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()].Normalize();
										}
										break;
									case 1:
										for(int q = 0; q < 4; ++q){
											renderingPoints[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = nextPos + currentAxis*halfwidth*pow(-1.0, q/2) 
												+ curnormal*halfthickness*(q%3 == 0 ? -1 : 1);
											renderingNormals[(x*NUM_SEGMENTS + sect) + q*interpPoints.size()] = curnormal*(q%3 == 0 ? -1 : 1);
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

					if(aHelices[i].selected == true){
						glPushAttrib(GL_LIGHTING_BIT);

						if(featureVecs.size() > 0){
							OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
							DrawSphere(featureVecs[i].get<0>(), 1.0);
							OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
							DrawSphere(featureVecs[i].get<1>(), 1.0);
						}else{
							OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
							DrawSphere(atoms[aHelices[i].atomHashes[0]].GetPosition(), 1.0);
							OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
							DrawSphere(atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size()-1]].GetPosition(), 1.0);
						}

						glPopAttrib();
						Vector3DFloat pos1 = atoms[aHelices[i].atomHashes[0]].GetPosition();
						Vector3DFloat pos2 = atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size()-1]].GetPosition();
						//printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i+1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

						fflush(stdout);
					}

					for(unsigned int j = 0; j < PDBIndices.size(); ++j){
						if(PDBIndices[j] == i){
							glPushAttrib(GL_LIGHTING_BIT);
							OpenGLUtils::SetColor(1.0, 0.0, 0.0, 1.0);
							DrawSphere(atoms[aHelices[i].atomHashes[0]].GetPosition(), 1.0);
							OpenGLUtils::SetColor(0.0, 0.0, 1.0, 1.0);
							DrawSphere(atoms[aHelices[i].atomHashes[aHelices[i].atomHashes.size()-1]].GetPosition(), 1.0);
							glPopAttrib();								
						}
					}
				}
				//}
				break;
			case 1: // Strands
				for (unsigned int i = 0; i < bStrands.size(); ++i){
					if(selectEnabled){
						glLoadName(i);
					}

					int atom_counter = 0;
					PDBAtom lastEnd;

					HermiteCurve curve;
					Vector3DFloat m0, m1, dir1, dir2;

					glPushAttrib(GL_LIGHTING_BIT);

					Secel currentSecel = bStrands[i];

					if(currentSecel.selected == true){
						glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
						glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
					}

					if(currentSecel.atomHashes.size() > 0){
						PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
						PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size()-1])->second;
						Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
						Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

						vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
						int num_interp_points = (points.size() - 1)*NUM_SEGMENTS + 1;
						int num_rendering_points = num_interp_points*4;
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
						vector<Vector3DFloat> renderingPoints(num_rendering_points);
						vector<Vector3DFloat> normals = CreateStrandNormals(points, preSecelAtomPos, postSecelAtomPos);
						vector<Vector3DFloat> renderingNormals(renderingPoints);
						double arrowhead_factor = 1.0;
						//vector<Vector3DFloat> boxpositions(8);
						//vector<Vector3DFloat> boxnormals(8);

						bool LAPLACIAN_SMOOTHING = true;
						int SMOOTHING_STEPS = 1;
						if(LAPLACIAN_SMOOTHING){
							points = LaplacianSmoothing(points, SMOOTHING_STEPS);
						}

						for(unsigned int i = 0; i < points.size()-1; ++i){
							if(i == 0){
								m0 = points[i+1] - points[i];
							} else {
								if(i + 2 < points.size()){
									m0 = points[i+2] - points[i];
								} else {
									m0 = postSecelAtomPos - points[i];
								}
							}

							if(i + 3 >= points.size()){
								m1 = postSecelAtomPos - points[i];
							} else {
								m1 = points[i+3] - points[i+1];
							}

							m0 = m0*STRAND_HERMITE_FACTOR;
							m1 = m1*STRAND_HERMITE_FACTOR;

							if(i == 0){
								dir2 = points[1] - points[0];
								dir2.Normalize();
							}

							dir1 = dir2;
							if(i + 2 < points.size()){
								dir2 = points[i+2] - points[i];
							} else {
								dir2 = postSecelAtomPos - points[i];
							}

							dir2.Normalize();

							curve.setCurve(points[i], points[i+1], m0, m1);

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
								
								double tsect = ((double)sect)/((double)NUM_SEGMENTS);
								if(i > points.size() - 3){
									arrowhead_factor = 2.0*(1-tsect);
								}
								direction = dir1*(1.0 - tsect) + dir2*tsect;
								currentNormal = normals[i]*(1.0 - tsect) + normals[i+1]*(tsect);
								side = currentNormal^direction;
								side.Normalize();
								Vector3DFloat nextPos = curve.getPos(tsect);

								switch(renderingType){
										case 0:
											for(int q = 0; q < 2; ++q){
												renderingPoints[i*NUM_SEGMENTS + sect + q*num_interp_points] = nextPos + side*(0.5*WIDTH*arrowhead_factor + LOOP_RADIUS/2.0)*(q%3 == 0 ? -1 : 1);
												renderingNormals[i*NUM_SEGMENTS + sect + q*num_interp_points] = currentNormal;//*pow(-1.0, q/2);
											}
											break;
										case 1:
											for(int q = 0; q < 4; ++q){
												renderingPoints[i*NUM_SEGMENTS + sect + q*num_interp_points] = nextPos + side*(0.5*WIDTH*arrowhead_factor + LOOP_RADIUS/2.0)*(q%3 == 0 ? -1 : 1) 
													+ currentNormal*0.5*THICKNESS*pow(-1.0, q/2);
												renderingNormals[i*NUM_SEGMENTS + sect + q*num_interp_points] = currentNormal*pow(-1.0, q/2);
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

					if(currentSecel.selected == true){
						glPushAttrib(GL_LIGHTING_BIT);

						Vector3DFloat pos1 = atoms[currentSecel.atomHashes[0]].GetPosition();
						Vector3DFloat pos2 = atoms[currentSecel.atomHashes[currentSecel.atomHashes.size()-1]].GetPosition();
						//printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i+1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

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
					if(selectEnabled){
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

					if(currentSecel.selected == true){
						glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
						glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
					}

					if(currentSecel.atomHashes.size() > 1){
						PDBAtom firstAtom = atoms.find(currentSecel.atomHashes[0])->second;
						PDBAtom lastAtom = atoms.find(currentSecel.atomHashes[currentSecel.atomHashes.size()-1])->second;
						Vector3DFloat preSecelAtomPos = atoms.find(firstAtom.GetPrevCAHash())->second.GetPosition();
						Vector3DFloat postSecelAtomPos = atoms.find(lastAtom.GetNextCAHash())->second.GetPosition();

						vector<Vector3DFloat> points = CreatePointVector(firstAtom, lastAtom);
						vector<Vector3DFloat> normals(points.size());
						normals = CreateStrandNormals(points, preSecelAtomPos, postSecelAtomPos);

						// generate smoothed interpolated points
						vector<Vector3DFloat> interpolatedPoints = InterpolateLoopPoints(points, preSecelAtomPos, postSecelAtomPos, NUM_SEGMENTS);
						int ptsize = interpolatedPoints.size();

						// create vectors to hold the vertices and normals of our loop polygon
						int renderingsize = (ptsize)*NUM_SLICES;
						switch (renderingType){
							case 0:
								renderingsize = ptsize*2;
								break;
							case 1:
								renderingsize = (ptsize)*NUM_SLICES;
								break;
							default:
								break;
						}
						vector<Vector3DFloat> renderingPoints(renderingsize);
						vector<Vector3DFloat> renderingNormals(renderingPoints.size());

						Vector3DFloat nextVector = interpolatedPoints[1]-interpolatedPoints[0];
						Vector3DFloat previousVector = preSecelAtomPos - interpolatedPoints[0];
						if(previousVector.Length() < .001){
							if(2 < interpolatedPoints.size()){
								previousVector = interpolatedPoints[2] - interpolatedPoints[0];
							} else {
								previousVector = postSecelAtomPos - interpolatedPoints[0];
							}
						}
						Vector3DFloat curAxis = nextVector^previousVector;
						Vector3DFloat curNormal = nextVector^curAxis;
						curNormal.Normalize();
						Vector3DFloat curPos = interpolatedPoints[0];
						nextVector.Normalize();

						// generate first stack of points
						Vector3DFloat outward = normals[0] ^ (interpolatedPoints[1]-interpolatedPoints[0]);
						outward.Normalize();
						switch (renderingType){
							case 0:
								for(unsigned int k = 0; k < 2; ++k){
									renderingPoints[k*ptsize] = curPos + outward*LOOP_RADIUS*(k%2 == 0 ? 1 : -1);
									renderingNormals[k*ptsize] = normals[0];
								}
								break;
							case 1:
								for(unsigned int k = 0; k < NUM_SLICES; ++k){
									Vector3DFloat outwardNormal = (curNormal.Rotate(nextVector, ((double)k)*2*PI/((double)NUM_SLICES)));
									outwardNormal.Normalize();
									renderingPoints[k*(ptsize)] = curPos+outwardNormal*LOOP_RADIUS;
									renderingNormals[k*(ptsize)] = outwardNormal; //renderingPoints[j*(ptsize)]-curPos;
								}
								break;
							default:
								break;
						}

						for(unsigned int j = 1; j < ptsize; ++j){
							Vector3DFloat nextVector, previousVector;
							curPos = interpolatedPoints[j];
							if(j == ptsize - 1){
								nextVector = postSecelAtomPos - interpolatedPoints[j];
							} else {
								nextVector = interpolatedPoints[j+1]-interpolatedPoints[j];
							}

							if(j == 0){
								previousVector = interpolatedPoints[0] - preSecelAtomPos;
							} else {
								previousVector = interpolatedPoints[j-1] - interpolatedPoints[j];
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

								newAxis = (interpolatedPoints[nix] - interpolatedPoints[j])^(interpolatedPoints[pix] - interpolatedPoints[j]);
							}
							curAxis = newAxis;
							double alpha = asin(curAxis.Length()/(nextVector.Length()*previousVector.Length()));
							curAxis.Normalize();
							curNormal.Normalize();
							if (nextVector.Length() > .001 && previousVector.Length() > .001){
								curNormal = curNormal.Rotate(curAxis, -1*alpha);
							}

							nextVector.Normalize();

							Vector3DFloat dirtemp;
							if (j + 1 < ptsize){
								dirtemp = (interpolatedPoints[j+1] - interpolatedPoints[j]);
							} else {
								dirtemp = interpolatedPoints[j] - interpolatedPoints[j-1];
							}
							float tsect = ((float)(j%NUM_SEGMENTS))/((float)NUM_SEGMENTS);
							Vector3DFloat outward = (normals[j/NUM_SEGMENTS]*(1.0 - tsect) + normals[j/NUM_SEGMENTS + 1]*(tsect)) ^ dirtemp;
							outward.Normalize();
							switch(renderingType){
								case 0:
									for(unsigned int k = 0; k < 2; ++k){
										renderingPoints[j+k*(ptsize)] = curPos+outward*LOOP_RADIUS*(k%2 == 0 ? 1 : -1);
										renderingNormals[j+k*(ptsize)] = normals[j/NUM_SEGMENTS];
									}
									break;
								case 1:
									for(unsigned int k = 0; k < NUM_SLICES; ++k){
										Vector3DFloat outwardNormal = (curNormal.Rotate(nextVector, ((double)k*2*PI)/NUM_SLICES));
										outwardNormal.Normalize();
										renderingPoints[j+k*(ptsize)] = curPos+outwardNormal*LOOP_RADIUS;
										renderingNormals[j+k*(ptsize)] = outwardNormal;
									}
									break;
								default:
									cout << "shouldn't've gotten to the loop rendering default case" << endl;
									break;
							}


						}
						float r,g,b,a;
						if(currentSecel.selected == true){
							OpenGLUtils::GetColor(r,g,b,a);
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
						if(currentSecel.selected == true){
							OpenGLUtils::SetColor(r,g,b,a);
						}
					}

					glPopAttrib();

					if(currentSecel.selected == true){
						glPushAttrib(GL_LIGHTING_BIT);

						Vector3DFloat pos1 = atoms[currentSecel.atomHashes[0]].GetPosition();
						Vector3DFloat pos2 = atoms[currentSecel.atomHashes[currentSecel.atomHashes.size()-1]].GetPosition();
						//printf("Drawing PDB Spheres at PDB ID %d with end #1 [%f, %f, %f] and #2 [%f, %f, %f]\n", i+1, pos1.X(), pos1.Y(), pos1.Z(), pos2.X(), pos2.Y(), pos2.Z());

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

			if(selectEnabled) {
				glPopName();
				glPopName();
			}
		}

		void CAlphaRenderer::DrawSideChainModel(int subSceneIndex, bool selectEnabled) {
			GLfloat emissionColor[4] = {1.0, 1.0, 1.0, 1.0};
			float r,g,b,a;

			if(subSceneIndex == 0) { // Drawing Atoms				
				if(selectEnabled) {
					atomHashKeys.clear();
					glPushName(0);
					glPushName(0);
				}
				for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
					glPushAttrib(GL_LIGHTING_BIT);					
					if(i->second.GetSelected()) {
						glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
						glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
					} else {
						i->second.GetColor(r, g, b, a);
						OpenGLUtils::SetColor(r,g,b,a);
					}					

					if(selectEnabled){
						//TODO: possibly implement mouse picking using ray intersection
						atomHashKeys.push_back(i->first); // adding the hash key
						glLoadName(static_cast<GLuint>(atomHashKeys.size() - 1)); // using the index of the element just added
					}
					if(i->second.GetVisible()) {
						DrawSphere(i->second.GetPosition(), i->second.GetAtomRadius() * 0.3);
					}

					glPopAttrib();

				}
				if(selectEnabled) {
					glPopName();
					glPopName();
				}
			} else if(subSceneIndex == 1) { // Drawing Bonds
				if(selectEnabled) {
					glPushName(1);
					glPushName(0);
				}
				Vector3DFloat v1, vc, v2;


				for(int i=0; i < (int)sidechainBonds.size(); i++) {
					glPushAttrib(GL_LIGHTING_BIT);
					if(selectEnabled){
						glLoadName(i);
					}

					if(atoms[sidechainBonds[i].GetAtom0Ix()].GetVisible() && atoms[sidechainBonds[i].GetAtom1Ix()].GetVisible()) {
						v1 = atoms[sidechainBonds[i].GetAtom0Ix()].GetPosition();
						v2 = atoms[sidechainBonds[i].GetAtom1Ix()].GetPosition();
						if(sidechainBonds[i].GetSelected()) {
							glMaterialfv(GL_FRONT, GL_EMISSION, emissionColor);
							glMaterialfv(GL_BACK, GL_EMISSION, emissionColor);
							DrawCylinder(v1, v2, 0.1, 6, 2);
						} else {						
							vc = (v1 + v2) * 0.5;
							atoms[sidechainBonds[i].GetAtom0Ix()].GetColor(r, g, b, a);
							OpenGLUtils::SetColor(r,g,b,a);
							DrawCylinder(v1, vc, 0.1, 6, 2);
							atoms[sidechainBonds[i].GetAtom1Ix()].GetColor(r, g, b, a);
							OpenGLUtils::SetColor(r,g,b,a);
							DrawCylinder(vc, v2, 0.1, 6, 2);
						}
					}
					glPopAttrib();
				}
				if(selectEnabled) {
					glPopName();
					glPopName();
				}
			} else if(subSceneIndex == 2) { // Drawing spheres to cover up the cylinder edges					
				for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
					glPushAttrib(GL_LIGHTING_BIT);					
					i->second.GetColor(r, g, b, a);
					OpenGLUtils::SetColor(r,g,b,a);
					DrawSphere(i->second.GetPosition(), 0.1);										
					glPopAttrib();
				}
			}
		}

		void CAlphaRenderer::Draw(int subSceneIndex, bool selectEnabled) {
			switch(displayStyle) {
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
			if(subsceneIndex == 0) {
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
			vector<unsigned long long> eraseKeys;
			eraseKeys.clear();

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {				
				if(i->second.GetName().compare("CA") != 0) {					
					eraseKeys.push_back(i->first);
				} 
			}

			for(unsigned int i = 0; i < eraseKeys.size(); i++) {
				atoms.erase(atoms.find(eraseKeys[i]));
			}

			eraseKeys.clear();

			std::list<SerialAndHashType> sortedSerials;
			SerialAndHashType elem;
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				elem.hashKey = i->first;
				elem.serial = i->second.GetSerial();

				sortedSerials.push_back(elem);				
			}
			sortedSerials.sort(SerialAndHashTypePredicate());


			std::list<SerialAndHashType>::iterator oldAtom = sortedSerials.begin();
			std::list<SerialAndHashType>::iterator startAtom = sortedSerials.begin();

			startAtom++;
			for(std::list<SerialAndHashType>::iterator i = startAtom; i != sortedSerials.end(); i++) {
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

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				tempFactor = i->second.GetTempFactor();
				if(tempFactor > maxTempFactor) {
					maxTempFactor = tempFactor;
				}
				if(tempFactor < minTempFactor) {
					minTempFactor = tempFactor;
				}
			}
			float r, g, b;

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				i->second.SetAtomRadius(3.0);
				tempFactor = i->second.GetTempFactor();
				if(tempFactor < 0) {
					tempFactor = (tempFactor / minTempFactor);
					r = 1.0f - tempFactor;
					g = 1.0f - tempFactor;
					b = 1.0f;
				} else {
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
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				i->second.SetTempFactor( i->second.GetTotalScore(correlationCoeff, skeletonCoeff, geometryCoeff) );
			}
			ColorSSEHunterAtoms();
		}

		void CAlphaRenderer::ColorSSEHunterAtoms() {
			float maxTempFactor = -10000.0f, minTempFactor = 10000.0f;
			float tempFactor;

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				tempFactor = i->second.GetTempFactor();
				if(tempFactor > maxTempFactor) {
					maxTempFactor = tempFactor;
				}
				if(tempFactor < minTempFactor) {
					minTempFactor = tempFactor;
				}
			}
			float r, g, b;

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				i->second.SetAtomRadius(3.0);
				tempFactor = i->second.GetTempFactor();
				if(tempFactor < 0) {
					tempFactor = (tempFactor / minTempFactor);
					r = 1.0f - tempFactor;
					g = 1.0f - tempFactor;
					b = 1.0f;
				} else {
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
			for(unsigned int i = 0; i < bonds.size(); i++) {
				if(bonds[i].GetSelected()) {
					count++;
				}
			}
			return count;
		}

		int CAlphaRenderer::SelectionAtomCount(){
			int count = 0;
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {					
				if(i->second.GetSelected()) {
					count++;
				}
			}
			return count;
		}


		Vector3DFloat CAlphaRenderer::SelectionCenterOfMass() {
			int count = 0;
			Vector3DFloat centerOfMass = Vector3DFloat(0,0,0);
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {					
				if(i->second.GetSelected()) {
					count++;
					centerOfMass = centerOfMass + i->second.GetPosition();
				}
			}

			for(unsigned int i = 0; i < bonds.size(); i++) {
				if(bonds[i].GetSelected()) {
					count++;
					centerOfMass = centerOfMass + (atoms[bonds[i].GetAtom0Ix()].GetPosition() + atoms[bonds[i].GetAtom1Ix()].GetPosition()) * 0.5;
				}
			}
			if(count == 0) {
				centerOfMass = Renderer::SelectionCenterOfMass();
			} else {
				centerOfMass = centerOfMass * (1.0f/(float)count);
			}
			return centerOfMass;
		}

		bool CAlphaRenderer::SelectionRotate(Vector3DFloat centerOfMass, Vector3DFloat rotationAxis, float angle) {
			bool rotated = false;
			Vector3 centerOfMassP3 = Vector3(centerOfMass.X(), centerOfMass.Y(), centerOfMass.Z());
			Vector3 rotationV3 = Vector3(rotationAxis.X(), rotationAxis.Y(), rotationAxis.Z());

			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {
				if(i->second.GetSelected()) {
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
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {					
				if(i->second.GetSelected()) {
					i->second.SetPosition(i->second.GetPosition() + moveDirection);
					i->second.SetFlag(1);
					moved = true;
				} else {
					i->second.SetFlag(0);
				}
			}

			for(unsigned int i = 0; i < bonds.size(); i++) {
				if(bonds[i].GetSelected()) {
					PDBAtom a = atoms[bonds[i].GetAtom0Ix()];
					if(a.GetFlag() == 0) {
						a.SetPosition(a.GetPosition() + moveDirection);
						a.SetFlag(1);
						moved = true;
					}

					a = atoms[bonds[i].GetAtom1Ix()];
					if(a.GetFlag() == 0) {
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
			if(Renderer::SelectionClear()) {					
				for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {					
					i->second.SetSelected(false);
				}

				for(unsigned int i = 0; i < bonds.size(); i++) {
					bonds[i].SetSelected(false);
				}
				for(unsigned int i = 0; i < aHelices.size(); i++) {
					aHelices[i].selected = false;
				}
				for(unsigned int i = 0; i < bStrands.size(); ++i) {
					bStrands[i].selected = false;
				}
				for(unsigned int i = 0; i < loops.size(); ++i) {
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
				switch(displayStyle) {
					case CALPHA_DISPLAY_STYLE_BACKBONE:
					case CALPHA_DISPLAY_STYLE_SIDE_CHAIN:
						it = atoms.find(atomHashKeys.at(ix0));
						if (it != atoms.end()) {
							a = &(it->second);
							a->SetSelected( forceTrue || !a->GetSelected() );
						}
						break;
					case CALPHA_DISPLAY_STYLE_RIBBON:

						if(aHelices[ix0].selected == true && !forceTrue){
							aHelices[ix0].selected = false;
						}
						else{
							cout << "Updating selectedHelix" << " ix0=" << ix0 << " forceTrue=" << forceTrue << endl;
							aHelices[ix0].selected = true;
							selectedHelixIndices.push_back(ix0);
						}
						break;
				} 
			} else if((subsceneIndex == 1) && (ix0 >= 0) && (ix0 <= (int)bonds.size())) {
				switch(displayStyle) {
					case CALPHA_DISPLAY_STYLE_BACKBONE:
						bonds[ix0].SetSelected(forceTrue || !bonds[ix0].GetSelected());
						break;
					case CALPHA_DISPLAY_STYLE_RIBBON:
						if(bStrands[ix0].selected == true && !forceTrue){
							bStrands[ix0].selected = false;
						} else {
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
			} else if((subsceneIndex == 2) && (ix0 != NULL)) {
				switch(displayStyle) {
					case CALPHA_DISPLAY_STYLE_BACKBONE:
						break;
					case CALPHA_DISPLAY_STYLE_RIBBON:
						if(loops[ix0].selected == true && !forceTrue){
							loops[ix0].selected = false;
						} else {
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
			if(atoms.size() > 0) {
				for(int i = 0; i < 3; i++) {
					minPts[i] = atoms.begin()->second.GetPosition().values[i];
					maxPts[i] = atoms.begin()->second.GetPosition().values[i];
				}

				for(AtomMapType::iterator j = atoms.begin(); j != atoms.end(); j++) {
					for(int i = 0; i < 3; i++) {
						minPts[i] = min(minPts[i], j->second.GetPosition().values[i]);
						maxPts[i] = max(maxPts[i], j->second.GetPosition().values[i]);
					}
				}
			} else {
				for(int i = 0; i < 3; i++) {
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
			for(AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
				if(it->second.GetSelected()) {
					if(count == selectionId) {
						return &it->second;
					}
					count++;
				}
			}
			return NULL;
		}

		vector<unsigned long long> CAlphaRenderer::GetAtomHashes() {
			vector<unsigned long long> atomHashes;
			for (AtomMapType::iterator it = atoms.begin(); it != atoms.end(); it++) {
				atomHashes.push_back(it->first);
			}
			return atomHashes;
		}

		int CAlphaRenderer::GetBondIndex(unsigned long long atom0, unsigned long long atom1) {
			for(unsigned int i = 0; i < bonds.size(); i++) {
				if(((bonds[i].GetAtom0Ix() == atom0) && (bonds[i].GetAtom1Ix() == atom1)) ||
					((bonds[i].GetAtom0Ix() == atom1) && (bonds[i].GetAtom1Ix() == atom0))) {
						return i;
				}
			}
			return -1;
		}

		int CAlphaRenderer::GetSideChainBondIndex(unsigned long long atom0, unsigned long long atom1) {
			for(unsigned int i = 0; i < sidechainBonds.size(); i++) {
				if(((sidechainBonds[i].GetAtom0Ix() == atom0) && (sidechainBonds[i].GetAtom1Ix() == atom1)) ||
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
			switch(subsceneIndex) {
				case(0):
					if((ix0 >= 0) && (ix0 <= (int)atoms.size())) {
						PDBAtom * a = (PDBAtom*)ix0;
						position = a->GetPosition();
					}
					break;
				case(1):
					if((ix0 >= 0) && (ix0 <= (int)bonds.size())) {
						position = (atoms[bonds[ix0].GetAtom0Ix()].GetPosition() + atoms[bonds[ix0].GetAtom1Ix()].GetPosition()) * 0.5;
					}
					break;
				default:
					position = Vector3DFloat(0,0,0);
					break;
			}
			return position;
		}

		void CAlphaRenderer::TransformAllAtomLocations(MatrixFloat transform) {
			for(AtomMapType::iterator i = atoms.begin(); i != atoms.end(); i++) {					
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

		void CAlphaRenderer::SetHelixCorrs(  vector < int > flatCorrespondences){
			if(flatCorrespondences.size() %2 != 0)
				return;
			else
				corrs.clear();
			for(int i=0; i < flatCorrespondences.size(); i = i+2){
				corrs.push_back(boost::tuple<int, int>(flatCorrespondences[i], flatCorrespondences[i + 1]));
			}
		}

		void CAlphaRenderer::SetFeatureVecs(vector<Vector3DFloat> flatFeatureVecs){
			if(flatFeatureVecs.size() %2 != 0)
				return;
			else
				featureVecs.clear();
			for(int i=0; i < flatFeatureVecs.size(); i = i+2){
				featureVecs.push_back(boost::tuple<Vector3DFloat, Vector3DFloat>(flatFeatureVecs[i], flatFeatureVecs[i + 1]));
			}

		}
		void CAlphaRenderer::SetSelectedSSEHelices(vector<int> indices){
			selectedSSEHelices.clear();
			selectedSSEHelices = indices;
		}

		void CAlphaRenderer::SetHelixColor(int helixNum, float r, float g, float b){
			//cout << "setting helix color " << helixNum << " to (" << r << ", " << g << ", " << b << ")" <<endl;
			helixColors.erase(helixNum);
			helixColors.insert(pair<int, boost::tuple<float, float, float> >(helixNum, boost::tuple<float, float, float>(r, g, b)));
		}

        void CAlphaRenderer::ClearHelixColors() {
            helixColors.clear();
        }

		// creates a vector of Vector3DFloats that represents the locations of all the PDBAtoms
		// starting with start and ending with end; it does not error check, so incorrectly
		// ordered points will break this method.  there are more efficient ways to handle this
		// functionality, but this seems simple and flexible enough
		vector<Vector3DFloat> CAlphaRenderer::CreatePointVector(PDBAtom start, PDBAtom end){
			vector<Vector3DFloat> points;

			PDBAtom current = start;
			while(current.GetHashKey() != end.GetHashKey()){
				points.push_back(current.GetPosition());
				if(current.GetHashKey() == current.GetNextCAHash()){
					break;
				}
				current = atoms.find(current.GetNextCAHash())->second;
			}

			points.push_back(end.GetPosition());
			return points;
		}

		// implementation of Laplacian smoothing for a vector of Vector3DFloats (treats them like points)
		// creating copies of "points" twice seems unnecessary, but I am unsure about the performance cost,
		// so I am leaving it for simplicity of implementation
		vector<Vector3DFloat> CAlphaRenderer::LaplacianSmoothing(vector<Vector3DFloat> points, int steps){
			vector<Vector3DFloat> pointsTemp(points);
			vector<Vector3DFloat> smoothedPoints(points);

			for(int i = 0; i < steps; ++i){
				for(int j = 1; j < points.size()-1; ++j){
					smoothedPoints[j] = (pointsTemp[j-1] + pointsTemp[j+1])*.5;
					smoothedPoints[j] = (smoothedPoints[j] + pointsTemp[j])*.5;
				}
				pointsTemp = smoothedPoints;
			}
			return pointsTemp;
		}

		// unsure of what behavior should be if points.size() < 3; in molscript the strand is skipped in this case
		vector<Vector3DFloat> CAlphaRenderer::CreateStrandNormals(vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next){
			vector<Vector3DFloat> normals(points);
			int ptsSize = points.size();

			for(int i = 1, length = ptsSize - 1; i < length; ++i){
				Vector3DFloat newPos = (points[i-1] + points[i+1])*.5;
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
			for(int j = 0, size = ptsSize - 1; j < size; ++j){
				if(normals[j]*normals[j+1] < 0){
					normals[j+1] = normals[j+1]*-1;
				}
			}

			// "smooth normals, one iteration" - molscript/graphics.c
			vector<Vector3DFloat> smoothedNormals(normals);

			for(int k = 1, size = ptsSize - 1; k < size; ++k){
				smoothedNormals[k] = normals[k-1] + normals[k] + normals[k+1];
				smoothedNormals[k].Normalize();
			}

			// "normals exactly perpendicular to strand" - molscript/graphics.c
			Vector3DFloat direction = points[1] - points[0];
			Vector3DFloat side = direction^smoothedNormals[0];
			smoothedNormals[0] = side ^ direction;
			smoothedNormals[0].Normalize();

			for(int i = 1, size = ptsSize - 1; i < size; ++i){
				direction = points[i+1] - points[i-1];
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

		void CAlphaRenderer::CreateHelixAxesTangentsAndPoints(vector<Vector3DFloat>& axes, vector<Vector3DFloat>& tangents, vector<Vector3DFloat>& interpPoints, std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, double HELIX_ALPHA, double HELIX_BETA, double HELIX_HERMITE_FACTOR){
			if(points.size() > 2){

				for(int i = 0; i < points.size() - 1; ++i){
					
					if(i > 0){
						Vector3DFloat cvec = points[i+1] - points[i-1];
						cvec.Normalize();
	
						Vector3DFloat rvec = (points[i]-points[i-1])^(points[i+1]-points[i]);
						rvec.Normalize();
	
						axes[i] = rvec*sin(HELIX_ALPHA) + cvec*cos(HELIX_ALPHA);
						tangents[i] = rvec*sin(HELIX_BETA) + cvec*cos(HELIX_BETA);
						tangents[i] = tangents[i]*HELIX_HERMITE_FACTOR;
					}
				}
				axes[0] = axes[1];
				axes[axes.size()-1] = axes[axes.size()-2];

				tangents[0] = previous - points[1];
				tangents[0].Normalize();
				tangents[0] = tangents[0]*HELIX_HERMITE_FACTOR;
				tangents[tangents.size()-1] = next - points[points.size()-2];
				tangents[tangents.size()-1].Normalize();
				tangents[tangents.size()-1] = tangents[tangents.size()-1]*HELIX_HERMITE_FACTOR;
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
					glColorMaterial ( GL_FRONT_AND_BACK, GL_SPECULAR | GL_AMBIENT_AND_DIFFUSE );
					glEnable ( GL_COLOR_MATERIAL );
					break;
				default:
					break;
			}
			for(int i = 0, runlength = points.size(); i < runlength; ++i){
				if (i%(stacks+1) == 0){
					glBegin(GL_TRIANGLE_STRIP);
				}

				int nextSliceIx = (i+stacks+1)%runlength;

				glNormal3f(normals[i].X(), normals[i].Y(), normals[i].Z());
				glVertex3f(points[i].X(), points[i].Y(), points[i].Z());

				glNormal3f(normals[nextSliceIx].X(), normals[nextSliceIx].Y(), normals[nextSliceIx].Z());
				glVertex3f(points[nextSliceIx].X(), points[nextSliceIx].Y(), points[nextSliceIx].Z());

				if((i+1)%(stacks+1) == 0){
					glEnd();
				}
			}
		}

		vector<Vector3DFloat> CAlphaRenderer::InterpolateLoopPoints(std::vector<Vector3DFloat> points, Vector3DFloat previous, Vector3DFloat next, int NUM_SECTIONS){
			HermiteCurve curve;
			Vector3DFloat m0, m1;
			vector<Vector3DFloat> pointstemp(points);
			bool LAPLACIAN_SMOOTHING = true;
			int SMOOTHING_STEPS = 1;
			double HERMITE_FACTOR = 0.5;
			int LOOP_SLICES = 10;
			if(LAPLACIAN_SMOOTHING){
				pointstemp = LaplacianSmoothing(points, SMOOTHING_STEPS);
			}

			vector<Vector3DFloat> interpolatedPoints((pointstemp.size()-1)*(NUM_SEGMENTS));

			for(unsigned int i = 0; i < points.size()-1; ++i){
				if(i == 0){
					m0 = pointstemp[i+1] - previous;
				} else {
					m0 = pointstemp[i+1] - pointstemp[i-1];
					m0 = m0*HERMITE_FACTOR;
				}

				if(i + 2 > pointstemp.size() - 1){
					m1 = next - pointstemp[i];
				} else {
					m1 = pointstemp[i+2] - pointstemp[i];
					m1 = m1*HERMITE_FACTOR;
				}

				curve.setCurve(pointstemp[i], pointstemp[i+1], m0, m1);
				interpolatedPoints[i*(NUM_SEGMENTS)] = pointstemp[i];
				for (int sect = 1; sect < NUM_SEGMENTS; ++sect){
					double tsect = ((double)sect)/((double)NUM_SEGMENTS);
					interpolatedPoints[i*(NUM_SEGMENTS) + sect] = curve.getPos(tsect);
				}
			}
			interpolatedPoints[interpolatedPoints.size()-1] = points[points.size() -1];
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
			hlt_r = ((double)col)/1000.0;
			//thinRibbThickness = hlt_r;
			cout << "hlt_r: " << hlt_r << endl;
		}

		void CAlphaRenderer::SetHltGValue(int col){
			hlt_g = ((double)col)/100.0;
			cout << "hlt_g: " << hlt_g << endl;
		}

		void CAlphaRenderer::SetHltBValue(int col){
			hlt_b = ((double)col)/100.0;
			cout << "hlt_b: " << hlt_b << endl;
		}

		void CAlphaRenderer::SetHltAValue(int col){
			hlt_a = ((double)col)/100.0;
			cout << "hlt_a: " << hlt_a << endl;
		}
	}
}

#endif
