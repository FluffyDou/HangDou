// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu), Refactored by Sasakthi Abeysinghe (sasakthi.abeysinghe@wustl.edu)
// Description:   Volumetric data definition


#ifndef SKELETON_MAKER_VOLUME_H
#define SKELETON_MAKER_VOLUME_H
#ifndef BOOST_MATH_TOOLS_SERIES_INCLUDED
#define BOOST_MATH_TOOLS_SERIES_INCLUDED

#define MAX_SHEETS 100000
#define MAX_QUEUELEN 5000000
#define MAX_ERODE 1000

#include "VolumeData.h"
#include "ThinningTemplate.h"
#include "GridQueue.h"
#include "GridQueue2.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "PriorityQueue.h"
#include <vector>
#include <MathTools/Vector3D.h>
#include <sstream>
#include <string>
#include <fstream>
#include <MathTools/LinearSolver.h>
#include <MathTools/DataStructures.h>
#include <MathTools/MathLib.h>
#include <MathTools/MatlabWrapper.h>
#include <boost/python.hpp>
#include <tbb/task_scheduler_init.h>
#include <map>
//#include <MathTools/VectorLib.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range3d.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Graph;

using namespace std;
using namespace wustl_mm::MathTools;

using namespace tbb;
using namespace Eigen;
using namespace std;
using namespace wustl_mm::MathTools;
namespace py = boost::python;
//using namespace tbb;

namespace wustl_mm {
    namespace SkeletonMaker {

        const int neighbor6[6][3] = { { 0, 0, 1 }, { 0, 0, -1 }, { 0, 1, 0 }, { 0, -1, 0 }, { 1, 0, 0 }, { -1, 0, 0 } };
        const int neighbor4[4][2] = { { 0, 1 }, { 0, -1 }, { 1, 0 }, { -1, 0 } };
        const int neighbor64[6][4][3] = {
            { { 0, 1, 0 }, { 0, -1, 0 }, { 1, 0, 0 }, { -1, 0, 0 } },
            { { 0, 1, 0 }, { 0, -1, 0 }, { 1, 0, 0 }, { -1, 0, 0 } },
            { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 } },
            { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 } },
            { { 0, 0, 1 }, { 0, 0, -1 }, { 0, 1, 0 }, { 0, -1, 0 } },
            { { 0, 0, 1 }, { 0, 0, -1 }, { 0, 1, 0 }, { 0, -1, 0 } } };

        const int sheetNeighbor[12][4][3] = {
            { { 0, -1, -1 }, { 0, -1, 0 }, { 0, 0, -1 }, { 0, 0, 0 } },
            { { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, 0 }, { 0, 0, 1 } },
            { { 0, 0, -1 }, { 0, 0, 0 }, { 0, 1, -1 }, { 0, 1, 0 } },
            { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, { 0, 1, 1 } },

            { { -1, 0, -1 }, { -1, 0, 0 }, { 0, 0, -1 }, { 0, 0, 0 } },
            { { -1, 0, 0 }, { -1, 0, 1 }, { 0, 0, 0 }, { 0, 0, 1 } },
            { { 0, 0, -1 }, { 0, 0, 0 }, { 1, 0, -1 }, { 1, 0, 0 } },
            { { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 1 } },

            { { -1, -1, 0 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, 0 } },
            { { -1, 0, 0 }, { -1, 1, 0 }, { 0, 0, 0 }, { 0, 1, 0 } },
            { { 0, -1, 0 }, { 0, 0, 0 }, { 1, -1, 0 }, { 1, 0, 0 } },
            { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 1, 1, 0 } }
        };

        const int faceCells[12][2] = { { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 } };

        const int cubeFaces[6][4] =
        { { 1, 5, 7, 3 }, { 0, 2, 6, 4 }, { 2, 3, 7, 6 }, { 0, 4, 5, 1 }, { 5, 4, 6, 7 }, { 0, 1, 3, 2 } };

        const int faceEdges[12][2] = { { 3, 1 }, { 3, 0 }, { 2, 1 }, { 2, 0 },
        { 5, 1 }, { 5, 0 }, { 4, 1 }, { 4, 0 },
        { 5, 3 }, { 5, 2 }, { 4, 3 }, { 4, 2 } };

        const int edgeFaces[6][4] = { { 1, 3, 5, 7 }, { 0, 2, 4, 6 }, { 2, 3, 9, 11 }, { 0, 1, 8, 10 }, { 6, 7, 10, 11 }, { 4, 5, 8, 9 } };

        struct gridPoint
        {
            int x, y, z;
        };
        /**
        struct Edge
        {
        float edgePoint[3];
        float a;
        float b;
        float w;
        int edgeTag = 0;	// 0 for non-extremal, 1 for max, 2 for min
        int f;
        bool extremal = false;
        bool valid = false;
        };
        **/
        class Edge
        {
        public:
            Edge(void);
            ~Edge(void);

            float edgePoint[3];
            float a;
            float b;
            float w;
            int edgeTag;	// 0 for non-extremal, 1 for max, 2 for min
            int f;
            bool extremal;
            bool valid;
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;

        };

        Edge::Edge()
        {
            extremal = false;
            valid = false;
            edgeTag = 0;
        }

        Edge::~Edge()
        {
        }

        class Cell
        {
        public:
            Cell(void);
            ~Cell(void);

            float centroid[3], extremalPoint[3];
            //float area;
            /*float relativeSaliencies[3];
            float localIntensity;
            float firstEigenvalue;*/
            int projType, extrPType, vertIdx;
            bool extremal, extrP, valid;
            float eigenValue1, eigenValue2, eigenValue3;
            std::vector<float> eigenVectors;
        };

        Cell::Cell(void)
        {
            extremal = false;
            extrP = false;
            valid = false;
            projType = 0;
            vertIdx = -1;
        }

        Cell::~Cell(void)
        {
        }

        class Face
        {
        public:
            Face(void);
            ~Face(void);

            /*float subdVerts[subdNum][subdNum][3];
            float gradients[subdNum][subdNum][3];
            float areas[(subdNum-1)*(subdNum-1)*2];
            float centers[(subdNum-1)*(subdNum-1)*2][3];*/
            float spherArea;
            float facePoint[3];
            //float h[4];
            /*float relativeSaliencies[3];
            float localIntensity;
            float firstEigenvalue;*/
            int faceTag;	// 0 for non-extremal, 1 for max, 2 for min, 3 for saddle
            bool extremal, valid;
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;
        };

        Face::Face(void)
        {
            extremal = false;
            valid = false;
            faceTag = 0;
        }

        Face::~Face(void)
        {
        }

        class Quad
        {
        public:
            Quad(void);
            ~Quad(void);

            //float vertices[4][3];
            int vertIdxs[4];
            float relativeSaliencies[3];
            float localIntensity;
            float firstEigenvalue;
            int type;
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;
        };

        Quad::Quad(void)
        {
        }

        Quad::~Quad(void)
        {
        }

        class Point
        {
        public:
            Point(void);
            ~Point(void);

            //float position[3];
            int vertIdx;
            float relativeSaliencies[3];
            float localIntensity;
            float firstEigenvalue;
            int type;
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;
        };

        Point::Point(void)
        {
        }

        Point::~Point(void)
        {
        }

        class Segment
        {
        public:
            Segment(void);
            ~Segment(void);

            //float vertices[2][3];
            int vertIdxs[2];
            float relativeSaliencies[3];
            float localIntensity;
            float firstEigenvalue;
            int type;
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;
        };

        bool operator< (const Segment &c1, const Segment &c2)
        {
            return c1.vertIdxs[0] < c2.vertIdxs[0];
        }

        Segment::Segment(void)
        {
        }

        Segment::~Segment(void)
        {
        }

        class Vertex
        {
        public:
            Vertex(void);
            ~Vertex(void);

            float position[3];
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;

        };

        Vertex::Vertex(void)
        {
        }

        Vertex::~Vertex(void)
        {
        }

        class DispCell
        {
        public:
            DispCell(void);
            ~DispCell(void);

            float center[3];
            float corners[8][3];
            float gradients[8][3];
            float v1s[8][3];
            float v3s[8][3];
            Cell cell;
            Face faces[6];
            float eigenvalues[6][3];
            float sumAB[6];
            Edge edges[12];
            std::vector<float> sampleP[12], sampleG[12], sampleV1[12], sampleV3[12], sampleX[12], sampleY[12], sampleProjG[12];
            std::vector<float> eigenVectors;
            float eigenValue1, eigenValue2, eigenValue3;

        };

        DispCell::DispCell(void)
        {
        }

        DispCell::~DispCell(void)
        {
        }



        class CGNSParallelFace
        {
            int gridx, gridy, gridz;
            Edge *edges;
            Face *faces;
            int *totalFacePoints;
            int dataType;
            float change;
            static const int subdNum = 3;
            float *gridPoints;
            float *coef;

            void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, int index1, int index2, int index3, int index4,
                int axis, int i, int j, int k, Face *lfaces, int *ltotalFacePoints) const;
            int getFaceIndex(int axis, int x, int y, int z) const;
            int getEdgeIndex(int axis, int x, int y, int z) const;
            void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int axis) const;
            void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h, int *ltotalFacePoints) const;
            void getCurveType(Face *face) const;
            void get3DWindingNum(int index1, int index3, Face *face, int axis) const;
            float smallAbsA(float A) const;
            void getGridPointPosD(int index, float v[3]) const;
            void rotate2D(float p[2], float theta) const;
            float signedArea3D(float v1[3], float v2[3], float v3[3]) const;
            void getGridPointPos(int index, float v[3]) const;
            float getDihedral(float b1[3], float b2[3], float b3[3]) const;



        public:
            bool *storeGrid_2DGradMag;
            float *grid_2DGradMag;

            CGNSParallelFace(int lgridx, int lgridy, int lgridz, int *ltotalFacePoints, Edge *ledges, Face *lfaces, int ldataType,
                float *lgridPoints, float lchange, bool *lstoreGrid_2DGradMag, float *lgrid_2DGradMag, float *lcoef)
            {
                gridx = lgridx;
                gridy = lgridy;
                gridz = lgridz;
                totalFacePoints = ltotalFacePoints;
                edges = ledges;
                faces = lfaces;
                dataType = ldataType;
                gridPoints = lgridPoints;
                change = lchange;
                storeGrid_2DGradMag = lstoreGrid_2DGradMag;
                grid_2DGradMag = lgrid_2DGradMag;
                coef = lcoef;
            };
        };

        class General_Data
        {
        public:
            General_Data();
            ~General_Data(void);
            virtual void dataGeneration(char *filename) {};
            virtual void buildGrid() {};
            void edgePhase(float* scalars, float* tensors, float* gradients, float *edgeTable);
            float *getGridShowPositions();
            float *getGridScalars();
            float *getGridV1s();
            float *getGridV2s();
            float *getGridV3s();
            float *getGridGradients();
            float *getEdgePoints();
            float *getFacePoints();
            float *getCellPoints();
            std::vector<Segment> *getSegments();
            std::vector<Quad> *getQuads();
            std::vector<Point> *getPoints();
            std::vector<Vertex> *getVertices();
            int getMaxGridIndex();
            float getMaxGradientNorm();
            float getMaxGrid();
            void facePhase(float* scalars, float* tensors, float* gradients, float *faceTable);
            void cellPhase(float* scalars, float* tensors, float* gradients);
            void edgePhase(float *coef);
            void facePhase(float *coef);
            void cellPhase(float *coef);
            void buildCurve();
            void buildSurface();
            void moveCell(int axis, bool inc);
            DispCell getDispCell();
            void getShowPos(float position[3], float showPos[3]);
            void reverseShowPos(float position[3], float showPos[3]);
            int getCellX();
            int getCellY();
            int getCellZ();
            void setCellX(int value);
            void setCellY(int value);
            void setCellZ(int value);
            void saveDisplay(bool resultDebug);
            void loadDisplay(bool resultDebug);
            void generateCubicVolume();
            void saveOFF(ofstream &ofs);
            void loadOFF(ifstream &ifs);
            void saveVPT(ofstream &ofs);
            void loadVPT(ifstream &ifs);
            Face* getFaces();
            void setVolumeData(VolumeData * vData);
            //void hideSegments();

            char dataName[256];
            int *totalEdgePoints, *totalFacePoints, *totalCellPoints;
            float maxIntensity, minIntensity, maxScalar, minScalar, *maxEigenvalue, *minEigenvalue;
            int kernelSize, height, width, slices, dataType;
            float reverse, thickness, resolution;
            bool allCubic, resultDebug;
            float halfDataSize, rightBackTop[3], isovalue;
            int gridx, gridy, gridz;
            int sizex, sizey, sizez, halfSize;
            static const int subdNum = 3;
            VolumeData * volData;
            concurrent_vector<Point> *conpoints;
            concurrent_vector<Vertex> *convertices;
            EigenVectorsAndValues3D tensorMVectorsAndValues;

        protected:
            virtual void getScalarGP(int index, float *scalar) {};
            virtual void getTensorGP(int index, float tensor[6]) {};
            virtual void getGradientGP(int index, float gradient[3]) {};
            virtual void getScalar(float position[3], float *scalar) {};
            virtual void getTensor(float position[3], float tensor[6]) {};
            virtual void getGradient(float position[3], float gradient[3]) {};
            virtual void getScalarCubic(float position[3], float *scalar) {};
            virtual void getTensorCubic(float position[3], float tensor[6]) {};
            virtual void getGradientCubic(float position[3], float gradient[3]) {};
            virtual void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6]) {};
            virtual void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3]) {};
            virtual void getGradientBicubicTable(int axis, int i, int j, int k, float sGradients[subdNum][subdNum][3]) {};

            int maxGrid, maxGridIndex;
            float *gridPoints;
            float gridSize;
            float maxGradientNorm;

        private:

            bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            bool getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            bool getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            bool getEigensolver2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2);
            bool getEigensolverCubic2(float p1[3], float p2[3], SelfAdjointEigenSolver<Matrix3f> *e1, SelfAdjointEigenSolver<Matrix3f> *e2);
            void computeGridShowPositions();
            void computeGridScalars();
            void computeGridVectors();
            void computeGridGradients();
            void computeEdgePoints();
            void computeFacePoints();
            void computeCellPoints();
            void getGridPointPos(int index, float v[3]);
            void getGridPointPosD(int index, float v[3]);
            void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
            void getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
            void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]);
            float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3], float eigenValue1, float eigenValue2, std::vector<float> eigenVectors);
            std::vector<float> getEigenVectors(SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            std::vector<float> getFirstTwoEigenvalues(SelfAdjointEigenSolver<Matrix3f> *eigensolver);
            void generateCubicMRC();
            void generateCubicDTI();


            void edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k);
            int getEdgeIndex(int axis, int x, int y, int z);
            void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC);
            float getDihedral(float b1[3], float b2[3], float b3[3]);


            bool adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
                float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG);
            /*bool adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
            float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG);*/
            bool adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
                float v3_1[3], float v3_2[3], std::list<float> *GL, std::list<float> *V1L, std::list<float> *V3L, int depth);
            /*bool adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
            float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth);*/
            void extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum);
            void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum);
            void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum);
            void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum);
            void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum);


            void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4,
                int index1, int index2, int index3, int index4, int axis, int i, int j, int k);
            int getFaceIndex(int axis, int x, int y, int z);
            void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int index1, int index2, int index3, int index4, int axis);
            void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h);
            void rotate2D(float p[2], float theta);
            void getCurveType(Face *face);
            void get3DWindingNum(Face *face, int ii, int jj, int kk, int axis);
            float signedArea3D(float v1[3], float v2[3], float v3[3]);


            void cellPhaseIteration(int i, int j, int k);
            bool getCentroid(int i, int j, int k, Cell *cell);
            void centroidProjection(int i, int j, int k, Cell *cell);
            void getProjDirection(int type, float position[3], float q[3]);
            float smallAbsA(float A);
            void buildCurveIteration(Face *face, Cell *cell1, Cell *cell2);
            void buildSurfaceIteration(Edge *edge, Cell *cell1, Cell *cell2, Cell *cell3, Cell *cell4);
            void buildPointIteration(Cell *cell);


            bool dispAdaptiveSampling(float p1[3], float p2[3], std::vector<float> *sampleP,
                std::vector<float> *sampleG, std::vector<float> *sampleV1, std::vector<float> *sampleV3);
            void dispAdaptiveRecursion(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3], float v3_1[3],
                float v3_2[3], std::vector<float> *PL, std::vector<float> *GL, std::vector<float> *V1L, std::vector<float> *V3L, int depth);
            void dispOrientV3(std::vector<float> *sampleV3);
            void dispPropagateXY(std::vector<float> *sampleV3, std::vector<float> *sampleX, std::vector<float> *sampleY);
            void dispCombineAB(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, float *sumAB);


            float *gridScalars;
            float *gridShowPositions, *gridV1s, *gridV2s, *gridV3s, *gridGradients,
                *edgePoints, *facePoints, *cellPoints;
            bool storeGridShowPos, storeGridVector, storeGridGradient, storeGridScalar,
                storeEdgePoint, storeFacePoint, storeCellPoint, storeDispCell;
            int maxEdgeIndex, maxFaceIndex, maxCellIndex;
            Edge *edges;
            Face *faces;
            Cell *cells;
            float change;
            float globalVec[3];
            float *grid_2DGradMag;
            bool *storeGrid_2DGradMag;
            std::vector<Segment> *segments;
            std::vector<Quad> *quads;
            std::vector<Point> *points;
            std::vector<Vertex> *vertices;
            float thresh;
            int dispCellX, dispCellY, dispCellZ;
            DispCell dispCell;

            float pointRatio, curveRatio, surfaceRatio, localIntensityThreshMinG, localIntensityThreshMaxG, eigenvalueThresh;
            bool showGridScalar, showGridV1, showGridV2, showGridV3, showGridGradient, showEdgePoint, showFacePoint, showCellPoint,
                saddleCurve, saddlePoint, hideCurve, hideSurface, hidePoint, minCurve, minSurface, minPoint,
                maxCurve, maxSurface, maxPoint, showCell, showSample, showFaceSampleGradient,
                showIsosurface, hideSaliencyRatio, hideLocalIntensity, hideEigenvalue,
                isosurfaceFace, cylinderCurve, shadingFace, shadingWireframe;
        };

        void General_Data::getGridPointPosD(int index, float v[3])
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        static int getIndex(int type, int axis, int x, int y, int z, int sizex, int sizey, int sizez) {
            return type * (x * sizey * sizez + y * sizez + z) + axis;
        }

        static void tricubic_3f_table(int sizex, int sizey, int sizez, float *table, int axis, int ii, int jj, int kk, int index, float *data, float value[3])
        {
            float validData[4][3];
            int count = 0;
            switch (axis)
            {
            case 1:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(3, 0, ii + i - 1, jj, kk, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    count++;
                }
                break;
            case 2:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(3, 0, ii, jj + i - 1, kk, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    count++;
                }
                break;
            case 3:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(3, 0, ii, jj, kk + i - 1, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    count++;
                }
                break;
            }
            value[0] = 0;
            value[1] = 0;
            value[2] = 0;
            int tableIndex = 4 * index;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    value[j] += table[tableIndex + i] * validData[i][j];
                }
            }
        }

        static float getNorm2(float v[2]) {
            return sqrt(v[0] * v[0] + v[1] * v[1]);
        }

        static float getNorm6(float v[6]) {
            return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);
        }
        /**
        static float getDotP(float v1[3], float v2[3]) {
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        }
        **/
        static void getCrossP(float v1[3], float v2[3], float v[3]) {
            v[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }

        static float getNorm(float v[3]) {
            return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }

        static void getTriangleNorm(float v1[3], float v2[3], float v3[3], float n[3]) {
            float p1[3], p2[3];
            p1[0] = v1[0] - v2[0];
            p1[1] = v1[1] - v2[1];
            p1[2] = v1[2] - v2[2];
            p2[0] = v1[0] - v3[0];
            p2[1] = v1[1] - v3[1];
            p2[2] = v1[2] - v3[2];
            getCrossP(p1, p2, n);
            float normN = getNorm(n);
            if (normN == 0) return;
            n[0] /= normN;
            n[1] /= normN;
            n[2] /= normN;
        }

        static void getQuadNorm(float v1[3], float v2[3], float v3[3], float v4[3], float n[3]) {
            /*float n1[3], n2[3];
            getTriangleNorm(v1, v2, v3, n1);
            getTriangleNorm(v3, v4, v1, n2);
            n[0] = n1[0] + n2[0];
            n[1] = n1[1] + n2[1];
            n[2] = n1[2] + n2[2];*/
            float p1[3], p2[3];

            /*p1[0] = v2[0] - v3[0];
            p1[1] = v2[1] - v3[1];
            p1[2] = v2[2] - v3[2];
            p2[0] = v1[0] - v4[0];
            p2[1] = v1[1] - v4[1];
            p2[2] = v1[2] - v4[2];*/

            p1[0] = v1[0] - v3[0];
            p1[1] = v1[1] - v3[1];
            p1[2] = v1[2] - v3[2];
            p2[0] = v2[0] - v4[0];
            p2[1] = v2[1] - v4[1];
            p2[2] = v2[2] - v4[2];
            getCrossP(p1, p2, n);

            float normN = getNorm(n);
            if (normN == 0) return;
            n[0] /= normN;
            n[1] /= normN;
            n[2] /= normN;

            /**
            float normal[3];
            normal[0] = 0.0;
            normal[1] = 0.0;
            normal[2] = 0.0;

            std::vector<float *> quadVertices;
            quadVertices.push_back(v1);
            quadVertices.push_back(v2);
            quadVertices.push_back(v3);
            quadVertices.push_back(v4);
            for(int i = 0; i < 4; i++)
            {
            int j = (i+1)%4;
            float * currentV = quadVertices[i];
            float * nextV = quadVertices[j];
            normal[0] += (currentV[1] - nextV[1]) * (currentV[2] + nextV[2]);
            normal[1] += (currentV[2] - nextV[2]) * (currentV[0] + nextV[0]);
            normal[2] += (currentV[0] - nextV[0]) * (currentV[1] + nextV[1]);
            }
            float normN = getNorm(normal);
            if (normN == 0) {
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 0.0;
            }
            else {
            normal[0] /= normN;
            normal[1] /= normN;
            normal[2] /= normN;
            n[0] = normal[0];
            n[1] = normal[1];
            n[2] = normal[2];
            }**/
        }



        static void linearInterpolate(float x1, float x2, float x, float y1[3], float y2[3], float y[3])
        {
            float mu = (x - x1) / (x2 - x1);
            y[0] = y1[0] * (1 - mu) + y2[0] * mu;
            y[1] = y1[1] * (1 - mu) + y2[1] * mu;
            y[2] = y1[2] * (1 - mu) + y2[2] * mu;
        }

        static float trilinear(int type, int axis, float x, float y, float z, int iIndx, int iIndy, int iIndz,
            int sizex, int sizey, int sizez, float *data)
        {
            float V000 = data[getIndex(type, axis, iIndx, iIndy, iIndz, sizex, sizey, sizez)];
            float V = V000 * (1 - x) * (1 - y) * (1 - z);
            float V100 = data[getIndex(type, axis, iIndx + 1, iIndy, iIndz, sizex, sizey, sizez)];
            V += V100 * x * (1 - y) * (1 - z);
            float V010 = data[getIndex(type, axis, iIndx, iIndy + 1, iIndz, sizex, sizey, sizez)];
            V += V010 * (1 - x) * y * (1 - z);
            float V001 = data[getIndex(type, axis, iIndx, iIndy, iIndz + 1, sizex, sizey, sizez)];
            V += V001 * (1 - x) * (1 - y) * z;
            float V101 = data[getIndex(type, axis, iIndx + 1, iIndy, iIndz + 1, sizex, sizey, sizez)];
            V += V101 * x * (1 - y) * z;
            float V011 = data[getIndex(type, axis, iIndx, iIndy + 1, iIndz + 1, sizex, sizey, sizez)];
            V += V011 * (1 - x) * y * z;
            float V110 = data[getIndex(type, axis, iIndx + 1, iIndy + 1, iIndz, sizex, sizey, sizez)];
            V += V110 * x * y * (1 - z);
            float V111 = data[getIndex(type, axis, iIndx + 1, iIndy + 1, iIndz + 1, sizex, sizey, sizez)];
            V += V111 * x * y * z;
            return V;
        }

        static void trilinear_f(int sizex, int sizey, int sizez, float *data, float position[3], float *V)
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            *V = (float)trilinear(1, 0, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
        }

        static void trilinear_3f(int sizex, int sizey, int sizez, float *data, float position[3], float V[3])
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            V[0] = (float)trilinear(3, 0, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[1] = (float)trilinear(3, 1, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[2] = (float)trilinear(3, 2, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
        }

        static void trilinear_6f(int sizex, int sizey, int sizez, float *data, float position[3], float V[6])
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            V[0] = trilinear(6, 0, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[1] = trilinear(6, 1, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[2] = trilinear(6, 2, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[3] = trilinear(6, 3, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[4] = trilinear(6, 4, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
            V[5] = trilinear(6, 5, x, y, z, iIndx, iIndy, iIndz, sizex, sizey, sizez, data);
        }

        static float cubic(float data[4], float x)
        {
            return (float)(x * x * x * (-0.5 * data[0] + 1.5 * data[1] - 1.5 * data[2] + 0.5 * data[3]) +
                x * x * (data[0] - 2.5 * data[1] + 2 * data[2] - 0.5 * data[3]) +
                x * (-0.5 * data[0] + 0.5 * data[2]) +
                data[1]);
        }

        static float bicubic(float data[4][4], float x, float y)
        {
            float interData[4];
            interData[0] = cubic(data[0], y);
            interData[1] = cubic(data[1], y);
            interData[2] = cubic(data[2], y);
            interData[3] = cubic(data[3], y);
            return cubic(interData, x);
        }

        static float tricubic(float data[4][4][4], float x, float y, float z)
        {
            float interData[4];
            interData[0] = bicubic(data[0], y, z);
            interData[1] = bicubic(data[1], y, z);
            interData[2] = bicubic(data[2], y, z);
            interData[3] = bicubic(data[3], y, z);
            return cubic(interData, x);
        }

        static void tricubic_f(int sizex, int sizey, int sizez, float *data, float position[3], float *value)
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            iIndx--;
            iIndy--;
            iIndz--;
            float validData0[4][4][4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++)
                    {
                        int index0 = getIndex(1, 0, iIndx + i, iIndy + j, iIndz + k, sizex, sizey, sizez);
                        validData0[i][j][k] = data[index0];
                    }
                    *value = (float)tricubic(validData0, x, y, z);
                }
            }
        }

        static void tricubic_3f(int sizex, int sizey, int sizez, float *data, float position[3], float value[3])
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            iIndx--;
            iIndy--;
            iIndz--;
            float validData0[4][4][4];
            float validData1[4][4][4];
            float validData2[4][4][4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++)
                    {
                        int index = getIndex(3, 0, iIndx + i, iIndy + j, iIndz + k, sizex, sizey, sizez);
                        validData0[i][j][k] = data[index];
                        validData1[i][j][k] = data[index + 1];
                        validData2[i][j][k] = data[index + 2];
                    }
                }
            }
            value[0] = (float)tricubic(validData0, x, y, z);
            value[1] = (float)tricubic(validData1, x, y, z);
            value[2] = (float)tricubic(validData2, x, y, z);
        }

        static void tricubic_6f(int sizex, int sizey, int sizez, float *data, float position[3], float value[6])
        {
            float indx, indy, indz, x, y, z;
            x = modf(position[0], &indx);
            y = modf(position[1], &indy);
            z = modf(position[2], &indz);
            int iIndx = (int)indx;
            int iIndy = (int)indy;
            int iIndz = (int)indz;
            iIndx--;
            iIndy--;
            iIndz--;
            float validData0[4][4][4];
            float validData1[4][4][4];
            float validData2[4][4][4];
            float validData3[4][4][4];
            float validData4[4][4][4];
            float validData5[4][4][4];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++)
                    {
                        int index = getIndex(6, 0, iIndx + i, iIndy + j, iIndz + k, sizex, sizey, sizez);
                        validData0[i][j][k] = data[index];
                        validData1[i][j][k] = data[index + 1];
                        validData2[i][j][k] = data[index + 2];
                        validData3[i][j][k] = data[index + 3];
                        validData4[i][j][k] = data[index + 4];
                        validData5[i][j][k] = data[index + 5];
                    }
                }
            }
            value[0] = tricubic(validData0, x, y, z);
            value[1] = tricubic(validData1, x, y, z);
            value[2] = tricubic(validData2, x, y, z);
            value[3] = tricubic(validData3, x, y, z);
            value[4] = tricubic(validData4, x, y, z);
            value[5] = tricubic(validData5, x, y, z);
        }





        static void cubicCoef(float coef[4], float x)
        {
            float x2 = x * x;
            float x3 = x2 * x;
            coef[0] = -0.5 * x3 + x2 - 0.5 * x;
            coef[1] = 1.5 * x3 - 2.5 * x2 + 1;
            coef[2] = -1.5 * x3 + 2 * x2 + 0.5 * x;
            coef[3] = 0.5 * x3 - 0.5 * x2;
        }

        static void buildEdgeTable(float *table, int level, int slot, float left, float right)
        {
            if (level > 10) return;
            float mid = (left + right) / 2;
            cubicCoef(&(table[slot * 4]), mid);
            float* currentSlotVal = &(table[slot * 4]);
            //cout << "slot " << slot << " " << currentSlotVal[0] << " " << currentSlotVal[1] << " " << currentSlotVal[2] << endl;
            buildEdgeTable(table, level + 1, slot * 2 - 1, left, mid);
            buildEdgeTable(table, level + 1, slot * 2, mid, right);
        }

        static void bicubicCoef(float coef[16], float x, float y)
        {
            float x2 = x * x;
            float x3 = x2 * x;
            float cx[4];
            cx[0] = -0.5 * x3 + x2 - 0.5 * x;
            cx[1] = 1.5 * x3 - 2.5 * x2 + 1;
            cx[2] = -1.5 * x3 + 2 * x2 + 0.5 * x;
            cx[3] = 0.5 * x3 - 0.5 * x2;
            float y2 = y * y;
            float y3 = y2 * y;
            float cy[4];
            cy[0] = -0.5 * y3 + y2 - 0.5 * y;
            cy[1] = 1.5 * y3 - 2.5 * y2 + 1;
            cy[2] = -1.5 * y3 + 2 * y2 + 0.5 * y;
            cy[3] = 0.5 * y3 - 0.5 * y2;
            int count = 0;
            for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                coef[count] = cx[i] * cy[j];
                count++;
            }
        }

        static void tricubic_6f_table(int sizex, int sizey, int sizez, float *table, int axis, int ii, int jj, int kk, int index, float *data, float value[6])
        {
            float validData[4][6];
            int count = 0;
            switch (axis)
            {
            case 1:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(6, 0, ii + i - 1, jj, kk, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    validData[count][3] = data[idx + 3];
                    validData[count][4] = data[idx + 4];
                    validData[count][5] = data[idx + 5];
                    count++;
                }
                break;
            case 2:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(6, 0, ii, jj + i - 1, kk, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    validData[count][3] = data[idx + 3];
                    validData[count][4] = data[idx + 4];
                    validData[count][5] = data[idx + 5];
                    count++;
                }
                break;
            case 3:
                for (int i = 0; i < 4; i++)
                {
                    int idx = getIndex(6, 0, ii, jj, kk + i - 1, sizex, sizey, sizez);
                    validData[count][0] = data[idx];
                    validData[count][1] = data[idx + 1];
                    validData[count][2] = data[idx + 2];
                    validData[count][3] = data[idx + 3];
                    validData[count][4] = data[idx + 4];
                    validData[count][5] = data[idx + 5];
                    count++;
                }
                break;
            }
            value[0] = 0;
            value[1] = 0;
            value[2] = 0;
            value[3] = 0;
            value[4] = 0;
            value[5] = 0;
            int tableIndex = 4 * index;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    value[j] += table[tableIndex + i] * validData[i][j];
                }
            }
        }

        class ParallelEdge
        {
        private:

            float thresh;
            int gridx, gridy, gridz;
            int *totalEdgePoints;
            int dataType;
            float* globalVec;
            Edge *edges;
            float *gridPoints;

            float change;
            int sizex, sizey, sizez, halfSize;
            float*scalars, *tensors, *gradients;
            float *edgeTable;
            float v1, v2, v3;


            void getGridPointPosD(int index, float v[3]) const;

            bool adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
                float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG, float *edgeTable) const;
            bool adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
                float v3_1[3], float v3_2[3], std::list<float> *GL, std::list<float> *V1L, std::list<float> *V3L, int depth, float *edgeTable) const;
            void extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints) const;
            void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum) const;
            void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum) const;
            void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum) const;
            void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum) const;

            void edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints, float *edgeTable) const;


            bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            bool getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            bool getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver, float *edgeTable) const;
            void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
            void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;
            void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC) const;
            float getDihedral(float b1[3], float b2[3], float b3[3]) const;
            int getEdgeIndex(int axis, int x, int y, int z) const;

        public:
            bool allCubic;
            bool *storeGrid_2DGradMag;
            float *grid_2DGradMag;
            ParallelEdge(int lgridx, int lgridy, int lgridz, int *ltotalEdgePoints, Edge *ledges, bool lallCubic, float lglobalVec[3], int ldataType, float *lgridPoints, float lchange, bool *lstoreGrid_2DGradMag, float *lgrid_2DGradMag, int lsizex, int lsizey, int lsizez, int lhalfSize, float* lscalars, float* ltensors, float* lgradients, float *ledgeTable)
            {
                gridx = lgridx;
                gridy = lgridy;
                gridz = lgridz;
                totalEdgePoints = ltotalEdgePoints;
                edges = ledges;
                allCubic = lallCubic;
                globalVec = lglobalVec;
                dataType = ldataType;
                gridPoints = lgridPoints;
                change = lchange;
                storeGrid_2DGradMag = lstoreGrid_2DGradMag;
                grid_2DGradMag = lgrid_2DGradMag;
                thresh = cos(20 * (float)M_PI / 180);
                sizex = lsizex;
                sizey = lsizey;
                sizez = lsizez;
                halfSize = lhalfSize;
                scalars = lscalars;
                tensors = ltensors;
                gradients = lgradients;
                edgeTable = ledgeTable;


            };

            void getScalarGP(int index, float *scalar) const
            {
                *scalar = scalars[index];
            }

            void getGradientGP(int index, float gradient[3]) const
            {
                gradient[0] = gradients[index];
                gradient[1] = gradients[index + 1];
                gradient[2] = gradients[index + 2];
            }

            void getTensorGP(int index, float tensor[6]) const
            {
                tensor[0] = tensors[index];
                tensor[1] = tensors[index + 1];
                tensor[2] = tensors[index + 2];
                tensor[3] = tensors[index + 3];
                tensor[4] = tensors[index + 4];
                tensor[5] = tensors[index + 5];
            }

            void getScalar(float position[3], float *scalar) const
            {
                trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradient(float position[3], float gradient[3]) const
            {
                trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensor(float position[3], float tensor[6]) const
            {
                trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
            };

            void getScalarCubic(float position[3], float *scalar) const
            {
                tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradientCubic(float position[3], float gradient[3]) const
            {
                tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensorCubic(float position[3], float tensor[6]) const
            {
                tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
            };

            void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3], float *edgeTable) const
            {
                tricubic_3f_table(sizex, sizey, sizez, edgeTable, axis, i + 2, j + 2, k + 2, index, gradients, gradient);
            };

            void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6], float *edgeTable) const
            {
                tricubic_6f_table(sizex, sizey, sizez, edgeTable, axis, i + 2, j + 2, k + 2, index, tensors, tensor);
            };
            // for tbb
            void operator()(const blocked_range3d<int>& r)	const {

                int *ltotalEdgePoints = totalEdgePoints;
                Edge *ledges = edges;

                for (int i = r.cols().begin(); i != r.cols().end(); ++i){
                    for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
                        for (int k = r.pages().begin(); k != r.pages().end(); ++k){
                            if (i + 1 == gridx)
                                continue;
                            int sIndex1 = getIndex(1, 0, i + 2, j + 2, k + 2, sizex, sizey, sizez);
                            int sIndex2 = getIndex(1, 0, i + 3, j + 2, k + 2, sizex, sizey, sizez);
                            int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                            int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
                            edgePhaseIteration(index1, index2, sIndex1, sIndex2, 1, i, j, k, ledges, ltotalEdgePoints, edgeTable);
                        }
                    }
                }
                for (int i = r.cols().begin(); i != r.cols().end(); ++i){
                    for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
                        for (int k = r.pages().begin(); k != r.pages().end(); ++k){
                            if (j + 1 == gridy)
                                continue;
                            int sIndex1 = getIndex(1, 0, i + 2, j + 2, k + 2, sizex, sizey, sizez);
                            int sIndex2 = getIndex(1, 0, i + 2, j + 3, k + 2, sizex, sizey, sizez);
                            int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                            int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
                            edgePhaseIteration(index1, index2, sIndex1, sIndex2, 2, i, j, k, ledges, ltotalEdgePoints, edgeTable);
                        }
                    }
                }
                for (int i = r.cols().begin(); i != r.cols().end(); ++i){
                    for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
                        for (int k = r.pages().begin(); k != r.pages().end(); ++k){
                            if (k + 1 == gridz)
                                continue;
                            int sIndex1 = getIndex(1, 0, i + 2, j + 2, k + 2, sizex, sizey, sizez);
                            int sIndex2 = getIndex(1, 0, i + 2, j + 2, k + 3, sizex, sizey, sizez);
                            int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                            int index2 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
                            edgePhaseIteration(index1, index2, sIndex1, sIndex2, 3, i, j, k, ledges, ltotalEdgePoints, edgeTable);
                        }
                    }
                }


            };


        };



        static float getDotP(float v1[3], float v2[3]) {
            return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        }

        static int sign(float val) {
            return (0 < val) - (val < 0);
        }

        void ParallelEdge::getGridPointPosD(int index, float v[3]) const
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        void ParallelEdge::edgePhaseIteration(int index1, int index2, int si1, int si2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints, float *edgeTable) const
        {
            float *sampleG, *sampleV1, *sampleV3, *X, *Y, *sampleProjG;
            int edgeSampleNum;
            if (!adaptiveSamplingArrayTable(si1, si2, axis, i, j, k, &edgeSampleNum, &sampleG, &sampleV1, &sampleV3, &X, &Y, &sampleProjG, edgeTable)) return;
            Edge *edge = &(ledges[getEdgeIndex(axis, i, j, k)]);
            edge->valid = true;
            extremalEdgeArray(index1, index2, edge, si1, si2, sampleV1, edgeSampleNum, ltotalEdgePoints);
            orientV3Array(edge, sampleV3, edgeSampleNum);
            propagateXYArray(edge, sampleV3, X, Y, edgeSampleNum);
            projectGArray(sampleG, X, Y, sampleProjG, edgeSampleNum);
            getWindingNumArray(edge, sampleProjG, index1, index2, edgeSampleNum);
            delete[] sampleG;
            delete[] sampleV1;
            delete[] sampleV3;
            delete[] X;
            delete[] Y;
            delete[] sampleProjG;
        }

        bool ParallelEdge::adaptiveSamplingArrayTable(int si1, int si2, int axis, int ii, int jj, int kk, int *edgeSampleNum, float **sampleG,
            float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG, float *edgeTable) const
        {
            float epsilon = 1e-12;
            int depth = 0;
            std::list<float> sampleGList, sampleV1List, sampleV3List;
            std::list<float> *sampleGLP, *sampleV1LP, *sampleV3LP;
            sampleGLP = &sampleGList;
            sampleV1LP = &sampleV1List;
            sampleV3LP = &sampleV3List;
            SelfAdjointEigenSolver<Matrix3f> eigensolver1, eigensolver2;
            float g1[3], g2[3];
            int gi1 = si1 * 3;
            int gi2 = si2 * 3;
            int ti1 = si1 * 6;
            int ti2 = si2 * 6;
            getGradientGP(gi1, g1);
            getGradientGP(gi2, g2);
            if (getNorm(g1) < epsilon || getNorm(g2) < epsilon) return false;
            if (!getEigensolverGP(ti1, &eigensolver1)) return false;
            if (!getEigensolverGP(ti2, &eigensolver2)) return false;
            float v1_1[3], v3_1[3], v1_2[3], v3_2[3];
            getV1(&eigensolver1, v1_1);
            getV3(&eigensolver1, v3_1);
            getV1(&eigensolver2, v1_2);
            getV3(&eigensolver2, v3_2);

            for (int j = 0; j < 3; j++)
            {
                sampleGLP->push_back(g1[j]);
                sampleV1LP->push_back(v1_1[j]);
                sampleV3LP->push_back(v3_1[j]);
            }
            if (!adaptiveRecursionArrayTable(axis, ii, jj, kk, 2, g1, g2, v1_1, v1_2, v3_1, v3_2, sampleGLP, sampleV1LP, sampleV3LP, depth, edgeTable)) return false;
            for (int j = 0; j < 3; j++)
            {
                sampleGLP->push_back(g2[j]);
                sampleV1LP->push_back(v1_2[j]);
                sampleV3LP->push_back(v3_2[j]);
            }
            //for (std::list<float>::iterator it=sampleGLP->begin(); it != sampleGLP->end(); ++it) {
            //std::cout << ' ' << *it;
            //}

            *edgeSampleNum = (int)(sampleGLP->size() / 3);
            *sampleG = new float[3 * (*edgeSampleNum)];
            *sampleV1 = new float[3 * (*edgeSampleNum)];
            *sampleV3 = new float[3 * (*edgeSampleNum)];
            *X = new float[3 * (*edgeSampleNum)];
            *Y = new float[3 * (*edgeSampleNum)];
            *sampleProjG = new float[3 * (*edgeSampleNum)];
            std::list<float>::iterator ig = sampleGLP->begin();
            std::list<float>::iterator iv1 = sampleV1LP->begin();
            std::list<float>::iterator iv3 = sampleV3LP->begin();
            for (int i = 0; ig != sampleGLP->end(); ig++, iv1++, iv3++, i++)
            {
                (*sampleG)[i] = *ig;
                (*sampleV1)[i] = *iv1;
                (*sampleV3)[i] = *iv3;
            }
            return true;
        }

        void ParallelEdge::extremalEdgeArray(int index1, int index2, Edge *edge, int si1, int si2, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints) const
        {
            for (int i = 1; i < edgeSampleNum; i++)
            {
                float fSign = (float)sign(getDotP(sampleV1 + 3 * i, sampleV1 + 3 * (i - 1)));
                if (fSign == 0) fSign = 1;
                sampleV1[3 * i] = fSign*sampleV1[3 * i];
                sampleV1[3 * i + 1] = fSign*sampleV1[3 * i + 1];
                sampleV1[3 * i + 2] = fSign*sampleV1[3 * i + 2];
            }
            float p1[3], p2[3];
            getGridPointPosD(index1, p1);
            getGridPointPosD(index2, p2);
            float g1[3], g2[3];
            getGradientGP(si1 * 3, g1);
            getGradientGP(si2 * 3, g2);
            float d1 = getDotP(sampleV1, g1);
            float d2 = getDotP(sampleV1 + 3 * (edgeSampleNum - 1), g2);
            if (sign(d1*d2) == -1)
            {
                (*ltotalEdgePoints)++;
                float p[3];
                linearInterpolate(d1, d2, 0, p1, p2, p);
                edge->edgePoint[0] = p[0];
                edge->edgePoint[1] = p[1];
                edge->edgePoint[2] = p[2];
                SelfAdjointEigenSolver<Matrix3f> eigensolver;
                getEigensolverCubic(p, &eigensolver);

                float v1[3];
                float p10[3], p11[3], p13[3], p14[3];
                getV1(&eigensolver, v1);
                p10[0] = p[0] - 2 * change*v1[0];
                p10[1] = p[1] - 2 * change*v1[1];
                p10[2] = p[2] - 2 * change*v1[2];
                p11[0] = p[0] - change*v1[0];
                p11[1] = p[1] - change*v1[1];
                p11[2] = p[2] - change*v1[2];
                p13[0] = p[0] + change*v1[0];
                p13[1] = p[1] + change*v1[1];
                p13[2] = p[2] + change*v1[2];
                p14[0] = p[0] + 2 * change*v1[0];
                p14[1] = p[1] + 2 * change*v1[1];
                p14[2] = p[2] + 2 * change*v1[2];

                float f10, f11, f13, f14, f;
                getScalarCubic(p10, &f10);
                getScalarCubic(p11, &f11);
                getScalarCubic(p13, &f13);
                getScalarCubic(p14, &f14);
                getScalarCubic(p, &f);
                edge->extremal = true;
                float secondDiff = -f10 + 16 * f11 - 30 * f + 16 * f13 - f14;
                if (secondDiff < 0) edge->edgeTag = 1;
                else edge->edgeTag = 2;
            }
        }

        void ParallelEdge::orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum) const
        {
            edge->f = 0;
            float oldEnd[3];
            oldEnd[0] = sampleV3[3 * (edgeSampleNum - 1)];
            oldEnd[1] = sampleV3[3 * (edgeSampleNum - 1) + 1];
            oldEnd[2] = sampleV3[3 * (edgeSampleNum - 1) + 2];
            for (int i = 1; i < edgeSampleNum; i++)
            {
                float fSign = (float)sign(getDotP(sampleV3 + 3 * i, sampleV3 + 3 * (i - 1)));
                if (fSign == 0) fSign = 1;
                sampleV3[3 * i + 0] = fSign*sampleV3[3 * i];
                sampleV3[3 * i + 1] = fSign*sampleV3[3 * i + 1];
                sampleV3[3 * i + 2] = fSign*sampleV3[3 * i + 2];
            }
            if (sign(getDotP(oldEnd, sampleV3 + 3 * (edgeSampleNum - 1))) == -1) edge->f = 1;
        }
        void ParallelEdge::propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum) const
        {
            float epsilon = 1e-12;
            float ox[3], oz[3];
            ox[0] = 1;
            ox[1] = 0;
            ox[2] = 0;
            oz[0] = 0;
            oz[1] = 0;
            oz[2] = 1;
            for (int i = 0; i < edgeSampleNum; i++)
            {
                float *v = sampleV3 + 3 * i;
                float axle[3];
                getCrossP(oz, v, axle);
                float axleNorm = getNorm(axle);
                if (axleNorm > epsilon)
                {
                    axle[0] = axle[0] / axleNorm;
                    axle[1] = axle[1] / axleNorm;
                    axle[2] = axle[2] / axleNorm;
                    float angleCos = getDotP(oz, v);
                    if (angleCos > 1) angleCos = 1;
                    if (angleCos < -1) angleCos = -1;
                    float angle = acos(angleCos);
                    float a = getDotP(ox, axle)*(1 - angleCos);
                    float crossP[3];
                    getCrossP(axle, ox, crossP);
                    float b = sin(angle);
                    ox[0] = ox[0] * angleCos + axle[0] * a + crossP[0] * b;
                    ox[1] = ox[1] * angleCos + axle[1] * a + crossP[1] * b;
                    ox[2] = ox[2] * angleCos + axle[2] * a + crossP[2] * b;
                    oz[0] = v[0];
                    oz[1] = v[1];
                    oz[2] = v[2];
                }
                X[3 * i] = ox[0];
                X[3 * i + 1] = ox[1];
                X[3 * i + 2] = ox[2];
                getCrossP(v, ox, Y + 3 * i);
            }
            float area = 0;
            float phiC = 0;
            for (int i = 0; i < edgeSampleNum - 1; i++)
            {
                float *p1 = sampleV3 + 3 * i;
                float *p2 = sampleV3 + 3 * (i + 1);
                float incArea, incPhiC;
                signedSphericalTriangleArea(p1, p2, &incArea, &incPhiC);
                area += incArea;
                phiC += incPhiC;
            }
            edge->a = area;
            edge->b = phiC;
        }

        void ParallelEdge::projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum) const
        {
            for (int i = 0; i < edgeSampleNum; i++)
            {
                float *g = sampleG + 3 * i;
                float *x = X + 3 * i;
                float *y = Y + 3 * i;
                sampleProjG[3 * i] = getDotP(g, x);
                sampleProjG[3 * i + 1] = getDotP(g, y);
                sampleProjG[3 * i + 2] = 0;
            }
        }

        void ParallelEdge::getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum) const
        {
            float totalAngle = 0;
            float epsilon = 1e-12;
            float *end1 = sampleProjG;
            for (int i = 1; i < edgeSampleNum; i++)
            {
                float *end2 = sampleProjG + 3 * i;
                /*float ang1, ang2;
                if (end1[0] != 0)
                {
                ang1 = atan(end1[1]/end1[0]);
                }
                else
                {
                ang1 = (float)M_PI/2;
                }
                if (end2[0] != 0)
                {
                ang2 = atan(end2[1]/end2[0]);
                }
                else
                {
                ang2 = (float)M_PI/2;
                }
                if (end1[0] < 0) ang1 = ang1-(float)M_PI;
                if (end2[0] < 0) ang2 = ang2+(float)M_PI;
                float localAngle = ang2-ang1;
                if (localAngle > M_PI) localAngle = localAngle-2*(float)M_PI;
                if (localAngle < -M_PI) localAngle = localAngle+2*(float)M_PI;*/
                float end1Norm = getNorm2(end1);
                float end2Norm = getNorm2(end2);
                if (end1Norm > 0 && end2Norm > 0)
                {
                    float angleCos = getDotP(end1, end2) / (end1Norm*end2Norm);
                    if (angleCos > 1) angleCos = 1;
                    if (angleCos < -1) angleCos = -1;
                    float localAngle = acos(angleCos);
                    float crossP[3];
                    getCrossP(end1, end2, crossP);
                    localAngle *= sign(crossP[2]);
                    totalAngle += localAngle;
                }
                end1 = end2;
            }
            edge->w = totalAngle;
            if (!storeGrid_2DGradMag[index1])
            {
                storeGrid_2DGradMag[index1] = true;
                grid_2DGradMag[index1] = getNorm(sampleProjG);
            }
            if (!storeGrid_2DGradMag[index2])
            {
                storeGrid_2DGradMag[index2] = true;
                grid_2DGradMag[index2] = getNorm(sampleProjG + 3 * (edgeSampleNum - 1));
            }
        }

        bool ParallelEdge::getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorCubic(position, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }
        bool ParallelEdge::getEigensolver(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensor(position, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }
        bool ParallelEdge::getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorGP(index, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }
        bool ParallelEdge::getEigensolverCubicTable(int axis, int i, int j, int k, int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver, float *edgeTable) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorCubicTable(axis, i, j, k, index, tensorV, edgeTable);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }

        void ParallelEdge::getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev;
            switch (dataType)
            {
            case 1:
                ev = eigensolver->eigenvectors().col(2).normalized();
                break;
            case 2:
                ev = eigensolver->eigenvectors().col(0).normalized();
                break;
            }
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }
        void ParallelEdge::getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev;
            switch (dataType)
            {
            case 1:
                ev = eigensolver->eigenvectors().col(0).normalized();
                break;
            case 2:
                ev = eigensolver->eigenvectors().col(2).normalized();
                break;
            }
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }

        void ParallelEdge::signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC) const
        {
            *phiC = getDihedral(p1, globalVec, p2);
            float phiP1 = getDihedral(p2, p1, globalVec);
            float phiP2 = getDihedral(globalVec, p2, p1);
            *area = *phiC + phiP1 + phiP2 - (float)M_PI;
            float first[3], second[3], up[3];
            first[0] = p1[0] - globalVec[0];
            first[1] = p1[1] - globalVec[1];
            first[2] = p1[2] - globalVec[2];
            second[0] = p2[0] - p1[0];
            second[1] = p2[1] - p1[1];
            second[2] = p2[2] - p1[2];
            getCrossP(first, second, up);
            int aSign = sign(getDotP(up, globalVec));
            *area *= aSign;
            *phiC *= aSign;
        }
        float ParallelEdge::getDihedral(float b1[3], float b2[3], float b3[3]) const
        {
            float b21[3], b23[3];
            getCrossP(b2, b1, b21);
            getCrossP(b2, b3, b23);
            float norm = getNorm(b21);
            if (norm != 0)
            {
                b21[0] = b21[0] / norm;
                b21[1] = b21[1] / norm;
                b21[2] = b21[2] / norm;
            }
            norm = getNorm(b23);
            if (norm != 0)
            {
                b23[0] = b23[0] / norm;
                b23[1] = b23[1] / norm;
                b23[2] = b23[2] / norm;
            }
            float dotP = getDotP(b21, b23);
            if (dotP > 1) dotP = 1;
            if (dotP < -1) dotP = -1;
            float phi = acos(dotP);
            return phi;
        }
        int ParallelEdge::getEdgeIndex(int axis, int x, int y, int z) const {
            int index;
            switch (axis)
            {
            case 1:
                index = x * gridy * gridz + y * gridz + z;
                return index;
                break;
            case 2:
                index = (gridx - 1) * gridy * gridz +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            case 3:
                index = (gridx - 1) * gridy * gridz +
                    gridx * (gridy - 1) * gridz +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            }
            return 0;
        }

        bool ParallelEdge::adaptiveRecursionArrayTable(int axis, int ii, int jj, int kk, int index, float g1[3], float g2[3], float v1_1[3], float v1_2[3],
            float v3_1[3], float v3_2[3], std::list<float> *GL, std::list<float> *V1L, std::list<float> *V3L, int depth, float *edgeTable) const
        {
            depth++;
            if (depth > 10) return true;
            float cosG, cosV1, cosV3;
            float normG = getNorm(g1)*getNorm(g2);
            float normV1 = getNorm(v1_1)*getNorm(v1_2);
            float normV3 = getNorm(v3_1)*getNorm(v3_2);
            if (normG == 0) cosG = 1;
            else cosG = getDotP(g1, g2) / normG;
            if (normV1 == 0) cosV1 = 1;
            else cosV1 = getDotP(v1_1, v1_2) / normV1;
            if (normV3 == 0) cosV3 = 1;
            else cosV3 = getDotP(v3_1, v3_2) / normV3;
            if (cosG > thresh && fabs(cosV1) > thresh && fabs(cosV3) > thresh) return true;

            /*float p[3];
            p[0] = (p1[0]+p2[0])/2;
            p[1] = (p1[1]+p2[1])/2;
            p[2] = (p1[2]+p2[2])/2;*/
            float g[3], v1[3], v3[3];

            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getGradientCubicTable(axis, ii, jj, kk, index, g, edgeTable);
            if (!getEigensolverCubicTable(axis, ii, jj, kk, index, &eigensolver, edgeTable)) return false;
            /*if ( allCubic )
            {
            getGradientCubic(p, g);
            if ( !getEigensolverCubic(p, &eigensolver) ) return false;
            }
            else
            {
            getGradient(p, g);
            if ( !getEigensolver(p, &eigensolver) ) return false;
            }*/
            getV1(&eigensolver, v1);
            getV3(&eigensolver, v3);

            if (!adaptiveRecursionArrayTable(axis, ii, jj, kk, index * 2 - 1, g1, g, v1_1, v1, v3_1, v3, GL, V1L, V3L, depth, edgeTable)) return false;
            for (int j = 0; j < 3; j++)
            {
                GL->push_back(g[j]);
                V1L->push_back(v1[j]);
                V3L->push_back(v3[j]);
            }
            if (!adaptiveRecursionArrayTable(axis, ii, jj, kk, index * 2, g, g2, v1, v1_2, v3, v3_2, GL, V1L, V3L, depth, edgeTable)) return false;
            return true;
        }
        /**
        class CGNSParallelEdge
        {
        private:
        float thresh;
        int gridx,gridy,gridz;
        int *totalEdgePoints;
        int dataType;
        float* globalVec;
        Edge *edges;
        float *gridPoints;
        float change;
        float *coef;

        void edgePhaseBegin(int cols, int rows, int pages);
        void getGridPointPosD(int index, float v[3]);
        void edgePhaseIteration(int index1, int index2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints);
        bool adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
        float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG);
        bool adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
        float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth);
        void extremalEdgeArray(float p1[3], float p2[3], Edge *edge, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints);
        void orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum);
        void propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum);
        void projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum);
        void getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum);
        void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC);
        float getDihedral(float b1[3], float b2[3], float b3[3]);
        int getEdgeIndex(int axis, int x, int y, int z);

        public:
        bool getEigensolverPN(float position[3]);
        //void signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC);
        void getV1(float v[3]);
        void getV3(float v[3]);
        bool allCubic;
        bool *storeGrid_2DGradMag;
        float *grid_2DGradMag;
        EigenVectorsAndValues3D tensorMVectorsAndValues;

        CGNSParallelEdge(int lgridx,int lgridy,int lgridz,int *ltotalEdgePoints,Edge *ledges,bool lallCubic,float lglobalVec[3],
        int ldataType,float *lgridPoints,float lchange,bool *lstoreGrid_2DGradMag,float *lgrid_2DGradMag, float *lcoef)
        {
        gridx=lgridx;
        gridy=lgridy;
        gridz=lgridz;
        totalEdgePoints=ltotalEdgePoints;
        edges = ledges;
        allCubic=lallCubic;
        globalVec = lglobalVec;
        dataType = ldataType;
        gridPoints = lgridPoints;
        change=lchange;
        storeGrid_2DGradMag = lstoreGrid_2DGradMag;
        grid_2DGradMag = lgrid_2DGradMag;
        thresh = cos(20 * (float)M_PI / 180);
        coef = lcoef;
        };
        void getScalarPN(float position[3], float *scalar)
        {
        float x = position[0];
        float y = position[1];
        float z = position[2];
        *scalar = coef[56]*pow(x,6) + coef[63]*pow(x,5)*y + coef[57]*pow(x,5)*z + coef[35]*pow(x,5) + coef[69]*pow(x,4)*pow(y,2) + coef[64]*pow(x,4)*y*z +
        coef[41]*pow(x,4)*y + coef[58]*pow(x,4)*pow(z,2) + coef[36]*pow(x,4)*z + coef[20]*pow(x,4) + coef[74]*pow(x,3)*pow(y,3) + coef[70]*pow(x,3)*pow(y,2)*z +
        coef[46]*pow(x,3)*pow(y,2) + coef[65]*pow(x,3)*y*pow(z,2) + coef[42]*pow(x,3)*y*z + coef[25]*pow(x,3)*y + coef[59]*pow(x,3)*pow(z,3) + coef[37]*pow(x,3)*pow(z,2) +
        coef[21]*pow(x,3)*z + coef[10]*pow(x,3) + coef[78]*pow(x,2)*pow(y,4) + coef[75]*pow(x,2)*pow(y,3)*z + coef[50]*pow(x,2)*pow(y,3) + coef[71]*pow(x,2)*pow(y,2)*pow(z,2) +
        coef[47]*pow(x,2)*pow(y,2)*z + coef[29]*pow(x,2)*pow(y,2) + coef[66]*pow(x,2)*y*pow(z,3) + coef[43]*pow(x,2)*y*pow(z,2) + coef[26]*pow(x,2)*y*z +
        coef[14]*pow(x,2)*y + coef[60]*pow(x,2)*pow(z,4) + coef[38]*pow(x,2)*pow(z,3) + coef[22]*pow(x,2)*pow(z,2) + coef[11]*pow(x,2)*z + coef[4]*pow(x,2) +
        coef[81]*x*pow(y,5) + coef[79]*x*pow(y,4)*z + coef[53]*x*pow(y,4) + coef[76]*x*pow(y,3)*pow(z,2) + coef[51]*x*pow(y,3)*z + coef[32]*x*pow(y,3) +
        coef[72]*x*pow(y,2)*pow(z,3) + coef[48]*x*pow(y,2)*pow(z,2) + coef[30]*x*pow(y,2)*z + coef[17]*x*pow(y,2) + coef[67]*x*y*pow(z,4) + coef[44]*x*y*pow(z,3) +
        coef[27]*x*y*pow(z,2) + coef[15]*x*y*z + coef[7]*x*y + coef[61]*x*pow(z,5) + coef[39]*x*pow(z,4) + coef[23]*x*pow(z,3) +
        coef[12]*x*pow(z,2) + coef[5]*x*z + coef[1]*x + coef[83]*pow(y,6) + coef[82]*pow(y,5)*z + coef[55]*pow(y,5) + coef[80]*pow(y,4)*pow(z,2) +
        coef[54]*pow(y,4)*z + coef[34]*pow(y,4) + coef[77]*pow(y,3)*pow(z,3) + coef[52]*pow(y,3)*pow(z,2) + coef[33]*pow(y,3)*z + coef[19]*pow(y,3) +
        coef[73]*pow(y,2)*pow(z,4) + coef[49]*pow(y,2)*pow(z,3) + coef[31]*pow(y,2)*pow(z,2) + coef[18]*pow(y,2)*z + coef[9]*pow(y,2) + coef[68]*y*pow(z,5) +
        coef[45]*y*pow(z,4) + coef[28]*y*pow(z,3) + coef[16]*y*pow(z,2) + coef[8]*y*z + coef[3]*y + coef[62]*pow(z,6) + coef[40]*pow(z,5) +
        coef[24]*pow(z,4) + coef[13]*pow(z,3) + coef[6]*pow(z,2) + coef[2]*z + coef[0];
        };

        void getGradientPN(float position[3], float gradient[3])
        {
        float x = position[0];
        float y = position[1];
        float z = position[2];
        gradient[0] = 6*coef[56]*pow(x,5) + 5*coef[63]*pow(x,4)*y + 5*coef[57]*pow(x,4)*z + 5*coef[35]*pow(x,4) + 4*coef[69]*pow(x,3)*pow(y,2) +
        4*coef[64]*pow(x,3)*y*z + 4*coef[41]*pow(x,3)*y + 4*coef[58]*pow(x,3)*pow(z,2) + 4*coef[36]*pow(x,3)*z + 4*coef[20]*pow(x,3) + 3*coef[74]*pow(x,2)*pow(y,3) +
        3*coef[70]*pow(x,2)*pow(y,2)*z + 3*coef[46]*pow(x,2)*pow(y,2) + 3*coef[65]*pow(x,2)*y*pow(z,2) + 3*coef[42]*pow(x,2)*y*z + 3*coef[25]*pow(x,2)*y +
        3*coef[59]*pow(x,2)*pow(z,3) + 3*coef[37]*pow(x,2)*pow(z,2) + 3*coef[21]*pow(x,2)*z + 3*coef[10]*pow(x,2) + 2*coef[78]*x*pow(y,4) + 2*coef[75]*x*pow(y,3)*z +
        2*coef[50]*x*pow(y,3) + 2*coef[71]*x*pow(y,2)*pow(z,2) + 2*coef[47]*x*pow(y,2)*z + 2*coef[29]*x*pow(y,2) + 2*coef[66]*x*y*pow(z,3) +
        2*coef[43]*x*y*pow(z,2) + 2*coef[26]*x*y*z + 2*coef[14]*x*y + 2*coef[60]*x*pow(z,4) + 2*coef[38]*x*pow(z,3) + 2*coef[22]*x*pow(z,2) +
        2*coef[11]*x*z + 2*coef[4]*x + coef[81]*pow(y,5) + coef[79]*pow(y,4)*z + coef[53]*pow(y,4) + coef[76]*pow(y,3)*pow(z,2) + coef[51]*pow(y,3)*z +
        coef[32]*pow(y,3) + coef[72]*pow(y,2)*pow(z,3) + coef[48]*pow(y,2)*pow(z,2) + coef[30]*pow(y,2)*z + coef[17]*pow(y,2) + coef[67]*y*pow(z,4) + coef[44]*y*pow(z,3) +
        coef[27]*y*pow(z,2) + coef[15]*y*z + coef[7]*y + coef[61]*pow(z,5) + coef[39]*pow(z,4) + coef[23]*pow(z,3) + coef[12]*pow(z,2) + coef[5]*z + coef[1];
        gradient[1] = coef[63]*pow(x,5) + 2*coef[69]*pow(x,4)*y + coef[64]*pow(x,4)*z + coef[41]*pow(x,4) + 3*coef[74]*pow(x,3)*pow(y,2) + 2*coef[70]*pow(x,3)*y*z +
        2*coef[46]*pow(x,3)*y + coef[65]*pow(x,3)*pow(z,2) + coef[42]*pow(x,3)*z + coef[25]*pow(x,3) + 4*coef[78]*pow(x,2)*pow(y,3) + 3*coef[75]*pow(x,2)*pow(y,2)*z +
        3*coef[50]*pow(x,2)*pow(y,2) + 2*coef[71]*pow(x,2)*y*pow(z,2) + 2*coef[47]*pow(x,2)*y*z + 2*coef[29]*pow(x,2)*y + coef[66]*pow(x,2)*pow(z,3) + coef[43]*pow(x,2)*pow(z,2) +
        coef[26]*pow(x,2)*z + coef[14]*pow(x,2) + 5*coef[81]*x*pow(y,4) + 4*coef[79]*x*pow(y,3)*z + 4*coef[53]*x*pow(y,3) + 3*coef[76]*x*pow(y,2)*pow(z,2) +
        3*coef[51]*x*pow(y,2)*z + 3*coef[32]*x*pow(y,2) + 2*coef[72]*x*y*pow(z,3) + 2*coef[48]*x*y*pow(z,2) + 2*coef[30]*x*y*z + 2*coef[17]*x*y +
        coef[67]*x*pow(z,4) + coef[44]*x*pow(z,3) + coef[27]*x*pow(z,2) + coef[15]*x*z + coef[7]*x + 6*coef[83]*pow(y,5) + 5*coef[82]*pow(y,4)*z +
        5*coef[55]*pow(y,4) + 4*coef[80]*pow(y,3)*pow(z,2) + 4*coef[54]*pow(y,3)*z + 4*coef[34]*pow(y,3) + 3*coef[77]*pow(y,2)*pow(z,3) + 3*coef[52]*pow(y,2)*pow(z,2) +
        3*coef[33]*pow(y,2)*z + 3*coef[19]*pow(y,2) + 2*coef[73]*y*pow(z,4) + 2*coef[49]*y*pow(z,3) + 2*coef[31]*y*pow(z,2) + 2*coef[18]*y*z +
        2*coef[9]*y + coef[68]*pow(z,5) + coef[45]*pow(z,4) + coef[28]*pow(z,3) + coef[16]*pow(z,2) + coef[8]*z + coef[3];
        gradient[2] = coef[57]*pow(x,5) + coef[64]*pow(x,4)*y + 2*coef[58]*pow(x,4)*z + coef[36]*pow(x,4) + coef[70]*pow(x,3)*pow(y,2) + 2*coef[65]*pow(x,3)*y*z +
        coef[42]*pow(x,3)*y + 3*coef[59]*pow(x,3)*pow(z,2) + 2*coef[37]*pow(x,3)*z + coef[21]*pow(x,3) + coef[75]*pow(x,2)*pow(y,3) + 2*coef[71]*pow(x,2)*pow(y,2)*z +
        coef[47]*pow(x,2)*pow(y,2) + 3*coef[66]*pow(x,2)*y*pow(z,2) + 2*coef[43]*pow(x,2)*y*z + coef[26]*pow(x,2)*y + 4*coef[60]*pow(x,2)*pow(z,3) + 3*coef[38]*pow(x,2)*pow(z,2) +
        2*coef[22]*pow(x,2)*z + coef[11]*pow(x,2) + coef[79]*x*pow(y,4) + 2*coef[76]*x*pow(y,3)*z + coef[51]*x*pow(y,3) + 3*coef[72]*x*pow(y,2)*pow(z,2) +
        2*coef[48]*x*pow(y,2)*z + coef[30]*x*pow(y,2) + 4*coef[67]*x*y*pow(z,3) + 3*coef[44]*x*y*pow(z,2) + 2*coef[27]*x*y*z + coef[15]*x*y +
        5*coef[61]*x*pow(z,4) + 4*coef[39]*x*pow(z,3) + 3*coef[23]*x*pow(z,2) + 2*coef[12]*x*z + coef[5]*x + coef[82]*pow(y,5) + 2*coef[80]*pow(y,4)*z +
        coef[54]*pow(y,4) + 3*coef[77]*pow(y,3)*pow(z,2) + 2*coef[52]*pow(y,3)*z + coef[33]*pow(y,3) + 4*coef[73]*pow(y,2)*pow(z,3) + 3*coef[49]*pow(y,2)*pow(z,2) +
        2*coef[31]*pow(y,2)*z + coef[18]*pow(y,2) + 5*coef[68]*y*pow(z,4) + 4*coef[45]*y*pow(z,3) + 3*coef[28]*y*pow(z,2) + 2*coef[16]*y*z + coef[8]*y +
        6*coef[62]*pow(z,5) + 5*coef[40]*pow(z,4) + 4*coef[24]*pow(z,3) + 3*coef[13]*pow(z,2) + 2*coef[6]*z + coef[2];
        };

        void getTensorPN(float position[3], float tensor[6])
        {
        float x = position[0];
        float y = position[1];
        float z = position[2];
        tensor[0] = 30*coef[56]*pow(x,4) + 20*coef[63]*pow(x,3)*y + 20*coef[57]*pow(x,3)*z + 20*coef[35]*pow(x,3) + 12*coef[69]*pow(x,2)*pow(y,2) +
        12*coef[64]*pow(x,2)*y*z + 12*coef[41]*pow(x,2)*y + 12*coef[58]*pow(x,2)*pow(z,2) + 12*coef[36]*pow(x,2)*z + 12*coef[20]*pow(x,2) + 6*coef[74]*x*pow(y,3) +
        6*coef[70]*x*pow(y,2)*z + 6*coef[46]*x*pow(y,2) + 6*coef[65]*x*y*pow(z,2) + 6*coef[42]*x*y*z + 6*coef[25]*x*y + 6*coef[59]*x*pow(z,3) +
        6*coef[37]*x*pow(z,2) + 6*coef[21]*x*z + 6*coef[10]*x + 2*coef[78]*pow(y,4) + 2*coef[75]*pow(y,3)*z + 2*coef[50]*pow(y,3) +
        2*coef[71]*pow(y,2)*pow(z,2) + 2*coef[47]*pow(y,2)*z + 2*coef[29]*pow(y,2) + 2*coef[66]*y*pow(z,3) + 2*coef[43]*y*pow(z,2) + 2*coef[26]*y*z +
        2*coef[14]*y + 2*coef[60]*pow(z,4) + 2*coef[38]*pow(z,3) + 2*coef[22]*pow(z,2) + 2*coef[11]*z + 2*coef[4];
        tensor[1] = 5*coef[63]*pow(x,4) + 8*coef[69]*pow(x,3)*y + 4*coef[64]*pow(x,3)*z + 4*coef[41]*pow(x,3) + 9*coef[74]*pow(x,2)*pow(y,2) +
        6*coef[70]*pow(x,2)*y*z + 6*coef[46]*pow(x,2)*y + 3*coef[65]*pow(x,2)*pow(z,2) + 3*coef[42]*pow(x,2)*z + 3*coef[25]*pow(x,2) + 8*coef[78]*x*pow(y,3) +
        6*coef[75]*x*pow(y,2)*z + 6*coef[50]*x*pow(y,2) + 4*coef[71]*x*y*pow(z,2) + 4*coef[47]*x*y*z + 4*coef[29]*x*y + 2*coef[66]*x*pow(z,3) +
        2*coef[43]*x*pow(z,2) + 2*coef[26]*x*z + 2*coef[14]*x + 5*coef[81]*pow(y,4) + 4*coef[79]*pow(y,3)*z + 4*coef[53]*pow(y,3) +
        3*coef[76]*pow(y,2)*pow(z,2) + 3*coef[51]*pow(y,2)*z + 3*coef[32]*pow(y,2) + 2*coef[72]*y*pow(z,3) + 2*coef[48]*y*pow(z,2) + 2*coef[30]*y*z +
        2*coef[17]*y + coef[67]*pow(z,4) + coef[44]*pow(z,3) + coef[27]*pow(z,2) + coef[15]*z + coef[7];
        tensor[2] = 5*coef[57]*pow(x,4) + 4*coef[64]*pow(x,3)*y + 8*coef[58]*pow(x,3)*z + 4*coef[36]*pow(x,3) + 3*coef[70]*pow(x,2)*pow(y,2) +
        6*coef[65]*pow(x,2)*y*z + 3*coef[42]*pow(x,2)*y + 9*coef[59]*pow(x,2)*pow(z,2) + 6*coef[37]*pow(x,2)*z + 3*coef[21]*pow(x,2) + 2*coef[75]*x*pow(y,3) +
        4*coef[71]*x*pow(y,2)*z + 2*coef[47]*x*pow(y,2) + 6*coef[66]*x*y*pow(z,2) + 4*coef[43]*x*y*z + 2*coef[26]*x*y + 8*coef[60]*x*pow(z,3) +
        6*coef[38]*x*pow(z,2) + 4*coef[22]*x*z + 2*coef[11]*x + coef[79]*pow(y,4) + 2*coef[76]*pow(y,3)*z + coef[51]*pow(y,3) +
        3*coef[72]*pow(y,2)*pow(z,2) + 2*coef[48]*pow(y,2)*z + coef[30]*pow(y,2) + 4*coef[67]*y*pow(z,3) + 3*coef[44]*y*pow(z,2) + 2*coef[27]*y*z +
        coef[15]*y + 5*coef[61]*pow(z,4) + 4*coef[39]*pow(z,3) + 3*coef[23]*pow(z,2) + 2*coef[12]*z + coef[5];
        tensor[3] = 2*coef[69]*pow(x,4) + 6*coef[74]*pow(x,3)*y + 2*coef[70]*pow(x,3)*z + 2*coef[46]*pow(x,3) + 12*coef[78]*pow(x,2)*pow(y,2) +
        6*coef[75]*pow(x,2)*y*z + 6*coef[50]*pow(x,2)*y + 2*coef[71]*pow(x,2)*pow(z,2) + 2*coef[47]*pow(x,2)*z + 2*coef[29]*pow(x,2) + 20*coef[81]*x*pow(y,3) +
        12*coef[79]*x*pow(y,2)*z + 12*coef[53]*x*pow(y,2) + 6*coef[76]*x*y*pow(z,2) + 6*coef[51]*x*y*z + 6*coef[32]*x*y + 2*coef[72]*x*pow(z,3) +
        2*coef[48]*x*pow(z,2) + 2*coef[30]*x*z + 2*coef[17]*x + 30*coef[83]*pow(y,4) + 20*coef[82]*pow(y,3)*z + 20*coef[55]*pow(y,3) +
        12*coef[80]*pow(y,2)*pow(z,2) + 12*coef[54]*pow(y,2)*z + 12*coef[34]*pow(y,2) + 6*coef[77]*y*pow(z,3) + 6*coef[52]*y*pow(z,2) + 6*coef[33]*y*z +
        6*coef[19]*y + 2*coef[73]*pow(z,4) + 2*coef[49]*pow(z,3) + 2*coef[31]*pow(z,2) + 2*coef[18]*z + 2*coef[9];
        tensor[4] = coef[64]*pow(x,4) + 2*coef[70]*pow(x,3)*y + 2*coef[65]*pow(x,3)*z + coef[42]*pow(x,3) + 3*coef[75]*pow(x,2)*pow(y,2) + 4*coef[71]*pow(x,2)*y*z +
        2*coef[47]*pow(x,2)*y + 3*coef[66]*pow(x,2)*pow(z,2) + 2*coef[43]*pow(x,2)*z + coef[26]*pow(x,2) + 4*coef[79]*x*pow(y,3) + 6*coef[76]*x*pow(y,2)*z +
        3*coef[51]*x*pow(y,2) + 6*coef[72]*x*y*pow(z,2) + 4*coef[48]*x*y*z + 2*coef[30]*x*y + 4*coef[67]*x*pow(z,3) + 3*coef[44]*x*pow(z,2) +
        2*coef[27]*x*z + coef[15]*x + 5*coef[82]*pow(y,4) + 8*coef[80]*pow(y,3)*z + 4*coef[54]*pow(y,3) + 9*coef[77]*pow(y,2)*pow(z,2) +
        6*coef[52]*pow(y,2)*z + 3*coef[33]*pow(y,2) + 8*coef[73]*y*pow(z,3) + 6*coef[49]*y*pow(z,2) + 4*coef[31]*y*z + 2*coef[18]*y +
        5*coef[68]*pow(z,4) + 4*coef[45]*pow(z,3) + 3*coef[28]*pow(z,2) + 2*coef[16]*z + coef[8];
        tensor[5] = 2*coef[58]*pow(x,4) + 2*coef[65]*pow(x,3)*y + 6*coef[59]*pow(x,3)*z + 2*coef[37]*pow(x,3) + 2*coef[71]*pow(x,2)*pow(y,2) +
        6*coef[66]*pow(x,2)*y*z + 2*coef[43]*pow(x,2)*y + 12*coef[60]*pow(x,2)*pow(z,2) + 6*coef[38]*pow(x,2)*z + 2*coef[22]*pow(x,2) + 2*coef[76]*x*pow(y,3) +
        6*coef[72]*x*pow(y,2)*z + 2*coef[48]*x*pow(y,2) + 12*coef[67]*x*y*pow(z,2) + 6*coef[44]*x*y*z + 2*coef[27]*x*y + 20*coef[61]*x*pow(z,3) +
        12*coef[39]*x*pow(z,2) + 6*coef[23]*x*z + 2*coef[12]*x + 2*coef[80]*pow(y,4) + 6*coef[77]*pow(y,3)*z + 2*coef[52]*pow(y,3) +
        12*coef[73]*pow(y,2)*pow(z,2) + 6*coef[49]*pow(y,2)*z + 2*coef[31]*pow(y,2) + 20*coef[68]*y*pow(z,3) + 12*coef[45]*y*pow(z,2) + 6*coef[28]*y*z +
        2*coef[16]*y + 30*coef[62]*pow(z,4) + 20*coef[40]*pow(z,3) + 12*coef[24]*pow(z,2) + 6*coef[13]*z + 2*coef[6];
        };
        };



        void CGNSParallelEdge::getGridPointPosD(int index, float v[3])
        {
        v[0] = gridPoints[3*index];
        v[1] = gridPoints[3*index+1];
        v[2] = gridPoints[3*index+2];
        }



        int CGNSParallelEdge::getEdgeIndex(int axis, int x, int y, int z) {
        int index;
        switch (axis)
        {
        case 1:
        index = x * gridy * gridz + y * gridz + z;
        return index;
        break;
        case 2:
        index = (gridx - 1) * gridy * gridz +
        x * (gridy - 1) * gridz + y * gridz + z;
        return index;
        break;
        case 3:
        index = (gridx - 1) * gridy * gridz +
        gridx * (gridy - 1) * gridz +
        x * gridy * (gridz - 1) + y * (gridz - 1) + z;
        return index;
        break;
        }
        return 0;
        }



        bool CGNSParallelEdge::getEigensolverPN(float position[3])
        {
        int count = 0;
        float epsilon = 1e-12;
        float tensorV[6];
        getTensorPN(position, tensorV);
        for ( int i = 0; i < 6; i ++ )
        {
        if ( fabs(tensorV[i]) < epsilon )
        {
        count ++;
        tensorV[i] = epsilon;
        }
        }
        float tensorM[3][3];
        tensorM[0][0] = tensorV[0];
        tensorM[0][1] = tensorV[1];
        tensorM[0][2] = tensorV[2];
        tensorM[1][0] = tensorV[1];
        tensorM[1][1] = tensorV[3];
        tensorM[1][2] = tensorV[4];
        tensorM[2][0] = tensorV[2];
        tensorM[2][1] = tensorV[4];
        tensorM[2][2] = tensorV[5];
        //EigenVectorsAndValues3D eigenData;
        for(int r = 0; r < 3; r++) {
        for(int c = 0; c < 3; c++) {
        tensorMVectorsAndValues.structureTensor[r][c] = 0;

        }
        }
        for(int r = 0; r < 3; r++) {
        for(int c = 0; c < 3; c++) {
        tensorMVectorsAndValues.structureTensor[r][c] = tensorM[r][c];
        }
        }
        MathLib * math = new MathLib();
        math->EigenAnalysis(tensorMVectorsAndValues);
        if ( count > 1 ) return false;
        return true;
        }

        void CGNSParallelEdge::getV1(float v[3])
        {


        Vector3DFloat ev;
        Vector3DFloat currentEdge1(tensorMVectorsAndValues.structureTensor[2][0], tensorMVectorsAndValues.structureTensor[2][1], tensorMVectorsAndValues.structureTensor[2][2]);
        currentEdge1 = normalize(currentEdge1);
        ev = currentEdge1;


        v[0] = (float)ev.X();
        v[1] = (float)ev.Y();
        v[2] = (float)ev.Z();
        }

        void CGNSParallelEdge::getV3(float v[3])
        {


        Vector3DFloat ev;
        Vector3DFloat currentEdge1(tensorMVectorsAndValues.structureTensor[0][0], tensorMVectorsAndValues.structureTensor[0][1], tensorMVectorsAndValues.structureTensor[0][2]);
        currentEdge1 = normalize(currentEdge1);
        ev = currentEdge1;


        v[0] = (float)ev.X();
        v[1] = (float)ev.Y();
        v[2] = (float)ev.Z();
        }

        bool CGNSParallelEdge::adaptiveRecursionArray(float p1[3], float p2[3], float g1[3], float g2[3], float v1_1[3], float v1_2[3],
        float v3_1[3], float v3_2[3], list<float> *GL, list<float> *V1L, list<float> *V3L, int depth)
        {
        depth ++;
        if ( depth > 10 ) return true;
        float cosG, cosV1, cosV3;
        float normG = getNorm(g1)*getNorm(g2);
        float normV1 = getNorm(v1_1)*getNorm(v1_2);
        float normV3 = getNorm(v3_1)*getNorm(v3_2);
        if (normG == 0) cosG = 1;
        else cosG = getDotP(g1, g2)/normG;
        if (normV1 == 0) cosV1 = 1;
        else cosV1 = getDotP(v1_1, v1_2)/normV1;
        if (normV3 == 0) cosV3 = 1;
        else cosV3 = getDotP(v3_1, v3_2)/normV3;
        if ( cosG > thresh && fabs(cosV1) > thresh && fabs(cosV3) > thresh ) return true;

        float p[3];
        p[0] = (p1[0]+p2[0])/2;
        p[1] = (p1[1]+p2[1])/2;
        p[2] = (p1[2]+p2[2])/2;
        float g[3], v1[3], v3[3];


        //SelfAdjointEigenSolver<Matrix3f> eigensolver;
        getGradientPN(p, g);
        if ( !getEigensolverPN(p) ) return false;
        getV1(v1);
        getV3(v3);


        //if ( !getEigensolverPN(p, &eigensolver) ) return false;
        //getV1(&eigensolver, v1);
        //getV3(&eigensolver, v3);

        if ( !adaptiveRecursionArray(p1, p, g1, g, v1_1, v1, v3_1, v3, GL, V1L, V3L, depth) ) return false;
        for ( int j = 0; j < 3; j ++ )
        {
        GL->push_back(g[j]);
        V1L->push_back(v1[j]);
        V3L->push_back(v3[j]);
        }
        if ( !adaptiveRecursionArray(p, p2, g, g2, v1, v1_2, v3, v3_2, GL, V1L, V3L, depth) ) return false;

        return true;
        }

        bool CGNSParallelEdge::adaptiveSamplingArray(float p1[3], float p2[3], int *edgeSampleNum, float **sampleG,
        float **sampleV1, float **sampleV3, float **X, float **Y, float **sampleProjG)
        {
        float epsilon = 1e-12;
        int depth = 0;
        list<float> sampleGList, sampleV1List, sampleV3List;
        list<float> *sampleGLP, *sampleV1LP, *sampleV3LP;
        sampleGLP = &sampleGList;
        sampleV1LP = &sampleV1List;
        sampleV3LP = &sampleV3List;
        //SelfAdjointEigenSolver<Matrix3f> eigensolver1, eigensolver2;
        float g1[3], g2[3];
        getGradientPN(p1, g1);
        getGradientPN(p2, g2);
        if ( getNorm(g1) < epsilon || getNorm(g2) < epsilon ) return false;
        if ( !getEigensolverPN(p1) ) return false;
        if ( !getEigensolverPN(p2) ) return false;
        //if ( !getEigensolverPN(p1, &eigensolver1) ) return false;
        //	if ( !getEigensolverPN(p2, &eigensolver2) ) return false;
        float v1_1[3], v3_1[3], v1_2[3], v3_2[3];
        getV1(v1_1);
        getV3(v3_1);
        getV1(v1_2);
        getV3(v3_2);

        //getV1(&eigensolver1, v1_1);
        //getV3(&eigensolver1, v3_1);
        //getV1(&eigensolver2, v1_2);
        //getV3(&eigensolver2, v3_2);


        for ( int j = 0; j < 3; j ++ )
        {
        sampleGLP->push_back(g1[j]);
        sampleV1LP->push_back(v1_1[j]);
        sampleV3LP->push_back(v3_1[j]);
        }
        if ( !adaptiveRecursionArray(p1, p2, g1, g2, v1_1, v1_2, v3_1, v3_2, sampleGLP, sampleV1LP, sampleV3LP, depth) ) return false;
        for ( int j = 0; j < 3; j ++ )
        {
        sampleGLP->push_back(g2[j]);
        sampleV1LP->push_back(v1_2[j]);
        sampleV3LP->push_back(v3_2[j]);
        }

        *edgeSampleNum = (int)(sampleGLP->size()/3);
        *sampleG = new float[3*(*edgeSampleNum)];
        *sampleV1 = new float[3*(*edgeSampleNum)];
        *sampleV3 = new float[3*(*edgeSampleNum)];
        *X = new float[3*(*edgeSampleNum)];
        *Y = new float[3*(*edgeSampleNum)];
        *sampleProjG = new float[3*(*edgeSampleNum)];
        list<float>::iterator ig = sampleGLP->begin();
        list<float>::iterator iv1 = sampleV1LP->begin();
        list<float>::iterator iv3 = sampleV3LP->begin();
        for ( int i = 0 ; ig != sampleGLP->end() ; ig ++, iv1 ++, iv3 ++, i ++ )
        {
        (*sampleG)[i] = *ig;
        (*sampleV1)[i] = *iv1;
        (*sampleV3)[i] = *iv3;
        }

        return true;
        }

        void CGNSParallelEdge::extremalEdgeArray(float p1[3], float p2[3], Edge *edge, float *sampleV1, int edgeSampleNum, int *ltotalEdgePoints)
        {
        for ( int i = 1 ; i < edgeSampleNum ; i ++ )
        {
        float fSign = (float)sign(getDotP(sampleV1+3*i, sampleV1+3*(i-1)));
        if (fSign == 0) fSign = 1;
        sampleV1[3*i] = fSign*sampleV1[3*i];
        sampleV1[3*i+1] = fSign*sampleV1[3*i+1];
        sampleV1[3*i+2] = fSign*sampleV1[3*i+2];
        }
        float g1[3], g2[3];
        getGradientPN(p1, g1);
        getGradientPN(p2, g2);
        float d1 = getDotP(sampleV1, g1);
        float d2 = getDotP(sampleV1+3*(edgeSampleNum-1), g2);
        if (sign(d1*d2) == -1)
        {
        (*ltotalEdgePoints)++;
        float p[3];
        linearInterpolate(d1, d2, 0, p1, p2, p);
        edge->edgePoint[0] = p[0];
        edge->edgePoint[1] = p[1];
        edge->edgePoint[2] = p[2];
        //SelfAdjointEigenSolver<Matrix3f> eigensolver;
        //getEigensolverPN(p, &eigensolver);
        getEigensolverPN(p);

        float v1[3];
        float p10[3], p11[3], p13[3], p14[3];
        //getV1(&eigensolver, v1);
        getV1(v1);
        p10[0] = p[0]-2*change*v1[0];
        p10[1] = p[1]-2*change*v1[1];
        p10[2] = p[2]-2*change*v1[2];
        p11[0] = p[0]-change*v1[0];
        p11[1] = p[1]-change*v1[1];
        p11[2] = p[2]-change*v1[2];
        p13[0] = p[0]+change*v1[0];
        p13[1] = p[1]+change*v1[1];
        p13[2] = p[2]+change*v1[2];
        p14[0] = p[0]+2*change*v1[0];
        p14[1] = p[1]+2*change*v1[1];
        p14[2] = p[2]+2*change*v1[2];

        float f10, f11, f13, f14, f;
        getScalarPN(p10, &f10);
        getScalarPN(p11, &f11);
        getScalarPN(p13, &f13);
        getScalarPN(p14, &f14);
        getScalarPN(p, &f);
        edge->extremal = true;
        float secondDiff = -f10+16*f11-30*f+16*f13-f14;
        if (secondDiff < 0) edge->edgeTag = 1;
        else edge->edgeTag = 2;
        }
        }


        void CGNSParallelEdge::signedSphericalTriangleArea(float p1[3], float p2[3], float *area, float *phiC)
        {
        //*phiC = getDihedral(p1, globalVec, p2);
        //float phiP1 = getDihedral(p2, p1, globalVec);
        //float phiP2 = getDihedral(globalVec, p2, p1);
        float phiP1 = 0.0;
        float phiP2 = 0.0;
        *area = *phiC+phiP1+phiP2-(float)M_PI;
        float first[3], second[3], up[3];
        first[0] = p1[0]-globalVec[0];
        first[1] = p1[1]-globalVec[1];
        first[2] = p1[2]-globalVec[2];
        second[0] = p2[0]-p1[0];
        second[1] = p2[1]-p1[1];
        second[2] = p2[2]-p1[2];
        getCrossP(first, second, up);
        int aSign = sign(getDotP(up, globalVec));
        *area *= aSign;
        *phiC *= aSign;
        }

        void CGNSParallelEdge::orientV3Array(Edge *edge, float *sampleV3, int edgeSampleNum)
        {
        edge->f = 0;
        float oldEnd[3];
        oldEnd[0] = sampleV3[3*(edgeSampleNum-1)];
        oldEnd[1] = sampleV3[3*(edgeSampleNum-1)+1];
        oldEnd[2] = sampleV3[3*(edgeSampleNum-1)+2];
        for ( int i = 1 ; i < edgeSampleNum ; i ++ )
        {
        float fSign = (float)sign(getDotP(sampleV3+3*i, sampleV3+3*(i-1)));
        if (fSign == 0) fSign = 1;
        sampleV3[3*i+0] = fSign*sampleV3[3*i];
        sampleV3[3*i+1] = fSign*sampleV3[3*i+1];
        sampleV3[3*i+2] = fSign*sampleV3[3*i+2];
        }
        if (sign(getDotP(oldEnd, sampleV3+3*(edgeSampleNum-1))) == -1) edge->f = 1;
        }

        void CGNSParallelEdge::propagateXYArray(Edge *edge, float *sampleV3, float *X, float *Y, int edgeSampleNum)
        {
        float epsilon = 1e-12;
        float ox[3], oz[3];
        ox[0] = 1;
        ox[1] = 0;
        ox[2] = 0;
        oz[0] = 0;
        oz[1] = 0;
        oz[2] = 1;
        for ( int i = 0 ; i < edgeSampleNum ; i ++ )
        {
        float *v = sampleV3+3*i;
        float axle[3];
        getCrossP(oz, v, axle);
        float axleNorm = getNorm(axle);
        if (axleNorm > epsilon)
        {
        axle[0] = axle[0] / axleNorm;
        axle[1] = axle[1] / axleNorm;
        axle[2] = axle[2] / axleNorm;
        float angleCos = getDotP(oz, v);
        if (angleCos > 1) angleCos = 1;
        if (angleCos < -1) angleCos = -1;
        float angle = acos(angleCos);
        float a = getDotP(ox, axle)*(1-angleCos);
        float crossP[3];
        getCrossP(axle, ox, crossP);
        float b = sin(angle);
        ox[0] = ox[0]*angleCos + axle[0]*a + crossP[0]*b;
        ox[1] = ox[1]*angleCos + axle[1]*a + crossP[1]*b;
        ox[2] = ox[2]*angleCos + axle[2]*a + crossP[2]*b;
        oz[0] = v[0];
        oz[1] = v[1];
        oz[2] = v[2];
        }
        X[3*i] = ox[0];
        X[3*i+1] = ox[1];
        X[3*i+2] = ox[2];
        getCrossP(v, ox, Y+3*i);
        }
        float area = 0;
        float phiC = 0;
        for ( int i = 0 ; i < edgeSampleNum-1 ; i ++ )
        {
        float *p1 = sampleV3+3*i;
        float *p2 = sampleV3+3*(i+1);
        float incArea, incPhiC;
        signedSphericalTriangleArea(p1, p2, &incArea, &incPhiC);
        area += incArea;
        phiC += incPhiC;
        }
        edge->a = area;
        edge->b = phiC;
        }

        void CGNSParallelEdge::projectGArray(float *sampleG, float *X, float *Y, float *sampleProjG, int edgeSampleNum)
        {
        for ( int i = 0 ; i < edgeSampleNum ; i ++ )
        {
        float *g = sampleG+3*i;
        float *x = X+3*i;
        float *y = Y+3*i;
        sampleProjG[3*i] = getDotP(g, x);
        sampleProjG[3*i+1] = getDotP(g, y);
        sampleProjG[3*i+2] = 0;
        }
        }

        void CGNSParallelEdge::getWindingNumArray(Edge *edge, float *sampleProjG, int index1, int index2, int edgeSampleNum)
        {
        float totalAngle = 0;
        float epsilon = 1e-12;
        float *end1 = sampleProjG;
        for ( int i = 1 ; i < edgeSampleNum ; i ++ )
        {
        float *end2 = sampleProjG+3*i;
        float end1Norm = getNorm2(end1);
        float end2Norm = getNorm2(end2);
        if (end1Norm > 0 && end2Norm > 0)
        {
        float angleCos = getDotP(end1, end2) / (end1Norm*end2Norm);
        if (angleCos > 1) angleCos = 1;
        if (angleCos < -1) angleCos = -1;
        float localAngle = acos( angleCos );
        float crossP[3];
        getCrossP(end1, end2, crossP);
        localAngle *= sign(crossP[2]);
        totalAngle += localAngle;
        }
        end1 = end2;
        }
        edge->w = totalAngle;
        if (!storeGrid_2DGradMag[index1])
        {
        storeGrid_2DGradMag[index1] = true;
        grid_2DGradMag[index1] = getNorm(sampleProjG);
        }
        if (!storeGrid_2DGradMag[index2])
        {
        storeGrid_2DGradMag[index2] = true;
        grid_2DGradMag[index2] = getNorm(sampleProjG+3*(edgeSampleNum-1));
        }
        }

        void CGNSParallelEdge::edgePhaseIteration(int index1, int index2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints) {
        float *sampleG, *sampleV1, *sampleV3, *X, *Y, *sampleProjG;
        int edgeSampleNum;
        float p1[3], p2[3];
        getGridPointPosD(index1, p1);
        getGridPointPosD(index2, p2);

        if (!adaptiveSamplingArray(p1, p2, &edgeSampleNum, &sampleG, &sampleV1, &sampleV3, &X, &Y, &sampleProjG)) return;
        Edge *edge = &(ledges[getEdgeIndex(axis, i, j, k)]);
        edge->valid = true;

        extremalEdgeArray(p1, p2, edge, sampleV1, edgeSampleNum, ltotalEdgePoints);

        orientV3Array(edge, sampleV3, edgeSampleNum);
        propagateXYArray(edge, sampleV3, X, Y, edgeSampleNum);
        projectGArray(sampleG, X, Y, sampleProjG, edgeSampleNum);
        getWindingNumArray(edge, sampleProjG, index1, index2, edgeSampleNum);
        delete[] sampleG;
        delete[] sampleV1;
        delete[] sampleV3;
        delete[] X;
        delete[] Y;
        delete[] sampleProjG;

        }


        /**
        void operator()(const tbb::blocked_range3d<int>& r)	const {
        for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
        for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
        for( int k=r.pages().begin();k!=r.pages().end();++k){
        if (i+1==gridx)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 1, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
        for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
        for( int k=r.pages().begin();k!=r.pages().end();++k){
        if(j+1==gridy)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 2, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        for( int i=r.cols().begin(); i!=r.cols().end(); ++i ){
        for( int j=r.rows().begin(); j!=r.rows().end(); ++j ) {
        for( int k=r.pages().begin();k!=r.pages().end();++k){
        if(k+1==gridz)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 3, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        };
        };
        **/
        /**
        void CGNSParallelEdge::getGridPointPosD(int index, float v[3]) const
        {
        v[0] = gridPoints[3*index];
        v[1] = gridPoints[3*index+1];
        v[2] = gridPoints[3*index+2];
        }



        void CGNSParallelEdge::edgePhaseIteration(int index1, int index2, int axis, int i, int j, int k, Edge *ledges, int *ltotalEdgePoints) const
        {
        float *sampleG, *sampleV1, *sampleV3, *X, *Y, *sampleProjG;
        int edgeSampleNum;
        float p1[3], p2[3];
        getGridPointPosD(index1, p1);
        getGridPointPosD(index2, p2);
        if (!adaptiveSamplingArray(p1, p2, &edgeSampleNum, &sampleG, &sampleV1, &sampleV3, &X, &Y, &sampleProjG)) return;
        Edge *edge = &(ledges[getEdgeIndex(axis, i, j, k)]);
        edge->valid = true;
        extremalEdgeArray(p1, p2, edge, sampleV1, edgeSampleNum, ltotalEdgePoints);
        orientV3Array(edge, sampleV3, edgeSampleNum);
        propagateXYArray(edge, sampleV3, X, Y, edgeSampleNum);
        projectGArray(sampleG, X, Y, sampleProjG, edgeSampleNum);
        getWindingNumArray(edge, sampleProjG, index1, index2, edgeSampleNum);
        delete[] sampleG;
        delete[] sampleV1;
        delete[] sampleV3;
        delete[] X;
        delete[] Y;
        delete[] sampleProjG;

        }


        void CGNSParallelEdge::edgePhaseBegin(int cols, int rows, int pages) {
        for( int i=0; i!=cols; ++i ){
        for( int j=0; j!=rows; ++j ) {
        for( int k=0;k!=pages;++k){
        if (i+1==gridx)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 1, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        for( int i=0; i!=cols; ++i ){
        for( int j=0; j!=rows; ++j ) {
        for( int k=0;k!=pages;++k){
        if(j+1==gridy)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 2, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        for( int i=0; i!=cols; ++i ){
        for( int j=0; j!=rows; ++j ) {
        for( int k=0;k!=pages;++k){
        if(k+1==gridz)
        continue;
        int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
        int index2 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);
        edgePhaseIteration(index1, index2, 3, i, j, k, edges, totalEdgePoints);
        }
        }
        }
        }
        **/
        void General_Data::edgePhase(float* scalars, float* tensors, float* gradients, float *edgeTable)
        {
            /*t_sampling = t_sampling_interp = t_sampling_other = t_extrEdge = t_extrEdge_orient = t_extrEdge_interp = t_extrEdge_other
            = t_orientV3 = t_propagate = t_propagate_getXY = t_propagate_getEdgeAB = t_projectG = t_getEdgeW
            = t_facePhase_other = t_facePhase_extremalPoint = t_cellPhase_other = t_cellPhase_extremalPoint = 0;*/
            //t_sampling_interp = t_extrEdge_interp = t_faceSampling_interp = t_faceOther_interp = t_cellPhase_interp = 0;

            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Edge Phase Begin:" << endl;
            change = 0.017f * gridSize;
            maxEdgeIndex = (gridx - 1) * gridy * gridz + gridx * (gridy - 1) * gridz + gridx * gridy * (gridz - 1);
            edges = new Edge[maxEdgeIndex];
            grid_2DGradMag = new float[maxGridIndex];
            storeGrid_2DGradMag = new bool[maxGridIndex];
            for (int index = 0; index < maxGridIndex; index++)
            {
                storeGrid_2DGradMag[index] = false;
            }
            totalEdgePoints = new int(0);

            tbb::task_scheduler_init init(task_scheduler_init::automatic);
            ParallelEdge parallel_edge(gridx, gridy, gridz, totalEdgePoints, edges, allCubic, globalVec, dataType, gridPoints, change, storeGrid_2DGradMag, grid_2DGradMag, sizex, sizey, sizez, halfSize, scalars, tensors, gradients, edgeTable);
            parallel_for(blocked_range3d<int>(0, gridz, 0, gridy, 0, gridx), parallel_edge, auto_partitioner());

            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            /*cout << "t_sampling: " << t_sampling << endl;
            cout << "\tt_sampling_interp: " << t_sampling_interp << endl;
            cout << "\tt_sampling_other: " << t_sampling_other << endl;
            cout << "t_extrEdge: " << t_extrEdge << endl;
            cout << "\tt_extrEdge_orient: " << t_extrEdge_orient << endl;
            cout << "\tt_extrEdge_interp: " << t_extrEdge_interp << endl;
            cout << "\tt_extrEdge_other: " << t_extrEdge_other << endl;
            cout << "t_orientV3: " << t_orientV3 << endl;
            cout << "t_propagate: " << t_propagate << endl;
            cout << "\tt_propagate_getXY: " << t_propagate_getXY << endl;
            cout << "\tt_propagate_getEdgeAB: " << t_propagate_getEdgeAB << endl;
            cout << "t_projectG: " << t_projectG << endl;
            cout << "t_getEdgeW: " << t_getEdgeW << endl;*/
            /*cout << "\tt_sampling_interp: " << t_sampling_interp/CLOCKS_PER_SEC << endl;
            cout << "\tt_extrEdge_interp: " << t_extrEdge_interp/CLOCKS_PER_SEC << endl;*/
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }

        /**
        void General_Data::edgePhase(float* scalars,float* tensors,float* gradients, float *edgeTable)
        {
        clock_t timeStart, timeEnd;
        timeStart = clock();
        cout << "Edge Phase Begin:" << endl;
        change = 0.017f * gridSize;
        maxEdgeIndex = (gridx-1) * gridy * gridz + gridx * (gridy-1) * gridz + gridx * gridy * (gridz-1);
        edges = new Edge[ maxEdgeIndex ] ;
        grid_2DGradMag = new float[ maxGridIndex ] ;
        storeGrid_2DGradMag = new bool[ maxGridIndex ] ;
        for ( int index = 0 ; index < maxGridIndex ; index ++ )
        {
        storeGrid_2DGradMag[index] = false;
        }
        totalEdgePoints = new int(0);

        //tbb::task_scheduler_init init(task_scheduler_init::automatic);
        ParallelEdge parallel_edge(gridx,gridy,gridz,totalEdgePoints,edges,allCubic,globalVec,dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,sizex,sizey,sizez,halfSize,scalars,tensors,gradients,edgeTable);
        //parallel_edge.edgePhaseBegin(gridx, gridy, gridz);
        parallel_edge.edgePhaseBegin(gridz, gridy, gridx);
        //parallel_for( blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_edge,auto_partitioner());
        //CGNSParallelEdge parallel_edge(gridx,gridy,gridz,totalEdgePoints,edges,allCubic,globalVec,
        //dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,coef);
        //parallel_edge.edgePhaseBegin(gridz, gridy, gridx);

        timeEnd = clock();
        cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
        cout << "Done!" << endl;
        cout << "****************************************************************" << endl;
        }
        **/


        void General_Data::setVolumeData(VolumeData * vData) {
            volData = vData;
        }
        /**
        void General_Data::edgePhase(float *coef)
        {
        clock_t timeStart, timeEnd;
        timeStart = clock();
        cout << "Edge Phase Begin:" << endl;
        change = 0.017f * gridSize;
        maxEdgeIndex = (gridx-1) * gridy * gridz + gridx * (gridy-1) * gridz + gridx * gridy * (gridz-1);
        edges = new Edge[ maxEdgeIndex ] ;
        grid_2DGradMag = new float[ maxGridIndex ] ;
        storeGrid_2DGradMag = new bool[ maxGridIndex ] ;
        for ( int index = 0 ; index < maxGridIndex ; index ++ )
        {
        storeGrid_2DGradMag[index] = false;
        }
        totalEdgePoints = new int(0);

        //	tbb::task_scheduler_init init(task_scheduler_init::automatic);
        CGNSParallelEdge parallel_edge(gridx,gridy,gridz,totalEdgePoints,edges,allCubic,globalVec,
        dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,coef);
        //tbb::parallel_for( tbb::blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_edge,tbb::auto_partitioner());

        timeEnd = clock();
        cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
        cout << "Done!" << endl;
        cout << "****************************************************************" << endl;
        }

        void General_Data::facePhase(float *coef)
        {
        clock_t timeStart, timeEnd;
        timeStart = clock();
        cout << "Face Phase Begin:" << endl;
        change = 0.017f * gridSize;
        maxFaceIndex = gridx * (gridy-1) * (gridz-1) + (gridx-1) * gridy * (gridz-1) + (gridx-1) * (gridy-1) * gridz;
        faces = new Face[ maxFaceIndex ] ;
        totalFacePoints = new int(0);

        //tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
        CGNSParallelFace parallel_face(gridx,gridy,gridz,totalFacePoints,edges,faces,dataType,
        gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,coef);
        //tbb::parallel_for( tbb::blocked_range3d<int>(0, gridz, 0, gridy,0, gridx),parallel_face,tbb::auto_partitioner());

        delete[] grid_2DGradMag;
        delete[] storeGrid_2DGradMag;
        timeEnd = clock();
        cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;
        cout << "Done!" << endl;
        cout << "****************************************************************" << endl;
        }
        **/







        static void buildFaceTable(float *table, int subdNum)
        {
            float step = 1.0f / (subdNum - 1);
            int count = 0;
            for (int i = 0; i < subdNum; i++)
            for (int j = 0; j < subdNum; j++)
            {
                bicubicCoef(&(table[count * 16]), i*step, j*step);
                count++;
            }
        }


        General_Data::General_Data()
        {
            allCubic = true;
            height = 256;
            width = 256;
            slices = 50;
            thickness = 1.00;
            resolution = 40;
            dataType = 1;

            storeGridShowPos = storeGridScalar = storeGridVector = storeGridGradient =
                storeEdgePoint = storeFacePoint = storeCellPoint = storeDispCell = false;

            maxEigenvalue = new float;
            minEigenvalue = new float;
            maxIntensity = *maxEigenvalue = -1e12;
            minIntensity = *minEigenvalue = 1e12;

            globalVec[0] = exp(1.0f);
            globalVec[1] = (float)M_PI;
            globalVec[2] = exp(2.0f);
            float norm = getNorm(globalVec);
            globalVec[0] = globalVec[0] / norm;
            globalVec[1] = globalVec[1] / norm;
            globalVec[2] = globalVec[2] / norm;
            thresh = cos(20 * (float)M_PI / 180);

            gridPoints = NULL;
            gridScalars = NULL;
            gridShowPositions = NULL;
            gridV1s = NULL;
            gridV2s = NULL;
            gridV3s = NULL;
            gridGradients = NULL;
            edgePoints = NULL;
            facePoints = NULL;
            cellPoints = NULL;
            edges = NULL;
            faces = NULL;
            cells = NULL;
            grid_2DGradMag = NULL;
            storeGrid_2DGradMag = NULL;
            segments = NULL;
            quads = NULL;
            points = NULL;
            vertices = NULL;
        }

        General_Data::~General_Data(void)
        {
            delete[] gridPoints;
            delete[] gridScalars;
            delete[] gridShowPositions;
            delete[] gridV1s;
            delete[] gridV2s;
            delete[] gridV3s;
            delete[] gridGradients;
            delete[] edgePoints;
            delete[] facePoints;
            delete[] cellPoints;
            delete[] edges;
            delete[] faces;
            delete[] cells;
            delete[] grid_2DGradMag;
            delete[] storeGrid_2DGradMag;
            delete segments;
            delete quads;
            delete points;
            delete vertices;
        }



        class ParallelFace
        {
            int gridx, gridy, gridz;
            Edge *edges;
            Face *faces;
            int *totalFacePoints;
            int dataType;
            float change;
            static const int subdNum = 3;

            float *gridPoints;

            int sizex, sizey, sizez, halfSize;
            float*scalars, *tensors, *gradients;
            float *faceTable;
            float eigenValue1, eigenValue2, eigenValue3;
            std::vector<float> eigenVectors;

            int getEdgeIndex(int axis, int x, int y, int z) const;
            int getFaceIndex(int axis, int x, int y, int z) const;
            void combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int axis) const;
            void getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h, int *ltotalFacePoints) const;
            void getCurveType(Face *face) const;
            void get3DWindingNum(Face *face, int ii, int jj, int kk, int axis, float *faceTable) const;
            float smallAbsA(float A) const;
            void getGridPointPos(int index, float v[3]) const;
            void getGridPointPosD(int index, float v[3]) const;
            void rotate2D(float p[2], float theta) const;
            bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
            void getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
            float signedArea3D(float v1[3], float v2[3], float v3[3]) const;
            float getDihedral(float b1[3], float b2[3], float b3[3]) const;
            void facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, int index1, int index2, int index3, int index4, int axis, int i, int j, int k, Face *lfaces, int *ltotalFacePoints, float *faceTable) const;

            void getScalarGP(int index, float *scalar) const
            {
                *scalar = scalars[index];
            }

            void getGradientGP(int index, float gradient[3]) const
            {
                gradient[0] = gradients[index];
                gradient[1] = gradients[index + 1];
                gradient[2] = gradients[index + 2];
            }

            void getTensorGP(int index, float tensor[6]) const
            {
                tensor[0] = tensors[index];
                tensor[1] = tensors[index + 1];
                tensor[2] = tensors[index + 2];
                tensor[3] = tensors[index + 3];
                tensor[4] = tensors[index + 4];
                tensor[5] = tensors[index + 5];
            }

            void getScalar(float position[3], float *scalar) const
            {
                trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradient(float position[3], float gradient[3]) const
            {
                trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensor(float position[3], float tensor[6]) const
            {
                trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
            };

            void getScalarCubic(float position[3], float *scalar) const
            {
                tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradientCubic(float position[3], float gradient[3]) const
            {
                tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensorCubic(float position[3], float tensor[6]) const
            {
                tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
            };

            void getGradientBicubicTable(int axis, int ii, int jj, int kk, float sGradients[subdNum][subdNum][3], float *faceTable) const
            {
                float validData[16][3];
                int count = 0;
                switch (axis)
                {
                case 1:
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            int idx = getIndex(3, 0, ii + 2, jj + i + 1, kk + j + 1, sizex, sizey, sizez);
                            validData[count][0] = gradients[idx];
                            validData[count][1] = gradients[idx + 1];
                            validData[count][2] = gradients[idx + 2];
                            count++;
                        }
                    }
                    break;
                case 2:
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            int idx = getIndex(3, 0, ii + j + 1, jj + 2, kk + i + 1, sizex, sizey, sizez);
                            validData[count][0] = gradients[idx];
                            validData[count][1] = gradients[idx + 1];
                            validData[count][2] = gradients[idx + 2];
                            count++;
                        }
                    }
                    break;
                case 3:
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            int idx = getIndex(3, 0, ii + i + 1, jj + j + 1, kk + 2, sizex, sizey, sizez);
                            validData[count][0] = gradients[idx];
                            validData[count][1] = gradients[idx + 1];
                            validData[count][2] = gradients[idx + 2];
                            count++;
                        }
                    }
                    break;
                }
                count = 0;
                for (int si = 0; si < subdNum; si++)
                {
                    for (int sj = 0; sj < subdNum; sj++)
                    {
                        int tableIndex = 16 * count;
                        count++;
                        sGradients[si][sj][0] = 0;
                        sGradients[si][sj][1] = 0;
                        sGradients[si][sj][2] = 0;
                        for (int i = 0; i < 16; i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                sGradients[si][sj][j] += faceTable[tableIndex + i] * validData[i][j];
                            }
                        }
                    }
                }
            }

        public:
            bool *storeGrid_2DGradMag;
            float *grid_2DGradMag;

            ParallelFace(int lgridx, int lgridy, int lgridz, int *ltotalFacePoints, Edge *ledges, Face *lfaces, int ldataType, float *lgridPoints, float lchange, bool *lstoreGrid_2DGradMag, float *lgrid_2DGradMag, int lsizex, int lsizey, int lsizez, int lhalfSize, float* lscalars, float* ltensors, float* lgradients, float *lfaceTable)
            {
                gridx = lgridx;
                gridy = lgridy;
                gridz = lgridz;
                totalFacePoints = ltotalFacePoints;
                edges = ledges;
                faces = lfaces;
                dataType = ldataType;
                gridPoints = lgridPoints;
                change = lchange;
                storeGrid_2DGradMag = lstoreGrid_2DGradMag;
                grid_2DGradMag = lgrid_2DGradMag;
                sizex = lsizex;
                sizey = lsizey;
                sizez = lsizez;
                halfSize = lhalfSize;
                scalars = lscalars;
                tensors = ltensors;
                gradients = lgradients;
                faceTable = lfaceTable;

            };
            void operator()(const blocked_range3d<int>& r)	const {
                int *ltotalFacePoints = totalFacePoints;
                //Edge *ledges = edges;
                Face *lfaces = faces;

                for (int i = r.cols().begin(); i != r.cols().end(); ++i)
                for (int j = r.rows().begin(); j != r.rows().end(); ++j)
                for (int k = r.pages().begin(); k != r.pages().end(); ++k)
                {
                    if (j == gridy - 1)
                        continue;
                    if (k == gridz - 1)
                        continue;
                    int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                    int index2 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);
                    int index3 = getIndex(1, 0, i, j + 1, k + 1, gridx, gridy, gridz);
                    int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);

                    Edge *edge1 = &(edges[getEdgeIndex(2, i, j, k)]);
                    if (!edge1->valid) continue;
                    Edge *edge2 = &(edges[getEdgeIndex(3, i, j + 1, k)]);
                    if (!edge2->valid) continue;
                    Edge *edge3 = &(edges[getEdgeIndex(2, i, j, k + 1)]);
                    if (!edge3->valid) continue;
                    Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
                    if (!edge4->valid) continue;


                    facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 1, i, j, k, lfaces, ltotalFacePoints, faceTable);
                }
                for (int i = r.cols().begin(); i != r.cols().end(); ++i)
                for (int j = r.rows().begin(); j != r.rows().end(); ++j)
                for (int k = r.pages().begin(); k != r.pages().end(); ++k)
                {
                    if (i == gridx - 1)
                        continue;
                    if (k == gridz - 1)
                        continue;
                    int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                    int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
                    int index3 = getIndex(1, 0, i + 1, j, k + 1, gridx, gridy, gridz);
                    int index4 = getIndex(1, 0, i, j, k + 1, gridx, gridy, gridz);


                    Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
                    if (!edge1->valid) continue;
                    Edge *edge2 = &(edges[getEdgeIndex(3, i + 1, j, k)]);
                    if (!edge2->valid) continue;
                    Edge *edge3 = &(edges[getEdgeIndex(1, i, j, k + 1)]);
                    if (!edge3->valid) continue;
                    Edge *edge4 = &(edges[getEdgeIndex(3, i, j, k)]);
                    if (!edge4->valid) continue;


                    facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 2, i, j, k, lfaces, ltotalFacePoints, faceTable);
                }
                for (int i = r.cols().begin(); i != r.cols().end(); ++i)
                for (int j = r.rows().begin(); j != r.rows().end(); ++j)
                for (int k = r.pages().begin(); k != r.pages().end(); ++k)
                {
                    if (i == gridx - 1)
                        continue;
                    if (j == gridy - 1)
                        continue;
                    int index1 = getIndex(1, 0, i, j, k, gridx, gridy, gridz);
                    int index2 = getIndex(1, 0, i + 1, j, k, gridx, gridy, gridz);
                    int index3 = getIndex(1, 0, i + 1, j + 1, k, gridx, gridy, gridz);
                    int index4 = getIndex(1, 0, i, j + 1, k, gridx, gridy, gridz);

                    Edge *edge1 = &(edges[getEdgeIndex(1, i, j, k)]);
                    if (!edge1->valid) continue;
                    Edge *edge2 = &(edges[getEdgeIndex(2, i + 1, j, k)]);
                    if (!edge2->valid) continue;
                    Edge *edge3 = &(edges[getEdgeIndex(1, i, j + 1, k)]);
                    if (!edge3->valid) continue;
                    Edge *edge4 = &(edges[getEdgeIndex(2, i, j, k)]);
                    if (!edge4->valid) continue;


                    facePhaseIteration(edge1, edge2, edge3, edge4, index1, index2, index3, index4, 3, i, j, k, lfaces, ltotalFacePoints, faceTable);
                }
            };
        };

        void ParallelFace::facePhaseIteration(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, int index1, int index2, int index3, int index4, int axis, int i, int j, int k, Face *lfaces, int *ltotalFacePoints, float *faceTable) const
        {
            Face *face = &(faces[getFaceIndex(axis, i, j, k)]);
            face->valid = true;
            float h[4];
            combineWindingNum(edge1, edge2, edge3, edge4, face, h, axis);
            getFaceIntersection(index1, index2, index3, index4, face, h, ltotalFacePoints);
            getCurveType(face);
            get3DWindingNum(face, i, j, k, axis, faceTable);
        }

        int ParallelFace::getFaceIndex(int axis, int x, int y, int z) const {
            int index;
            switch (axis)
            {
            case 1:
                index = x * (gridy - 1) * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 2:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 3:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    (gridx - 1) * gridy * (gridz - 1) +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            }
            return 0;
        }

        int ParallelFace::getEdgeIndex(int axis, int x, int y, int z) const {
            int index;
            switch (axis)
            {
            case 1:
                index = x * gridy * gridz + y * gridz + z;
                return index;
                break;
            case 2:
                index = (gridx - 1) * gridy * gridz +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            case 3:
                index = (gridx - 1) * gridy * gridz +
                    gridx * (gridy - 1) * gridz +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            }
            return 0;
        }

        void ParallelFace::combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int axis) const
        {
            int f[4];
            float a[4], b[4], w[4];
            f[0] = edge1->f;
            f[1] = edge2->f;
            f[2] = edge3->f;
            f[3] = edge4->f;
            a[0] = edge1->a;
            a[1] = edge2->a;
            a[2] = -edge3->a;
            a[3] = -edge4->a;
            b[0] = edge1->b;
            b[1] = edge2->b;
            b[2] = -edge3->b;
            b[3] = -edge4->b;
            w[0] = edge1->w;
            w[1] = edge2->w;
            w[2] = -edge3->w;
            w[3] = -edge4->w;
            int d[4] = { 0, 0, 1, 1 };
            int F = f[0] + f[1] + f[2] + f[3];
            if (F % 2 != 0) return;
            F = 0;
            float A = 0;
            float W = 0;
            h[0] = 0;
            for (int i = 0; i < 4; i++)
            {
                if ((F + d[i] * f[i]) % 2 != 0)
                {
                    a[i] = 2 * b[i] - a[i];
                    w[i] = -w[i];
                }
                F += f[i];
                A += a[i];
                W += w[i];
                if (i < 3)
                {
                    h[i + 1] = W;
                }
            }
            A = smallAbsA(A);
            float wn = (W + A) / (2 * (float)M_PI);

            wn = floor(wn + 0.5f);
            if ((int)wn % 2 != 0)
            {
                face->extremal = true;
                for (int i = 1; i < 4; i++)
                {
                    h[i] += (float)i / 4 * A;
                }
            }
        }

        void ParallelFace::getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h, int *ltotalFacePoints) const
        {
            if (!face->extremal) return;
            (*ltotalFacePoints++);
            float vertices[4][3];
            getGridPointPos(index1, vertices[0]);
            getGridPointPos(index2, vertices[1]);
            getGridPointPos(index3, vertices[2]);
            getGridPointPos(index4, vertices[3]);
            float projVertices[4][2];
            projVertices[0][0] = grid_2DGradMag[index1];
            projVertices[1][0] = grid_2DGradMag[index2];
            projVertices[2][0] = grid_2DGradMag[index3];
            projVertices[3][0] = grid_2DGradMag[index4];
            projVertices[0][1] = 0;
            projVertices[1][1] = 0;
            projVertices[2][1] = 0;
            projVertices[3][1] = 0;
            for (int i = 0; i < 4; i++)
            {
                rotate2D(projVertices[i], h[i]);
            }
            int facePointNum = 4;
            float spoke[4], rim[4], alpha[4], omega[4], lambda[4];
            for (int i = 0; i < facePointNum; i++)
            {
                spoke[i] = getNorm2(projVertices[i]);
            }
            float diff[2];
            for (int i = 0; i < facePointNum - 1; i++)
            {
                diff[0] = projVertices[i + 1][0] - projVertices[i][0];
                diff[1] = projVertices[i + 1][1] - projVertices[i][1];
                rim[i] = getNorm2(diff);
            }
            diff[0] = projVertices[0][0] - projVertices[facePointNum - 1][0];
            diff[1] = projVertices[0][1] - projVertices[facePointNum - 1][1];
            rim[facePointNum - 1] = getNorm2(diff);
            for (int i = 0; i < facePointNum - 1; i++)
            {
                if (2 * spoke[i] * spoke[i + 1] == 0)
                {
                    alpha[i] = 0;
                }
                else
                {
                    float angleCos = (pow(spoke[i], 2) + pow(spoke[i + 1], 2) - pow(rim[i], 2)) / (2 * spoke[i] * spoke[i + 1]);
                    if (angleCos > 1) angleCos = 1;
                    if (angleCos < -1) angleCos = -1;
                    alpha[i] = acos(angleCos);
                }
            }
            if (2 * spoke[facePointNum - 1] * spoke[0] == 0)
            {
                alpha[facePointNum - 1] = 0;
            }
            else
            {
                float angleCos = (pow(spoke[facePointNum - 1], 2) + pow(spoke[0], 2) - pow(rim[facePointNum - 1], 2))
                    / (2 * spoke[facePointNum - 1] * spoke[0]);
                if (angleCos > 1) angleCos = 1;
                if (angleCos < -1) angleCos = -1;
                alpha[facePointNum - 1] = acos(angleCos);
            }
            for (int i = 1; i < facePointNum; i++)
            {
                if (spoke[i] == 0) omega[i] = 0;
                else omega[i] = (tan(alpha[i - 1] / 2) + tan(alpha[i] / 2)) / spoke[i];
            }
            if (spoke[0] == 0) omega[0] = 0;
            else omega[0] = (tan(alpha[facePointNum - 1] / 2) + tan(alpha[0] / 2)) / spoke[0];
            float omegaSum = 0;
            for (int i = 0; i < facePointNum; i++)
            {
                omegaSum += omega[i];
            }
            if (omegaSum == 0)
            {
                for (int i = 0; i < facePointNum; i++)
                {
                    lambda[i] = 0.25;
                }
            }
            else
            {
                for (int i = 0; i < facePointNum; i++)
                {
                    lambda[i] = omega[i] / omegaSum;
                }
            }
            face->facePoint[0] = 0;
            face->facePoint[1] = 0;
            face->facePoint[2] = 0;
            for (int i = 0; i < facePointNum; i++)
            {
                face->facePoint[0] += lambda[i] * vertices[i][0];
                face->facePoint[1] += lambda[i] * vertices[i][1];
                face->facePoint[2] += lambda[i] * vertices[i][2];
            }
        }

        void ParallelFace::getCurveType(Face *face) const
        {
            if (!face->extremal) return;

            float *p = face->facePoint;
            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getEigensolverCubic(p, &eigensolver);

            float v1[3];
            getV1(&eigensolver, v1);
            float p10[3], p11[3], p13[3], p14[3];
            p10[0] = p[0] - 2 * change*v1[0];
            p10[1] = p[1] - 2 * change*v1[1];
            p10[2] = p[2] - 2 * change*v1[2];
            p11[0] = p[0] - change*v1[0];
            p11[1] = p[1] - change*v1[1];
            p11[2] = p[2] - change*v1[2];
            p13[0] = p[0] + change*v1[0];
            p13[1] = p[1] + change*v1[1];
            p13[2] = p[2] + change*v1[2];
            p14[0] = p[0] + 2 * change*v1[0];
            p14[1] = p[1] + 2 * change*v1[1];
            p14[2] = p[2] + 2 * change*v1[2];

            float v2[3];
            getV2(&eigensolver, v2);
            float p20[3], p21[3], p23[3], p24[3];
            p20[0] = p[0] - 2 * change*v2[0];
            p20[1] = p[1] - 2 * change*v2[1];
            p20[2] = p[2] - 2 * change*v2[2];
            p21[0] = p[0] - change*v2[0];
            p21[1] = p[1] - change*v2[1];
            p21[2] = p[2] - change*v2[2];
            p23[0] = p[0] + change*v2[0];
            p23[1] = p[1] + change*v2[1];
            p23[2] = p[2] + change*v2[2];
            p24[0] = p[0] + 2 * change*v2[0];
            p24[1] = p[1] + 2 * change*v2[1];
            p24[2] = p[2] + 2 * change*v2[2];

            float f, f10, f11, f13, f14, f20, f21, f23, f24;
            getScalarCubic(p, &f);
            getScalarCubic(p10, &f10);
            getScalarCubic(p11, &f11);
            getScalarCubic(p13, &f13);
            getScalarCubic(p14, &f14);
            getScalarCubic(p20, &f20);
            getScalarCubic(p21, &f21);
            getScalarCubic(p23, &f23);
            getScalarCubic(p24, &f24);
            float secondDiff1 = -f10 + 16 * f11 - 30 * f + 16 * f13 - f14;
            float secondDiff2 = -f20 + 16 * f21 - 30 * f + 16 * f23 - f24;
            if (secondDiff1 * secondDiff2 < 0) face->faceTag = 3;
            else if (secondDiff1 > 0 && secondDiff2 > 0) face->faceTag = 2;
            else face->faceTag = 1;
        }

        void ParallelFace::get3DWindingNum(Face *face, int ii, int jj, int kk, int axis, float *faceTable) const
        {
            float sGradients[subdNum][subdNum][3];
            getGradientBicubicTable(axis, ii, jj, kk, sGradients, faceTable);
            float normGrads[subdNum][subdNum][3];
            for (int i = 0; i < subdNum; i++)
            {
                for (int j = 0; j < subdNum; j++)
                {
                    //getShowPos(subdVerts[i][j], face->subdVerts[i][j]);
                    //getGradientCubic(subdVerts[i][j], gradients[i][j]);
                    float norm = getNorm(sGradients[i][j]);
                    if (norm != 0)
                    {
                        normGrads[i][j][0] = sGradients[i][j][0] / norm;
                        normGrads[i][j][1] = sGradients[i][j][1] / norm;
                        normGrads[i][j][2] = sGradients[i][j][2] / norm;
                    }
                }
            }
            float area = 0, currentArea;
            int count = 0;
            for (int i = 0; i < subdNum - 1; i++)
            {
                for (int j = 0; j < subdNum - 1; j++)
                {
                    currentArea = signedArea3D(normGrads[i][j], normGrads[i][j + 1], normGrads[i + 1][j]);
                    /*face->areas[count] = currentArea;
                    face->centers[count][0] = (face->subdVerts[i][j][0]+face->subdVerts[i][j+1][0]+face->subdVerts[i+1][j][0])/3;
                    face->centers[count][1] = (face->subdVerts[i][j][1]+face->subdVerts[i][j+1][1]+face->subdVerts[i+1][j][1])/3;
                    face->centers[count][2] = (face->subdVerts[i][j][2]+face->subdVerts[i][j+1][2]+face->subdVerts[i+1][j][2])/3;
                    count ++;*/
                    area += currentArea;
                    currentArea = signedArea3D(normGrads[i][j + 1], normGrads[i + 1][j + 1], normGrads[i + 1][j]);
                    /*face->areas[count] = currentArea;
                    face->centers[count][0] = (face->subdVerts[i][j+1][0]+face->subdVerts[i+1][j+1][0]+face->subdVerts[i+1][j][0])/3;
                    face->centers[count][1] = (face->subdVerts[i][j+1][1]+face->subdVerts[i+1][j+1][1]+face->subdVerts[i+1][j][1])/3;
                    face->centers[count][2] = (face->subdVerts[i][j+1][2]+face->subdVerts[i+1][j+1][2]+face->subdVerts[i+1][j][2])/3;
                    count ++;*/
                    area += currentArea;
                }
            }
            face->spherArea = area;
        }

        float ParallelFace::smallAbsA(float A) const
        {
            float A1 = fmod(A, 4 * (float)M_PI);
            float A2;
            if (A1 > 0) A2 = A1 - 4 * (float)M_PI;
            else A2 = A1 + 4 * (float)M_PI;
            if (fabs(A1) > fabs(A2)) return A2;
            else return A1;
        }

        void ParallelFace::getGridPointPosD(int index, float v[3]) const
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        void ParallelFace::rotate2D(float p[2], float theta) const
        {
            float ox = p[0];
            float oy = p[1];
            float x = ox*cos(theta) - oy*sin(theta);
            float y = ox*sin(theta) + oy*cos(theta);
            p[0] = x;
            p[1] = y;
        }

        bool ParallelFace::getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorCubic(position, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }
        bool ParallelFace::getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorGP(index, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }

        void ParallelFace::getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev;
            switch (dataType)
            {
            case 1:
                ev = eigensolver->eigenvectors().col(2).normalized();
                break;
            case 2:
                ev = eigensolver->eigenvectors().col(0).normalized();
                break;
            }
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }

        void ParallelFace::getV2(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev = eigensolver->eigenvectors().col(1).normalized();
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }

        float ParallelFace::signedArea3D(float v1[3], float v2[3], float v3[3]) const
        {
            float phiC = getDihedral(v1, v3, v2);
            float phiP1 = getDihedral(v2, v1, v3);
            float phiP2 = getDihedral(v3, v2, v1);
            float area = phiC + phiP1 + phiP2 - (float)M_PI;
            float first[3], second[3], up[3];
            first[0] = v1[0] - v3[0];
            first[1] = v1[1] - v3[1];
            first[2] = v1[2] - v3[2];
            second[0] = v2[0] - v1[0];
            second[1] = v2[1] - v1[1];
            second[2] = v2[2] - v1[2];
            getCrossP(first, second, up);
            int aSign = sign(getDotP(up, v3));
            area *= aSign;
            return area;
        }

        void ParallelFace::getGridPointPos(int index, float v[3]) const
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        float ParallelFace::getDihedral(float b1[3], float b2[3], float b3[3]) const
        {
            float b21[3], b23[3];
            getCrossP(b2, b1, b21);
            getCrossP(b2, b3, b23);
            float norm = getNorm(b21);
            if (norm != 0)
            {
                b21[0] = b21[0] / norm;
                b21[1] = b21[1] / norm;
                b21[2] = b21[2] / norm;
            }
            norm = getNorm(b23);
            if (norm != 0)
            {
                b23[0] = b23[0] / norm;
                b23[1] = b23[1] / norm;
                b23[2] = b23[2] / norm;
            }
            float dotP = getDotP(b21, b23);
            if (dotP > 1) dotP = 1;
            if (dotP < -1) dotP = -1;
            float phi = acos(dotP);
            return phi;
        }

        int General_Data::getFaceIndex(int axis, int x, int y, int z) {
            int index;
            switch (axis)
            {
            case 1:
                index = x * (gridy - 1) * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 2:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 3:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    (gridx - 1) * gridy * (gridz - 1) +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            }
            return 0;
        }

        void General_Data::getGridPointPos(int index, float v[3])
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        void General_Data::combineWindingNum(Edge *edge1, Edge *edge2, Edge *edge3, Edge *edge4, Face *face, float *h, int index1, int index2, int index3, int index4, int axis)
        {
            int f[4];
            float a[4], b[4], w[4];
            f[0] = edge1->f;
            f[1] = edge2->f;
            f[2] = edge3->f;
            f[3] = edge4->f;
            a[0] = edge1->a;
            a[1] = edge2->a;
            a[2] = -edge3->a;
            a[3] = -edge4->a;
            b[0] = edge1->b;
            b[1] = edge2->b;
            b[2] = -edge3->b;
            b[3] = -edge4->b;
            w[0] = edge1->w;
            w[1] = edge2->w;
            w[2] = -edge3->w;
            w[3] = -edge4->w;
            int d[4] = { 0, 0, 1, 1 };
            int F = f[0] + f[1] + f[2] + f[3];
            if (F % 2 != 0) return;
            F = 0;
            float A = 0;
            float W = 0;
            h[0] = 0;
            for (int i = 0; i < 4; i++)
            {
                if ((F + d[i] * f[i]) % 2 != 0)
                {
                    a[i] = 2 * b[i] - a[i];
                    w[i] = -w[i];
                }
                F += f[i];
                A += a[i];
                W += w[i];
                if (i < 3)
                {
                    h[i + 1] = W;
                }
            }
            A = smallAbsA(A);
            float wn = (W + A) / (2 * (float)M_PI);

            /*float vertices[4][3];
            float gradients[4][3];
            clock_t timeStart_i, timeEnd_i;
            for ( int i = 0; i < 3; i ++ )
            {
            getGridPointPosD(index1, vertices[0]);
            getGridPointPosD(index2, vertices[1]);
            getGridPointPosD(index3, vertices[2]);
            getGridPointPosD(index4, vertices[3]);
            timeStart_i = clock();
            getGradientCubic(vertices[0], gradients[0]);
            getGradientCubic(vertices[1], gradients[1]);
            getGradientCubic(vertices[2], gradients[2]);
            getGradientCubic(vertices[3], gradients[3]);
            timeEnd_i = clock();
            t_faceOther_interp = t_faceOther_interp + timeEnd_i - timeStart_i;
            }*/

            wn = floor(wn + 0.5f);
            if ((int)wn % 2 != 0)
            {
                face->extremal = true;
                for (int i = 1; i < 4; i++)
                {
                    h[i] += (float)i / 4 * A;
                }
            }
        }

        float General_Data::smallAbsA(float A)
        {
            float A1 = fmod(A, 4 * (float)M_PI);
            float A2;
            if (A1 > 0) A2 = A1 - 4 * (float)M_PI;
            else A2 = A1 + 4 * (float)M_PI;
            if (fabs(A1) > fabs(A2)) return A2;
            else return A1;
        }

        void General_Data::rotate2D(float p[2], float theta)
        {
            float ox = p[0];
            float oy = p[1];
            float x = ox*cos(theta) - oy*sin(theta);
            float y = ox*sin(theta) + oy*cos(theta);
            p[0] = x;
            p[1] = y;
        }


        void General_Data::getFaceIntersection(int index1, int index2, int index3, int index4, Face *face, float *h)
        {
            if (!face->extremal) return;
            (*totalFacePoints++);
            float vertices[4][3];
            getGridPointPos(index1, vertices[0]);
            getGridPointPos(index2, vertices[1]);
            getGridPointPos(index3, vertices[2]);
            getGridPointPos(index4, vertices[3]);
            float projVertices[4][2];
            projVertices[0][0] = grid_2DGradMag[index1];
            projVertices[1][0] = grid_2DGradMag[index2];
            projVertices[2][0] = grid_2DGradMag[index3];
            projVertices[3][0] = grid_2DGradMag[index4];
            projVertices[0][1] = 0;
            projVertices[1][1] = 0;
            projVertices[2][1] = 0;
            projVertices[3][1] = 0;
            for (int i = 0; i < 4; i++)
            {
                rotate2D(projVertices[i], h[i]);
            }
            int facePointNum = 4;
            float spoke[4], rim[4], alpha[4], omega[4], lambda[4];
            for (int i = 0; i < facePointNum; i++)
            {
                spoke[i] = getNorm2(projVertices[i]);
            }
            float diff[2];
            for (int i = 0; i < facePointNum - 1; i++)
            {
                diff[0] = projVertices[i + 1][0] - projVertices[i][0];
                diff[1] = projVertices[i + 1][1] - projVertices[i][1];
                rim[i] = getNorm2(diff);
            }
            diff[0] = projVertices[0][0] - projVertices[facePointNum - 1][0];
            diff[1] = projVertices[0][1] - projVertices[facePointNum - 1][1];
            rim[facePointNum - 1] = getNorm2(diff);
            for (int i = 0; i < facePointNum - 1; i++)
            {
                if (2 * spoke[i] * spoke[i + 1] == 0)
                {
                    alpha[i] = 0;
                }
                else
                {
                    float angleCos = (pow(spoke[i], 2) + pow(spoke[i + 1], 2) - pow(rim[i], 2)) / (2 * spoke[i] * spoke[i + 1]);
                    if (angleCos > 1) angleCos = 1;
                    if (angleCos < -1) angleCos = -1;
                    alpha[i] = acos(angleCos);
                }
            }
            if (2 * spoke[facePointNum - 1] * spoke[0] == 0)
            {
                alpha[facePointNum - 1] = 0;
            }
            else
            {
                float angleCos = (pow(spoke[facePointNum - 1], 2) + pow(spoke[0], 2) - pow(rim[facePointNum - 1], 2))
                    / (2 * spoke[facePointNum - 1] * spoke[0]);
                if (angleCos > 1) angleCos = 1;
                if (angleCos < -1) angleCos = -1;
                alpha[facePointNum - 1] = acos(angleCos);
            }
            for (int i = 1; i < facePointNum; i++)
            {
                if (spoke[i] == 0) omega[i] = 0;
                else omega[i] = (tan(alpha[i - 1] / 2) + tan(alpha[i] / 2)) / spoke[i];
            }
            if (spoke[0] == 0) omega[0] = 0;
            else omega[0] = (tan(alpha[facePointNum - 1] / 2) + tan(alpha[0] / 2)) / spoke[0];
            float omegaSum = 0;
            for (int i = 0; i < facePointNum; i++)
            {
                omegaSum += omega[i];
            }
            if (omegaSum == 0)
            {
                for (int i = 0; i < facePointNum; i++)
                {
                    lambda[i] = 0.25;
                }
            }
            else
            {
                for (int i = 0; i < facePointNum; i++)
                {
                    lambda[i] = omega[i] / omegaSum;
                }
            }
            face->facePoint[0] = 0;
            face->facePoint[1] = 0;
            face->facePoint[2] = 0;
            for (int i = 0; i < facePointNum; i++)
            {
                face->facePoint[0] += lambda[i] * vertices[i][0];
                face->facePoint[1] += lambda[i] * vertices[i][1];
                face->facePoint[2] += lambda[i] * vertices[i][2];
            }
        }



        Face *General_Data::getFaces() {
            return faces;
        }


        /**
        void General_Data::facePhase(float* scalars,float* tensors,float* gradients, float *faceTable)
        {
        clock_t timeStart, timeEnd;
        timeStart = clock();
        cout << "Face Phase Begin:" << endl;
        change = 0.017f * gridSize;
        maxFaceIndex = gridx * (gridy-1) * (gridz-1) + (gridx-1) * gridy * (gridz-1) + (gridx-1) * (gridy-1) * gridz;
        faces = new Face[ maxFaceIndex ] ;
        totalFacePoints = new int(0);



        ParallelFace parallel_face(gridx,gridy,gridz,totalFacePoints,edges,faces,dataType,gridPoints,change,storeGrid_2DGradMag,grid_2DGradMag,sizex,sizey,sizez,halfSize,scalars,tensors,gradients,faceTable);
        parallel_face.facePhaseBegin(gridz,gridy,gridx);

        delete[] grid_2DGradMag;
        delete[] storeGrid_2DGradMag;
        timeEnd = clock();
        //t_facePhase_other = timeEnd-timeStart - t_facePhase_extremalPoint;
        cout << "Spend Time: " << (timeEnd-timeStart)/CLOCKS_PER_SEC << endl;

        cout << "Done!" << endl;
        cout << "****************************************************************" << endl;
        }
        **/
        void General_Data::facePhase(float* scalars, float* tensors, float* gradients, float *faceTable)
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Face Phase Begin:" << endl;
            change = 0.017f * gridSize;
            maxFaceIndex = gridx * (gridy - 1) * (gridz - 1) + (gridx - 1) * gridy * (gridz - 1) + (gridx - 1) * (gridy - 1) * gridz;
            faces = new Face[maxFaceIndex];
            totalFacePoints = new int(0);

            tbb::task_scheduler_init init(task_scheduler_init::automatic);
            ParallelFace parallel_face(gridx, gridy, gridz, totalFacePoints, edges, faces, dataType, gridPoints, change, storeGrid_2DGradMag, grid_2DGradMag, sizex, sizey, sizez, halfSize, scalars, tensors, gradients, faceTable);
            parallel_for(blocked_range3d<int>(0, gridz, 0, gridy, 0, gridx), parallel_face, auto_partitioner());

            delete[] grid_2DGradMag;
            delete[] storeGrid_2DGradMag;
            timeEnd = clock();
            //t_facePhase_other = timeEnd-timeStart - t_facePhase_extremalPoint;
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            /*cout << "\tt_facePhase_other: " << t_facePhase_other << endl;
            cout << "\tt_facePhase_extremalPoint: " << t_facePhase_extremalPoint << endl;*/
            /*cout << "\tt_faceSampling_interp: " << t_faceSampling_interp/CLOCKS_PER_SEC << endl;
            cout << "\tt_faceOther_interp: " << t_faceOther_interp/CLOCKS_PER_SEC << endl;*/
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }

        class MVC
        {
            // Mesh
            float ** vs;
            int ** ts;
            int nv, nt;
            float ** vert;
            float * mags;

        public:
            // Constructor
            MVC(int numv, int numt, float ** verts, int ** tris)
            {
                vs = verts;
                ts = tris;
                nv = numv;
                nt = numt;

                vert = new float*[numv];
                mags = new float[numv];
                for (int i = 0; i < numv; i++)
                {
                    vert[i] = new float[3];
                }
            };

            // Compute weights
            int getWeights(float x[3], float* weights)
            {
                int i, j;
                float* v[3];
                float a[3], w[3], totalW = 0, totalWF = 0, totalA, vol;
                float l[3], t[3], s, sinS, sinT[3], cs[3], sinG[3];
                float EPSILN = 0.001;
                float mag;

                for (i = 0; i < nv; i++)
                {
                    weights[i] = 0;
                }

                for (i = 0; i < nv; i++)
                {
                    vert[i][0] = vs[i][0] - x[0];
                    vert[i][1] = vs[i][1] - x[1];
                    vert[i][2] = vs[i][2] - x[2];

                    mag = vert[i][0] * vert[i][0] + vert[i][1] * vert[i][1] + vert[i][2] * vert[i][2];
                    if (mag < EPSILN * EPSILN)
                    {
                        weights[i] = 1;
                        return 0;
                    }
                    mag = (float)sqrt(mag);
                    vert[i][0] /= mag;
                    vert[i][1] /= mag;
                    vert[i][2] /= mag;
                    mags[i] = mag;
                }

                for (i = 0; i < nt; i++)
                {
                    v[0] = vert[ts[i][0]];
                    v[1] = vert[ts[i][1]];
                    v[2] = vert[ts[i][2]];

                    l[0] = l[1] = l[2] = 0;
                    for (j = 0; j < 3; j++)
                    {
                        l[0] += (v[2][j] - v[1][j]) * (v[2][j] - v[1][j]);
                        l[1] += (v[0][j] - v[2][j]) * (v[0][j] - v[2][j]);
                        l[2] += (v[1][j] - v[0][j]) * (v[1][j] - v[0][j]);
                    }
                    l[0] = sqrt(l[0]);
                    l[1] = sqrt(l[1]);
                    l[2] = sqrt(l[2]);
                    assert(_finite(l[0]) && _finite(l[1]) && _finite(l[2]));

                    t[0] = 2 * (float)asin(l[0] / 2);
                    t[1] = 2 * (float)asin(l[1] / 2);
                    t[2] = 2 * (float)asin(l[2] / 2);
                    assert(_finite(t[0]) && _finite(t[1]) && _finite(t[2]));

                    s = 0.5f * (t[0] + t[1] + t[2]);
                    assert(_finite(s));

                    if (!(3.141592654 - s < EPSILN))
                    {
                        sinS = (float)sin(s);
                        assert(_finite(sinS));
                        sinT[0] = (float)sin(t[0]);
                        sinT[1] = (float)sin(t[1]);
                        sinT[2] = (float)sin(t[2]);
                        assert(_finite(sinT[0]) && _finite(sinT[1]) && _finite(sinT[2]));

                        cs[0] = (float)((2 * sinS * sin(s - t[0])) / (sinT[2] * sinT[1]) - 1);
                        cs[1] = (float)((2 * sinS * sin(s - t[1])) / (sinT[0] * sinT[2]) - 1);
                        cs[2] = (float)((2 * sinS * sin(s - t[2])) / (sinT[1] * sinT[0]) - 1);
                        //assert ( _finite ( cs [ 0 ] ) && _finite ( cs [ 1 ] ) && _finite ( cs [ 2 ] ) );

                        //			if ( fabs ( vol ) > EPSILN * EPSILN )
                        if ((1 - cs[0]) > EPSILN && (1 - cs[1]) > EPSILN && (1 - cs[2]) > EPSILN)
                        {
                            float cv[3];
                            cv[0] = v[1][1] * v[2][2] - v[1][2] * v[2][1];
                            cv[1] = v[1][2] * v[2][0] - v[1][0] * v[2][2];
                            cv[2] = v[1][0] * v[2][1] - v[1][1] * v[2][0];

                            vol = v[0][0] * cv[0] + v[0][1] * cv[1] + v[0][2] * cv[2];

                            if (vol < 0)
                            {
                                sinG[0] = -1 * sqrt(1 - cs[0] * cs[0]);
                                sinG[1] = -1 * sqrt(1 - cs[1] * cs[1]);
                                sinG[2] = -1 * sqrt(1 - cs[2] * cs[2]);
                            }
                            else
                            {
                                sinG[0] = sqrt(1 - cs[0] * cs[0]);
                                sinG[1] = sqrt(1 - cs[1] * cs[1]);
                                sinG[2] = sqrt(1 - cs[2] * cs[2]);
                            }

                            w[0] = (float)((t[0] - cs[1] * t[2] - cs[2] * t[1]) / (sinT[1] * sinG[2]));
                            w[1] = (float)((t[1] - cs[2] * t[0] - cs[0] * t[2]) / (sinT[2] * sinG[0]));
                            w[2] = (float)((t[2] - cs[0] * t[1] - cs[1] * t[0]) / (sinT[0] * sinG[1]));
                            // Added Begin
                            /*if ( !( _finite ( w [ 0 ] ) && _finite ( w [ 1 ] ) && _finite ( w [ 2 ] ) ) )
                            {
                            int breakpoint = 0;
                            }*/
                            // Added End
                            //assert ( _finite ( w [ 0 ] ) && _finite ( w [ 1 ] ) && _finite ( w [ 2 ] ) );

                            w[0] /= mags[ts[i][0]];
                            w[1] /= mags[ts[i][1]];
                            w[2] /= mags[ts[i][2]];

                            weights[ts[i][0]] += w[0];
                            weights[ts[i][1]] += w[1];
                            weights[ts[i][2]] += w[2];

                            totalW += (w[0] + w[1] + w[2]);
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        a[0] = (float)sin(t[0]) * mags[ts[i][2]] * mags[ts[i][1]];
                        a[1] = (float)sin(t[1]) * mags[ts[i][0]] * mags[ts[i][2]];
                        a[2] = (float)sin(t[2]) * mags[ts[i][1]] * mags[ts[i][0]];

                        totalA = a[0] + a[1] + a[2];

                        // we're on a face
                        a[0] /= totalA;
                        a[1] /= totalA;
                        a[2] /= totalA;

                        for (j = 0; j < nv; j++)
                        {
                            weights[j] = 0;
                        }

                        weights[ts[i][0]] = a[0];
                        weights[ts[i][1]] = a[1];
                        weights[ts[i][2]] = a[2];
                        return 0;
                    }
                }

                if (totalW > 0)
                {
                    for (j = 0; j < nv; j++)
                    {
                        weights[j] /= totalW;
                    }
                }

                if (totalW > 0)
                {
                    return 1;
                }
                else
                {
                    return -1;
                }
            };

            int _finite(float x)
            {
                if (fabs(x) < 1000)
                {
                    return 1;
                }
                return 0;
            }
        };

        class ParallelCell
        {


            int gridx, gridy, gridz;
            Edge *edges;
            Face *faces;
            Cell *cells;
            std::vector<Vertex> *vertices;
            std::vector<Point> *points;
            concurrent_vector<Point> *conpoints;
            concurrent_vector<Vertex> *convertices;

            float *maxEigenvalue, *minEigenvalue;
            float halfDataSize, *rightBackTop;
            float thickness;
            int *totalCellPoints;
            float gridSize;

            float eigenValue1, eigenValue2, eigenValue3;
            std::vector<float> eigenVectors;
            int dataType;
            float *gridPoints;
            float change;
            int sizex, sizey, sizez, halfSize;
            float*scalars, *tensors, *gradients;

            int getEdgeIndex(int axis, int x, int y, int z) const;
            int getFaceIndex(int axis, int x, int y, int z) const;

            void getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3])const;
            void getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const;
            void buildPointIteration(Cell *cell, concurrent_vector<Point> *points, concurrent_vector<Vertex> *vertices, float* maxEigenvalue, float* minEigenvalue) const;
            void getGridPointPosD(int index, float v[3]) const;
            void centroidProjection(int i, int j, int k, Cell *cell) const;
            bool getCentroid(int i, int j, int k, Cell *cell, int *ltotalCellPoints) const;
            void getProjDirection(int type, float position[3], float q[3]) const;
            bool getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            bool getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            //float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3],float* maxEigenvalue, float* minEigenvalue) const;
            float getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3], float* maxEigenvalue, float* minEigenvalue, float eigenValue1, float eigenValue2, std::vector<float> eigenVectors) const;
            std::vector<float> getEigenVectors(SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            std::vector<float> getFirstTwoEigenvalues(SelfAdjointEigenSolver<Matrix3f> *eigensolver) const;
            void cellPhaseIteration(int i, int j, int k, Cell *cells, concurrent_vector<Point> *lpoints, concurrent_vector<Vertex> *lvertices, int *ltotalCellPoints, float* lmaxEigenvalue, float* lminEigenvalue) const;

            void getShowPos(float position[3], float showPos[3]) const;

            void getScalarGP(int index, float *scalar) const
            {
                *scalar = scalars[index];
            }

            void getGradientGP(int index, float gradient[3]) const
            {
                gradient[0] = gradients[index];
                gradient[1] = gradients[index + 1];
                gradient[2] = gradients[index + 2];
            }

            void getTensorGP(int index, float tensor[6]) const
            {
                tensor[0] = tensors[index];
                tensor[1] = tensors[index + 1];
                tensor[2] = tensors[index + 2];
                tensor[3] = tensors[index + 3];
                tensor[4] = tensors[index + 4];
                tensor[5] = tensors[index + 5];
            }

            void getScalar(float position[3], float *scalar) const
            {
                trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradient(float position[3], float gradient[3]) const
            {
                trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensor(float position[3], float tensor[6]) const
            {
                trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
            };

            void getScalarCubic(float position[3], float *scalar) const
            {
                tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
            };

            void getGradientCubic(float position[3], float gradient[3]) const
            {
                tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
            };

            void getTensorCubic(float position[3], float tensor[6]) const
            {
                tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
            };
        public:

            ParallelCell(int lgridx, int lgridy, int lgridz, int *ltotalCellPoints, Edge *ledges, Face *lfaces, Cell *lcells, concurrent_vector<Point> *lpoints, int ldataType, float *lgridPoints, float lchange, int lsizex, int lsizey, int lsizez, int lhalfSize, float* lscalars, float* ltensors, float* lgradients, float lgridSize, float* lmaxEigenvalue, float* lminEigenvalue, concurrent_vector<Vertex> *lvertices, float lhalfDataSize, float *lrightBackTop, float lthickness)
            {
                gridx = lgridx;
                gridy = lgridy;
                gridz = lgridz;
                totalCellPoints = ltotalCellPoints;
                edges = ledges;
                faces = lfaces;
                cells = lcells;
                conpoints = lpoints;
                convertices = lvertices;
                dataType = ldataType;
                gridPoints = lgridPoints;


                change = lchange;
                sizex = lsizex;
                sizey = lsizey;
                sizez = lsizez;
                halfSize = lhalfSize;
                scalars = lscalars;
                tensors = ltensors;
                gradients = lgradients;
                gridSize = lgridSize;
                //
                maxEigenvalue = lmaxEigenvalue;
                minEigenvalue = lminEigenvalue;

                halfDataSize = lhalfDataSize;
                rightBackTop = lrightBackTop;
                thickness = lthickness;

            };

            void operator()(const blocked_range3d<int>& r)	const {
                int *ltotalCellPoints = totalCellPoints;
                //Edge *ledges = edges;
                //Face *lfaces = faces;
                Cell *lcells = cells;
                concurrent_vector<Point> *lpoints = conpoints;
                float *lmaxEigenvalue = maxEigenvalue;
                float *lminEigenvalue = minEigenvalue;
                concurrent_vector<Vertex> *lvertices = convertices;

                for (int i = r.cols().begin(); i != r.cols().end(); ++i)
                for (int j = r.rows().begin(); j != r.rows().end(); ++j)
                for (int k = r.pages().begin(); k != r.pages().end(); ++k)
                {
                    if (i == gridx - 1 || j == gridy - 1 || k == gridz - 1)
                        continue;
                    cellPhaseIteration(i, j, k, lcells, lpoints, lvertices, ltotalCellPoints, lmaxEigenvalue, lminEigenvalue);
                }
            };

        };


        void ParallelCell::cellPhaseIteration(int i, int j, int k, Cell *cells, concurrent_vector<Point> *points, concurrent_vector<Vertex> *vertices, int *ltotalCellPoints, float* lmaxEigenvalue, float* lminEigenvalue) const
        {
            Cell *cell = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
            if (!getCentroid(i, j, k, cell, ltotalCellPoints)) return;
            cell->valid = true;
            centroidProjection(i, j, k, cell);
            if (cell->extrP) buildPointIteration(cell, points, vertices, lmaxEigenvalue, lminEigenvalue);
        }

        bool ParallelCell::getCentroid(int i, int j, int k, Cell *cell, int *ltotalCellPoints) const
        {
            Face *cellFaces[6];
            cellFaces[0] = &(faces[getFaceIndex(1, i, j, k)]);
            if (!cellFaces[0]->valid) return false;
            cellFaces[1] = &(faces[getFaceIndex(1, i + 1, j, k)]);
            if (!cellFaces[1]->valid) return false;
            cellFaces[2] = &(faces[getFaceIndex(2, i, j, k)]);
            if (!cellFaces[2]->valid) return false;
            cellFaces[3] = &(faces[getFaceIndex(2, i, j + 1, k)]);
            if (!cellFaces[3]->valid) return false;
            cellFaces[4] = &(faces[getFaceIndex(3, i, j, k)]);
            if (!cellFaces[4]->valid) return false;
            cellFaces[5] = &(faces[getFaceIndex(3, i, j, k + 1)]);
            if (!cellFaces[5]->valid) return false;
            Edge *cellEdges[12];
            cellEdges[0] = &(edges[getEdgeIndex(1, i, j, k)]);
            cellEdges[1] = &(edges[getEdgeIndex(1, i, j + 1, k)]);
            cellEdges[2] = &(edges[getEdgeIndex(1, i, j + 1, k + 1)]);
            cellEdges[3] = &(edges[getEdgeIndex(1, i, j, k + 1)]);
            cellEdges[4] = &(edges[getEdgeIndex(2, i, j, k)]);
            cellEdges[5] = &(edges[getEdgeIndex(2, i + 1, j, k)]);
            cellEdges[6] = &(edges[getEdgeIndex(2, i + 1, j, k + 1)]);
            cellEdges[7] = &(edges[getEdgeIndex(2, i, j, k + 1)]);
            cellEdges[8] = &(edges[getEdgeIndex(3, i, j, k)]);
            cellEdges[9] = &(edges[getEdgeIndex(3, i + 1, j, k)]);
            cellEdges[10] = &(edges[getEdgeIndex(3, i + 1, j + 1, k)]);
            cellEdges[11] = &(edges[getEdgeIndex(3, i, j + 1, k)]);

            int ci = i;
            int cj = j;
            int ck = k;

            float area = 0;
            for (int i = 0; i < 3; i++)
            {
                area -= cellFaces[2 * i]->spherArea;
                area += cellFaces[2 * i + 1]->spherArea;
            }
            //cell->area = area;
            float wn = area / ((float)(4 * M_PI));
            wn = floor(wn + 0.5f);
            if ((int)wn % 2 != 0)
            {
                cell->extrP = true;

                float **grads;
                grads = new float*[8];
                for (int i = 0; i < 8; i++)
                {
                    grads[i] = new float[3];
                }
                float verts[8][3];
                for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                {
                    int index = getIndex(1, 0, ci + i, cj + j, ck + k, gridx, gridy, gridz);
                    getGridPointPosD(index, verts[i + 2 * j + 4 * k]);
                    index = getIndex(3, 0, ci + i + 2, cj + j + 2, ck + k + 2, sizex, sizey, sizez);
                    getGradientGP(index, grads[i + 2 * j + 4 * k]);
                }
                int **tris;
                tris = new int*[12];
                for (int i = 0; i < 12; i++)
                {
                    tris[i] = new int[3];
                }
                int static_tris[12][3] =
                { { 0, 2, 3 }, { 0, 3, 1 },
                { 5, 7, 6 }, { 5, 6, 4 },
                { 4, 6, 2 }, { 4, 2, 0 },
                { 1, 3, 7 }, { 1, 7, 5 },
                { 4, 0, 1 }, { 4, 1, 5 },
                { 2, 6, 7 }, { 2, 7, 3 } };
                for (int i = 0; i < 12; i++)
                for (int j = 0; j < 3; j++)
                    tris[i][j] = static_tris[i][j];

                MVC mvc(8, 12, grads, tris);
                float weights[8];
                float origin[3] = {};
                mvc.getWeights(origin, weights);
                float p[3];
                p[0] = 0;
                p[1] = 0;
                p[2] = 0;
                for (int i = 0; i < 8; i++)
                for (int j = 0; j < 3; j++)
                    p[j] += verts[i][j] * weights[i];
                if (p[0] < verts[0][0])
                {
                    p[0] = verts[0][0];
                }
                if (p[0] > verts[7][0])
                {
                    p[0] = verts[7][0];
                }
                if (p[1] < verts[0][1])
                {
                    p[1] = verts[0][1];
                }
                if (p[1] > verts[7][1])
                {
                    p[1] = verts[7][1];
                }
                if (p[2] < verts[0][2])
                {
                    p[2] = verts[0][2];
                }
                if (p[2] > verts[7][2])
                {
                    p[2] = verts[7][2];
                }
                for (int i = 0; i < 3; i++)
                    cell->extremalPoint[i] = p[i];

                float p10[3], p11[3], p13[3], p14[3];
                p10[0] = p[0] - 2 * change;
                p10[1] = p[1];
                p10[2] = p[2];
                p11[0] = p[0] - change;
                p11[1] = p[1];
                p11[2] = p[2];
                p13[0] = p[0] + change;
                p13[1] = p[1];
                p13[2] = p[2];
                p14[0] = p[0] + 2 * change;
                p14[1] = p[1];
                p14[2] = p[2];

                float p20[3], p21[3], p23[3], p24[3];
                p20[0] = p[0];
                p20[1] = p[1] - 2 * change;
                p20[2] = p[2];
                p21[0] = p[0];
                p21[1] = p[1] - change;
                p21[2] = p[2];
                p23[0] = p[0];
                p23[1] = p[1] + change;
                p23[2] = p[2];
                p24[0] = p[0];
                p24[1] = p[1] + 2 * change;
                p24[2] = p[2];

                float p30[3], p31[3], p33[3], p34[3];
                p30[0] = p[0];
                p30[1] = p[1];
                p30[2] = p[2] - 2 * change;
                p31[0] = p[0];
                p31[1] = p[1];
                p31[2] = p[2] - change;
                p33[0] = p[0];
                p33[1] = p[1];
                p33[2] = p[2] + change;
                p34[0] = p[0];
                p34[1] = p[1];
                p34[2] = p[2] + 2 * change;

                float f, f10, f11, f13, f14, f20, f21, f23, f24, f30, f31, f33, f34;

                getScalarCubic(p, &f);
                getScalarCubic(p10, &f10);
                getScalarCubic(p11, &f11);
                getScalarCubic(p13, &f13);
                getScalarCubic(p14, &f14);
                getScalarCubic(p20, &f20);
                getScalarCubic(p21, &f21);
                getScalarCubic(p23, &f23);
                getScalarCubic(p24, &f24);
                getScalarCubic(p30, &f30);
                getScalarCubic(p31, &f31);
                getScalarCubic(p33, &f33);
                getScalarCubic(p34, &f34);

                float secondDiff1 = -f10 + 16 * f11 - 30 * f + 16 * f13 - f14;
                float secondDiff2 = -f20 + 16 * f21 - 30 * f + 16 * f23 - f24;
                float secondDiff3 = -f30 + 16 * f31 - 30 * f + 16 * f33 - f34;
                if (secondDiff1<0 && secondDiff2<0 && secondDiff3<0) cell->extrPType = 1;
                else if (secondDiff1>0 && secondDiff2>0 && secondDiff3>0) cell->extrPType = 2;
                else cell->extrPType = 3;
            }


            cell->centroid[0] = 0;
            cell->centroid[1] = 0;
            cell->centroid[2] = 0;
            int count = 0;
            bool minEP, maxEP, minFP, maxFP, other;
            minEP = maxEP = minFP = maxFP = other = false;

            for (int i = 0; i < 12; i++)
            {
                Edge *edge = cellEdges[i];
                if (edge->extremal)
                {
                    cell->extremal = true;
                    if (edge->edgeTag == 1) maxEP = true;
                    else minEP = true;
                    cell->centroid[0] += edge->edgePoint[0];
                    cell->centroid[1] += edge->edgePoint[1];
                    cell->centroid[2] += edge->edgePoint[2];
                    count++;
                }
            }

            for (int i = 0; i < 6; i++)
            {
                Face *face = cellFaces[i];
                if (face->extremal)
                {
                    cell->extremal = true;
                    if (face->faceTag == 1) maxFP = true;
                    else if (face->faceTag == 2) minFP = true;
                    else other = true;
                    cell->centroid[0] += face->facePoint[0];
                    cell->centroid[1] += face->facePoint[1];
                    cell->centroid[2] += face->facePoint[2];
                    count++;
                }
            }

            if (cell->extremal)
            {
                cell->centroid[0] /= (float)count;
                cell->centroid[1] /= (float)count;
                cell->centroid[2] /= (float)count;
                (*ltotalCellPoints)++;
            }
            if (!other && !((minEP || minFP) && (maxEP || maxFP)))
            {
                if (minFP) cell->projType = 3;
                else if (maxFP) cell->projType = 4;
                else if (minEP) cell->projType = 1;
                else cell->projType = 2;
            }
            return true;
        }

        void ParallelCell::centroidProjection(int i, int j, int k, Cell *cell) const
        {
            if (cell->extremal && cell->projType > 0)
            {
                float border1[3], border2[3];
                getGridPointPosD(getIndex(1, 0, i, j, k, gridx, gridy, gridz), border1);
                getGridPointPosD(getIndex(1, 0, i + 1, j + 1, k + 1, gridx, gridy, gridz), border2);
                float projThresh = 0.001f * gridSize;
                float stepSize = 0.1f * gridSize;
                float p[3];
                p[0] = cell->centroid[0];
                p[1] = cell->centroid[1];
                p[2] = cell->centroid[2];
                int type = cell->projType;
                float q[3];
                getProjDirection(type, p, q);
                float qNorm = getNorm(q);
                int count = 0;
                float p1[3];
                float q1[3];
                float q1Norm;
                while (qNorm > projThresh && count < 10)
                {
                    count++;
                    p1[0] = p[0] - stepSize * (q[0] / qNorm);
                    p1[1] = p[1] - stepSize * (q[1] / qNorm);
                    p1[2] = p[2] - stepSize * (q[2] / qNorm);
                    int dCount = 0;
                    if (p1[0] < border1[0])
                    {
                        p1[0] = border1[0];
                        dCount++;
                    }
                    if (p1[0] > border2[0])
                    {
                        p1[0] = border2[0];
                        dCount++;
                    }
                    if (p1[1] < border1[1])
                    {
                        p1[1] = border1[1];
                        dCount++;
                    }
                    if (p1[1] > border2[1])
                    {
                        p1[1] = border2[1];
                        dCount++;
                    }
                    if (p1[2] < border1[2])
                    {
                        p1[2] = border1[2];
                        dCount++;
                    }
                    if (p1[2] > border2[2])
                    {
                        p1[2] = border2[2];
                        dCount++;
                    }
                    getProjDirection(type, p1, q1);
                    q1Norm = getNorm(q1);
                    if (q1Norm > qNorm)
                    {
                        stepSize /= 2;
                        continue;
                    }
                    p[0] = p1[0];
                    p[1] = p1[1];
                    p[2] = p1[2];
                    if (dCount == 3) break;
                    q[0] = q1[0];
                    q[1] = q1[1];
                    q[2] = q1[2];
                    qNorm = q1Norm;
                }
                cell->centroid[0] = p[0];
                cell->centroid[1] = p[1];
                cell->centroid[2] = p[2];
            }
            if (cell->extrP && (cell->extrPType == 1 || cell->extrPType == 2))
            {
                float border1[3], border2[3];
                getGridPointPosD(getIndex(1, 0, i, j, k, gridx, gridy, gridz), border1);
                getGridPointPosD(getIndex(1, 0, i + 1, j + 1, k + 1, gridx, gridy, gridz), border2);
                float projThresh = 0.001f * gridSize;
                float stepSize = 0.01f * gridSize;
                if (cell->extrPType == 1) stepSize = -stepSize;
                float p[3];
                p[0] = cell->extremalPoint[0];
                p[1] = cell->extremalPoint[1];
                p[2] = cell->extremalPoint[2];
                float q[3];
                getGradientCubic(p, q);
                float qNorm = getNorm(q);
                int count = 0;
                float p1[3];
                float q1[3];
                float q1Norm;
                while (qNorm > projThresh && count < 10)
                {
                    count++;
                    p1[0] = p[0] - stepSize * (q[0] / qNorm);
                    p1[1] = p[1] - stepSize * (q[1] / qNorm);
                    p1[2] = p[2] - stepSize * (q[2] / qNorm);
                    int dCount = 0;
                    if (p1[0] < border1[0])
                    {
                        p1[0] = border1[0];
                        dCount++;
                    }
                    if (p1[0] > border2[0])
                    {
                        p1[0] = border2[0];
                        dCount++;
                    }
                    if (p1[1] < border1[1])
                    {
                        p1[1] = border1[1];
                        dCount++;
                    }
                    if (p1[1] > border2[1])
                    {
                        p1[1] = border2[1];
                        dCount++;
                    }
                    if (p1[2] < border1[2])
                    {
                        p1[2] = border1[2];
                        dCount++;
                    }
                    if (p1[2] > border2[2])
                    {
                        p1[2] = border2[2];
                        dCount++;
                    }
                    getGradientCubic(p1, q1);
                    q1Norm = getNorm(q1);
                    if (q1Norm > qNorm)
                    {
                        stepSize /= 10;
                        continue;
                    }
                    p[0] = p1[0];
                    p[1] = p1[1];
                    p[2] = p1[2];
                    if (dCount == 3) break;
                    q[0] = q1[0];
                    q[1] = q1[1];
                    q[2] = q1[2];
                    qNorm = q1Norm;
                }
                cell->extremalPoint[0] = p[0];
                cell->extremalPoint[1] = p[1];
                cell->extremalPoint[2] = p[2];
            }
        }

        void ParallelCell::buildPointIteration(Cell *cell, concurrent_vector<Point> *points, concurrent_vector<Vertex> *vertices, float* lmaxEigenvalue, float* lminEigenvalue) const
        {
            Vertex vertex;
            Point point;
            getShowPos(cell->extremalPoint, vertex.position);
            point.vertIdx = vertices->size();
            vertices->push_back(vertex);
            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getEigensolverCubic(cell->extremalPoint, &eigensolver);
            point.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, point.relativeSaliencies, lmaxEigenvalue, lminEigenvalue, point.eigenValue1, point.eigenValue2, point.eigenVectors);
            point.eigenVectors = getEigenVectors(&eigensolver);
            std::vector<float> firstTwoEigenvalues = getFirstTwoEigenvalues(&eigensolver);
            point.eigenValue1 = firstTwoEigenvalues[0];
            point.eigenValue2 = firstTwoEigenvalues[1];
            getScalarCubic(cell->extremalPoint, &(point.localIntensity));
            point.type = cell->extrPType;
            points->push_back(point);
        }

        void ParallelCell::getProjDirection(int type, float position[3], float q[3]) const
        {
            float g[3];
            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getGradientCubic(position, g);
            getEigensolverCubic(position, &eigensolver);
            float v1[3], v3[3];
            getV1(&eigensolver, v1);
            getV3(&eigensolver, v3);
            float dotP, crossP[3];
            switch (type)
            {
            case 1:
            case 2:
                dotP = getDotP(g, v1);
                q[0] = dotP*v1[0];
                q[1] = dotP*v1[1];
                q[2] = dotP*v1[2];
                break;
            case 3:
            case 4:
                getCrossP(v3, g, crossP);
                getCrossP(crossP, v3, q);
                break;
            }
            if (type == 2 || type == 4)
            {
                q[0] = -q[0];
                q[1] = -q[1];
                q[2] = -q[2];
            }
        }

        int ParallelCell::getEdgeIndex(int axis, int x, int y, int z) const {
            int index;
            switch (axis)
            {
            case 1:
                index = x * gridy * gridz + y * gridz + z;
                return index;
                break;
            case 2:
                index = (gridx - 1) * gridy * gridz +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            case 3:
                index = (gridx - 1) * gridy * gridz +
                    gridx * (gridy - 1) * gridz +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            }
            return 0;
        }

        int ParallelCell::getFaceIndex(int axis, int x, int y, int z) const {
            int index;
            switch (axis)
            {
            case 1:
                index = x * (gridy - 1) * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 2:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            case 3:
                index = gridx * (gridy - 1) * (gridz - 1) +
                    (gridx - 1) * gridy * (gridz - 1) +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            }
            return 0;
        }

        void ParallelCell::getGridPointPosD(int index, float v[3]) const
        {
            v[0] = gridPoints[3 * index];
            v[1] = gridPoints[3 * index + 1];
            v[2] = gridPoints[3 * index + 2];
        }

        void ParallelCell::getShowPos(float position[3], float showPos[3]) const
        {
            showPos[0] = ((float)position[0] - 2 - rightBackTop[0]) / halfDataSize;
            showPos[1] = ((float)position[1] - 2 - rightBackTop[1]) / halfDataSize;
            showPos[2] = ((float)position[2] - 2 - rightBackTop[2]) / halfDataSize * thickness;
        }

        bool ParallelCell::getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorCubic(position, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }
        bool ParallelCell::getEigensolverGP(int index, SelfAdjointEigenSolver<Matrix3f> *eigensolver) const
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorGP(index, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }

        std::vector<float> ParallelCell::getEigenVectors(SelfAdjointEigenSolver<Matrix3f> *eigensolver) const {
            Vector3f eVectors0 = eigensolver->eigenvectors().col(0).normalized();
            Vector3f eVectors1 = eigensolver->eigenvectors().col(1).normalized();
            Vector3f eVectors2 = eigensolver->eigenvectors().col(2).normalized();
            std::vector<float> eVectors;
            eVectors.push_back(eVectors0(0));
            eVectors.push_back(eVectors0(1));
            eVectors.push_back(eVectors0(2));
            eVectors.push_back(eVectors1(0));
            eVectors.push_back(eVectors1(1));
            eVectors.push_back(eVectors1(2));
            eVectors.push_back(eVectors2(0));
            eVectors.push_back(eVectors2(1));
            eVectors.push_back(eVectors2(2));
            return eVectors;
        }

        std::vector<float> ParallelCell::getFirstTwoEigenvalues(SelfAdjointEigenSolver<Matrix3f> *eigensolver) const {
            Vector3f ev = eigensolver->eigenvalues();
            std::vector<float> twoEigenvalues;
            twoEigenvalues.push_back(ev(0));
            twoEigenvalues.push_back(ev(1));
            return twoEigenvalues;
        }

        float ParallelCell::getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3], float* maxEigenvalue, float* minEigenvalue, float eigenValue1, float eigenValue2, std::vector<float> eigenVectors) const
        {
            Vector3f ev = eigensolver->eigenvalues();
            if (ev(2) > 0)
            {
                relativeSaliencies[0] = float((ev(2) - ev(1)) / ev(2));
                relativeSaliencies[1] = float((ev(1) - ev(0)) / ev(2));
                relativeSaliencies[2] = float(ev(0) / ev(2));
                if ((*minEigenvalue) > ev(2)) (*minEigenvalue) = ev(2);
                if ((*maxEigenvalue) < ev(2)) (*maxEigenvalue) = ev(2);
            }
            else
            {
                relativeSaliencies[0] = 0;
                relativeSaliencies[1] = 0;
                relativeSaliencies[2] = 0;
            }
            eigenValue1 = ev(0);
            eigenValue2 = ev(1);
            Vector3f eVectors0 = eigensolver->eigenvectors().col(0).normalized();
            Vector3f eVectors1 = eigensolver->eigenvectors().col(1).normalized();
            Vector3f eVectors2 = eigensolver->eigenvectors().col(2).normalized();
            eigenVectors.push_back(eVectors0(0));
            eigenVectors.push_back(eVectors0(1));
            eigenVectors.push_back(eVectors0(2));
            eigenVectors.push_back(eVectors1(0));
            eigenVectors.push_back(eVectors1(1));
            eigenVectors.push_back(eVectors1(2));
            eigenVectors.push_back(eVectors2(0));
            eigenVectors.push_back(eVectors2(1));
            eigenVectors.push_back(eVectors2(2));


            return float(ev(2));
        }

        void ParallelCell::getV1(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev;
            switch (dataType)
            {
            case 1:
                ev = eigensolver->eigenvectors().col(2).normalized();
                break;
            case 2:
                ev = eigensolver->eigenvectors().col(0).normalized();
                break;
            }
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }
        void ParallelCell::getV3(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float v[3]) const
        {
            Vector3f ev;
            switch (dataType)
            {
            case 1:
                ev = eigensolver->eigenvectors().col(0).normalized();
                break;
            case 2:
                ev = eigensolver->eigenvectors().col(2).normalized();
                break;
            }
            v[0] = (float)ev(0);
            v[1] = (float)ev(1);
            v[2] = (float)ev(2);
        }

        void General_Data::cellPhase(float* scalars, float* tensors, float* gradients)
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Cell Phase Begin:" << endl;
            points = new std::vector<Point>;
            vertices = new std::vector<Vertex>;
            conpoints = new concurrent_vector<Point>;
            convertices = new concurrent_vector<Vertex>;
            change = 0.017f * gridSize;
            dispCellX = (int)floor(gridx / 2.0) - 1;
            dispCellY = (int)floor(gridy / 2.0) - 1;
            dispCellZ = (int)floor(gridz / 2.0) - 1;
            maxCellIndex = (gridx - 1) * (gridy - 1) * (gridz - 1);
            cells = new Cell[maxCellIndex];
            totalCellPoints = new int(0);

            tbb::task_scheduler_init init(task_scheduler_init::automatic);
            ParallelCell parallel_cell(gridx, gridy, gridz, totalCellPoints, edges, faces, cells, conpoints, dataType, gridPoints, change, sizex, sizey, sizez, halfSize, scalars, tensors, gradients, gridSize, maxEigenvalue, minEigenvalue, convertices, halfDataSize, rightBackTop, thickness);
            parallel_for(blocked_range3d<int>(0, gridz, 0, gridy, 0, gridx), parallel_cell, auto_partitioner());

            concurrent_vector<Point>::iterator iter = conpoints->begin();
            while (iter != conpoints->end())
            {
                Point p;
                p.vertIdx = (*iter).vertIdx;
                p.relativeSaliencies[0] = (*iter).relativeSaliencies[0];
                p.relativeSaliencies[1] = (*iter).relativeSaliencies[1];
                p.relativeSaliencies[2] = (*iter).relativeSaliencies[2];
                p.localIntensity = (*iter).localIntensity;
                p.firstEigenvalue = (*iter).firstEigenvalue;
                p.type = (*iter).type;
                points->push_back(p);
                iter++;
            }

            concurrent_vector<Vertex>::iterator iterv = convertices->begin();
            while (iterv != convertices->end())
            {
                Vertex p;
                p.position[0] = (*iterv).position[0];
                p.position[1] = (*iterv).position[1];
                p.position[2] = (*iterv).position[2];
                vertices->push_back(p);
                iterv++;
            }
            timeEnd = clock();
            //t_cellPhase_other = timeEnd-timeStart - t_cellPhase_extremalPoint;
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            /*cout << "\tt_cellPhase_other: " << t_cellPhase_other << endl;
            cout << "\tt_cellPhase_extremalPoint: " << t_cellPhase_extremalPoint << endl;*/
            //cout << "\tt_cellPhase_interp: " << t_cellPhase_interp/CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }

        void General_Data::getShowPos(float position[3], float showPos[3])
        {
            showPos[0] = ((float)position[0] - 2 - rightBackTop[0]) / halfDataSize;
            showPos[1] = ((float)position[1] - 2 - rightBackTop[1]) / halfDataSize;
            showPos[2] = ((float)position[2] - 2 - rightBackTop[2]) / halfDataSize;
        }

        void General_Data::reverseShowPos(float position[3], float showPos[3])
        {
            //sizex = sizex-2-2*halfSize;
            //				sizey = sizey-2-2*halfSize;
            //			sizez = sizez-2-2*halfSize;
            //		rightBackTop[0] = ((float)sizex-5) / 2;
            //	rightBackTop[1] = ((float)sizey-5) / 2;
            //rightBackTop[2] = ((float)sizez-5) / 2;

            float xScale = ((float)volData->GetSizeX()) / ((float)volData->GetSizeX() - 8.0);
            float yScale = ((float)volData->GetSizeY()) / ((float)volData->GetSizeY() - 8.0);
            float zScale = ((float)volData->GetSizeZ()) / ((float)volData->GetSizeZ() - 8.0);

            float spacingTransX = (volData->GetSpacingX() - 1.0) * ((float)volData->GetSizeX() - 8.0);
            float spacingTransY = (volData->GetSpacingY() - 1.0) * ((float)volData->GetSizeY() - 8.0);
            float spacingTransZ = (volData->GetSpacingZ() - 1.0) * ((float)volData->GetSizeZ() - 8.0);
            // + rightBackTop[0] + 2.0 + kernelSize
            position[0] = volData->GetSpacingX() * ((showPos[0] * xScale * (halfDataSize - 2.0)) + rightBackTop[0] + (kernelSize / 2.0) + 3.0);
            position[1] = volData->GetSpacingY() * ((showPos[1] * yScale * (halfDataSize - 2.0)) + rightBackTop[1] + (kernelSize / 2.0) + 3.0);
            position[2] = volData->GetSpacingZ() * ((showPos[2] * zScale * (halfDataSize - 2.0)) + rightBackTop[2] + (kernelSize / 2.0) + 3.0);
        }



        bool General_Data::getEigensolverCubic(float position[3], SelfAdjointEigenSolver<Matrix3f> *eigensolver)
        {
            int count = 0;
            float epsilon = 1e-12;
            float tensorV[6];
            getTensorCubic(position, tensorV);
            for (int i = 0; i < 6; i++)
            {
                if (fabs(tensorV[i]) < epsilon)
                {
                    count++;
                    tensorV[i] = epsilon;
                }
            }
            Matrix3f tensorM;
            tensorM << tensorV[0], tensorV[1], tensorV[2],
                tensorV[1], tensorV[3], tensorV[4],
                tensorV[2], tensorV[4], tensorV[5];
            eigensolver->compute(tensorM);
            if (count > 1) return false;
            return true;
        }

        std::vector<float> General_Data::getEigenVectors(SelfAdjointEigenSolver<Matrix3f> *eigensolver) {
            Vector3f eVectors0 = eigensolver->eigenvectors().col(0).normalized();
            Vector3f eVectors1 = eigensolver->eigenvectors().col(1).normalized();
            Vector3f eVectors2 = eigensolver->eigenvectors().col(2).normalized();
            std::vector<float> eVectors;
            eVectors.push_back(eVectors0(0));
            eVectors.push_back(eVectors0(1));
            eVectors.push_back(eVectors0(2));
            eVectors.push_back(eVectors1(0));
            eVectors.push_back(eVectors1(1));
            eVectors.push_back(eVectors1(2));
            eVectors.push_back(eVectors2(0));
            eVectors.push_back(eVectors2(1));
            eVectors.push_back(eVectors2(2));
            return eVectors;
        }

        std::vector<float> General_Data::getFirstTwoEigenvalues(SelfAdjointEigenSolver<Matrix3f> *eigensolver) {
            Vector3f ev = eigensolver->eigenvalues();
            std::vector<float> first2Eigenvalues;
            first2Eigenvalues.push_back(ev(0));
            first2Eigenvalues.push_back(ev(0));
            return first2Eigenvalues;
        }

        float General_Data::getFirstEigenvalueAndRelativeSaliencies(SelfAdjointEigenSolver<Matrix3f> *eigensolver, float relativeSaliencies[3], float eigenValue1, float eigenValue2, std::vector<float> eigenVectors)
        {
            Vector3f ev = eigensolver->eigenvalues();
            if (ev(2) > 0)
            {
                relativeSaliencies[0] = float((ev(2) - ev(1)) / ev(2));
                relativeSaliencies[1] = float((ev(1) - ev(0)) / ev(2));
                relativeSaliencies[2] = float(ev(0) / ev(2));
                if ((*minEigenvalue) > ev(2)) (*minEigenvalue) = ev(2);
                if ((*maxEigenvalue) < ev(2)) (*maxEigenvalue) = ev(2);
            }
            else
            {
                relativeSaliencies[0] = 0;
                relativeSaliencies[1] = 0;
                relativeSaliencies[2] = 0;
            }
            eigenValue1 = ev(0);
            eigenValue2 = ev(1);
            /**
            Vector3f eVectors0 = eigensolver->eigenvectors().col(0).normalized();
            Vector3f eVectors1 = eigensolver->eigenvectors().col(1).normalized();
            Vector3f eVectors2 = eigensolver->eigenvectors().col(2).normalized();

            eigenVectors.push_back(eVectors0(0));
            eigenVectors.push_back(eVectors0(1));
            eigenVectors.push_back(eVectors0(2));
            eigenVectors.push_back(eVectors1(0));
            eigenVectors.push_back(eVectors1(1));
            eigenVectors.push_back(eVectors1(2));
            eigenVectors.push_back(eVectors2(0));
            eigenVectors.push_back(eVectors2(1));
            eigenVectors.push_back(eVectors2(2));
            **/
            return float(ev(2));
        }
        /**
        float General_Data::getFirstEigenvalueSaliencies(float* eigenvalues, float relativeSaliencies[3])
        {
        Vector3DFloat ev = Vector3DFloat(eigenvalues[0], eigenvalues[1], eigenvalues[2]);
        if ( ev.Z() > 0 )
        {
        relativeSaliencies[0] = float((ev.Z() - ev.Y() ) / ev.Z() );
        relativeSaliencies[1] = float((ev.Y()  - ev.X() ) / ev.Z() );
        relativeSaliencies[2] = float(ev.X() / ev.Z() );
        if ((*minEigenvalue) > ev.Z() ) (*minEigenvalue) = ev.Z() ;
        if ((*maxEigenvalue) < ev.Z() ) (*maxEigenvalue) = ev.Z() ;
        }
        else
        {
        relativeSaliencies[0] = 0;
        relativeSaliencies[1] = 0;
        relativeSaliencies[2] = 0;
        }
        return float(ev.Z() );
        }
        **/

        void General_Data::buildCurveIteration(Face *face, Cell *cell1, Cell *cell2)
        {
            Segment segment;
            if (cell1->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell1->centroid, vertex.position);
                segment.vertIdxs[0] = vertices->size();
                cell1->vertIdx = segment.vertIdxs[0];
                vertices->push_back(vertex);
            }
            else
            {
                segment.vertIdxs[0] = cell1->vertIdx;
            }
            if (cell2->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell2->centroid, vertex.position);
                segment.vertIdxs[1] = vertices->size();
                cell2->vertIdx = segment.vertIdxs[1];
                vertices->push_back(vertex);
            }
            else
            {
                segment.vertIdxs[1] = cell2->vertIdx;
            }
            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getEigensolverCubic(face->facePoint, &eigensolver);
            segment.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, segment.relativeSaliencies, segment.eigenValue1, segment.eigenValue2, segment.eigenVectors);
            segment.eigenVectors = getEigenVectors(&eigensolver);
            std::vector<float> first2Eigenvalues = getFirstTwoEigenvalues(&eigensolver);
            segment.eigenValue1 = first2Eigenvalues[0];
            segment.eigenValue2 = first2Eigenvalues[1];
            float localI;
            getScalarCubic(face->facePoint, &localI);
            if (minIntensity > localI) minIntensity = localI;
            if (maxIntensity < localI) maxIntensity = localI;
            segment.localIntensity = localI;
            segment.type = face->faceTag;
            segments->push_back(segment);
        }


        void General_Data::buildCurve()
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Build Curve Begin:" << endl;
            segments = new std::vector<Segment>;
            for (int i = 0; i < gridx; i++)
            for (int j = 0; j < gridy - 1; j++)
            for (int k = 0; k < gridz - 1; k++)
            {
                Face *face = &(faces[getFaceIndex(1, i, j, k)]);
                if (!face->valid) continue;
                if (!face->extremal || i == 0 || i == gridx - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i - 1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                buildCurveIteration(face, cell1, cell2);
            }
            for (int i = 0; i < gridx - 1; i++)
            for (int j = 0; j < gridy; j++)
            for (int k = 0; k < gridz - 1; k++)
            {
                Face *face = &(faces[getFaceIndex(2, i, j, k)]);
                if (!face->valid) continue;
                if (!face->extremal || j == 0 || j == gridy - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i, j - 1, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                buildCurveIteration(face, cell1, cell2);
            }
            for (int i = 0; i < gridx - 1; i++)
            for (int j = 0; j < gridy - 1; j++)
            for (int k = 0; k < gridz; k++)
            {
                Face *face = &(faces[getFaceIndex(3, i, j, k)]);
                if (!face->valid) continue;
                if (!face->extremal || k == 0 || k == gridz - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i, j, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                buildCurveIteration(face, cell1, cell2);
            }
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }


        int General_Data::getEdgeIndex(int axis, int x, int y, int z) {
            int index;
            switch (axis)
            {
            case 1:
                index = x * gridy * gridz + y * gridz + z;
                return index;
                break;
            case 2:
                index = (gridx - 1) * gridy * gridz +
                    x * (gridy - 1) * gridz + y * gridz + z;
                return index;
                break;
            case 3:
                index = (gridx - 1) * gridy * gridz +
                    gridx * (gridy - 1) * gridz +
                    x * gridy * (gridz - 1) + y * (gridz - 1) + z;
                return index;
                break;
            }
            return 0;
        }


        std::vector<Point> * General_Data::getPoints() {
            return points;
        }

        void General_Data::buildSurfaceIteration(Edge *edge, Cell *cell1, Cell *cell2, Cell *cell3, Cell *cell4)
        {
            Quad quad;
            if (cell1->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell1->centroid, vertex.position);
                quad.vertIdxs[0] = vertices->size();
                cell1->vertIdx = quad.vertIdxs[0];
                vertices->push_back(vertex);
            }
            else
            {
                quad.vertIdxs[0] = cell1->vertIdx;
            }
            if (cell2->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell2->centroid, vertex.position);
                quad.vertIdxs[1] = vertices->size();
                cell2->vertIdx = quad.vertIdxs[1];
                vertices->push_back(vertex);
            }
            else
            {
                quad.vertIdxs[1] = cell2->vertIdx;
            }
            if (cell3->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell3->centroid, vertex.position);
                quad.vertIdxs[2] = vertices->size();
                cell3->vertIdx = quad.vertIdxs[2];
                vertices->push_back(vertex);
            }
            else
            {
                quad.vertIdxs[2] = cell3->vertIdx;
            }
            if (cell4->vertIdx < 0)
            {
                Vertex vertex;
                getShowPos(cell4->centroid, vertex.position);
                quad.vertIdxs[3] = vertices->size();
                cell4->vertIdx = quad.vertIdxs[3];
                vertices->push_back(vertex);
            }
            else
            {
                quad.vertIdxs[3] = cell4->vertIdx;
            }
            SelfAdjointEigenSolver<Matrix3f> eigensolver;
            getEigensolverCubic(edge->edgePoint, &eigensolver);
            quad.firstEigenvalue = getFirstEigenvalueAndRelativeSaliencies(&eigensolver, quad.relativeSaliencies, quad.eigenValue1, quad.eigenValue2, quad.eigenVectors);
            quad.eigenVectors = getEigenVectors(&eigensolver);
            std::vector<float> first2Eigenvalues = getFirstTwoEigenvalues(&eigensolver);

            quad.eigenValue1 = first2Eigenvalues[0];
            quad.eigenValue2 = first2Eigenvalues[1];
            getScalarCubic(edge->edgePoint, &(quad.localIntensity));
            quad.type = edge->edgeTag;
            quads->push_back(quad);
        }

        void General_Data::buildSurface()
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Build Surface Begin:" << endl;
            quads = new std::vector<Quad>;
            for (int i = 0; i < gridx - 1; i++)
            for (int j = 0; j < gridy; j++)
            for (int k = 0; k < gridz; k++)
            {
                Edge *edge = &(edges[getEdgeIndex(1, i, j, k)]);
                if (!edge->valid) continue;
                if (!edge->extremal || j == 0 || j == gridy - 1 || k == 0 || k == gridz - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i, j - 1, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                Cell *cell3 = &(cells[getIndex(1, 0, i, j - 1, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell3->valid) continue;
                Cell *cell4 = &(cells[getIndex(1, 0, i, j, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell4->valid) continue;
                buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
            }
            for (int i = 0; i < gridx; i++)
            for (int j = 0; j < gridy - 1; j++)
            for (int k = 0; k < gridz; k++)
            {
                Edge *edge = &(edges[getEdgeIndex(2, i, j, k)]);
                if (!edge->valid) continue;
                if (!edge->extremal || i == 0 || i == gridx - 1 || k == 0 || k == gridz - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i - 1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                Cell *cell3 = &(cells[getIndex(1, 0, i - 1, j, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell3->valid) continue;
                Cell *cell4 = &(cells[getIndex(1, 0, i, j, k - 1, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell4->valid) continue;
                buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
            }
            for (int i = 0; i < gridx; i++)
            for (int j = 0; j < gridy; j++)
            for (int k = 0; k < gridz - 1; k++)
            {
                Edge *edge = &(edges[getEdgeIndex(3, i, j, k)]);
                if (!edge->valid) continue;
                if (!edge->extremal || i == 0 || i == gridx - 1 || j == 0 || j == gridy - 1) continue;
                Cell *cell1 = &(cells[getIndex(1, 0, i, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell1->valid) continue;
                Cell *cell2 = &(cells[getIndex(1, 0, i, j - 1, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell2->valid) continue;
                Cell *cell3 = &(cells[getIndex(1, 0, i - 1, j - 1, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell3->valid) continue;
                Cell *cell4 = &(cells[getIndex(1, 0, i - 1, j, k, gridx - 1, gridy - 1, gridz - 1)]);
                if (!cell4->valid) continue;
                buildSurfaceIteration(edge, cell1, cell2, cell3, cell4);
            }
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }


        class MRC_Data :
            public General_Data
        {
        public:
            MRC_Data();
            ~MRC_Data(void);
            void getScalarGP(int index, float *scalar);
            void getTensorGP(int index, float tensor[6]);
            void getGradientGP(int index, float gradient[3]);
            void getScalar(float position[3], float *scalar);
            void getTensor(float position[3], float tensor[6]);
            void getGradient(float position[3], float gradient[3]);
            void getScalarCubic(float position[3], float *scalar);
            void getTensorCubic(float position[3], float tensor[6]);
            void getGradientCubic(float position[3], float gradient[3]);
            void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6]);
            void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3]);
            void dataGeneration(char *filename);
            void buildGrid();
            void getGradientBicubicTable(int axis, int ii, int jj, int kk, float sGradients[subdNum][subdNum][3]);
        private:
            void MRCGeneration(char *filename);
            //void DTIGeneration(char *filename);


        public:
            float *scalars;
            float *tensors;
            float *gradients;
            int iteration;
            float *edgeTable;
            float *faceTable;
        };



        MRC_Data::MRC_Data() :General_Data()
        {
            iteration = 4;
            scalars = NULL;
            tensors = NULL;
            gradients = NULL;
        }

        MRC_Data::~MRC_Data(void)
        {
            delete[] scalars;
            delete[] tensors;
            delete[] gradients;
            delete[] edgeTable;
            delete[] faceTable;
        }

        void MRC_Data::dataGeneration(char *filename)
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Read Data Begin:" << endl;
            switch (dataType)
            {
            case 1:
                MRCGeneration(filename);
                break;
            case 2:
                MRCGeneration(filename);
                //DTIGeneration(filename);
                break;
            }
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }

        void gaussianiir3d(float *volume, long width, long height, long depth,
            float sigma, int numsteps)
        {
            const long plane = width*height;
            const long numel = plane*depth;
            double lambda, dnu;
            float nu, boundaryscale, postscale;
            float *ptr;
            long i, x, y, z;
            int step;

            if (sigma <= 0 || numsteps < 0)
                return;

            lambda = (sigma*sigma) / (2.0*numsteps);
            dnu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda)) / (2.0*lambda);
            nu = (float)dnu;
            boundaryscale = (float)(1.0 / (1.0 - dnu));
            postscale = (float)(pow(dnu / lambda, 3 * numsteps));

            /* Filter horizontally along each row */
            for (z = 0; z < depth; z++)
            {
                for (y = 0; y < height; y++)
                {
                    for (step = 0; step < numsteps; step++)
                    {
                        ptr = volume + width*(y + height*z);
                        ptr[0] *= boundaryscale;

                        /* Filter rightwards */
                        for (x = 1; x < width; x++)
                            ptr[x] += nu*ptr[x - 1];

                        ptr[x = width - 1] *= boundaryscale;

                        /* Filter leftwards */
                        for (; x > 0; x--)
                            ptr[x - 1] += nu*ptr[x];
                    }
                }
            }

            /* Filter vertically along each column */
            for (z = 0; z < depth; z++)
            {
                for (x = 0; x < width; x++)
                {
                    for (step = 0; step < numsteps; step++)
                    {
                        ptr = volume + x + plane*z;
                        ptr[0] *= boundaryscale;

                        /* Filter downwards */
                        for (i = width; i < plane; i += width)
                            ptr[i] += nu*ptr[i - width];

                        ptr[i = plane - width] *= boundaryscale;

                        /* Filter upwards */
                        for (; i > 0; i -= width)
                            ptr[i - width] += nu*ptr[i];
                    }
                }
            }

            /* Filter along z-dimension */
            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    for (step = 0; step < numsteps; step++)
                    {
                        ptr = volume + x + width*y;
                        ptr[0] *= boundaryscale;

                        for (i = plane; i < numel; i += plane)
                            ptr[i] += nu*ptr[i - plane];

                        ptr[i = numel - plane] *= boundaryscale;

                        for (; i > 0; i -= plane)
                            ptr[i - plane] += nu*ptr[i];
                    }
                }
            }

            for (i = 0; i < numel; i++)
                volume[i] *= postscale;

            return;
        }

        void MRC_Data::MRCGeneration(char *filename)
        {
            //strcpy(dataName,filename);

            float sigma = ((float)kernelSize + 1.0) / 6.0;
            halfSize = (int)ceil(((float)kernelSize - 1.0) / 2.0);
            cout << "MRC Format" << endl;
            sizex = volData->GetSizeX();
            sizey = volData->GetSizeY();
            sizez = volData->GetSizeZ();

            float* rawData = new float[sizex * sizey * sizez];
            //ofstream myfile;
            //myfile.open("d.txt"); 
            bool first = true;
            for (int i = 0; i < sizez; i++) {
                for (int j = 0; j < sizey; j++) {
                    for (int k = 0; k < sizex; k++) {
                        float d = volData->GetDataAt(k, j, i);
                        //myfile << " " << d << endl;
                        rawData[getIndex(1, 0, k, j, i, sizex, sizey, sizez)] = d;
                        if (first)
                        {
                            first = false;
                            minScalar = maxScalar = d;
                        }
                        if (minScalar > d) {
                            minScalar = d;
                        }
                        if (maxScalar < d){
                            maxScalar = d;
                        }


                    }
                }
            }
            //myfile.close();
            float *Ix2 = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            float *Iy2 = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            float *Iz2 = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            float *IxIy = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            float *IyIz = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            float *IxIz = new float[(sizex - 2) * (sizey - 2) * (sizez - 2)];
            for (int i = 0; i < sizex - 2; i++) {
                for (int j = 0; j < sizey - 2; j++) {
                    for (int k = 0; k < sizez - 2; k++)
                    {
                        int xForward = getIndex(1, 0, i + 2, j + 1, k + 1, sizex, sizey, sizez);
                        int xBackward = getIndex(1, 0, i, j + 1, k + 1, sizex, sizey, sizez);
                        int yForward = getIndex(1, 0, i + 1, j + 2, k + 1, sizex, sizey, sizez);
                        int yBackward = getIndex(1, 0, i + 1, j, k + 1, sizex, sizey, sizez);
                        int zForward = getIndex(1, 0, i + 1, j + 1, k + 2, sizex, sizey, sizez);
                        int zBackward = getIndex(1, 0, i + 1, j + 1, k, sizex, sizey, sizez);
                        float Ix = (rawData[xForward] - rawData[xBackward]) / 2;
                        float Iy = (rawData[yForward] - rawData[yBackward]) / 2;
                        float Iz = (rawData[zForward] - rawData[zBackward]) / 2;
                        int index = getIndex(1, 0, i, j, k, sizex - 2, sizey - 2, sizez - 2);
                        Ix2[index] = Ix * Ix;
                        Iy2[index] = Iy * Iy;
                        Iz2[index] = Iz * Iz;
                        IxIy[index] = Ix * Iy;
                        IyIz[index] = Iy * Iz;
                        IxIz[index] = Ix * Iz;

                    }
                }
            }

            gaussianiir3d(Ix2, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);
            gaussianiir3d(Iy2, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);
            gaussianiir3d(Iz2, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);
            gaussianiir3d(IxIy, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);
            gaussianiir3d(IyIz, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);
            gaussianiir3d(IxIz, sizez - 2, sizey - 2, sizex - 2, sigma, iteration);

            tensors = new float[6 * (sizex - 2 - 2 * halfSize) * (sizey - 2 - 2 * halfSize) * (sizez - 2 - 2 * halfSize)];
            //myfile.open("tensorsscalars.txt");
            for (int i = 0; i < sizex - 2 - 2 * halfSize; i++) {
                for (int j = 0; j < sizey - 2 - 2 * halfSize; j++) {
                    for (int k = 0; k < sizez - 2 - 2 * halfSize; k++)
                    {
                        int oldIndex = getIndex(1, 0, i + halfSize, j + halfSize, k + halfSize, sizex - 2, sizey - 2, sizez - 2);
                        int index = getIndex(6, 0, i, j, k, sizex - 2 - 2 * halfSize, sizey - 2 - 2 * halfSize, sizez - 2 - 2 * halfSize);
                        tensors[index] = Ix2[oldIndex];
                        tensors[index + 1] = IxIy[oldIndex];
                        tensors[index + 2] = IxIz[oldIndex];
                        tensors[index + 3] = Iy2[oldIndex];
                        tensors[index + 4] = IyIz[oldIndex];
                        tensors[index + 5] = Iz2[oldIndex];
                        //myfile << std::to_string(tensors[index]) << " " << std::to_string(tensors[index+1]) << " " << std::to_string(tensors[index+2]) << " " <<  std::to_string(tensors[index+3]) << " " << std::to_string(tensors[index+4]) << " " << std::to_string(tensors[index+5]) << endl;
                    }
                }
            }
            gaussianiir3d(rawData, sizez, sizey, sizex, sigma, iteration);
            gradients = new float[3 * (sizex - 2 - 2 * halfSize) * (sizey - 2 - 2 * halfSize) * (sizez - 2 - 2 * halfSize)];
            maxGradientNorm = 0;
            for (int i = 0; i < sizex - 2 - 2 * halfSize; i++) {
                for (int j = 0; j < sizey - 2 - 2 * halfSize; j++) {
                    for (int k = 0; k < sizez - 2 - 2 * halfSize; k++)
                    {
                        int xForward = getIndex(1, 0, i + 2 + halfSize, j + 1 + halfSize, k + 1 + halfSize, sizex, sizey, sizez);
                        int xBackward = getIndex(1, 0, i + halfSize, j + 1 + halfSize, k + 1 + halfSize, sizex, sizey, sizez);
                        int yForward = getIndex(1, 0, i + 1 + halfSize, j + 2 + halfSize, k + 1 + halfSize, sizex, sizey, sizez);
                        int yBackward = getIndex(1, 0, i + 1 + halfSize, j + halfSize, k + 1 + halfSize, sizex, sizey, sizez);
                        int zForward = getIndex(1, 0, i + 1 + halfSize, j + 1 + halfSize, k + 2 + halfSize, sizex, sizey, sizez);
                        int zBackward = getIndex(1, 0, i + 1 + halfSize, j + 1 + halfSize, k + halfSize, sizex, sizey, sizez);
                        float Ix = (rawData[xForward] - rawData[xBackward]) / 2;
                        float Iy = (rawData[yForward] - rawData[yBackward]) / 2;
                        float Iz = (rawData[zForward] - rawData[zBackward]) / 2;
                        int index = getIndex(3, 0, i, j, k, sizex - 2 - 2 * halfSize, sizey - 2 - 2 * halfSize, sizez - 2 - 2 * halfSize);
                        gradients[index] = Ix;
                        gradients[index + 1] = Iy;
                        gradients[index + 2] = Iz;
                        //myfile << std::to_string(gradients[index]) << " " << std::to_string(gradients[index+1]) << " " << std::to_string(gradients[index+2]) << endl;
                        if (maxGradientNorm < getNorm(gradients + index)) maxGradientNorm = getNorm(gradients + index);
                    }
                }
            }
            scalars = new float[(sizex - 2 - 2 * halfSize) * (sizey - 2 - 2 * halfSize) * (sizez - 2 - 2 * halfSize)];
            for (int i = 0; i < sizex - 2 - 2 * halfSize; i++) {
                for (int j = 0; j < sizey - 2 - 2 * halfSize; j++) {
                    for (int k = 0; k < sizez - 2 - 2 * halfSize; k++)
                    {
                        int oldIndex = getIndex(1, 0, i + 1 + halfSize, j + 1 + halfSize, k + 1 + halfSize, sizex, sizey, sizez);
                        int index = getIndex(1, 0, i, j, k, sizex - 2 - 2 * halfSize, sizey - 2 - 2 * halfSize, sizez - 2 - 2 * halfSize);
                        scalars[index] = rawData[oldIndex];
                        //myfile << std::to_string(scalars[index]) << endl;
                    }
                }
            }
            //myfile.close();
            sizex = sizex - 2 - 2 * halfSize;
            sizey = sizey - 2 - 2 * halfSize;
            sizez = sizez - 2 - 2 * halfSize;
            rightBackTop[0] = ((float)sizex - 5) / 2;
            rightBackTop[1] = ((float)sizey - 5) / 2;
            rightBackTop[2] = ((float)sizez - 5) / 2;

            delete[] rawData;
            delete[] Ix2;
            delete[] Iy2;
            delete[] Iz2;
            delete[] IxIy;
            delete[] IyIz;
            delete[] IxIz;
            //cout << "sizex " << sizex << endl;
            cout << "Data Size: " << sizex << " * " << sizey << " * " << sizez << endl;
        }




        void MRC_Data::getScalarGP(int index, float *scalar)
        {
            *scalar = scalars[index];
        }

        void MRC_Data::getGradientGP(int index, float gradient[3])
        {
            gradient[0] = gradients[index];
            gradient[1] = gradients[index + 1];
            gradient[2] = gradients[index + 2];
        }

        void MRC_Data::getTensorGP(int index, float tensor[6])
        {
            tensor[0] = tensors[index];
            tensor[1] = tensors[index + 1];
            tensor[2] = tensors[index + 2];
            tensor[3] = tensors[index + 3];
            tensor[4] = tensors[index + 4];
            tensor[5] = tensors[index + 5];
        }

        void MRC_Data::getScalar(float position[3], float *scalar)
        {
            trilinear_f(sizex, sizey, sizez, scalars, position, scalar);
        }

        void MRC_Data::getGradient(float position[3], float gradient[3])
        {
            trilinear_3f(sizex, sizey, sizez, gradients, position, gradient);
        }

        void MRC_Data::getTensor(float position[3], float tensor[6])
        {
            trilinear_6f(sizex, sizey, sizez, tensors, position, tensor);
        }

        void MRC_Data::getScalarCubic(float position[3], float *scalar)
        {
            tricubic_f(sizex, sizey, sizez, scalars, position, scalar);
        }

        void MRC_Data::getGradientCubic(float position[3], float gradient[3])
        {
            tricubic_3f(sizex, sizey, sizez, gradients, position, gradient);
        }

        void MRC_Data::getTensorCubic(float position[3], float tensor[6])
        {
            tricubic_6f(sizex, sizey, sizez, tensors, position, tensor);
        }

        void MRC_Data::getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3])
        {
            tricubic_3f_table(sizex, sizey, sizez, edgeTable, axis, i + 2, j + 2, k + 2, index, gradients, gradient);
        }

        void MRC_Data::getGradientBicubicTable(int axis, int ii, int jj, int kk, float sGradients[subdNum][subdNum][3])
        {
            float validData[16][3];
            int count = 0;
            switch (axis)
            {
            case 1:
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        int idx = getIndex(3, 0, ii + 2, jj + i + 1, kk + j + 1, sizex, sizey, sizez);
                        validData[count][0] = gradients[idx];
                        validData[count][1] = gradients[idx + 1];
                        validData[count][2] = gradients[idx + 2];
                        count++;
                    }
                }
                break;
            case 2:
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        int idx = getIndex(3, 0, ii + j + 1, jj + 2, kk + i + 1, sizex, sizey, sizez);
                        validData[count][0] = gradients[idx];
                        validData[count][1] = gradients[idx + 1];
                        validData[count][2] = gradients[idx + 2];
                        count++;
                    }
                }
                break;
            case 3:
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        int idx = getIndex(3, 0, ii + i + 1, jj + j + 1, kk + 2, sizex, sizey, sizez);
                        validData[count][0] = gradients[idx];
                        validData[count][1] = gradients[idx + 1];
                        validData[count][2] = gradients[idx + 2];
                        count++;
                    }
                }
                break;
            }
            count = 0;
            for (int si = 0; si < subdNum; si++)
            {
                for (int sj = 0; sj < subdNum; sj++)
                {
                    int tableIndex = 16 * count;
                    count++;
                    sGradients[si][sj][0] = 0;
                    sGradients[si][sj][1] = 0;
                    sGradients[si][sj][2] = 0;
                    for (int i = 0; i < 16; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            sGradients[si][sj][j] += faceTable[tableIndex + i] * validData[i][j];
                        }
                    }
                }
            }
        }

        void MRC_Data::getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6])
        {
            tricubic_6f_table(sizex, sizey, sizez, edgeTable, axis, i + 2, j + 2, k + 2, index, tensors, tensor);
        }

        void MRC_Data::buildGrid()
        {
            clock_t timeStart, timeEnd;
            timeStart = clock();
            cout << "Build Grid Begin:" << endl;
            int maxSize = sizex;
            halfDataSize = rightBackTop[0];
            if (maxSize < sizey)
            {
                maxSize = sizey;
                halfDataSize = rightBackTop[1];
            }
            if (maxSize < sizez)
            {
                maxSize = sizez;
                halfDataSize = rightBackTop[2];
            }

            gridx = sizex - 4;
            gridy = sizey - 4;
            gridz = sizez - 4;
            gridSize = 1;
            gridPoints = new float[3 * gridx * gridy * gridz];
            for (int i = 0; i < gridx; i++) {
                for (int j = 0; j < gridy; j++) {
                    for (int k = 0; k < gridz; k++)
                    {
                        int index = getIndex(3, 0, i, j, k, gridx, gridy, gridz);
                        gridPoints[index] = i + 2.0f;
                        gridPoints[index + 1] = j + 2.0f;
                        gridPoints[index + 2] = k + 2.0f;
                    }
                }
            }

            maxGrid = gridx;
            if (maxGrid < gridy) maxGrid = gridy;
            if (maxGrid < gridz) maxGrid = gridz;
            maxGridIndex = gridx * gridy * gridz;
            cout << "Grid Size: " << gridx << " * " << gridy << " * " << gridz << endl;
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;

            timeStart = clock();
            cout << "Calculate Edge Lookup Table Begin:" << endl;
            edgeTable = new float[4100];
            buildEdgeTable(edgeTable, 1, 2, 0, 1);
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;

            timeStart = clock();
            cout << "Calculate Face Lookup Table Begin:" << endl;
            faceTable = new float[16 * subdNum*subdNum];
            buildFaceTable(faceTable, subdNum);
            timeEnd = clock();
            cout << "Spend Time: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;
            cout << "Done!" << endl;
            cout << "****************************************************************" << endl;
        }

        class Volume {
        public:
            Volume(int x, int y, int z);
            Volume(int x, int y, int z, float val);
            Volume(int x, int y, int z, int offx, int offy, int offz, Volume * vol);
            Volume(const Volume& obj);
            ~Volume();

            py::list getExtremalBonds1();
            py::list getExtremalBonds2();
            py::list getMinCurveBonds1();
            py::list getMinCurveBonds2();
            py::list getSaddleCurveBonds1();
            py::list getSaddleCurveBonds2();
            py::list hideSegments();
            py::list getDisplayQuads();
            py::list getExtremalNormals();
            py::list getQuadTypes();
            py::list hideSurfaces();
            py::list getVertexPos();

            py::list getMinCurvePos();
            py::list getMaxCurvePos();
            py::list getSaddleCurvePos();
            py::list getMinSurfacePos();
            py::list getMaxSurfacePos();

            py::list getQuadSaliencies();
            py::list getQuadIntensities();
            py::list getQuadEigenvalues();

            py::list getMaxPointSaliencies();
            py::list getMaxPointIntensities();
            py::list getMaxPointEigenvalues();
            py::list getMinPointSaliencies();
            py::list getMinPointIntensities();
            py::list getMinPointEigenvalues();

            py::list getSaddlePointSaliencies();
            py::list getSaddlePointIntensities();
            py::list getSaddlePointEigenvalues();

            py::list getMaxCurveEigenvalues0();

            py::list getMaxCurveEigenvalues1();

            py::list getMinCurveEigenvalues0();

            py::list getMinCurveEigenvalues1();

            py::list getSaddleCurveEigenvalues0();

            py::list getSaddleCurveEigenvalues1();

            py::list getMaxPointEigenvalues0();

            py::list getMaxPointEigenvalues1();

            py::list getMinPointEigenvalues0();

            py::list getMinPointEigenvalues1();

            py::list getSaddlePointEigenvalues0();

            py::list getSaddlePointEigenvalues1();

            py::list getQuadEigenvalues0();

            py::list getQuadEigenvalues1();
            py::list getMaxCurveEigenvectors();
            py::list getMinCurveEigenvectors();
            py::list getSaddleCurveEigenvectors();

            py::list getQuadEigenvectors();
            py::list getMaxPointEigenvectors();
            py::list getMinPointEigenvectors();
            py::list getSaddlePointEigenvectors();

            float getSpacingX();
            float getSpacingY();
            float getSpacingZ();
            float getOriginX();
            float getOriginY();
            float getOriginZ();
            int getSizeX();
            int getSizeY();
            int getSizeZ();
            int getIndex(int x, int y, int z);
            double getDataAt(int x, int y, int z);
            double getDataAt(int index);
            void setSpacing(float spx, float spy, float spz);
            void setOrigin(float orgX, float orgY, float orgZ);
            void setDataAt(int x, int y, int z, double d);
            void setDataAt(int index, double d);

            void setMinMaxIntensitiesEigenvalues(float minI, float maxI, float minE, float maxE);
            float getMinI();
            float getMaxI();
            float getMinE();
            float getMaxE();

            Volume * getPseudoDensity();
            Volume * getDistanceField(int rad, float randf);
            int getNonZeroVoxelCount();
            void print();
            void subtract(Volume * vol);
            void pad(int padBy, double padValue);
            void applyMask(Volume * maskVol, double maskValue, bool keepMaskValue);
            void setExtremalWidth(int width);
            void setExtremalHeight(int height);
            void setExtremalSlices(int slices);
            void setExtremalGridResolution(int resolution);
            void setExtremalSlicesThickness(int thickness);
            void setExtremalIsovalue(double isovalue);
            void setMaxPointOn(bool maxOn);
            void setMinPointOn(bool minOn);
            void setSaddlePointOn(bool saddleOn);



            double getMin();
            double getMax();
            double getMaxValuePosition(int& maxX, int& maxY, int& maxZ);
            double getLocalMax(int x, int y, int z, int radius);
            double getLocalMin(int x, int y, int z, int radius);

            float getMean(); // Returns the mean value of all the voxels
            float getEdgeMean(); // Returns the mean value of all the surface voxels but no interior voxels
            float getStdDev(); // Returns the population standard deviation of the values at all the voxels
            Vector3DFloat getCenterOfMass(); // Returns the center of mass of the image in pixels (not angstroms)
            float* getArrayCopy(int padX = 0, int padY = 0, int padZ = 0, float padValue = 0);

            void fill(double val);
            int isBertrandBorder(int ox, int oy, int oz, int dir);
            int isBertrandEndPoint(int ox, int oy, int oz);
            int isHelix(int ox, int oy, int oz);
            int isSheet(int ox, int oy, int oz);
            Volume * getSheets(int minSize);
            Volume * getHelices(int minSize);
            int isEndPoint(int ox, int oy, int oz);
            int getNumNeighbor6(int ox, int oy, int oz);
            int testIsSheetEnd(int ox, int oy, int oz);
            int isNoiseSheetEnd(int ox, int oy, int oz);
            int isInternal(int ox, int oy, int oz);
            int isInternal2(int ox, int oy, int oz);
            int hasIsolatedFace(int ox, int oy, int oz);
            int hasIsolatedEdge(int ox, int oy, int oz);
            int countFace(int ox, int oy, int oz, int m);
            int hasCell(int ox, int oy, int oz);
            Volume * markCellFace();
            Volume * markFaceEdge();
            int hasCompleteSheet(int ox, int oy, int oz, Volume * fvol);
            int hasCompleteSheet(int ox, int oy, int oz);
            int hasCompleteSheetSlow(int ox, int oy, int oz);
            int hasCompleteHelix(int ox, int oy, int oz);
            int hasCompleteHelix(int ox, int oy, int oz, Volume * fvol);
            int isHelixEnd(int ox, int oy, int oz, Volume * nvol);
            int isFeature18(int ox, int oy, int oz);
            int isEdgeEnd(int ox, int oy, int oz);
            int isFaceEnd(int ox, int oy, int oz);
            int isNoise(int ox, int oy, int oz, int noise);
            int isNoiseHelixEnd(int ox, int oy, int oz);
            int isHelixEnd(int ox, int oy, int oz);
            int isSheetEnd(int ox, int oy, int oz, Volume * nvol);
            int getNumFaces(int ox, int oy, int oz);
            int getNumCells(int ox, int oy, int oz);
            int getNumIsolatedEdges(int ox, int oy, int oz);
            int getNumIsolatedFaces(int ox, int oy, int oz);
            int isFeatureFace2(int ox, int oy, int oz);
            int isFeatureFace(int ox, int oy, int oz);
            int hasFeatureFace(int ox, int oy, int oz);
            int isSheetEnd(int ox, int oy, int oz);
            int isSimple(int ox, int oy, int oz);
            int isPiercable(int ox, int oy, int oz);
            int isSimple2(int v[3][3][3]);
            int getNumPotComplex3(int ox, int oy, int oz);
            int getNumPotComplex4(int ox, int oy, int oz);
            int getNumPotComplex(int ox, int oy, int oz);
            int getNumPotComplex2(int ox, int oy, int oz);
            int getNumNeighbor(int ox, int oy, int oz);
            void setScoreNeighbor(GridQueue* queue);
            int components6(int vox[3][3][3]);
            int components26(int vox[3][3][3]);
            int countExt(double vox[3][3][3]);
            int countInt(double vox[3][3][3]);
            int countIntEuler(int ox, int oy, int oz);
            void erodeNoTopo(float thr, int wid);
            void erodeTopo(float thr, int wid);
            void erode2(float thr, int wid);
            void erodeShapeTopo(float thr, int wid);
            void erodeAtom(float thr, int wid, Volume* avol);
            void curveSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol);
            void curveSkeleton(float thr, Volume* svol);
            void curveSkeleton2D(float thr, Volume* svol);
            void skeleton(float thr, int off);
            void skeleton2(float thr, int off);
            void pointSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol, Volume* hvol);
            void skeleton(float thr, Volume* svol, Volume* hvol);
            void erodeHelix();
            void erodeHelix(int disthr);
            int erodeSheet();
            int erodeSheet(int disthr);
            void erodeSheetOld(int disthr);
            void addNoise(float thr, float pos);
            void sequentialSkeleton(float thr, int type, int noise);
            void dumbsurfaceSkeleton(float thr);
            void surfaceSkeleton(Volume* grayvol, float lowthr, float highthr);
            void surfaceSkeleton(float thr);
            void surfaceSkeleton(float thr, Volume* svol);
            void extremalCurveSkeleton(int kernelSize, Volume* svol, int width, int height, int slices, int resolution, float thickness, float isovalue, int kernelSizeGone);
            void surfaceSkeletonOld(float thr);
            void surfaceSkeletonPres(float thr, Volume * preserve);
            void bertrandSurfaceSkeleton2(float thr);
            void bertrandSurfaceSkeleton(float thr);
            void palagyiSurfaceSkeleton(float thr);
            void threshold(double thr);
            void threshold(double thr, int out, int in);
            void threshold(double thr, int out, int in, int boundary);
            void threshold(double thr, int out, int in, int boundary, bool markBoundary);
            void threshold2(double thr, int out, int in);
            void smooth(float alpha);
            void normalize(double min, double max);
            void normalize(double min, double max, double thresh, double ithresh);
            Volume * getDataRange(int x, int y, int z, int radius);
            double getInterpDataAt(double x, double y, double z);
            void rotateX(double a);
            void toMathematicaFile(char* fname);
            void toMathematicaFile(char* fname, int lx, int hx, int ly, int hy, int lz, int hz);
            void toOFFCells(char* fname);
            void toOFFCells2(char* fname);
            void toOFFCells2(char* fname, float thr);
            void toOFFCells(char* fname, float thr);
            void segment(float threshold, Volume* lowvol, Volume* highvol, char* mrcfile);
            void segment(float threshold, Volume* vol, int maxDis, char* mrcfile);
            void writeSegmentation(float threshold, Volume* segvol, char* txtfile, char* mrcfile);
            void floodFill(float thr);
            void reduceComponent(int size);
            void reduceComponent2(int num);
            void floodFillPQR(int offset);
            void writeDistances(char* fname, int maxDis);
            void toPQRFile(char* fname, float spc, float minx, float miny, float minz, int padding);
            void toMRCFile(char* fname);
            void buildHistogram(int binCount);
            int getHistogramBinValue(int binIx);
            int getMaxSize();

            void setOrigSize(int newSize);
            int getOrigSize();
            void setExtremalParams(float curve, float point, float surface, float minG, float maxG, float eigen);
            void setCurrentGenData(General_Data * genData);
            void setExtremalChecks(bool saliency, bool intensity, bool eigenvalue);
            py::list getExtremalMaxPoints();
            py::list getExtremalMinPoints();
            py::list getExtremalSaddlePoints();

            void setXMin(int min);
            void setYMin(int min);
            void setZMin(int min);
            void setXMax(int max);
            void setYMax(int max);
            void setZMax(int max);

            py::list getMaxSaliencies();
            py::list getMaxIntensities();
            py::list getMaxEigenvalues();

            py::list getMinSaliencies();
            py::list getMinIntensities();
            py::list getMinEigenvalues();

            py::list getSaddleSaliencies();
            py::list getSaddleIntensities();
            py::list getSaddleEigenvalues();

            py::list getMaxHashes();
            py::list getMinHashes();
            py::list getSaddleHashes();
            py::list getQuadIndices();

            std::vector<string> maxCurveBonds;
            std::vector<string> minCurveBonds;
            std::vector<string> saddleCurveBonds;


            std::vector<float> extremalPoints;
            std::vector<int> extremalBonds1;
            std::vector<int> extremalBonds2;

            std::vector<int> minCurve1;
            std::vector<int> minCurve2;

            std::vector<int> saddleCurve1;
            std::vector<int> saddleCurve2;

            std::vector<string> quadDisplays;
            //std::vector<int> quadDisplays;
            std::vector<float> quadNormals;
            std::vector<int> quadTypes;
            std::vector<float> quadSaliencies;
            std::vector<float> quadIntensities;
            std::vector<float> quadEigenvalues;
            std::vector<float> quadEigenvalues0;
            std::vector<float> quadEigenvalues1;

            std::vector<int> quadIndices;

            std::vector<float> maxCurveSaliencies;
            std::vector<float> maxCurveIntensities;
            std::vector<float> maxCurveEigenvalues;

            std::vector<float> maxCurveEigenvalues0;
            std::vector<float> maxCurveEigenvalues1;
            std::vector<float> maxCurveEigenvectors;

            std::vector<float> quadEigenvectors;
            std::vector<float> maxPointEigenvectors;
            std::vector<float> minPointEigenvectors;
            std::vector<float> saddlePointEigenvectors;

            std::vector<float> minCurveSaliencies;
            std::vector<float> minCurveIntensities;
            std::vector<float> minCurveEigenvalues;

            std::vector<float> minCurveEigenvectors;
            std::vector<float> saddleCurveEigenvectors;

            std::vector<float> minCurveEigenvalues0;
            std::vector<float> minCurveEigenvalues1;


            std::vector<float> saddleCurveSaliencies;
            std::vector<float> saddleCurveIntensities;
            std::vector<float> saddleCurveEigenvalues;

            std::vector<float> saddleCurveEigenvalues0;
            std::vector<float> saddleCurveEigenvalues1;

            std::vector<float> minCurvePos;
            std::vector<float> maxCurvePos;
            std::vector<float> saddleCurvePos;
            std::vector<float> minSurfacePos;
            std::vector<float> maxSurfacePos;
            std::vector<float> vertexPos;

            std::vector<float> minPointSaliencies;
            std::vector<float> minPointIntensities;
            std::vector<float> minPointEigenvalues;

            std::vector<float> minPointEigenvalues0;
            std::vector<float> minPointEigenvalues1;

            std::vector<float> maxPointSaliencies;
            std::vector<float> maxPointIntensities;
            std::vector<float> maxPointEigenvalues;

            std::vector<float> maxPointEigenvalues0;
            std::vector<float> maxPointEigenvalues1;

            std::vector<float> saddlePointSaliencies;
            std::vector<float> saddlePointIntensities;
            std::vector<float> saddlePointEigenvalues;

            std::vector<float> saddlePointEigenvalues0;
            std::vector<float> saddlePointEigenvalues1;

            int extremalWidth, extremalHeight, extremalSlices, extremalResolution, extremalThickness, extremalIsovalue;
            double thickness;

            float curveRatio, pointRatio, surfaceRatio, minGeo, maxGeo, eigenValue;
            bool saliencyCheck, intensityCheck, eigenvalueCheck;
            bool minOn, maxOn, saddleOn;
            General_Data * currentExtremalData;
            int xMin, yMin, zMin, xMax, yMax, zMax;
            float minEigen, maxEigen, minIntensity, maxIntensity;

        private:
            VolumeData * getVolumeData();
            std::vector<int> histogram;



            /**
            float change;
            int maxEdgeIndex;
            int *totalEdgePoints;
            float *grid_2DGradMag;
            bool *storeGrid_2DGradMag;
            //Edge *edges;
            float globalVec[3];
            float thresh;
            bool storeGridShowPos;
            bool storeGridScalar;
            bool storeGridVector;
            bool storeGridGradient;
            bool storeEdgePoint;
            bool storeFacePoint;
            bool storeCellPoint;
            bool storeDispCell;
            **/

        private:
            VolumeData * volData;
            int origSize;
            //View *view3D;
        };

        Volume::Volume(const Volume& obj) {
            volData = new VolumeData(*(obj.volData));
            histogram = obj.histogram;
        }

        Volume::Volume(int x, int y, int z) {
            volData = new VolumeData(x, y, z);
        }

        Volume::Volume(int x, int y, int z, float val) {
            volData = new VolumeData(x, y, z, val);
        }

        Volume::Volume(int x, int y, int z, int offx, int offy, int offz, Volume * vol) {
            volData = new VolumeData(x, y, z, offx, offy, offz, vol->getVolumeData());
        }


        Volume::~Volume()
        {
            delete volData;
        }

        int Volume::getSizeX() {
            return volData->GetSizeX();
        }

        int Volume::getSizeY() {
            return volData->GetSizeY();
        }

        int Volume::getSizeZ() {
            return volData->GetSizeZ();
        }

        int Volume::getOrigSize() {
            return origSize;
        }

        void Volume::setXMin(int min) {
            xMin = min;
        }

        void Volume::setYMin(int min) {
            yMin = min;
        }

        void Volume::setZMin(int min) {
            zMin = min;
        }

        void Volume::setXMax(int max) {
            xMax = max;
        }

        void Volume::setYMax(int max) {
            yMax = max;
        }



        void Volume::setZMax(int max){
            zMax = max;
        }

        void Volume::setMaxPointOn(bool pointMaxOn) {
            maxOn = pointMaxOn;
        }

        void Volume::setSaddlePointOn(bool pointSaddleOn) {
            saddleOn = pointSaddleOn;
        }

        void Volume::setMinMaxIntensitiesEigenvalues(float minI, float maxI, float minE, float maxE) {
            minIntensity = minI;
            maxIntensity = maxI;

            minEigen = minE;
            maxEigen = maxE;
        }

        float Volume::getMinI() {
            return minIntensity;
        }


        float Volume::getMaxI() {
            return maxIntensity;
        }

        float Volume::getMinE() {
            return minEigen;
        }

        float Volume::getMaxE() {
            return maxEigen;
        }



        void Volume::setMinPointOn(bool pointMinOn) {
            minOn = pointMinOn;
        }

        void Volume::setExtremalChecks(bool saliency, bool intensity, bool eigenvalue) {
            saliencyCheck = saliency;
            intensityCheck = intensity;
            eigenvalueCheck = eigenvalue;
        }

        void Volume::setExtremalWidth(int width) {
            extremalWidth = width;
        }

        void Volume::setExtremalHeight(int height) {
            extremalHeight = height;
        }

        void Volume::setExtremalSlices(int slices) {
            extremalSlices = slices;
        }

        void Volume::setExtremalGridResolution(int resolution) {
            extremalResolution = resolution;
        }

        void Volume::setExtremalSlicesThickness(int thickness) {
            extremalThickness = thickness;
        }

        void Volume::setExtremalIsovalue(double isovalue) {
            extremalIsovalue = isovalue;
        }

        void Volume::setCurrentGenData(General_Data * genData) {
            currentExtremalData = genData;
        }

        void Volume::setExtremalParams(float curve, float point, float surface, float minG, float maxG, float eigen) {
            curveRatio = curve;
            pointRatio = point;
            surfaceRatio = surface;
            minGeo = minG;
            maxGeo = maxG;
            eigenValue = eigen;
        }

        void Volume::setOrigSize(int newSize) {
            origSize = newSize;
        }

        int Volume::getIndex(int x, int y, int z) {
            return volData->GetIndex(x, y, z);
        }

        int Volume::getMaxSize() {
            int sizex = getSizeX();
            int sizey = getSizeY();
            int sizez = getSizeZ();
            int maxSize = sizex;
            if (sizey > maxSize) {
                maxSize = sizey;
            }
            if (sizez > maxSize) {
                maxSize = sizex;
            }
            return maxSize;
        }

        void Volume::setDataAt(int x, int y, int z, double d) {
            volData->SetDataAt(x, y, z, (float)d);
        }

        void Volume::setDataAt(int index, double d) {
            volData->SetDataAt(index, (float)d);
        }

        double Volume::getDataAt(int x, int y, int z)
        {
            return volData->GetDataAt(x, y, z);
        }

        double Volume::getDataAt(int index) {
            return volData->GetDataAt(index);
        }

        VolumeData * Volume::getVolumeData() {
            return volData;
        }

        void Volume::setSpacing(float spx, float spy, float spz) {
            volData->SetSpacing(spx, spy, spz);
        }

        void Volume::setOrigin(float orgX, float orgY, float orgZ) {
            volData->SetOrigin(orgX, orgY, orgZ);
        }

        float Volume::getSpacingX() {
            return volData->GetSpacingX();
        }

        float Volume::getSpacingY() {
            return volData->GetSpacingY();
        }

        float Volume::getSpacingZ() {
            return volData->GetSpacingZ();
        }

        float Volume::getOriginX() {
            return volData->GetOriginX();
        }

        float Volume::getOriginY() {
            return volData->GetOriginY();
        }

        float Volume::getOriginZ() {
            return volData->GetOriginZ();
        }

        Volume * Volume::getPseudoDensity() {
            // This function assumes the volume is binary (1/0), and builds a pseudo density volume from the 1 voxels
            // First: assign a density value at each point
            int i, j, k;
            Volume * res = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);
            int size = getSizeX() * getSizeY() * getSizeZ();
            srand(123);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if (res->getDataAt(i, j, k) > 0) {
                    int ct = 0;
                    for (int m = 0; m < 6; m++) {
                        if (res->getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) > 0) {
                            ct++;
                        }
                    }
                    res->setDataAt(i, j, k, (k / (float)getSizeZ())*(k / (float)getSizeZ()));
                    if (ct > 2) {
                        //res->setDataAt( i,j,k, rand() / (double) RAND_MAX / 2.0f ) ;
                    }
                    else {
                        //res->setDataAt( i,j,k, rand() / (double) RAND_MAX ) ;
                    }
                }
            }

            /* Next, smooth
            for ( i = 0 ; i < 20 ; i ++ )
            {
            printf("Smoothing round %d\n", i) ;
            res->smooth( 0.5f ) ;
            }
            */

            Volume * tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, res);
            float d, ad, ct, temp;
            for (int it = 0; it < 3; it++)
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if ((d = (float)tvol->getDataAt(i, j, k)) > 0) {
                    ad = 0; ct = 0;
                    for (int m = 0; m < 6; m++) {
                        if ((temp = (float)tvol->getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2])) > 0) {
                            ad += temp;
                            ct++;
                        }
                    }
                    if (ct > 0) {
                        res->setDataAt(i, j, k, (d + ad / ct) / 2);
                    }
                }
            }

            delete tvol;
            tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, res);
            for (i = 0; i < 40; i++)
            {
                printf("Smoothing round %d\n", i);
                res->smooth(0.5f);
                continue;

                //for ( j = 0 ; j < size ; j ++ )
                //{
                //	if ( tvol->getDataAt( j ) > 0 )
                //	{
                //		res->setDataAt( j, tvol->getDataAt( j ) ) ;
                //	}
                //	
                //}

            }


            return res;
        }


        Volume * Volume::getDistanceField(int rad, float randf) {
            // This function assumes the volume is binary (1/0), and builds a pseudo density volume from the 1 voxels
            // rad is the radius of each distance function (e.g., center pixel gets 1, pixels at rad from the center gets 0)
            // randf is how much noise you want to add. this means the center pixel will maximally have value 1+randf.

            // First: assign a density value at each point
            int i, j, k;
            Volume * res = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);
            srand(123);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if (getDataAt(i, j, k) > 0) {
                    float mag = 1 + randf * (float)rand() / (float)RAND_MAX;
                    int lx = max(0, i - rad);
                    int ly = max(0, j - rad);
                    int lz = max(0, k - rad);
                    int hx = min(getSizeX() - 1, i + rad);
                    int hy = min(getSizeY() - 1, j + rad);
                    int hz = min(getSizeZ() - 1, k + rad);
                    int x, y, z;
                    for (x = lx; x <= hx; x++)
                    for (y = ly; y <= hy; y++)
                    for (z = lz; z <= hz; z++) {
                        float val = 1 - (float)sqrt((double)((x - i)*(x - i) + (y - j)*(y - j) + (z - k)*(z - k))) / (float)rad;
                        val *= mag;
                        if (res->getDataAt(x, y, z) < val) {
                            res->setDataAt(x, y, z, val);
                        }
                    }
                }
            }

            /* Next, smooth */
            for (i = 0; i < 2; i++)
            {
                printf("Smoothing round %d\n", i);
                res->smooth(0.5f);
            }


            return res;
        }


        int Volume::getNonZeroVoxelCount() {
            int count = 0;
            for (int x = 0; x < getSizeX(); x++){
                for (int y = 0; y < getSizeY(); y++){
                    for (int z = 0; z < getSizeZ(); z++){
                        if (this->getDataAt(x, y, z) > 0.0) {
                            count++;
                        }
                    }
                }
            }
            return count;
        }
        void Volume::print() {
            for (int x = 0; x < getSizeX(); x++) {
                printf("{ ");
                for (int y = 0; y < getSizeY(); y++) {
                    printf("{ ");
                    for (int z = 0; z < getSizeZ(); z++) {
                        printf("%f, ", getDataAt(x, y, z));
                    }
                    printf("} ");
                }
                printf("} ");
            }
            printf("\n");
        }


        void Volume::subtract(Volume* vol) {
            int i, j, k;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if (getDataAt(i, j, k) > 0) {
                    if (vol->getDataAt(i, j, k) > 0) {
                        setDataAt(i, j, k, 0);
                    }
                }
            }

        }

        void Volume::pad(int padBy, double padValue) {
            volData->Pad(padBy, padValue);
        }

        void Volume::applyMask(Volume * maskVol, double maskValue, bool keepMaskValue) {
            for (int x = 0; x < maskVol->getSizeX(); x++) {
                for (int y = 0; y < maskVol->getSizeY(); y++) {
                    for (int z = 0; z < maskVol->getSizeZ(); z++) {
                        if (((maskVol->getDataAt(x, y, z) == maskValue) && !keepMaskValue) ||
                            ((maskVol->getDataAt(x, y, z) != maskValue) && keepMaskValue)) {
                            setDataAt(x, y, z, 0);
                        }
                    }
                }
            }
        }

        double Volume::getMin() {
            int size = volData->GetMaxIndex();
            double rvalue = volData->GetDataAt(0);
            for (int i = 1; i < size; i++) {
                float val = volData->GetDataAt(i);
                if (rvalue > val) {
                    rvalue = val;
                }
            }
            return rvalue;
        }

        double Volume::getMax() {
            int size = volData->GetMaxIndex();
            double rvalue = volData->GetDataAt(0);
            for (int i = 1; i < size; i++) {
                float val = volData->GetDataAt(i);
                if (rvalue < val) {
                    rvalue = val;
                }
            }
            return rvalue;
        }

        double Volume::getMaxValuePosition(int& maxX, int& maxY, int& maxZ) {
            double maxVal = getDataAt(0, 0, 0);
            maxX = 0; maxY = 0; maxZ = 0;
            double data;

            for (int x = 0; x < getSizeX(); x++) {
                for (int y = 0; y < getSizeY(); y++) {
                    for (int z = 0; z < getSizeZ(); z++) {
                        data = getDataAt(x, y, z);
                        if (data > maxVal) {
                            maxVal = data;
                            maxX = x; maxY = y; maxZ = z;
                        }
                    }
                }
            }
            return maxVal;
        }

        double Volume::getLocalMax(int x, int y, int z, int radius) {
            double mx = getDataAt(x, y, z);
            for (int xx = x - radius; xx <= x + radius; xx++) {
                for (int yy = y - radius; yy <= y + radius; yy++) {
                    for (int zz = z - radius; zz <= z + radius; zz++) {
                        mx = max(mx, getDataAt(xx, yy, zz));
                    }
                }
            }
            return mx;
        }

        double Volume::getLocalMin(int x, int y, int z, int radius) {
            double mn = getDataAt(x, y, z);
            for (int xx = x - radius; xx <= x + radius; xx++) {
                for (int yy = y - radius; yy <= y + radius; yy++) {
                    for (int zz = z - radius; zz <= z + radius; zz++) {
                        mn = min(mn, getDataAt(xx, yy, zz));
                    }
                }
            }
            return mn;
        }

        void Volume::fill(double val)
        {
            for (int x = 0; x < getSizeX(); x++) {
                for (int y = 0; y < getSizeY(); y++) {
                    for (int z = 0; z < getSizeZ(); z++) {
                        setDataAt(x, y, z, val);
                    }
                }
            }
        }

        int Volume::isBertrandBorder(int ox, int oy, int oz, int dir) {
            int nx = ox + neighbor6[dir][0];
            int ny = oy + neighbor6[dir][1];
            int nz = oz + neighbor6[dir][2];
            if (getDataAt(nx, ny, nz) < 0) {
                return 1;
            }
            return 0;
        }

        int Volume::isBertrandEndPoint(int ox, int oy, int oz) {
            double vox[3][3][3];

            int i, j, k;
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                double tval = getDataAt(ox + i, oy + j, oz + k);
                vox[i + 1][j + 1][k + 1] = tval;
            }

            // Count #X^x_6
            int xx6 = 0;
            for (i = 0; i < 6; i++) {
                if (vox[neighbor6[i][0] + 1][neighbor6[i][1] + 1][neighbor6[i][2] + 1] >= 0) {
                    xx6++;
                }
            }

            // Now, count A26 and /A26
            int a26 = 0, na26 = 0;
            for (i = 0; i < 3; i += 2)
            for (j = 0; j < 3; j += 2)
            for (k = 0; k < 3; k += 2) {
                if (vox[1][j][k] < 0) {
                    continue;
                }
                if (vox[i][1][k] < 0) {
                    continue;
                }
                if (vox[i][j][1] < 0) {
                    continue;
                }
                if (vox[i][1][1] < 0) {
                    continue;
                }
                if (vox[1][j][1] < 0) {
                    continue;
                }
                if (vox[1][1][k] < 0) {
                    continue;
                }
                if (vox[i][j][k] >= 0) {
                    a26++;
                }
                else {
                    na26++;
                }
            }

            if ((na26 == 0) && (a26 != 0) && (xx6 <= a26 + 2))  {
                return 0;
            }
            if (isFeatureFace(ox, oy, oz)) {
                //return 0 ;
            }
            return 1;
        }

        int Volume::isHelix(int ox, int oy, int oz) {
            int cn = 12;
            int nx, ny, nz;
            int i, j, k;

            double vox[3][3][3];
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                vox[i + 1][j + 1][k + 1] = getDataAt(ox + i, oy + j, oz + k);
            }

            for (i = 0; i < 12; i++) {
                for (j = 0; j < 4; j++) {
                    nx = sheetNeighbor[i][j][0] + 1;
                    ny = sheetNeighbor[i][j][1] + 1;
                    nz = sheetNeighbor[i][j][2] + 1;

                    if (vox[nx][ny][nz] <= 0) {
                        cn--;
                        break;
                    }
                }
            }

            if (cn >= 1) {
                return 0;
            }
            else {
                return 1;
            }
        }

        int Volume::isSheet(int ox, int oy, int oz) {
            int cn = 12;
            int nx, ny, nz;

            for (int i = 0; i < 12; i++) {
                for (int j = 0; j < 4; j++) {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) <= 0) {
                        cn--;
                        break;
                    }
                }
            }
            return (cn >= 3);
        }

        Volume * Volume::getSheets(int minSize) {
            int i, j, k;

            //Initialize volume
            printf("Initialize volume at %d %d %d\n", getSizeX(), getSizeY(), getSizeZ());
            Volume* svol = new Volume(getSizeX(), getSizeY(), getSizeZ());

            //Initialize cluster counters
            int sheets[MAX_SHEETS];
            for (i = 0; i < MAX_SHEETS; i++) {
                sheets[i] = 0;
            }
            int totSheets = 1;

            //Start clustering
            printf("Start clustering...\n");
            int ox, oy, oz;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if (getDataAt(i, j, k) <= 0 || svol->getDataAt(i, j, k) != 0) {
                    // Not a data point or has been visited
                    continue;
                }
                if (!isSheet(i, j, k)) {
                    // Not a sheet point
                    continue;
                }

                //Initialize queue
                int numNodes = 1;
                svol->setDataAt(i, j, k, totSheets);
                GridQueue* queue = new GridQueue();
                queue->pushQueue(i, j, k);
                while (queue->popQueue(ox, oy, oz)) {
                    // Test if neighbors satisfy sheet condition
                    if (isSheet(ox, oy, oz)) {
                        for (int m = 0; m < 6; m++) {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];

                            if (getDataAt(nx, ny, nz) > 0 && svol->getDataAt(nx, ny, nz) == 0) {
                                svol->setDataAt(nx, ny, nz, totSheets);
                                queue->pushQueue(nx, ny, nz);
                                numNodes++;
                            }
                        }
                    }
                }

                delete queue;
                if (numNodes > 0) {
                    //	printf("Sheet %d contain %d nodes.\n", totSheets, numNodes) ;
                    sheets[totSheets] = numNodes;
                    totSheets++;
                }
            }

            // Removing clusters less than minSize
            printf("Removing small clusters.\n");
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                int cnt = (int)svol->getDataAt(i, j, k);
                if (cnt > 0 && sheets[cnt] < minSize) {
                    svol->setDataAt(i, j, k, -1);
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            svol->threshold(0.1, 0, 1);

            return svol;
        }

        Volume * Volume::getHelices(int minSize) {
            printf("Segmenting helices from eroded volume.\n");
            int i, j, k;

            //Initialize volume
            printf("Initialize volume at %d %d %d\n", getSizeX(), getSizeY(), getSizeZ());
            Volume* svol = new Volume(getSizeX(), getSizeY(), getSizeZ());

            //Initialize cluster counters
            int helices[MAX_SHEETS];
            for (i = 0; i < MAX_SHEETS; i++) {
                helices[i] = 0;
            }
            int totHelices = 1;

            //Start clustering
            printf("Start clustering...\n");
            int ox, oy, oz;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                if (getDataAt(i, j, k) <= 0 || svol->getDataAt(i, j, k) != 0) {
                    // Not a data point or has been visited
                    continue;
                }
                if (!isHelix(i, j, k)) {
                    // Not a helix point
                    continue;
                }

                //Initialize queue
                int numNodes = 1;
                svol->setDataAt(i, j, k, totHelices);
                GridQueue* queue = new GridQueue();
                queue->pushQueue(i, j, k);
                while (queue->popQueue(ox, oy, oz))
                {
                    // Test if neighbors satisfy helix condition
                    if (isHelix(ox, oy, oz)) {
                        for (int m = 0; m < 6; m++) {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];

                            if (getDataAt(nx, ny, nz) > 0 && svol->getDataAt(nx, ny, nz) == 0) {
                                svol->setDataAt(nx, ny, nz, totHelices);
                                queue->pushQueue(nx, ny, nz);
                                numNodes++;
                            }
                        }
                    }
                }

                delete queue;
                if (numNodes > 0) {
                    //	printf("Helix %d contain %d nodes.\n", totHelices, numNodes) ;
                    helices[totHelices] = numNodes;
                    totHelices++;
                }
            }

            // Removing clusters less than minSize
            printf("Removing small clusters.\n");
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++) {
                int cnt = (int)svol->getDataAt(i, j, k);
                if (cnt > 0 && helices[cnt] < minSize) {
                    svol->setDataAt(i, j, k, -1);
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            svol->threshold(0.1, 0, 1);

            return svol;

        }

        int Volume::isEndPoint(int ox, int oy, int oz) {
            if (getDataAt(ox - 1, oy, oz) < 0 && getDataAt(ox + 1, oy, oz) < 0) {
                return 1;
            }
            if (getDataAt(ox, oy - 1, oz) < 0 && getDataAt(ox, oy + 1, oz) < 0) {
                return 1;
            }
            if (getDataAt(ox, oy, oz - 1) < 0 && getDataAt(ox, oy, oz + 1) < 0) {
                return 1;
            }
            return 0;
        }

        int Volume::getNumNeighbor6(int ox, int oy, int oz) {
            int rvalue = 0;
            for (int i = 0; i < 6; i++) {
                int nx = ox + neighbor6[i][0];
                int ny = oy + neighbor6[i][1];
                int nz = oz + neighbor6[i][2];
                if (getDataAt(nx, ny, nz) >= 0) {
                    rvalue++;
                }
            }

            return rvalue;
        }
        int Volume::testIsSheetEnd(int ox, int oy, int oz) {
            // Returns 1 if it lies on the edge of a non-noise sheet
            int i, j;
            int nx, ny, nz;

            int edge[6] = { 0, 0, 0, 0, 0, 0 };
            int faceflag[12];
            int hasNoiseFace = 0;
            int tot = 0;

            for (i = 0; i < 12; i++) {
                faceflag[i] = 1;
                int hasNoise = 0;

                for (j = 0; j < 4; j++) {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) == 0) {
                        hasNoise = 1;
                    }
                    if (getDataAt(nx, ny, nz) < 0) {
                        faceflag[i] = 0;
                        break;
                    }
                }
                if (faceflag[i] == 1 && hasNoise) {
                    hasNoiseFace++;
                    // return 0 ;
                }

                if (faceflag[i]) {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                    tot++;
                }
            }

            if (hasNoiseFace == tot) {
                return 0;
            }

            if (tot == 0) {
                return 0;
            }

            // Removing 1s
            int numones = 0;
            for (i = 0; i < 6; i++) {
                if (edge[i] == 1) {
                    numones++;
                }
            }
            while (numones > 0) {
                int e;
                for (i = 0; i < 6; i++) {
                    if (edge[i] == 1) {
                        e = i;
                        break;
                    }
                }

                int f, e2;
                for (j = 0; j < 4; j++) {
                    f = edgeFaces[e][j];
                    if (faceflag[f]) {
                        break;
                    }
                }

                if (faceEdges[f][0] == e) {
                    e2 = faceEdges[f][1];
                }
                else {
                    e2 = faceEdges[f][0];
                }

                edge[e] --;
                numones--;
                edge[e2] --;
                faceflag[f] = 0;
                tot--;

                if (edge[e2] == 1) {
                    numones++;
                }
                else if (edge[e2] == 0) {
                    numones--;
                }
            }

            if (tot > 0) {
                return 0;
            }
            return 1;
        }

        int Volume::isNoiseSheetEnd(int ox, int oy, int oz) {
            int i, j;
            int nx, ny, nz;

            int edge[6] = { 0, 0, 0, 0, 0, 0 };
            int faceflag[12];
            int hasNoiseFace = 0;
            int tot = 0;

            for (i = 0; i < 12; i++) {
                faceflag[i] = 1;
                int hasNoise = 0;

                for (j = 0; j < 4; j++) {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) == 0) {
                        hasNoise = 1;
                    }
                    if (getDataAt(nx, ny, nz) < 0) {
                        faceflag[i] = 0;
                        break;
                    }
                }
                if (faceflag[i] == 1 && hasNoise) {
                    hasNoiseFace++;
                    // return 0 ;
                }

                if (faceflag[i]) {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                    tot++;
                }
            }

            if (hasNoiseFace < tot) {
                return 0;
            }

            if (tot == 0) {
                return 0;
            }

            // Removing 1s
            int numones = 0;
            for (i = 0; i < 6; i++) {
                if (edge[i] == 1)
                {
                    numones++;
                }
            }
            while (numones > 0) {
                int e;
                for (i = 0; i < 6; i++) {
                    if (edge[i] == 1) {
                        e = i;
                        break;
                    }
                }

                int f, e2;
                for (j = 0; j < 4; j++) {
                    f = edgeFaces[e][j];
                    if (faceflag[f]) {
                        break;
                    }
                }

                if (faceEdges[f][0] == e) {
                    e2 = faceEdges[f][1];
                }
                else {
                    e2 = faceEdges[f][0];
                }

                edge[e] --;
                numones--;
                edge[e2] --;
                faceflag[f] = 0;
                tot--;

                if (edge[e2] == 1) {
                    numones++;
                }
                else if (edge[e2] == 0) {
                    numones--;
                }
            }

            if (tot > 0) {
                return 0;
            }
            return 1;
        }

        int Volume::isInternal(int ox, int oy, int oz) {
            // assuming it's 0/1 volume
            int i, j, k;

            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                if (getDataAt(ox + i, oy + j, oz + k) <= 0) {
                    return 0;
                }
            }
            return 1;
        }

        int Volume::isInternal2(int ox, int oy, int oz) {
            // assuming it's -1/0 volume
            int i, j, k;

            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                if (getDataAt(ox + i, oy + j, oz + k) < 0) {
                    return 0;
                }
            }

            return 1;
        }

        int Volume::hasIsolatedFace(int ox, int oy, int oz) {
            int i, j, k;
            int nx, ny, nz;

            double vox[3][3][3];
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                vox[i + 1][j + 1][k + 1] = getDataAt(ox + i, oy + j, oz + k);
            }

            int cells[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
            for (i = 0; i < 8; i++) {
                int x = ((i >> 2) & 1);
                int y = ((i >> 1) & 1);
                int z = (i & 1);
                for (j = 0; j < 8; j++) {
                    nx = x + ((j >> 2) & 1);
                    ny = y + ((j >> 1) & 1);
                    nz = z + (j & 1);

                    if (vox[nx][ny][nz] < 0) {
                        cells[i] = 0;
                        break;
                    }
                }
            }

            for (i = 0; i < 12; i++) {
                if (cells[faceCells[i][0]] == 1 || cells[faceCells[i][1]] == 1) {
                    continue;
                }
                int flag = 1;
                for (j = 0; j < 4; j++) {
                    nx = 1 + sheetNeighbor[i][j][0];
                    ny = 1 + sheetNeighbor[i][j][1];
                    nz = 1 + sheetNeighbor[i][j][2];

                    if (vox[nx][ny][nz] < 0) {
                        flag = 0;
                        break;
                    }
                }
                if (flag) {
                    return 1;
                }
            }

            return 0;
        }

        int Volume::hasIsolatedEdge(int ox, int oy, int oz) {
            int i, j, k;
            int nx, ny, nz;

            double vox[3][3][3];
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++) {
                vox[i + 1][j + 1][k + 1] = getDataAt(ox + i, oy + j, oz + k);
            }

            int edge[6] = { 0, 0, 0, 0, 0, 0 };

            for (i = 0; i < 12; i++) {
                int flag = 1;
                for (j = 0; j < 4; j++) {
                    nx = 1 + sheetNeighbor[i][j][0];
                    ny = 1 + sheetNeighbor[i][j][1];
                    nz = 1 + sheetNeighbor[i][j][2];

                    if (vox[nx][ny][nz] < 0) {
                        flag = 0;
                        break;
                    }
                }

                if (flag) {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                }
            }

            for (i = 0; i < 6; i++) {
                if (edge[i]) {
                    continue;
                }

                nx = 1 + neighbor6[i][0];
                ny = 1 + neighbor6[i][1];
                nz = 1 + neighbor6[i][2];

                if (vox[nx][ny][nz] >= 0) {
                    return 1;
                }

            }

            return 0;
        }

        int Volume::countFace(int ox, int oy, int oz, int m) {
            int facenum = 4;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int nx = ox + sheetNeighbor[edgeFaces[m][i]][j][0];
                    int ny = oy + sheetNeighbor[edgeFaces[m][i]][j][1];
                    int nz = oz + sheetNeighbor[edgeFaces[m][i]][j][2];

                    if (getDataAt(nx, ny, nz) < 0) {
                        facenum--;
                        break;
                    }
                }
            }

            return facenum;
        }
        int Volume::hasCell(int ox, int oy, int oz) {
            for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                if (getDataAt(ox + i, oy + j, oz + k) < 0) {
                    return 0;
                }
            }
            return 1;
        }

        Volume * Volume::markCellFace() {
            int i, j, k;
            Volume* fvol = new Volume(getSizeX(), getSizeY(), getSizeZ());

            //return fvol ;

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (hasCell(i, j, k))
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = i + neighbor6[m][0];
                            int ny = j + neighbor6[m][1];
                            int nz = k + neighbor6[m][2];
                            if (!hasCell(nx, ny, nz))
                            {
                                fvol->setDataAt(i, j, k, (double)(1 << m));
                                break;
                            }
                        }
                    }
                }
            }


            return fvol;
        }

        Volume * Volume::markFaceEdge() {
            int x, y, z, i, j;
            Volume* fvol = new Volume(getSizeX(), getSizeY(), getSizeZ());

            //return fvol ;

            for (x = 0; x < getSizeX(); x++)
            for (y = 0; y < getSizeY(); y++)
            for (z = 0; z < getSizeZ(); z++)
            {
                if (getDataAt(x, y, z) >= 0)
                {

                    for (i = 0; i < 3; i++)
                    {
                        int hasFace = 1;
                        for (j = 0; j < 4; j++)
                        {
                            int nx = x + sheetNeighbor[4 * i + 3][j][0];
                            int ny = y + sheetNeighbor[4 * i + 3][j][1];
                            int nz = z + sheetNeighbor[4 * i + 3][j][2];

                            if (getDataAt(nx, ny, nz) < 0)
                            {
                                hasFace = 0;
                                break;
                            }
                        }

                        if (hasFace)
                        {
                            // Look for open edges
                            switch (i)
                            {
                            case 0:
                                if (countFace(x, y, z, 0) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 0));
                                    break;
                                }
                                if (countFace(x, y, z, 2) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 2));
                                    break;
                                }
                                if (countFace(x, y + 1, z, 0) == 1)
                                {
                                    fvol->setDataAt(x, y + 1, z, (double)(1 << 0));
                                    break;
                                }
                                if (countFace(x, y, z + 1, 2) == 1)
                                {
                                    fvol->setDataAt(x, y, z + 1, (double)(1 << 2));
                                    break;
                                }
                                printf("Hmmm... a face with no open edges.\n");
                                break;
                            case 1:
                                if (countFace(x, y, z, 0) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 0));
                                    break;
                                }
                                if (countFace(x, y, z, 4) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 4));
                                    break;
                                }
                                if (countFace(x + 1, y, z, 0) == 1)
                                {
                                    fvol->setDataAt(x + 1, y, z, (double)(1 << 0));
                                    break;
                                }
                                if (countFace(x, y, z + 1, 4) == 1)
                                {
                                    fvol->setDataAt(x, y, z + 1, (double)(1 << 4));
                                    break;
                                }
                                printf("Hmmm... a face with no open edges.\n");
                                break;
                            case 2:
                                if (countFace(x, y, z, 2) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 2));
                                    break;
                                }
                                if (countFace(x, y, z, 4) == 1)
                                {
                                    fvol->setDataAt(x, y, z, (double)(1 << 4));
                                    break;
                                }
                                if (countFace(x + 1, y, z, 2) == 1)
                                {
                                    fvol->setDataAt(x + 1, y, z, (double)(1 << 2));
                                    break;
                                }
                                if (countFace(x, y + 1, z, 4) == 1)
                                {
                                    fvol->setDataAt(x, y + 1, z, (double)(1 << 4));
                                    break;
                                }
                                printf("Hmmm... a face with no open edges.\n");
                                break;
                            }
                        }
                    }

                }
            }


            return fvol;
        }

        int Volume::hasCompleteSheet(int ox, int oy, int oz, Volume* fvol) {
            int i, j, k;
            int nx, ny, nz;

            int edge[6] = { 0, 0, 0, 0, 0, 0 };
            int faceflag[12];
            int tot = 0;
            int cellflag[8];

            int ct = 0;
            for (i = -1; i < 1; i++)
            for (j = -1; j < 1; j++)
            for (k = -1; k < 1; k++)
            {
                if (hasCell(ox + i, oy + j, oz + k))
                {
                    cellflag[ct] = 1;
                }
                else
                {
                    cellflag[ct] = 0;
                }
                ct++;
            }

            for (i = 0; i < 12; i++)
            {
                faceflag[i] = 1;
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        faceflag[i] = 0;
                        break;
                    }
                }

                if (faceflag[i])
                {
                    if (cellflag[faceCells[i][0]] ^ cellflag[faceCells[i][1]])
                    {
                        int v1 = (int)(fvol->getDataAt(
                            ox - 1 + ((faceCells[i][0] >> 2) & 1),
                            oy - 1 + ((faceCells[i][0] >> 1) & 1),
                            oz - 1 + ((faceCells[i][0]) & 1)));
                        int v2 = (int)(fvol->getDataAt(
                            ox - 1 + ((faceCells[i][1] >> 2) & 1),
                            oy - 1 + ((faceCells[i][1] >> 1) & 1),
                            oz - 1 + ((faceCells[i][1]) & 1)));
                        if (((v1 >> (2 * (2 - i / 4))) & 1) ||
                            ((v2 >> (2 * (2 - i / 4) + 1)) & 1))
                        {
                            faceflag[i] = 0;
                        }
                    }
                }

                if (faceflag[i])
                {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                    tot++;
                }
            }

            // Removing 1s
            int numones = 0;
            for (i = 0; i < 6; i++)
            {
                if (edge[i] == 1)
                {
                    numones++;
                }
            }
            while (numones > 0)
            {
                int e;
                for (i = 0; i < 6; i++)
                {
                    if (edge[i] == 1)
                    {
                        e = i;
                        break;
                    }
                }
                /*
                if ( edge[ e ] != 1 )
                {
                printf("Wrong Again!********\n") ;
                }
                */

                int f, e2;
                for (j = 0; j < 4; j++)
                {
                    f = edgeFaces[e][j];
                    if (faceflag[f])
                    {
                        break;
                    }
                }

                /*
                if ( faceflag[ f ] == 0 )
                {
                printf("Wrong!********\n") ;
                }
                */

                if (faceEdges[f][0] == e)
                {
                    e2 = faceEdges[f][1];
                }
                else
                {
                    e2 = faceEdges[f][0];
                }


                edge[e] --;
                numones--;
                edge[e2] --;
                faceflag[f] = 0;
                tot--;

                if (edge[e2] == 1)
                {
                    numones++;
                }
                else if (edge[e2] == 0)
                {
                    numones--;
                }
            }

            if (tot > 0)
            {
                return 1;
            }

            return 0;
        }

        int Volume::hasCompleteSheet(int ox, int oy, int oz) {
            // Returns 1 if it lies in the middle of a sheet
            int temp = countIntEuler(ox, oy, oz);
            if (temp > 0)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        int Volume::hasCompleteSheetSlow(int ox, int oy, int oz) {
            int i, j;
            int nx, ny, nz;

            int edge[6] = { 0, 0, 0, 0, 0, 0 };
            int faceflag[12];
            int tot = 0;

            for (i = 0; i < 12; i++)
            {
                faceflag[i] = 1;
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        faceflag[i] = 0;
                        break;
                    }
                }

                if (faceflag[i])
                {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                    tot++;
                }
            }

            // Removing 1s
            int numones = 0;
            for (i = 0; i < 6; i++)
            {
                if (edge[i] == 1)
                {
                    numones++;
                }
            }
            while (numones > 0)
            {
                int e;
                for (i = 0; i < 6; i++)
                {
                    if (edge[i] == 1)
                    {
                        e = i;
                        break;
                    }
                }
                /*
                if ( edge[ e ] != 1 )
                {
                printf("Wrong Again!********\n") ;
                }
                */

                int f, e2;
                for (j = 0; j < 4; j++)
                {
                    f = edgeFaces[e][j];
                    if (faceflag[f])
                    {
                        break;
                    }
                }

                /*
                if ( faceflag[ f ] == 0 )
                {
                printf("Wrong!********\n") ;
                }
                */

                if (faceEdges[f][0] == e)
                {
                    e2 = faceEdges[f][1];
                }
                else
                {
                    e2 = faceEdges[f][0];
                }


                edge[e] --;
                numones--;
                edge[e2] --;
                faceflag[f] = 0;
                tot--;

                if (edge[e2] == 1)
                {
                    numones++;
                }
                else if (edge[e2] == 0)
                {
                    numones--;
                }
            }

            //		int temp = countIntEuler( ox, oy, oz ) ;
            if (tot > 0)
            {
                //			if ( temp <= 0 )
                {
                    //				printf("Counting wrong: %d\n", temp ) ;
                }
                return 1;
            }
            else
            {
                //			if ( temp > 0 )
                {
                    //				printf("Counting wrong: %d\n", temp ) ;
                }
                return 0;
            }
        }
        int Volume::hasCompleteHelix(int ox, int oy, int oz)
        {
            // Returns 1 if it has a complete helix
            int i;
            int c1 = 0;
            int nx, ny, nz;
            int j;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];
                if (getDataAt(nx, ny, nz) >= 0)
                {
                    c1++;
                    j = i;
                }

            }

            if (c1 > 1) // || c1 == 0 )
            {
                return 1;
            }

            return 0;

            /*
            ox = ox + neighbor6[j][0] ;
            oy = oy + neighbor6[j][1] ;
            oz = oz + neighbor6[j][2] ;
            c1 = 0 ;
            for ( i = 0 ; i < 6 ; i ++ )
            {
            nx = ox + neighbor6[i][0] ;
            ny = oy + neighbor6[i][1] ;
            nz = oz + neighbor6[i][2] ;
            if ( getDataAt( nx, ny, nz ) >= 0 )
            {
            c1 ++ ;
            }

            }

            if ( c1 > 1 )
            {
            return 0 ;
            }
            else
            {
            return 1 ;
            }
            */
        }

        int Volume::hasCompleteHelix(int ox, int oy, int oz, Volume* fvol)
        {

            int i;
            int c1 = 0;
            int nx, ny, nz;
            int j;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];
                if (getDataAt(nx, ny, nz) >= 0)
                {
                    if (i % 2 == 0)
                    {
                        nx = ox;
                        ny = oy;
                        nz = oz;
                    }

                    int val = (int)fvol->getDataAt(nx, ny, nz);
                    if (((val >> (2 * (i / 2))) & 1) == 0)
                    {
                        c1++;
                        j = i;
                    }
                }

            }

            if (c1 > 1)
            {
                return 1;
            }

            return 0;

            /*
            ox = ox + neighbor6[j][0] ;
            oy = oy + neighbor6[j][1] ;
            oz = oz + neighbor6[j][2] ;
            c1 = 0 ;
            for ( i = 0 ; i < 6 ; i ++ )
            {
            nx = ox + neighbor6[i][0] ;
            ny = oy + neighbor6[i][1] ;
            nz = oz + neighbor6[i][2] ;
            if ( getDataAt( nx, ny, nz ) >= 0 )
            {
            c1 ++ ;
            }

            }

            if ( c1 > 1 )
            {
            return 0 ;
            }
            else
            {
            return 1 ;
            }
            */
        }

        int Volume::isHelixEnd(int ox, int oy, int oz, Volume* nvol) {
            // Returns 1 if it is a curve endpoint				
            int i;
            int c1 = 0, c2 = 0;
            int nx, ny, nz;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                double val = getDataAt(nx, ny, nz);

                if (val >= 0)
                {
                    c1++;
                    if (val > 0 && val < MAX_ERODE && nvol->getDataAt(nx, ny, nz) == 0)
                    {
                        c2++;
                    }
                }

            }

            if (c1 == 1 && c2 == 1)
            {
                return 1;
            }

            return 0;
        }

        int Volume::isFeature18(int ox, int oy, int oz)
        {
            // Border: > 0
            // Interior: == 0
            // Outside: < 0 
            if (getDataAt(ox, oy, oz) <= 0)
            {
                return 0;
            }

            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                if (i == 1 || j == 1 || k == 1)
                {
                    if (getDataAt(ox - 1 + i, oy - 1 + j, oz - 1 + k) == 0)
                    {
                        return 0;
                    }
                }
            }
            return 1;
        }

        int Volume::isEdgeEnd(int ox, int oy, int oz)
        {
            int i;
            int c1 = 0;
            int nx, ny, nz;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                double val = getDataAt(nx, ny, nz);

                if (val >= 0)
                {
                    c1++;
                }

            }

            if (c1 == 1)
            {
                return 1;
            }

            return 0;
        }

        int Volume::isFaceEnd(int ox, int oy, int oz)
        {
            // return isSheetEnd(ox,oy,oz) ;

            int i, j, k;
            int nx, ny, nz;

            double vox[3][3][3];
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                vox[i + 1][j + 1][k + 1] = getDataAt(ox + i, oy + j, oz + k);
            }

            int edge[6] = { 4, 4, 4, 4, 4, 4 };
            int edge2[6] = { 4, 4, 4, 4, 4, 4 };

            for (i = 0; i < 12; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    nx = 1 + sheetNeighbor[i][j][0];
                    ny = 1 + sheetNeighbor[i][j][1];
                    nz = 1 + sheetNeighbor[i][j][2];

                    if (vox[nx][ny][nz] < 0)
                    {
                        edge[faceEdges[i][0]] --;
                        edge[faceEdges[i][1]] --;
                        break;
                    }
                }
                for (j = 0; j < 4; j++)
                {
                    nx = 1 + sheetNeighbor[i][j][0];
                    ny = 1 + sheetNeighbor[i][j][1];
                    nz = 1 + sheetNeighbor[i][j][2];

                    if (vox[nx][ny][nz] < 2)
                    {
                        edge2[faceEdges[i][0]] --;
                        edge2[faceEdges[i][1]] --;
                        break;
                    }
                }
            }

            for (i = 0; i < 6; i++)
            {
                if (edge[i] == 1 && edge2[i] == 1)
                {
                    return 1;
                }
            }

            return 0;
        }

        int Volume::isNoise(int ox, int oy, int oz, int noise)
        {
            if (getDataAt(ox, oy, oz) == 1)
            {
                return 1;
            }

            if (noise > 0)
            {
                for (int i = 0; i < 6; i++)
                {
                    int nx = ox + neighbor6[i][0];
                    int ny = oy + neighbor6[i][1];
                    int nz = oz + neighbor6[i][2];

                    if (getDataAt(nx, ny, nz) > 0)
                    {
                        if (isNoise(nx, ny, nz, noise - 1))
                        {
                            return 1;
                        }
                    }
                }
            }

            return 0;

        }

        int Volume::isNoiseHelixEnd(int ox, int oy, int oz)
        {

            int i;
            int c1 = 0, c2 = 0;
            int nx, ny, nz;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                double val = getDataAt(nx, ny, nz);

                if (val >= 0)
                {
                    c1++;
                    if (val > 0 && val < MAX_ERODE)
                    {
                        c2++;
                    }
                }

            }

            if (c1 == 1 && c2 == 0)
            {
                return 1;
            }

            return 0;
        }


        int Volume::isHelixEnd(int ox, int oy, int oz) {

            int i;
            int c1 = 0, c2 = 0;
            int nx, ny, nz;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                double val = getDataAt(nx, ny, nz);

                if (val >= 0)
                {
                    c1++;
                    if (getNumNeighbor6(nx, ny, nz) < 6) // if ( val > 0 && val < MAX_ERODE ) 
                    {
                        c2++;
                    }
                }

            }

            if (c1 == 1 && c2 == 1)
            {
                return 1;
            }

            return 0;
        }
        int Volume::isSheetEnd(int ox, int oy, int oz, Volume* nvol)
        {
            // Returns 1 if it contains a sheet boundary. Noise-resistant
            int i, j;
            int nx, ny, nz;

            int edge[6] = { 0, 0, 0, 0, 0, 0 };
            int faceflag[12];
            int hasFeatureFace = 0;
            int tot = 0;

            for (i = 0; i < 12; i++)
            {
                faceflag[i] = 1;
                int hasFeature = 1;

                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) == 0 || nvol->getDataAt(nx, ny, nz) == 1)
                    {
                        hasFeature = 0;
                    }
                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        faceflag[i] = 0;
                        break;
                    }
                }
                if (faceflag[i] == 1 && hasFeature)
                {
                    hasFeatureFace++;
                    // return 0 ;
                }

                if (faceflag[i])
                {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                    tot++;
                }
            }

            if (tot == 0 || hasFeatureFace == 0)
            {
                return 0;
            }

            // Removing 1s
            int numones = 0;
            for (i = 0; i < 6; i++)
            {
                if (edge[i] == 1)
                {
                    numones++;
                }
            }
            while (numones > 0)
            {
                int e;
                for (i = 0; i < 6; i++)
                {
                    if (edge[i] == 1)
                    {
                        e = i;
                        break;
                    }
                }
                /*
                if ( edge[ e ] != 1 )
                {
                printf("Wrong Again!********\n") ;
                }
                */

                int f, e2;
                for (j = 0; j < 4; j++)
                {
                    f = edgeFaces[e][j];
                    if (faceflag[f])
                    {
                        break;
                    }
                }

                /*
                if ( faceflag[ f ] == 0 )
                {
                printf("Wrong!********\n") ;
                }
                */

                if (faceEdges[f][0] == e)
                {
                    e2 = faceEdges[f][1];
                }
                else
                {
                    e2 = faceEdges[f][0];
                }


                edge[e] --;
                numones--;
                edge[e2] --;
                faceflag[f] = 0;
                tot--;

                if (edge[e2] == 1)
                {
                    numones++;
                }
                else if (edge[e2] == 0)
                {
                    numones--;
                }
            }

            if (tot > 0)
            {
                return 0;
            }

            return 1;
        }

        int Volume::getNumFaces(int ox, int oy, int oz)
        {
            int i, j;
            int nx, ny, nz;

            int faces = 12;
            for (i = 0; i < 12; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        faces--;
                        break;
                    }
                }
            }

            return faces;
        }

        int Volume::getNumCells(int ox, int oy, int oz)
        {
            int i, j, k;

            int cells = 0;

            for (i = -1; i < 1; i++)
            for (j = -1; j < 1; j++)
            for (k = -1; k < 1; k++)
            {
                if (hasCell(ox + i, oy + j, oz + k))
                {
                    cells++;
                }
            }

            // printf("Faces: %d\n", faces);
            return cells;
        }


        int Volume::getNumIsolatedEdges(int ox, int oy, int oz)
        {
            int i, j, k;
            int nx, ny, nz;

            double vox[3][3][3];
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                vox[i + 1][j + 1][k + 1] = getDataAt(ox + i, oy + j, oz + k);
            }


            int edge[6] = { 0, 0, 0, 0, 0, 0 };

            for (i = 0; i < 12; i++)
            {
                int flag = 1;
                for (j = 0; j < 4; j++)
                {
                    nx = 1 + sheetNeighbor[i][j][0];
                    ny = 1 + sheetNeighbor[i][j][1];
                    nz = 1 + sheetNeighbor[i][j][2];

                    if (vox[nx][ny][nz] < 0)
                    {
                        flag = 0;
                        break;
                    }
                }

                if (flag)
                {
                    int e0 = faceEdges[i][0], e1 = faceEdges[i][1];
                    edge[e0] ++;
                    edge[e1] ++;
                }
            }

            int edges = 0;
            for (i = 0; i < 6; i++)
            {
                if (edge[i])
                {
                    continue;
                }

                nx = 1 + neighbor6[i][0];
                ny = 1 + neighbor6[i][1];
                nz = 1 + neighbor6[i][2];

                if (vox[nx][ny][nz] >= 0)
                {
                    edges++;
                }

            }

            return edges;
        }

        int Volume::getNumIsolatedFaces(int ox, int oy, int oz)
        {
            int i, j, k;
            int nx, ny, nz;

            int faces = 0;
            int cellflag[8];

            int ct = 0;
            for (i = -1; i < 1; i++)
            for (j = -1; j < 1; j++)
            for (k = -1; k < 1; k++)
            {
                if (hasCell(ox + i, oy + j, oz + k))
                {
                    cellflag[ct] = 1;
                }
                else
                {
                    cellflag[ct] = 0;
                }
                ct++;
            }

            for (i = 0; i < 12; i++)
            {
                int hasFace = 1;
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        hasFace = 0;
                        break;
                    }
                }

                if (hasFace)
                {
                    if (cellflag[faceCells[i][0]] == 0 && cellflag[faceCells[i][1]] == 0)
                    {
                        faces++;
                    }
                }
            }

            // printf("Faces: %d\n", faces);
            return faces;
        }

        int Volume::isFeatureFace2(int ox, int oy, int oz)
        {
            int i;
            int nx, ny, nz;

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                double val = getDataAt(nx, ny, nz);

                if (val == 0)
                {
                    return 0;
                }

            }

            return 1;
        }

        int Volume::isFeatureFace(int ox, int oy, int oz)
        {
            // return 1 ;

            int i, j;
            int nx, ny, nz;

            int faces = 12;
            for (i = 0; i < 12; i++)
            {
                int ct = 0;
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        ct = -1;
                        break;
                    }
                    else if (getNumNeighbor6(nx, ny, nz) == 6)
                    {
                        ct = -1;
                        break;
                    }
                    //				else if ( getDataAt( nx, ny, nz ) == 0 )
                    //				{
                    //					ct ++ ;
                    //				}


                }
                if (ct == -1 || ct >= 1)
                {
                    faces--;
                }
            }

            if (faces > 0)
            {
                return 1;
            }
            return 0;

        }

        int Volume::hasFeatureFace(int ox, int oy, int oz)
        {
            int i, j, k;
            int nx, ny, nz;

            int faceflag;
            int cellflag[8];

            // Get cells
            int ct = 0;
            for (i = -1; i < 1; i++)
            for (j = -1; j < 1; j++)
            for (k = -1; k < 1; k++)
            {
                if (hasCell(ox + i, oy + j, oz + k))
                {
                    cellflag[ct] = 1;
                }
                else
                {
                    cellflag[ct] = 0;
                }
                ct++;
            }

            // Find isolated and boundary faces
            for (i = 0; i < 12; i++)
            {
                faceflag = 1;
                for (j = 0; j < 4; j++)
                {
                    nx = ox + sheetNeighbor[i][j][0];
                    ny = oy + sheetNeighbor[i][j][1];
                    nz = oz + sheetNeighbor[i][j][2];

                    if (getDataAt(nx, ny, nz) < 0)
                    {
                        faceflag = 0;
                        break;
                    }
                }

                if (faceflag)
                {
                    if (cellflag[faceCells[i][0]] == 0 && cellflag[faceCells[i][1]] == 0)
                    {
                        return 1;
                    }
                }
            }

            return 0;

        }

        int Volume::isSheetEnd(int ox, int oy, int oz)
        {
            return ((hasCompleteSheet(ox, oy, oz) == 0) && isFeatureFace(ox, oy, oz));

            //// return testIsSheetEnd( ox, oy, oz ) ;
            //
            //int i, j, k ;
            //int nx, ny, nz ;

            //double vox[3][3][3] ;
            //for ( i = -1 ; i < 2 ; i ++ )
            //	for ( j = -1 ; j < 2 ; j ++ )
            //		for ( k = -1 ; k < 2 ; k ++ )
            //		{
            //			vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
            //		}

            //int edge[6] = { 4,4,4,4,4,4 } ;
            //int edge2[6] = { 4,4,4,4,4,4 } ;

            //for ( i = 0 ; i < 12 ; i ++ )
            //{	
            //	for ( j = 0 ; j < 4 ; j ++ )
            //	{
            //		nx = 1 + sheetNeighbor[i][j][0] ;
            //		ny = 1 + sheetNeighbor[i][j][1] ;
            //		nz = 1 + sheetNeighbor[i][j][2] ;

            //		if ( vox[nx][ny][nz] < 0 )
            //		{
            //			edge[ faceEdges[ i ][ 0 ] ] -- ;
            //			edge[ faceEdges[ i ][ 1 ] ] -- ;
            //			break ;
            //		}
            //	}

            //	for ( j = 0 ; j < 4 ; j ++ )
            //	{
            //		nx = 1 + sheetNeighbor[i][j][0] ;
            //		ny = 1 + sheetNeighbor[i][j][1] ;
            //		nz = 1 + sheetNeighbor[i][j][2] ;

            //		if ( vox[nx][ny][nz] <= 0 )
            //		{
            //			edge2[ faceEdges[ i ][ 0 ] ] -- ;
            //			edge2[ faceEdges[ i ][ 1 ] ] -- ;
            //			break ;
            //		}
            //	}
            //}

            //
            ///*
            //for ( i = 0 ; i < 6 ; i ++ )
            //{
            //	nx = 1 + neighbor6[i][0] ;
            //	ny = 1 + neighbor6[i][1] ;
            //	nz = 1 + neighbor6[i][2] ;
            //	if ( edge[i] == 0 && vox[nx][ny][nz] >= 0 ) 
            //	{
            //		return 0 ;
            //	}
            //}
            //*/
            //


            //for ( i = 0 ; i < 6 ; i ++ )
            //{
            //	if ( edge[ i ] == 1 && edge2[ i ] == 1 )
            //	{
            //		return 1 ;
            //	}
            //}

            //return 0 ;
        }

        int Volume::isSimple(int ox, int oy, int oz)
        {
            /* Test if this is a simple voxel */
            // int flag = 0 ;
            double vox[3][3][3];

            int i, j, k;
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                double tval = getDataAt(ox + i, oy + j, oz + k);

                /*
                if ( tval == 0 || tval > (va )
                {
                flag = 1 ;
                }
                */
                /*
                if ( tval < 0 && tval == - curwid )
                {
                printf("Here %d", (int)-tval) ;
                vox[ i + 1 ][ j + 1 ][ k + 1 ] = - tval ;
                }
                else
                */
                {
                    vox[i + 1][j + 1][k + 1] = tval;
                }
            }

            /* Debugging
            printf("{") ;
            for ( i = 0 ; i < 3 ; i ++ )
            {
            if ( i ) printf(",") ;
            printf("{") ;
            for ( j = 0 ; j < 3 ; j ++ )
            {
            if ( j ) printf(",") ;
            printf("{") ;
            for ( k = 0 ; k < 3 ; k ++ )
            {
            if ( k ) printf(",") ;
            printf("%d", (vox[i][j][k] >=0 ? 1: 0));
            }
            printf("}") ;
            }
            printf("}") ;
            }
            printf("} Int: %d, Ext: %d\n", countInt( vox ), countExt( vox )) ;
            */

            if (countInt(vox) == 1 && countExt(vox) == 1)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }


        int Volume::isPiercable(int ox, int oy, int oz)
        {
            /* Test if this is a simple voxel */
            // int flag = 0 ;
            double vox[3][3][3];

            int i, j, k;
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                double tval = getDataAt(ox + i, oy + j, oz + k);

                /*
                if ( tval == 0 || tval > (va )
                {
                flag = 1 ;
                }
                */
                /*
                if ( tval < 0 && tval == - curwid )
                {
                printf("Here %d", (int)-tval) ;
                vox[ i + 1 ][ j + 1 ][ k + 1 ] = - tval ;
                }
                else
                */
                {
                    vox[i + 1][j + 1][k + 1] = tval;
                }
            }

            /* Debugging
            printf("{") ;
            for ( i = 0 ; i < 3 ; i ++ )
            {
            if ( i ) printf(",") ;
            printf("{") ;
            for ( j = 0 ; j < 3 ; j ++ )
            {
            if ( j ) printf(",") ;
            printf("{") ;
            for ( k = 0 ; k < 3 ; k ++ )
            {
            if ( k ) printf(",") ;
            printf("%d", (vox[i][j][k] >=0 ? 1: 0));
            }
            printf("}") ;
            }
            printf("}") ;
            }
            printf("} Int: %d, Ext: %d\n", countInt( vox ), countExt( vox )) ;
            */

            if (countInt(vox) == 1 && countExt(vox) != 1)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }


        int Volume::isSimple2(int v[3][3][3])
        {
            // int flag = 0 ;
            double vox[3][3][3];

            int i, j, k;
            for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
            {
                if (v[i][j][k] == 0)
                {
                    vox[i][j][k] = 1;
                }
                else
                {
                    vox[i][j][k] = -1;
                }
            }
            if (countInt(vox) == 1 && countExt(vox) == 1)
            {
                return 1;
            }
            else
            {
                printf("Int: %d Ext: %d\n", countInt(vox), countExt(vox));
                return 0;
            }
        }

        int Volume::getNumPotComplex3(int ox, int oy, int oz)
        {
            // return 0 ;


            int i, j, k;
            if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz))
            {
                return 60;
            }
            double val = getDataAt(ox, oy, oz);

            int nx, ny, nz;


            int numSimple = 0;
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                nx = ox + i;
                ny = oy + j;
                nz = oz + k;
                if (getDataAt(nx, ny, nz) > 0)
                {
                    if (isSimple(nx, ny, nz) && !isSheetEnd(nx, ny, nz))
                    {
                        numSimple++;
                    }
                }
            }
            numSimple--;

            setDataAt(ox, oy, oz, -val);

            int numPotSimple = 0;
            for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
            for (k = -1; k < 2; k++)
            {
                nx = ox + i;
                ny = oy + j;
                nz = oz + k;
                if (getDataAt(nx, ny, nz) > 0)
                {
                    if (isSimple(nx, ny, nz) && !isSheetEnd(nx, ny, nz))
                    {
                        numPotSimple++;
                    }
                }
            }


            setDataAt(ox, oy, oz, val);




            return 30 - (numPotSimple - numSimple);
        }

        int Volume::getNumPotComplex4(int ox, int oy, int oz)
        {

            int cells = getNumCells(ox, oy, oz);
            //int ifaces = getNumIsolatedFaces(ox,oy,oz) ;
            int faces = getNumFaces(ox, oy, oz);
            //int iedges = getNumIsolatedEdges(ox,oy,oz) ;

            return (cells * 6 - 2 * faces); // + 2 * ( faces * 4 - 4 * iedges ) ;
        }

        int Volume::getNumPotComplex(int ox, int oy, int oz)
        {
            //return 0 ;

            int i;
            double val = getDataAt(ox, oy, oz);
            if (val <= 0)
            {
                //		return 70 ;
            }

            // return getNumNeighbor6( ox, oy, oz ) ;

            // int v = ((getNumNeighbor6( ox, oy, oz ) & 255) << 24) ;
            //int v = 0  ;

            int rvalue = 0, nx, ny, nz;
            setDataAt(ox, oy, oz, -val);

            /*
            for ( i = -1 ; i < 2 ; i ++ )
            for ( j = -1 ; j < 2 ; j ++ )
            for ( k = -1 ; k < 2 ; k ++ )
            {
            nx = ox + i ;
            ny = oy + j ;
            nz = oz + k ;
            if ( getDataAt( nx, ny, nz ) == val )
            {
            if ( isSheetEnd( nx, ny, nz) || ! isSimple ( nx, ny, nz ) )
            {
            rvalue ++ ;
            }
            }
            }
            */

            for (i = 0; i < 6; i++)
            {
                nx = ox + neighbor6[i][0];
                ny = oy + neighbor6[i][1];
                nz = oz + neighbor6[i][2];

                if (getDataAt(nx, ny, nz) >= 0)
                {
                    int num = getNumNeighbor6(nx, ny, nz);
                    if (num > rvalue)
                    {
                        rvalue = num;
                    }
                }
            }


            setDataAt(ox, oy, oz, val);

            return rvalue + getNumNeighbor6(ox, oy, oz) * 10;
            /*
            int v = (((rvalue + getNumNeighbor6( ox, oy, oz ) * 10) & 255) << 24) ;
            v |= ( ( ox & 255 ) << 16 )  ;
            v |= ( ( oy & 255 ) << 8 ) ;
            v |= ( ( oz & 255 ) ) ;
            return v ;
            */

        }

        int Volume::getNumPotComplex2(int ox, int oy, int oz)
        {
            return getNumPotComplex(ox, oy, oz);

            //int i, j, k ;
            //double val = getDataAt( ox, oy, oz ) ;
            //if ( val <= 0 )
            //{
            //	return 0 ;
            //}

            //int rvalue = 0, nx, ny, nz ;
            //setDataAt( ox, oy, oz, -val ) ;

            //for ( i = -1 ; i < 2 ; i ++ )
            //	for ( j = -1 ; j < 2 ; j ++ )
            //		for ( k = -1 ; k < 2 ; k ++ )
            //		{
            //			nx = ox + i ;
            //			ny = oy + j ;
            //			nz = oz + k ;
            //			if ( getDataAt( nx, ny, nz ) == val )
            //			{
            //				if ( isHelixEnd( nx, ny, nz) || ! isSimple ( nx, ny, nz ) )
            //				{
            //					rvalue ++ ;
            //				}
            //			}
            //		}

            //setDataAt( ox, oy, oz, val ) ;

            //return rvalue + getNumNeighbor6( ox, oy, oz ) * 30 ;
        }

        int Volume::getNumNeighbor(int ox, int oy, int oz)
        {
            int rvalue = 0;
            double val = getDataAt(ox, oy, oz);
            for (int i = 0; i < 6; i++)
            {
                int nx = ox + neighbor6[i][0];
                int ny = oy + neighbor6[i][1];
                int nz = oz + neighbor6[i][2];

                if (getDataAt(nx, ny, nz) == val)
                {
                    rvalue++;
                }
            }
            /*
            for ( int i = -1 ; i < 2 ; i ++ )
            for ( int j = -1 ; j < 2 ; j ++ )
            for ( int k = -1 ; k < 2 ; k ++ )
            {
            int nx = ox + i ;
            int ny = oy + j ;
            int nz = oz + k ;

            if ( getDataAt( nx, ny, nz ) == val )
            {
            rvalue ++ ;
            }
            }
            */
            return rvalue;
        }


        void Volume::setScoreNeighbor(GridQueue* queue)
        {
            //printf("Scoring each node with number of neighbors...\n") ;
            gridQueueEle* ele = queue->getHead();
            while (ele != NULL)
            {
                ele->score = getNumNeighbor(ele->x, ele->y, ele->z);
                ele = ele->next;
            }

            queue->sort(queue->getNumElements());
        }


        int Volume::components6(int vox[3][3][3])
        {
            // Stupid flood fill 
            int tot = 0;
            int queue[27][3];
            int vis[3][3][3];
            int head = 0, tail = 1;
            int i, j, k;
            for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
            {
                vis[i][j][k] = 0;
                if (vox[i][j][k])
                {
                    if (tot == 0)
                    {
                        queue[0][0] = i;
                        queue[0][1] = j;
                        queue[0][2] = k;
                        vis[i][j][k] = 1;
                    }
                    tot++;
                }
            }
            if (tot == 0)
            {
                return 0;
            }
            // printf("total: %d\n", tot) ;

            int ct = 1;
            while (head != tail)
            {
                int x = queue[head][0];
                int y = queue[head][1];
                int z = queue[head][2];
                head++;

                for (i = 0; i < 6; i++)
                {
                    int nx = x + neighbor6[i][0];
                    int ny = y + neighbor6[i][1];
                    int nz = z + neighbor6[i][2];
                    if (nx >= 0 && nx < 3 && ny >= 0 && ny < 3 && nz >= 0 && nz < 3)
                    {
                        if (vox[nx][ny][nz] && !vis[nx][ny][nz])
                        {
                            queue[tail][0] = nx;
                            queue[tail][1] = ny;
                            queue[tail][2] = nz;
                            tail++;
                            vis[nx][ny][nz] = 1;
                            ct++;
                        }
                    }
                }
            }

            if (ct == tot)
            {
                return 1;
            }
            else
            {
                return 2;
            }

        }
        int Volume::components26(int vox[3][3][3])
        {
            // Stupid flood fill 
            int tot = 0;
            int queue[27][3];
            int vis[3][3][3];
            int head = 0, tail = 1;
            int i, j, k;
            for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
            {
                vis[i][j][k] = 0;
                if (vox[i][j][k])
                {
                    if (tot == 0)
                    {
                        queue[0][0] = i;
                        queue[0][1] = j;
                        queue[0][2] = k;
                        vis[i][j][k] = 1;
                    }
                    tot++;
                }
            }
            if (tot == 0)
            {
                return 0;
            }

            int ct = 1;
            while (head != tail)
            {
                int x = queue[head][0];
                int y = queue[head][1];
                int z = queue[head][2];
                head++;

                for (i = -1; i < 2; i++)
                for (j = -1; j < 2; j++)
                for (k = -1; k < 2; k++)
                {
                    int nx = x + i;
                    int ny = y + j;
                    int nz = z + k;
                    if (nx >= 0 && nx < 3 && ny >= 0 && ny < 3 && nz >= 0 && nz < 3)
                    {
                        if (vox[nx][ny][nz] && !vis[nx][ny][nz])
                        {
                            queue[tail][0] = nx;
                            queue[tail][1] = ny;
                            queue[tail][2] = nz;
                            tail++;
                            vis[nx][ny][nz] = 1;
                            ct++;
                        }
                    }
                }
            }

            if (ct == tot)
            {
                return 1;
            }
            else
            {
                return 2;
            }

        }

        int Volume::countExt(double vox[3][3][3])
        {
            int tvox[3][3][3];

            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                if (vox[i][j][k] < 0)
                {
                    tvox[i][j][k] = 1;
                }
                else
                {
                    tvox[i][j][k] = 0;
                }
            }

            return components26(tvox);
        }

        int Volume::countInt(double vox[3][3][3])
        {
            int i, j, k;
            int tvox[3][3][3];

            for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
            {
                tvox[i][j][k] = 0;
            }

            for (i = 0; i < 6; i++)
            {
                int nx = 1 + neighbor6[i][0];
                int ny = 1 + neighbor6[i][1];
                int nz = 1 + neighbor6[i][2];
                if (vox[nx][ny][nz] >= 0)
                {
                    tvox[nx][ny][nz] = 1;
                    for (j = 0; j < 4; j++)
                    {
                        int nnx = neighbor64[i][j][0] + nx;
                        int nny = neighbor64[i][j][1] + ny;
                        int nnz = neighbor64[i][j][2] + nz;
                        if (vox[nnx][nny][nnz] >= 0)
                        {
                            tvox[nnx][nny][nnz] = 1;
                        }
                    }
                }
            }

            return components6(tvox);
        }

        int Volume::countIntEuler(int ox, int oy, int oz)
        {
            int nv = 0, ne = 0, nc = 0;

            int i, j, k;
            int tvox[3][3][3];
            double vox[3][3][3];

            for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
            {
                vox[i][j][k] = getDataAt(ox - 1 + i, oy - 1 + j, oz - 1 + k);
                tvox[i][j][k] = 0;
            }

            for (i = 0; i < 6; i++)
            {
                int nx = 1 + neighbor6[i][0];
                int ny = 1 + neighbor6[i][1];
                int nz = 1 + neighbor6[i][2];
                if (vox[nx][ny][nz] >= 0)
                {
                    tvox[nx][ny][nz] = 1;

                    nv++;

                    for (j = 0; j < 4; j++)
                    {
                        int nnx = neighbor64[i][j][0] + nx;
                        int nny = neighbor64[i][j][1] + ny;
                        int nnz = neighbor64[i][j][2] + nz;
                        if (vox[nnx][nny][nnz] >= 0)
                        {
                            if (tvox[nnx][nny][nnz] == 0)
                            {
                                tvox[nnx][nny][nnz] = 1;
                                nv++;
                            }

                            ne++;
                        }
                    }
                }
            }

            nc = components6(tvox);

            return (nc - (nv - ne));
        }

        void Volume::erodeNoTopo(float thr, int wid)
        {
            /* No topology check */
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue* queue = new GridQueue();

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            setDataAt(i, j, k, 1);
                            queue->pushQueue(i, j, k);
                            break;
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue->getNumElements());

            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif
            double val = 0;
            int ox, oy, oz;
            int curwid = 0;
            int total = 0, ct = 0;
            while (1)
            {
                if (ct == total)
                {
#ifdef VERBOSE
                    printf("Layer %d has %d nodes.\n", (int)curwid, total);
#endif
                    curwid++;
                    ct = 0;
                    total = queue->getNumElements();
                    if (total == 0)
                    {
                        break;
                    }
                }

                queue->popQueue(ox, oy, oz);
                val = getDataAt(ox, oy, oz);
                if (val > wid)
                {
                    break;
                }
                ct++;

                setDataAt(ox, oy, oz, -val);


                for (int m = 0; m < 6; m++)
                {
                    int nx = ox + neighbor6[m][0];
                    int ny = oy + neighbor6[m][1];
                    int nz = oz + neighbor6[m][2];
                    if (getDataAt(nx, ny, nz) == 0)
                    {
                        setDataAt(nx, ny, nz, val + 1);
                        queue->pushQueue(nx, ny, nz);
                    }
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }

        void Volume::erodeTopo(float thr, int wid)
        {
            /* Minimal topology check */
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue* queue = new GridQueue();

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            setDataAt(i, j, k, 1);
                            queue->pushQueue(i, j, k);
                            break;
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue->getNumElements());

            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif

            double val = 0;
            int ox, oy, oz;
            int curwid = 0;
            int total = 0, ct = 0;
            while (1)
            {
                if (ct == total)
                {
#ifdef VERBOSE
                    printf("Layer %d has %d nodes.\n", (int)curwid, total);
#endif
                    curwid++;
                    ct = 0;
                    total = queue->getNumElements();
                    if (total == 0)
                    {
                        break;
                    }
                }

                queue->popQueue(ox, oy, oz);
                val = getDataAt(ox, oy, oz);
                if (val > wid)
                {
                    break;
                }
                ct++;

                if (isSimple(ox, oy, oz))
                {
                    // Simple node, remove
                    setDataAt(ox, oy, oz, -val);
                }
                else
                {
                    // Preserve for topology
                    setDataAt(ox, oy, oz, val + 1);
                    queue->pushQueue(ox, oy, oz);
                }


                for (int m = 0; m < 6; m++)
                {
                    int nx = ox + neighbor6[m][0];
                    int ny = oy + neighbor6[m][1];
                    int nz = oz + neighbor6[m][2];
                    if (getDataAt(nx, ny, nz) == 0)
                    {
                        setDataAt(nx, ny, nz, val + 1);
                        queue->pushQueue(nx, ny, nz);
                    }
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }

        void Volume::erode2(float thr, int wid)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue = new GridQueue2();
            GridQueue2* queue2 = new GridQueue2();

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            setDataAt(i, j, k, 1);
                            queue->prepend(i, j, k);
                            break;
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue->getNumElements());

            // Perform erosion 
            // wid = MAX_ERODE ;
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            int ox, oy, oz;

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue->getNumElements(), curwid);
#endif

                /* set nodes for next layer
                while ( ( ele = queue->getNext() ) != NULL )
                {
                for ( int m = 0 ; m < 6 ; m ++ )
                {
                int nx = ele->x + neighbor6[m][0] ;
                int ny = ele->y + neighbor6[m][1] ;
                int nz = ele->z + neighbor6[m][2] ;
                if ( getDataAt( nx, ny, nz ) == 0 )
                {
                setDataAt( nx, ny, nz, curwid + 1 ) ;
                queue2->prepend( nx, ny, nz ) ;
                }
                }

                }
                */

                // erosion
                int seed[3] = { -1, -1, -1 };
                queue->reset();
                while (queue->getNumElements() > 0)
                {
                    if (seed[0] < 0) printf("After initial scoring...\n");
                    queue->reset();
                    ele = queue->getNext();

                    // Secure complex nodes
                    while (ele != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        // Check simple only if within the last modified range
                        if (seed[0] < 0 ||
                            (ox < seed[0] + 2 && ox > seed[0] - 2 &&
                            oy < seed[1] + 2 && oy > seed[1] - 2 &&
                            oz < seed[2] + 2 && oz > seed[2] - 2))
                        {
                            if (!isSimple(ox, oy, oz))
                            {
                                // Complex, set to next layer
                                setDataAt(ox, oy, oz, curwid + 1);
                                queue2->prepend(ox, oy, oz);
                                ele = queue->remove();

                                numComplex++;
                            }
                            else
                            {
                                ele = queue->getNext();
                            }
                        }
                        else
                        {
                            ele = queue->getNext();
                        }
                    }

                    // Remove the simple node with the most potential neighbors
                    queue->reset();
                    ele = queue->getNext();
                    int preScore = -1;
                    while (ele != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        // Update score only if within the last modified range
                        if (seed[0] < 0 ||
                            (ox < seed[0] + 3 && ox > seed[0] - 3 &&
                            oy < seed[1] + 3 && oy > seed[1] - 3 &&
                            oz < seed[2] + 3 && oz > seed[2] - 3))
                        {
                            ele->score = getNumPotComplex(ox, oy, oz);
                        }


                        if (ele->score < preScore)
                        {
                            // Swap
                            ele = queue->swap();
                        }
                        else
                        {
                            preScore = ele->score;
                        }

                        // At the end of the queue, remove this simple node
                        if (ele->next == NULL)
                        {
                            ox = ele->x;
                            oy = ele->y;
                            oz = ele->z;
                            setDataAt(ox, oy, oz, -1);
                            //						printf("%d %d %d\n", ox, oy, oz) ;
                            seed[0] = ox;
                            seed[1] = oy;
                            seed[2] = oz;
                            queue->remove();
                            // printf("Highest score: %d\n", preScore) ;


                            for (int m = 0; m < 6; m++)
                            {
                                int nx = ox + neighbor6[m][0];
                                int ny = oy + neighbor6[m][1];
                                int nz = oz + neighbor6[m][2];
                                if (getDataAt(nx, ny, nz) == 0)
                                {
                                    setDataAt(nx, ny, nz, curwid + 1);
                                    queue2->prepend(nx, ny, nz);
                                }
                            }


                            numSimple++;
                            ele = NULL;
                        }
                        else
                        {
                            ele = queue->getNext();
                        }
                    }

                }

                delete queue;
                queue = queue2;
                queue2 = new GridQueue2();
#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }

        void Volume::erodeShapeTopo(float thr, int wid)
        {
            /* Faster version of erode2 using priority queue */
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            setDataAt(i, j, k, 1);
                            queue2->prepend(i, j, k);
                            break;
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            // wid = MAX_ERODE ;
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                // First, 
                // check for complex nodes in queue2 
                // move them from queue2 to queue3
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue3->prepend(ox, oy, oz);
                        ele = queue2->remove();

                        numComplex++;
                    }
                    else
                    {
                        ele = queue2->getNext();
                    }
                }

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    queue->add(gp, -score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    score = -score;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    // Remove this simple node
                    setDataAt(ox, oy, oz, -1);
                    numSimple++;
                    // printf("Highest score: %d\n", score) ;

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for (i = -1; i < 2; i++)
                    for (j = -1; j < 2; j++)
                    for (k = -1; k < 2; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        // Check simple 
                        if (getDataAt(nx, ny, nz) == curwid && !isSimple(nx, ny, nz))
                        {
                            // Complex, set to next layer
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                            numComplex++;
                        }
                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue
                    /*
                    for ( i = -2 ; i < 3 ;i ++ )
                    for ( j = -2 ; j < 3 ; j ++ )
                    for ( k = -2 ; k < 3 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    if ( getDataAt( nx, ny, nz ) == curwid )
                    {
                    // Compute score
                    score = getNumPotComplex( nx, ny, nz ) ;

                    if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
                    {
                    // printf("Update\n") ;
                    scrvol->setDataAt( nx, ny, nz, score ) ;
                    // Push to queue
                    gp = new gridPoint ;
                    gp->x = nx ;
                    gp->y = ny ;
                    gp->z = nz ;
                    queue->add( gp, -score ) ;
                    }
                    }
                    }
                    */

                }


#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif
                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
            delete queue;

        }


        void Volume::erodeAtom(float thr, int wid, Volume* avol)
        {
            /* Erode to atoms */
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (avol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                setDataAt(i, j, k, 1);
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            // wid = MAX_ERODE ;
            printf("Start erosion to %d...\n", wid);
#endif

            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                // First, 
                // check for complex nodes in queue2 
                // move them from queue2 to queue3
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue3->prepend(ox, oy, oz);
                        ele = queue2->remove();

                        numComplex++;
                    }
                    else
                    {
                        ele = queue2->getNext();
                    }
                }

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    queue->add(gp, -score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    score = -score;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    // Remove this simple node
                    setDataAt(ox, oy, oz, -1);
                    numSimple++;
                    // printf("Highest score: %d\n", score) ;

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for (i = -1; i < 2; i++)
                    for (j = -1; j < 2; j++)
                    for (k = -1; k < 2; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        // Check simple 
                        if (getDataAt(nx, ny, nz) == curwid && !isSimple(nx, ny, nz))

                        {
                            // Complex, set to next layer
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                            numComplex++;
                        }
                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue
                    /*
                    for ( i = -2 ; i < 3 ;i ++ )
                    for ( j = -2 ; j < 3 ; j ++ )
                    for ( k = -2 ; k < 3 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    if ( getDataAt( nx, ny, nz ) == curwid )
                    {
                    // Compute score
                    score = getNumPotComplex( nx, ny, nz ) ;

                    if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
                    {
                    // printf("Update\n") ;
                    scrvol->setDataAt( nx, ny, nz, score ) ;
                    // Push to queue
                    gp = new gridPoint ;
                    gp->x = nx ;
                    gp->y = ny ;
                    gp->z = nz ;
                    queue->add( gp, -score ) ;
                    }
                    }
                    }
                    */

                }
#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif			

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
            delete queue;

        }

        void General_Data::computeGridShowPositions()
        {
            gridShowPositions = new float[3 * maxGridIndex];
            for (int index = 0; index < maxGridIndex; index++)
            {
                float position[3];
                float showPos[3];
                getGridPointPosD(index, position);
                getShowPos(position, showPos);
                gridShowPositions[3 * index] = showPos[0];
                gridShowPositions[3 * index + 1] = showPos[1];
                gridShowPositions[3 * index + 2] = showPos[2];
            }
        }

        float *General_Data::getGridShowPositions()
        {
            return gridShowPositions;
        }

        float *General_Data::getEdgePoints()
        {
            return edgePoints;
        }


        std::vector<Vertex> *General_Data::getVertices()
        {
            return vertices;
        }

        void padTo(std::string &str, const size_t num, const char paddingChar = ' ')
        {
            if (num > str.size()) {
                str.insert(0, num - str.size(), paddingChar);
            }
        }

        inline double doubleRound(double val)
        {
            if (val < 0) return ceil(val - 0.5);
            return floor(val + 0.5);
        }




        std::vector<Segment> *General_Data::getSegments()
        {
            return segments;
        }




        //template<class T>
        //py::list std_vector_to_py_list(const std::vector<T>& v)
        //{
        //    py::object get_iter = py::iterator<std::vector<T> >();
        //    py::object iter = get_iter(v);
        //    py::list l(iter);
        //    return l;
        //}

        template<class T>
        py::list std_vector_to_py_list(const std::vector<T>& v)
        {
            py::list temp;
            for (int i = 0; i < v.size(); i++) {
                temp.append(v[i]);
            }
            return temp;
        }

        py::list Volume::getVertexPos() {
            return std_vector_to_py_list(vertexPos);
        }

        py::list Volume::getMaxCurveEigenvectors() {
            return std_vector_to_py_list(maxCurveEigenvectors);
        }

        py::list Volume::getQuadEigenvectors() {
            return std_vector_to_py_list(quadEigenvectors);
        }

        py::list Volume::getMaxPointEigenvectors() {
            return std_vector_to_py_list(maxPointEigenvectors);
        }

        py::list Volume::getMinPointEigenvectors() {
            return std_vector_to_py_list(minPointEigenvectors);
        }

        py::list Volume::getSaddlePointEigenvectors() {
            return std_vector_to_py_list(saddlePointEigenvectors);
        }

        py::list Volume::getMinCurveEigenvectors() {
            return std_vector_to_py_list(minCurveEigenvectors);
        }

        py::list Volume::getSaddleCurveEigenvectors() {
            return std_vector_to_py_list(saddleCurveEigenvectors);
        }

        py::list Volume::getMinCurvePos() {
            return std_vector_to_py_list(minCurvePos);
        }

        py::list Volume::getMaxCurvePos() {
            return std_vector_to_py_list(maxCurvePos);
        }

        py::list Volume::getSaddleCurvePos() {
            return std_vector_to_py_list(saddleCurvePos);
        }

        py::list Volume::getMinSurfacePos() {
            return std_vector_to_py_list(minSurfacePos);
        }

        py::list Volume::getMaxSurfacePos() {
            return std_vector_to_py_list(maxSurfacePos);
        }

        py::list Volume::getMaxSaliencies() {
            return std_vector_to_py_list(maxCurveSaliencies);
        }

        py::list Volume::getMaxIntensities(){
            return std_vector_to_py_list(maxCurveIntensities);
        }

        py::list Volume::getMaxEigenvalues(){
            return std_vector_to_py_list(maxCurveEigenvalues);
        }

        py::list Volume::getMinSaliencies() {
            return std_vector_to_py_list(minCurveSaliencies);
        }

        py::list Volume::getMinIntensities(){
            return std_vector_to_py_list(minCurveIntensities);
        }

        py::list Volume::getMinEigenvalues(){
            return std_vector_to_py_list(minCurveEigenvalues);
        }

        py::list Volume::getSaddleSaliencies() {
            return std_vector_to_py_list(saddleCurveSaliencies);
        }

        py::list Volume::getSaddleIntensities(){
            return std_vector_to_py_list(saddleCurveIntensities);
        }

        py::list Volume::getSaddleEigenvalues(){
            return std_vector_to_py_list(saddleCurveEigenvalues);
        }

        py::list Volume::getQuadTypes() {
            return std_vector_to_py_list(quadTypes);
        }

        py::list Volume::hideSurfaces()
        {
            std::vector<int> hides;
            std::vector<Quad> quads;

            quads = (*currentExtremalData->getQuads());
            float pointR = pointRatio;
            float curveR = curveRatio;
            float surfaceR = surfaceRatio;
            float minG = minGeo;
            float maxG = maxGeo;
            float eigenV = eigenValue;

            for (int i = 0; i < quads.size(); i++) {
                Quad quad = quads[i];
                int type = quad.type;
                bool hide = false;
                if (saliencyCheck) {
                    float *saliencies = quad.relativeSaliencies;
                    if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])) hide = true;
                }
                if (intensityCheck) {
                    float maxI = currentExtremalData->maxIntensity;
                    float minI = currentExtremalData->minIntensity;
                    float localI = quad.localIntensity;
                    if ((type == 1) && ((localI - minI) / (maxI - minI) < maxG)) hide = true;
                    if ((type == 2) && ((localI - minI) / (maxI - minI) > minG)) hide = true;
                }

                if (eigenvalueCheck) {
                    float maxE = (*currentExtremalData->maxEigenvalue);
                    float minE = (*currentExtremalData->minEigenvalue);
                    float localE = quad.firstEigenvalue;
                    if ((localE - minE) / (maxE - minE) < eigenV) hide = true;
                }
                if (hide) {
                    hides.push_back(quad.vertIdxs[0]);
                    hides.push_back(quad.vertIdxs[1]);
                    hides.push_back(quad.vertIdxs[2]);
                    hides.push_back(quad.vertIdxs[3]);

                }
            }
            return std_vector_to_py_list(hides);
        }

        py::list Volume::hideSegments()
        {
            std::vector<int> hides;
            std::vector<Segment> segments;

            segments = (*currentExtremalData->getSegments());
            float pointR = pointRatio;
            float curveR = curveRatio;
            float surfaceR = surfaceRatio;
            float minG = minGeo;
            float maxG = maxGeo;
            float eigenV = eigenValue;

            for (int i = 0; i < segments.size(); i++) {
                Segment segment = segments[i];
                int type = segment.type;
                bool hide = false;
                if (saliencyCheck) {
                    float *saliencies = segment.relativeSaliencies;
                    if ((curveR * saliencies[1] < surfaceR * saliencies[0]) || (curveR * saliencies[1] < pointR * saliencies[2])) hide = true;
                }
                if (intensityCheck) {
                    float maxI = currentExtremalData->maxIntensity;
                    float minI = currentExtremalData->minIntensity;
                    float localI = segment.localIntensity;
                    if ((type == 1) && ((localI - minI) / (maxI - minI) < maxG)) hide = true;
                    if ((type == 2) && ((localI - minI) / (maxI - minI) > minG)) hide = true;
                }

                if (eigenvalueCheck) {
                    float maxE = (*currentExtremalData->maxEigenvalue);
                    float minE = (*currentExtremalData->minEigenvalue);
                    float localE = segment.firstEigenvalue;
                    if ((localE - minE) / (maxE - minE) < eigenV) hide = true;
                }
                if (hide) {
                    hides.push_back(segment.vertIdxs[0]);
                    hides.push_back(segment.vertIdxs[1]);
                }
            }
            return std_vector_to_py_list(hides);
        }

        py::list Volume::getExtremalBonds1() {
            return std_vector_to_py_list(extremalBonds1);
        }

        py::list Volume::getExtremalBonds2() {
            return std_vector_to_py_list(extremalBonds2);
        }

        py::list Volume::getMinCurveBonds1() {
            return std_vector_to_py_list(minCurve1);
        }

        py::list Volume::getMinCurveBonds2() {
            return std_vector_to_py_list(minCurve2);
        }

        py::list Volume::getSaddleCurveBonds1() {
            return std_vector_to_py_list(saddleCurve1);
        }

        py::list Volume::getSaddleCurveBonds2() {
            return std_vector_to_py_list(saddleCurve2);
        }

        py::list Volume::getDisplayQuads() {
            py::list quadHashes;
            for (int i = 0; i < quadDisplays.size(); i++) {
                quadHashes.append(quadDisplays[i]);
            }
            return quadHashes;
        }

        py::list Volume::getQuadSaliencies() {
            return std_vector_to_py_list(quadSaliencies);
        }

        py::list Volume::getQuadIntensities() {
            return std_vector_to_py_list(quadIntensities);
        }

        py::list Volume::getQuadEigenvalues() {
            return std_vector_to_py_list(quadEigenvalues);
        }

        py::list Volume::getMaxCurveEigenvalues0() {
            return std_vector_to_py_list(maxCurveEigenvalues0);
        }

        py::list Volume::getMaxCurveEigenvalues1() {
            return std_vector_to_py_list(maxCurveEigenvalues1);
        }

        py::list Volume::getMinCurveEigenvalues0() {
            return std_vector_to_py_list(minCurveEigenvalues0);
        }

        py::list Volume::getMinCurveEigenvalues1() {
            return std_vector_to_py_list(minCurveEigenvalues1);
        }

        py::list Volume::getSaddleCurveEigenvalues0() {
            return std_vector_to_py_list(saddleCurveEigenvalues0);
        }

        py::list Volume::getSaddleCurveEigenvalues1() {
            return std_vector_to_py_list(saddleCurveEigenvalues1);
        }

        py::list Volume::getMaxPointEigenvalues0() {
            return std_vector_to_py_list(maxPointEigenvalues0);
        }

        py::list Volume::getMaxPointEigenvalues1() {
            return std_vector_to_py_list(maxPointEigenvalues1);
        }

        py::list Volume::getMinPointEigenvalues0() {
            return std_vector_to_py_list(minPointEigenvalues0);
        }

        py::list Volume::getMinPointEigenvalues1() {
            return std_vector_to_py_list(minPointEigenvalues1);
        }

        py::list Volume::getSaddlePointEigenvalues0() {
            return std_vector_to_py_list(saddlePointEigenvalues0);
        }

        py::list Volume::getSaddlePointEigenvalues1() {
            return std_vector_to_py_list(saddlePointEigenvalues1);
        }

        py::list Volume::getQuadEigenvalues0() {
            return std_vector_to_py_list(quadEigenvalues0);
        }

        py::list Volume::getQuadEigenvalues1() {
            return std_vector_to_py_list(quadEigenvalues1);
        }

        py::list Volume::getExtremalNormals() {
            return std_vector_to_py_list(quadNormals);
        }

        py::list Volume::getQuadIndices() {
            return std_vector_to_py_list(quadIndices);
        }

        py::list Volume::getMaxHashes() {
            py::list maxHashes;
            for (int i = 0; i < maxCurveBonds.size(); i++) {
                maxHashes.append(maxCurveBonds[i]);
            }
            return maxHashes;
            //return std_vector_to_py_list(maxCurveBonds);
        }

        py::list Volume::getMinHashes() {
            py::list minHashes;
            for (int i = 0; i < minCurveBonds.size(); i++) {
                minHashes.append(minCurveBonds[i]);
            }
            return minHashes;
        }

        py::list Volume::getSaddleHashes() {
            py::list saddleHashes;
            for (int i = 0; i < saddleCurveBonds.size(); i++) {
                saddleHashes.append(saddleCurveBonds[i]);
            }
            return saddleHashes;
        }

        std::vector<Quad> *General_Data::getQuads()
        {
            return quads;
        }

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

        template <class T>
        inline std::string to_string(const T& t)
        {
            std::stringstream ss;
            ss << t;
            return ss.str();
        }

        void Volume::extremalCurveSkeleton(int kernelSize, Volume* svol, int width, int height, int slices, int resolution, float thickness, float isovalue, int kernelSizeGone)
        {
            General_Data* currentData = new MRC_Data();
            currentData->setVolumeData(svol->getVolumeData());
            currentData->height = height;
            currentData->width = width;
            currentData->slices = slices;
            currentData->resolution = resolution;
            currentData->thickness = thickness;
            currentData->isovalue = isovalue;
            currentData->kernelSize = kernelSize;
            //cout << "Arrrr origins " << svol->getOriginX() << " " << svol->getOriginY() << " " << svol->getOriginZ() << endl;
            /**
            for(int i = 0; i < svol->getSizeX(); i++) {
            for(int j = 0; j < svol->getSizeY(); j++) {
            for(int k = 0; k < svol->getSizeZ(); k++) {

            cout << svol->getDataAt(i, j, k) << endl;
            }
            }
            }
            **/
            float originX = svol->getOriginX();
            float originY = svol->getOriginY();
            float originZ = svol->getOriginZ();
            currentData->dataGeneration("");
            currentData->buildGrid();
            MRC_Data *derivedData = dynamic_cast<MRC_Data*>(currentData);
            float *scalars = derivedData->scalars;
            float *tensors = derivedData->tensors;
            float *gradients = derivedData->gradients;
            float *edgeTable = derivedData->edgeTable;
            float *faceTable = derivedData->faceTable;
            
            currentData->edgePhase(scalars, tensors, gradients, edgeTable);

            currentData->facePhase(scalars, tensors, gradients, faceTable);

            currentData->cellPhase(scalars, tensors, gradients);


            currentData->buildCurve();
            currentData->buildSurface();
            svol->setCurrentGenData(currentData);
            setMinMaxIntensitiesEigenvalues(currentData->minIntensity, currentData->maxIntensity, (*currentData->minEigenvalue), (*currentData->maxEigenvalue));
            std::vector<Segment> segments = (*currentData->getSegments());

            for (int i = 0; i < segments.size(); i++) {
                if (segments[i].type == 1) {

                    extremalBonds1.push_back(segments[i].vertIdxs[0]);
                    extremalBonds2.push_back(segments[i].vertIdxs[1]);
                    unsigned long long hash1 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[0], "CA");
                    unsigned long long hash2 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[1], "CA");
                    //cout << "hash 1 " << hash1 << endl;
                    maxCurveBonds.push_back(to_string(hash1));
                    maxCurveBonds.push_back(to_string(hash2));

                    maxCurveSaliencies.push_back(segments[i].relativeSaliencies[0]);
                    maxCurveSaliencies.push_back(segments[i].relativeSaliencies[1]);
                    maxCurveSaliencies.push_back(segments[i].relativeSaliencies[2]);
                    maxCurveIntensities.push_back(segments[i].localIntensity);
                    //float[3] origMaxEigenvalue;
                    //origMaxEigenvalue[0] = segments[i].eigenValue1;
                    //origMaxEigenvalue[1] = segments[i].eigenValue2;
                    //origMaxEigenvalue[2] = segments[i].firstEigenvalue;

                    //float[3] displayMaxEigenvalue;
                    //gorgonDisplayPos(displayMaxEigenvalue, origMaxEigenvalue);

                    //maxCurveEigenvalues.push_back(displayMaxEigenvalue[2]);
                    //maxCurveEigenvalues0.push_back(displayMaxEigenvalue[0]);
                    //maxCurveEigenvalues1.push_back(displayMaxEigenvalue[1]);
                    maxCurveEigenvalues.push_back(segments[i].firstEigenvalue);
                    maxCurveEigenvalues0.push_back(segments[i].eigenValue1);
                    maxCurveEigenvalues1.push_back(segments[i].eigenValue2);
                    //cout << "segment eigenvalues " << segments[i].eigenValue1 << " " << segments[i].eigenValue2 << " " << segments[i].firstEigenvalue << endl;
                    //int maxCurveSize = segments[i].eigenVectors.size() / 3;
                    for (int j = 0; j < segments[i].eigenVectors.size(); j++) {
                        //float origEigenvector[3];
                        //origEigenvector[0] = segments[i].eigenVectors[3*j];
                        //origEigenvector[1] = segments[i].eigenVectors[3*j+1];
                        //origEigenvector[2] = segments[i].eigenVectors[3*j+2];

                        //float displayEigenVector[3];
                        //currentData->reverseShowPos(displayEigenVector, origEigenvector);

                        //maxCurveEigenvectors.push_back(displayEigenVector[0]);
                        //maxCurveEigenvectors.push_back(displayEigenVector[1]);
                        //maxCurveEigenvectors.push_back(displayEigenVector[2]);

                        maxCurveEigenvectors.push_back(segments[i].eigenVectors[j]);
                    }

                }
                if (segments[i].type == 2) {
                    minCurve1.push_back(segments[i].vertIdxs[0]);
                    minCurve2.push_back(segments[i].vertIdxs[1]);
                    unsigned long long hash1 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[0], "CA");
                    unsigned long long hash2 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[1], "CA");
                    minCurveBonds.push_back(to_string(hash1));
                    minCurveBonds.push_back(to_string(hash2));
                    minCurveSaliencies.push_back(segments[i].relativeSaliencies[0]);
                    minCurveSaliencies.push_back(segments[i].relativeSaliencies[1]);
                    minCurveSaliencies.push_back(segments[i].relativeSaliencies[2]);
                    minCurveIntensities.push_back(segments[i].localIntensity);
                    minCurveEigenvalues.push_back(segments[i].firstEigenvalue);

                    float eigenValue0 = segments[i].relativeSaliencies[2] * segments[i].firstEigenvalue;
                    minCurveEigenvalues0.push_back(eigenValue0);
                    float eigenValue1 = (segments[i].relativeSaliencies[1] * segments[i].firstEigenvalue) + eigenValue0;
                    minCurveEigenvalues1.push_back(eigenValue1);
                    for (int j = 0; j < segments[i].eigenVectors.size(); j++) {
                        minCurveEigenvectors.push_back(segments[i].eigenVectors[j]);
                    }
                }
                if (segments[i].type == 3) {
                    saddleCurve1.push_back(segments[i].vertIdxs[0]);
                    saddleCurve2.push_back(segments[i].vertIdxs[1]);
                    unsigned long long hash1 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[0], "CA");
                    unsigned long long hash2 = ConstructHashKey("___3", 'A', (unsigned int)segments[i].vertIdxs[1], "CA");
                    saddleCurveBonds.push_back(to_string(hash1));
                    saddleCurveBonds.push_back(to_string(hash2));
                    saddleCurveSaliencies.push_back(segments[i].relativeSaliencies[0]);
                    saddleCurveSaliencies.push_back(segments[i].relativeSaliencies[1]);
                    saddleCurveSaliencies.push_back(segments[i].relativeSaliencies[2]);
                    saddleCurveIntensities.push_back(segments[i].localIntensity);
                    saddleCurveEigenvalues.push_back(segments[i].firstEigenvalue);

                    float eigenValue0 = segments[i].relativeSaliencies[2] * segments[i].firstEigenvalue;
                    saddleCurveEigenvalues0.push_back(eigenValue0);
                    float eigenValue1 = (segments[i].relativeSaliencies[1] * segments[i].firstEigenvalue) + eigenValue0;

                    saddleCurveEigenvalues1.push_back(eigenValue1);
                    for (int j = 0; j < segments[i].eigenVectors.size(); j++) {
                        saddleCurveEigenvectors.push_back(segments[i].eigenVectors[j]);
                    }
                }
            }


            std::vector<Vertex> vertices = (*currentData->getVertices());

            int sizex = currentData->gridx;
            int sizey = currentData->gridy;
            int sizez = currentData->gridz;

            float minPts[3];
            float maxPts[3];
            minPts[0] = minPts[1] = minPts[2] = 1000.0;
            maxPts[0] = maxPts[1] = maxPts[2] = -1000.0;

            for (int i = 0; i < vertices.size(); i++) {
                Vertex vertex = vertices[i];
                if (vertices[i].position[0] < minPts[0]) {
                    minPts[0] = vertices[i].position[0];
                }
                if (vertices[i].position[1] < minPts[1]) {
                    minPts[1] = vertices[i].position[1];
                }
                if (vertices[i].position[2] < minPts[2]) {
                    minPts[2] = vertices[i].position[2];
                }
                if (vertices[i].position[0] > maxPts[0]) {
                    maxPts[0] = vertices[i].position[0];
                }
                if (vertices[i].position[1] > maxPts[1]) {
                    maxPts[1] = vertices[i].position[1];
                }
                if (vertices[i].position[2] > maxPts[2]) {
                    maxPts[2] = vertices[i].position[2];
                }
            }

            float minPtsR[3];
            float maxPtsR[3];
            minPtsR[0] = minPtsR[1] = minPtsR[2] = 1000.0;
            maxPtsR[0] = maxPtsR[1] = maxPtsR[2] = -1000.0;
            for (int i = 0; i < vertices.size(); i++) {
                Vertex vertex = vertices[i];
                float gorgonDisplayPos[3];
                currentData->reverseShowPos(gorgonDisplayPos, vertex.position);
                if (gorgonDisplayPos[0] < minPtsR[0]) {
                    minPtsR[0] = gorgonDisplayPos[0];
                }
                if (gorgonDisplayPos[1] < minPtsR[1]) {
                    minPtsR[1] = gorgonDisplayPos[1];
                }
                if (gorgonDisplayPos[2] < minPtsR[2]) {
                    minPtsR[2] = gorgonDisplayPos[2];
                }
                if (gorgonDisplayPos[0] > maxPtsR[0]) {
                    maxPtsR[0] = gorgonDisplayPos[0];
                }
                if (gorgonDisplayPos[1] > maxPtsR[1]) {
                    maxPtsR[1] = gorgonDisplayPos[1];
                }
                if (gorgonDisplayPos[2] > maxPtsR[2]) {
                    maxPtsR[2] = gorgonDisplayPos[2];
                }
            }
            float centerX = (maxPtsR[0] - minPtsR[0]) / 2.0;
            float centerY = (maxPtsR[1] - minPtsR[1]) / 2.0;
            float centerZ = (maxPtsR[2] - minPtsR[2]) / 2.0;


            float xScale = (svol->getSizeX() / (maxPts[0] - minPts[0]));
            float yScale = (svol->getSizeY() / (maxPts[1] - minPts[1]));
            float zScale = (svol->getSizeZ() / (maxPts[2] - minPts[2]));

            float xScaleSpacing = ((float)svol->getSizeX()) / ((float)svol->getSizeX() - 8.0);
            float yScaleSpacing = ((float)svol->getSizeY()) / ((float)svol->getSizeY() - 8.0);
            float zScaleSpacing = ((float)svol->getSizeZ()) / ((float)svol->getSizeZ() - 8.0);

            ofstream myfile;
            myfile.open("extremal.pdb");

            for (int i = 0; i < vertices.size(); i++) {
                Vertex vertex = vertices[i];
                string vIndex = std::to_string(i);
                padTo(vIndex, 5);
                float gorgonDisplayPos[3];
                currentData->reverseShowPos(gorgonDisplayPos, vertex.position);
                float xPos = gorgonDisplayPos[0] + originX;
                float yPos = gorgonDisplayPos[1] + originY;
                float zPos = gorgonDisplayPos[2] + originZ;
                vertexPos.push_back(xPos);
                vertexPos.push_back(yPos);
                vertexPos.push_back(zPos);
                //vertexPos.push_back(svol->getSpacingX()*( svol->xMin + (svol->xMax - svol->xMin) * (vertex.position[0]-minPts[0])/(maxPts[0]-minPts[0])));
                //vertexPos.push_back(svol->getSpacingY()*( svol->yMin + (svol->yMax - svol->yMin) * (vertex.position[1]-minPts[1])/(maxPts[1]-minPts[1])));
                //vertexPos.push_back(svol->getSpacingZ()*( svol->zMin + (svol->zMax - svol->zMin) * (vertex.position[2]-minPts[2])/(maxPts[2]-minPts[2])));

                string posx = std::to_string(xPos);
                string posy = std::to_string(yPos);
                string posz = std::to_string(zPos);
                //string posx = std::to_string(svol->getSpacingX()*( svol->xMin + (svol->xMax - svol->xMin) * (vertex.position[0]-minPts[0])/(maxPts[0]-minPts[0]))  );

                //string posy = std::to_string(svol->getSpacingY()*( svol->yMin + (svol->yMax - svol->yMin)  * (vertex.position[1]-minPts[1])/(maxPts[1]-minPts[1])) );

                //string posz = std::to_string(svol->getSpacingZ()*( svol->zMin + (svol->zMax - svol->zMin) * (vertex.position[2]-minPts[2])/(maxPts[2]-minPts[2])) );

                myfile << "ATOM  " << vIndex << "  CA  ALA " << vIndex << "     " << posx.substr(0, 7) << " " << posy.substr(0, 7) << " " << posz.substr(0, 7) << "  1.00  1.00      S_00  0 " << endl;

            }
            myfile.close();
            std::vector<Quad> quads = (*currentData->getQuads());
            for (int i = 0; i < quads.size(); i++)
            {
                Quad quad = quads[i];
                float p1[3], p2[3], p3[3], p4[3];
                p1[0] = vertices[quad.vertIdxs[0]].position[0] * xScale;
                p1[1] = vertices[quad.vertIdxs[0]].position[1] * yScale;
                p1[2] = vertices[quad.vertIdxs[0]].position[2] * zScale;

                p2[0] = vertices[quad.vertIdxs[1]].position[0] * xScale;
                p2[1] = vertices[quad.vertIdxs[1]].position[1] * yScale;
                p2[2] = vertices[quad.vertIdxs[1]].position[2] * zScale;

                p3[0] = vertices[quad.vertIdxs[2]].position[0] * xScale;
                p3[1] = vertices[quad.vertIdxs[2]].position[1] * yScale;
                p3[2] = vertices[quad.vertIdxs[2]].position[2] * zScale;

                p4[0] = vertices[quad.vertIdxs[3]].position[0] * xScale;
                p4[1] = vertices[quad.vertIdxs[3]].position[1] * yScale;
                p4[2] = vertices[quad.vertIdxs[3]].position[2] * zScale;

                int type = quad.type;
                bool hide = false;
                float normal[3];
                getQuadNorm(p1, p2, p3, p4, normal);
                unsigned long long quadVert0 = ConstructHashKey("___3", 'A', (unsigned int)quad.vertIdxs[0], "CA");
                unsigned long long quadVert1 = ConstructHashKey("___3", 'A', (unsigned int)quad.vertIdxs[1], "CA");
                unsigned long long quadVert2 = ConstructHashKey("___3", 'A', (unsigned int)quad.vertIdxs[2], "CA");
                unsigned long long quadVert3 = ConstructHashKey("___3", 'A', (unsigned int)quad.vertIdxs[3], "CA");

                quadDisplays.push_back(to_string(quadVert0));
                quadDisplays.push_back(to_string(quadVert1));
                quadDisplays.push_back(to_string(quadVert2));
                quadDisplays.push_back(to_string(quadVert3));

                quadIndices.push_back(quad.vertIdxs[0]);
                quadIndices.push_back(quad.vertIdxs[1]);
                quadIndices.push_back(quad.vertIdxs[2]);
                quadIndices.push_back(quad.vertIdxs[3]);

                //quadDisplays.push_back(quad.vertIdxs[0]);
                //quadDisplays.push_back(quad.vertIdxs[1]);
                //quadDisplays.push_back(quad.vertIdxs[2]);
                //quadDisplays.push_back(quad.vertIdxs[3]);
                quadNormals.push_back(normal[0]);
                quadNormals.push_back(normal[1]);
                quadNormals.push_back(normal[2]);
                quadTypes.push_back(type);
                quadSaliencies.push_back(quad.relativeSaliencies[0]);
                quadSaliencies.push_back(quad.relativeSaliencies[1]);
                quadSaliencies.push_back(quad.relativeSaliencies[2]);
                quadIntensities.push_back(quad.localIntensity);
                quadEigenvalues.push_back(quad.firstEigenvalue);

                float eigenValue0 = quad.relativeSaliencies[2] * quad.firstEigenvalue;
                quadEigenvalues0.push_back(eigenValue0);
                float eigenValue1 = (quad.relativeSaliencies[1] * quad.firstEigenvalue) + eigenValue0;
                quadEigenvalues1.push_back(eigenValue1);

                for (int i = 0; i < quad.eigenVectors.size(); i++) {
                    quadEigenvectors.push_back(quad.eigenVectors[i]);
                }
            }
            std::vector<Point> points = (*currentData->getPoints());
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].eigenVectors.size(); j++) {
                    if (points[i].type == 1) {
                        maxPointEigenvectors.push_back(points[i].eigenVectors[j]);
                    }
                    if (points[i].type == 2) {
                        minPointEigenvectors.push_back(points[i].eigenVectors[j]);
                    }
                    if (points[i].type == 3) {
                        saddlePointEigenvectors.push_back(points[i].eigenVectors[j]);
                    }
                }
            }



        }

        py::list Volume::getMaxPointSaliencies() {
            return std_vector_to_py_list(maxPointSaliencies);
        }
        py::list Volume::getMaxPointIntensities() {
            return std_vector_to_py_list(maxPointIntensities);
        }
        py::list Volume::getMaxPointEigenvalues() {
            return std_vector_to_py_list(maxPointEigenvalues);
        }
        py::list Volume::getMinPointSaliencies() {
            return std_vector_to_py_list(minPointSaliencies);
        }
        py::list Volume::getMinPointIntensities() {
            return std_vector_to_py_list(minPointIntensities);
        }
        py::list Volume::getMinPointEigenvalues() {
            return std_vector_to_py_list(minPointEigenvalues);
        }

        py::list Volume::getSaddlePointSaliencies() {
            return std_vector_to_py_list(saddlePointSaliencies);
        }
        py::list Volume::getSaddlePointIntensities() {
            return std_vector_to_py_list(saddlePointIntensities);
        }
        py::list Volume::getSaddlePointEigenvalues() {
            return std_vector_to_py_list(saddlePointEigenvalues);
        }

        py::list Volume::getExtremalMinPoints() {
            std::vector<Point> points = (*currentExtremalData->getPoints());
            std::vector<Vertex> vertices = (*currentExtremalData->getVertices());
            std::vector<int> pointsIndices;
            for (int i = 0; i < points.size(); i++)
            {
                Point point = points[i];
                int type = point.type;
                float *saliencies = point.relativeSaliencies;
                float localI = point.localIntensity;
                float localE = point.firstEigenvalue;
                if (type == 2) {
                    minPointSaliencies.push_back(saliencies[0]);
                    minPointSaliencies.push_back(saliencies[1]);
                    minPointSaliencies.push_back(saliencies[2]);
                    minPointIntensities.push_back(localI);
                    minPointEigenvalues.push_back(localE);

                    float eigenValue0 = saliencies[2] * localE;
                    minPointEigenvalues0.push_back(eigenValue0);
                    float eigenValue1 = (saliencies[1] * localE) + eigenValue0;
                    minPointEigenvalues1.push_back(eigenValue1);
                    pointsIndices.push_back(point.vertIdx);
                    //cout << "min point eigenvalues " << eigenValue0 << " " << eigenValue1 << " " << localE << endl;

                }

            }
            return std_vector_to_py_list(pointsIndices);
        }


        py::list Volume::getExtremalMaxPoints() {
            std::vector<Point> points = (*currentExtremalData->getPoints());
            std::vector<Vertex> vertices = (*currentExtremalData->getVertices());
            std::vector<int> pointsIndices;
            for (int i = 0; i < points.size(); i++)
            {
                Point point = points[i];
                int type = point.type;
                float *saliencies = point.relativeSaliencies;
                float localI = point.localIntensity;
                float localE = point.firstEigenvalue;
                if (type == 1)
                {
                    maxPointSaliencies.push_back(saliencies[0]);
                    maxPointSaliencies.push_back(saliencies[1]);
                    maxPointSaliencies.push_back(saliencies[2]);
                    maxPointIntensities.push_back(localI);
                    maxPointEigenvalues.push_back(localE);
                    //pointsIndices.push_back(point.vertIdx);

                    float eigenValue0 = saliencies[2] * localE;
                    maxPointEigenvalues0.push_back(eigenValue0);
                    float eigenValue1 = (saliencies[1] * localE) + eigenValue0;
                    maxPointEigenvalues1.push_back(eigenValue1);
                    pointsIndices.push_back(point.vertIdx);
                    //cout << "max point eigenvalues " << eigenValue0 << " " << eigenValue1 << " " << localE << endl;
                }

            }
            return std_vector_to_py_list(pointsIndices);
        }

        py::list Volume::getExtremalSaddlePoints() {
            std::vector<Point> points = (*currentExtremalData->getPoints());
            std::vector<Vertex> vertices = (*currentExtremalData->getVertices());
            std::vector<int> pointsIndices;
            for (int i = 0; i < points.size(); i++)
            {
                Point point = points[i];
                int type = point.type;
                float *saliencies = point.relativeSaliencies;
                float localI = point.localIntensity;
                float localE = point.firstEigenvalue;
                if (type == 3)
                {
                    saddlePointSaliencies.push_back(saliencies[0]);
                    saddlePointSaliencies.push_back(saliencies[1]);
                    saddlePointSaliencies.push_back(saliencies[2]);
                    saddlePointIntensities.push_back(localI);
                    saddlePointEigenvalues.push_back(localE);
                    //pointsIndices.push_back(point.vertIdx);

                    float eigenValue0 = saliencies[2] * localE;
                    saddlePointEigenvalues0.push_back(eigenValue0);
                    float eigenValue1 = (saliencies[1] * localE) + eigenValue0;
                    saddlePointEigenvalues1.push_back(eigenValue1);
                    pointsIndices.push_back(point.vertIdx);
                    //cout << "saddle point eigenvalues " << eigenValue0 << " " << eigenValue1 << " " << localE << endl;

                }

            }
            return std_vector_to_py_list(pointsIndices);
        }


        /**
        bool hide = false;


        if(saliencyCheck) {
        float *saliencies = point.relativeSaliencies;
        if ( ( pointRatio * saliencies[2] < surfaceRatio * saliencies[0] ) ||
        ( pointRatio * saliencies[2] < curveRatio * saliencies[1] ) ) hide = true;
        }
        if(intensityCheck){
        float maxI = currentExtremalData ->maxIntensity;
        float minI = currentExtremalData ->minIntensity;
        float localI = point.localIntensity;
        if ( ( type == 1 ) && ( ( localI - minI ) / ( maxI - minI ) < maxGeo ) ) hide = true;
        if ( ( type == 2 ) && ( ( localI - minI ) / ( maxI - minI ) > minGeo ) ) hide = true;

        }
        if(eigenvalueCheck) {
        float maxE = *currentExtremalData ->maxEigenvalue;
        float minE = *currentExtremalData ->minEigenvalue;
        float localE = point.firstEigenvalue;
        if ( ( localE - minE ) / ( maxE - minE ) < eigenValue ) hide = true;

        }
        if (!hide && (type == 1)) {
        **/
        /* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
        *  Assuming the current volume has already been thresholded to 0/1
        */
        void Volume::curveSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(0.5f, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    float v = (float)grayvol->getDataAt(i, j, k);
                    if (v <= lowthr || v > highthr || svol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
            printf("Start erosion to %d...\n", wid);
#endif


            // Perform erosion 
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_HELIX
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_HELIX
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_HELIX <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                int cur, visited;
                for (cur = 1; cur < NOISE_DIS_HELIX; cur++)
                {
                    queue2->reset();
                    int count = 0;
                    visited = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            visited++;
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                visited++;
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }
                printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited);


#endif
                /* Commented out for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_HELIX
                if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented out for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_HELIX
                    // if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) ) 
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#else
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;


                        /* Adding ends */
                        // Move its neighboring unvisited node to queue2
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                // setDataAt( nx, ny, nz, curwid + 1 ) ;
                                queue2->prepend(nx, ny, nz);
                            }
                        }

                    }


                    /* Commented out for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_HELIX
                    ( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex2(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    if (queue2->getNumElements() > 0)
                    {
                        printf("*************************wierd here*************************\n");
                    }
                    break;
                }
            }

            // Remove all internal voxels (contained in manifold surfaces)
            queue2->reset();
            queue4->reset();
            ele = queue4->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;

                if (isPiercable(ox, oy, oz) == 1)  // hasCompleteSheet( ox, oy, oz ) == 1 ) //  
                {
                    queue2->prepend(ox, oy, oz);
                    //	setDataAt( ox, oy, oz, -1 ) ;
                }
                ele = queue4->remove();
            }

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0 && isPiercable(i, j, k)) //hasCompleteSheet(i,j,k) == 1) //  
                {
                    queue2->prepend(i, j, k);
                }
            }
            queue2->reset();
            ele = queue2->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;
                setDataAt(ox, oy, oz, -1);
                ele = queue2->remove();
            }


            // Finally, clean up
            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
            delete queue4;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }

        // Compute curve skeleton
        void Volume::curveSkeleton(float thr, Volume* svol)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (svol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }

            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
            printf("Start erosion to %d...\n", wid);
#endif


            // Perform erosion 
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_HELIX
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_HELIX
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_HELIX <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                int cur, visited;
                for (cur = 1; cur < NOISE_DIS_HELIX; cur++)
                {
                    queue2->reset();
                    int count = 0;
                    visited = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            visited++;
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                visited++;
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }
                printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited);


#endif
                /* Commented out for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_HELIX
                if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented out for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_HELIX
                    // if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) ) 
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#else
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented out for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_HELIX
                    ( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex2(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
            delete queue4;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }

        // Compute curve skeleton in 2D
        void Volume::curveSkeleton2D(float thr, Volume* svol)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (svol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 4; m++)
                        {
                            if (getDataAt(i + neighbor4[m][0], j + neighbor4[m][1], k) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
            printf("Start erosion to %d...\n", wid);
#endif


            // Perform erosion 
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_HELIX
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_HELIX
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_HELIX <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                int cur, visited;
                for (cur = 1; cur < NOISE_DIS_HELIX; cur++)
                {
                    queue2->reset();
                    int count = 0;
                    visited = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            visited++;
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                visited++;
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }
                printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited);


#endif
                /* Commented out for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_HELIX
                if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    //score = getNumNeighbor6( ox, oy, oz ) ;
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented out for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_HELIX
                    // if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) ) 
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#else
                    if (isHelixEnd(ox, oy, oz) || !isSimple(ox, oy, oz))
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 4; m++)
                    {
                        int nx = ox + neighbor4[m][0];
                        int ny = oy + neighbor4[m][1];
                        int nz = oz;
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented out for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_HELIX
                    ( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex2(nx, ny, nz);
                            //score = getNumNeighbor6( nx, ny, nz ) ;

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
            delete queue4;
        }

        // Compute minimal skeleton
        void Volume::skeleton(float thr, int off)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue = new GridQueue2();

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                setDataAt(i, j, k, 1);
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = 0;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
            printf("Start erosion to %d...\n", wid);
#endif

            // Perform erosion 
            gridQueueEle* ele;
            int ox, oy, oz;

            while (1)
            {
                wid++;

                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue is empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), wid);
#endif			

                // Rename queue2 to be queue, 
                // Clear queue2
                // From now on, queue2 holds nodes of next level
                delete queue;
                queue = queue2;
                queue2 = new GridQueue2();

                // Next, start queue iteration
                queue->reset();
                ele = queue->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    //				delete ele ;

                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        queue2->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    /*
                    else if ( ox == off || oy == off || oz == off ||
                    ox == getSizeX() - off - 1 || oy == getSizeY() - off - 1 || oz == getSizeZ() - off - 1 )
                    {
                    // Wall, don't erode, set to next layer
                    queue2->prepend( ox, oy, oz ) ;
                    numComplex ++ ;
                    }
                    */
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                setDataAt(nx, ny, nz, 1);
                                queue2->prepend(nx, ny, nz);
                            }
                        }

                    }

                    ele = queue->remove();
                }
#ifdef VERBOSE
                printf("Level %d: %d complex, %d simple\n", wid, numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
            delete queue;
            delete queue2;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }

        // Compute minimal skeleton
        void Volume::skeleton2(float thr, int off)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (i == off || j == off || k == off ||
                        i == getSizeX() - off - 1 || j == getSizeY() - off - 1 || k == getSizeZ() - off - 1)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                setDataAt(i, j, k, 1);
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif

            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;

            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }


            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Added for debugging */
                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue2->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue
                    /*
                    for ( i = -2 ; i < 3 ;i ++ )
                    for ( j = -2 ; j < 3 ; j ++ )
                    for ( k = -2 ; k < 3 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    if ( getDataAt( nx, ny, nz ) == curwid )
                    {
                    // Compute score
                    score = getNumPotComplex2( nx, ny, nz ) ;

                    if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
                    {
                    // printf("Update\n") ;
                    scrvol->setDataAt( nx, ny, nz, score ) ;
                    // Push to queue
                    gp = new gridPoint ;
                    gp->x = nx ;
                    gp->y = ny ;
                    gp->z = nz ;
                    // queue->add( gp, -score ) ;
                    queue->add( gp, score ) ;
                    }
                    }
                    }

                    */
                }
#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
            delete scrvol;
            delete queue;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }

        /* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
        *  Assuming the current volume has already been thresholded to 0/1
        */
        void Volume::pointSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol, Volume* hvol)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(0.5f, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    float v = (float)grayvol->getDataAt(i, j, k);
                    if (v <= lowthr || v > highthr || svol->getDataAt(i, j, k) > 0 || hvol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                setDataAt(i, j, k, 1);
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
#endif


            // Perform erosion 
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }


            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Added for debugging */
                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue2->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;
                        /* Adding ends */

                        // Move its neighboring unvisited node to queue2
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                setDataAt(nx, ny, nz, curwid + 1);
                                queue2->prepend(nx, ny, nz);
                            }
                        }

                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex2(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }
#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Remove all internal voxels (contained in manifold curves)
            queue2->reset();
            ele = queue2->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;

                if (hasCompleteHelix(ox, oy, oz) == 1)
                {
                    ele = queue2->getNext();
                }
                else
                {
                    ele = queue2->remove();
                }
            }

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0 && hasCompleteHelix(i, j, k) == 1)
                {
                    queue2->prepend(i, j, k);
                }
            }
            queue2->reset();
            ele = queue2->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;
                setDataAt(ox, oy, oz, -1);
                ele = queue2->remove();
            }

            // Finally, clean up
            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }


        // Compute minimal skeleton
        void Volume::skeleton(float thr, Volume* svol, Volume* hvol)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Next, apply convergent erosion 
            // by preserving: complex nodes, curve end-points, and sheet points

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (svol->getDataAt(i, j, k) > 0 || hvol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                setDataAt(i, j, k, 1);
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }


            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex2(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Added for debugging */
                    // Check simple 
                    if (!isSimple(ox, oy, oz))
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue2->prepend(ox, oy, oz);
                        numComplex++;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            setDataAt(nx, ny, nz, curwid + 1);
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex2(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif

                if (numSimple == 0)
                {
                    delete queue2;
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
            delete scrvol;
            delete queue;
            delete queue3;
        }


        // Apply helix erosion
        void Volume::erodeHelix()
        {
            erodeHelix(3);
        }

        //float Volume::getNorm(float v[3]) {
        //	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        //}

        void Volume::erodeHelix(int disthr)
        {
            int i, j, k;
            // First, threshold the volume
            //printf("Thresholding the volume to -1/0...\n") ;
            threshold(0.1f, -1, 0);

            /* Debug: remove faces */
            //Volume* facevol = markFaceEdge() ;
            /* End debugging */

            // Next, initialize the linked queue
            //printf("Initializing queue...\n") ;
            Volume * fvol = new Volume(getSizeX(), getSizeY(), getSizeZ());
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2** queues = new GridQueue2*[10000];

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (!hasCompleteHelix(i, j, k))
                        // if ( ! hasCompleteHelix( i, j, k, facevol ) )
                    {
                        queue2->prepend(i, j, k);
                        fvol->setDataAt(i, j, k, -1);
                    }
                }
            }
            //printf("Total %d nodes\n", queue2->getNumElements() ) ;

            // Start erosion
            gridQueueEle* ele;
            int dis = -1;
            while (queue2->getNumElements() > 0)
            {
                // First, set distance
                dis--;
                queues[-dis] = new GridQueue2();
                //printf("Distance transform to %d...", dis) ;
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, dis);
                    queues[-dis]->prepend(ele->x, ele->y, ele->z);
                }
                //printf("%d nodes\n", queues[-dis]->getNumElements()) ;

                // Next, find next layer
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ele->x + neighbor6[m][0];
                        int ny = ele->y + neighbor6[m][1];
                        int nz = ele->z + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            fvol->setDataAt(nx, ny, nz, dis);

                            if (!hasCompleteHelix(nx, ny, nz))
                                // if ( ! hasCompleteHelix( nx, ny, nz, facevol ) )
                            {
                                setDataAt(nx, ny, nz, 1);
                                queue3->prepend(nx, ny, nz);
                            }
                        }
                    }

                    ele = queue2->remove();
                }

                // Next, swap queues
                GridQueue2* temp = queue2;
                queue2 = queue3;
                queue3 = temp;
            }

            // Deal with closed rings
            dis--;
            queues[-dis] = new GridQueue2();
#ifdef VERBOSE
            printf("Detecting closed rings %d...", dis);
#endif
            int ftot = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0)
                {
                    setDataAt(i, j, k, dis);
                    queues[-dis]->prepend(i, j, k);
                    /*
                    int fval = (int) fvol->getDataAt( i, j, k ) ;
                    if ( fval == 0)
                    {
                    // queues[ -dis ]->prepend( i, j, k ) ;
                    }
                    else
                    {
                    setDataAt( i, j, k, fval - 1 ) ;
                    // queues[ -fval + 1 ]->prepend( i, j, k ) ;
                    }
                    */
                    ftot++;
                }
            }
#ifdef VERBOSE
            printf("%d nodes\n", ftot);
#endif


            // return ;

            /* Find local minimum: to help determine erosion level
            int cts[ 64 ] ;
            for ( i = 0 ; i <= - dis ; i ++ )
            {
            cts[ i ] = 0 ;
            }
            for ( i = 0 ; i < getSizeX() ; i ++ )
            for ( j = 0 ; j < getSizeY() ; j ++ )
            for ( k = 0 ; k < getSizeZ() ; k ++ )
            {
            double val = getDataAt( i, j, k ) ;
            if ( val < -1 )
            {
            int min = 1 ;
            for ( int m = 0 ; m < 6 ; m ++ )
            {
            int nx = i + neighbor6[m][0] ;
            int ny = j + neighbor6[m][1] ;
            int nz = k + neighbor6[m][2] ;
            if ( getDataAt( nx, ny, nz ) < val )
            {
            min = 0 ;
            break ;
            }
            }

            if ( min )
            {
            cts[ (int) - val ] ++ ;
            }
            }
            }

            for ( i = 2 ; i <= - dis ; i ++ )
            {
            printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
            }
            */

            // Dilate back
            // Starting from nodes with distance - 2 - disthr

            if (disthr + 2 > -dis)
            {
                disthr = -dis - 2;
            }
            int d;

            for (d = -dis; d > disthr + 1; d--)
            {
                queues[d]->reset();
                while ((ele = queues[d]->getNext()) != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, d);
                }
            }


            for (d = disthr + 1; d >= 2; d--)
            {
                //delete queue3 ;
                //queue3 = new GridQueue2( ) ;
                queues[d]->reset();
                ele = queues[d]->getNext();
                while (ele != NULL)
                {
                    int dilatable = 0;
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ele->x + neighbor6[m][0];
                        int ny = ele->y + neighbor6[m][1];
                        int nz = ele->z + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == d + 1)
                        {
                            dilatable = 1;
                            break;
                        }
                    }


                    if (!dilatable)
                    {
                        /*
                        setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
                        */

                        setDataAt(ele->x, ele->y, ele->z, -d + 1);
                        if (d > 2)
                        {
                            queues[d - 1]->prepend(ele->x, ele->y, ele->z);
                        }
                        ele = queues[d]->remove();
                    }
                    else
                    {
                        setDataAt(ele->x, ele->y, ele->z, d);
                        ele = queues[d]->getNext();
                    }

                }

                /* Debug: extract minimal set */
                while (1)
                {
                    int num = 0;
                    queues[d]->reset();
                    ele = queues[d]->getNext();
                    while (ele != NULL)
                    {
                        int critical = 0;
                        setDataAt(ele->x, ele->y, ele->z, -1);

                        for (i = -1; i < 2; i++)
                        {
                            for (j = -1; j < 2; j++)
                            {
                                for (k = -1; k < 2; k++)
                                {
                                    if (i != 0 && j != 0 && k != 0)
                                    {
                                        continue;
                                    }
                                    int nx = ele->x + i;
                                    int ny = ele->y + j;
                                    int nz = ele->z + k;
                                    if (getDataAt(nx, ny, nz) == d + 1 && !hasCompleteHelix(nx, ny, nz)) //, facevol ) )
                                    {
                                        critical = 1;
                                        break;
                                    }
                                }
                                if (critical)
                                {
                                    break;
                                }
                            }
                            if (critical)
                            {
                                break;
                            }
                        }

                        if (critical)
                        {
                            setDataAt(ele->x, ele->y, ele->z, d);
                            ele = queues[d]->getNext();
                        }
                        else
                        {
                            setDataAt(ele->x, ele->y, ele->z, -d + 1);
                            if (d > 2)
                            {
                                queues[d - 1]->prepend(ele->x, ele->y, ele->z);
                            }
                            ele = queues[d]->remove();
                            num++;
                        }

                    }

#ifdef VERBOSE
                    printf("Non-minimal: %d\n", num);
#endif

                    if (num == 0)
                    {
                        break;
                    }
                }


                /* End of debugging */

                /*
                queue3->reset() ;
                ele = queue3->getNext() ;
                while ( ele != NULL )
                {
                setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
                ele = queue3->remove() ;
                }
                */
            }

            // Finally, threshold the volume
#ifdef VERBOSE
            //printf("Thresholding the volume to 0/1...\n") ;
#endif
            //threshold( -1, 1, 0, 0 ) ;
            threshold(0, 0, 1);
            delete fvol;
            delete queue2;
            delete queue3;
            for (d = -dis; d >= 2; d--) {
                delete queues[d];
            }
            delete[] queues;

        }



        // Apply sheet erosion
        int Volume::erodeSheet()
        {
            return erodeSheet(3);
        }


        int Volume::erodeSheet(int disthr)
        {
            int i, j, k;
            // First, threshold the volume
            //printf("Thresholding the volume to -1/0...\n") ;
            threshold(0.1f, -1, 0);

            /* Debug: remove cells */
            Volume* facevol = markCellFace();
            /* End debugging */

            // Next, initialize the linked queue
            //printf("Initializing queue...\n") ;
            Volume * fvol = new Volume(getSizeX(), getSizeY(), getSizeZ());
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2** queues = new GridQueue2*[10000];

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (!hasCompleteSheet(i, j, k))
                        //if ( ! hasCompleteSheet( i, j, k, facevol ) )
                    {
                        queue2->prepend(i, j, k);
                        fvol->setDataAt(i, j, k, -1);
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
#endif

            // Start erosion
            gridQueueEle* ele;
            int dis = -1;
            while (queue2->getNumElements() > 0)
            {
                // First, set distance
                dis--;
                queues[-dis] = new GridQueue2();
                //printf("Distance transform to %d...", dis) ;
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, dis);
                    queues[-dis]->prepend(ele->x, ele->y, ele->z);
                }
                //printf("%d nodes\n", queues[-dis]->getNumElements()) ;

                // Next, find next layer
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    // for ( int m = 0 ; m < 6 ; m ++ )
                    for (int mx = -1; mx < 2; mx++)
                    for (int my = -1; my < 2; my++)
                    for (int mz = -1; mz < 2; mz++)
                    {
                        if (mx != 0 && my != 0 && mz != 0)
                        {
                            continue;
                        }
                        //int nx = ele->x + neighbor6[m][0] ;
                        //int ny = ele->y + neighbor6[m][1] ;
                        //int nz = ele->z + neighbor6[m][2] ;
                        int nx = ele->x + mx;
                        int ny = ele->y + my;
                        int nz = ele->z + mz;

                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            fvol->setDataAt(nx, ny, nz, dis);

                            if (!hasCompleteSheet(nx, ny, nz))
                                // if  ( ! hasCompleteSheet( nx, ny, nz, facevol ) )
                            {
                                setDataAt(nx, ny, nz, 1);
                                queue3->prepend(nx, ny, nz);
                            }
                        }
                    }

                    ele = queue2->remove();
                }

                // Next, swap queues
                GridQueue2* temp = queue2;
                queue2 = queue3;
                queue3 = temp;
            }

            /* Deal with closed rings */

            dis--;
            queues[-dis] = new GridQueue2();
#ifdef VERBOSE
            printf("Detecting closed rings %d...", dis);
#endif
            int ftot = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0)
                {
                    /*
                    int fval = (int) fvol->getDataAt( i, j, k ) ;
                    if ( fval == 0)
                    {
                    setDataAt( i, j, k, dis - 2 ) ;
                    // queues[ -dis ]->prepend( i, j, k ) ;
                    }
                    else
                    {
                    setDataAt( i, j, k, fval - 1 ) ;
                    queues[ -fval + 1 ]->prepend( i, j, k ) ;
                    }
                    */
                    setDataAt(i, j, k, dis);
                    queues[-dis]->prepend(i, j, k);

                    ftot++;
                }
            }
#ifdef VERBOSE
            printf("%d nodes\n", ftot);
#endif


            /* Find local minimum: to help determine erosion level
            int cts[ 64 ] ;
            for ( i = 0 ; i <= - dis ; i ++ )
            {
            cts[ i ] = 0 ;
            }
            for ( i = 0 ; i < getSizeX() ; i ++ )
            for ( j = 0 ; j < getSizeY() ; j ++ )
            for ( k = 0 ; k < getSizeZ() ; k ++ )
            {
            double val = getDataAt( i, j, k ) ;
            if ( val < -1 )
            {
            int min = 1 ;
            for ( int m = 0 ; m < 6 ; m ++ )
            {
            int nx = i + neighbor6[m][0] ;
            int ny = j + neighbor6[m][1] ;
            int nz = k + neighbor6[m][2] ;
            if ( getDataAt( nx, ny, nz ) < val )
            {
            min = 0 ;
            break ;
            }
            }

            if ( min )
            {
            cts[ (int) - val ] ++ ;
            }
            }
            }

            for ( i = 2 ; i <= - dis ; i ++ )
            {
            printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
            }
            */

            // return ;

            // Dilate back
            // Starting from nodes with distance - 2 - disthr
            int d;
            if (disthr + 2 > -dis)
            {
                disthr = -dis - 2;

            }
            for (d = -dis; d > disthr + 1; d--)
            {
                queues[d]->reset();
                while ((ele = queues[d]->getNext()) != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, d);
                }
            }

            for (d = disthr + 1; d >= 2; d--)
            {

                //delete queue3 ;
                //queue3 = new GridQueue2( ) ;
                queues[d]->reset();
                ele = queues[d]->getNext();
                while (ele != NULL)
                {
                    int dilatable = 0;
                    // for ( int m = 0 ; m < 6 ; m ++ )
                    /*
                    for ( int mx = -1 ; mx < 2 ; mx ++ )
                    for ( int my = -1 ; my < 2 ; my ++ )
                    for ( int mz = -1 ; mz < 2 ; mz ++ )
                    {
                    if ( mx == 0 || my == 0 || mz == 0 )
                    {
                    int nx = ele->x + mx ; // neighbor6[m][0] ;
                    int ny = ele->y + my ; // neighbor6[m][1] ;
                    int nz = ele->z + mz ; // neighbor6[m][2] ;
                    if ( getDataAt( nx, ny, nz ) == - d - 1 )
                    {
                    dilatable = 1 ;
                    break ;
                    }
                    }
                    }
                    */
                    for (i = 0; i < 12; i++)
                    {
                        int flag = 1, flag2 = 0;
                        for (j = 0; j < 4; j++)
                        {
                            int nx = ele->x + sheetNeighbor[i][j][0];
                            int ny = ele->y + sheetNeighbor[i][j][1];
                            int nz = ele->z + sheetNeighbor[i][j][2];

                            double val = getDataAt(nx, ny, nz);

                            if (val > -d && val < 0)
                            {
                                flag = 0;
                                break;
                            }
                            if (val == d + 1)
                            {
                                flag2++;
                            }
                        }

                        if (flag && flag2)
                        {
                            dilatable = 1;
                            break;
                        }
                    }

                    if (!dilatable)
                    {
                        // setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
                        // queue3->prepend( ele->x, ele->y, ele->z ) ;

                        setDataAt(ele->x, ele->y, ele->z, -d + 1);
                        if (d > 2)
                        {
                            queues[d - 1]->prepend(ele->x, ele->y, ele->z);
                        }
                        ele = queues[d]->remove();
                    }
                    else
                    {
                        setDataAt(ele->x, ele->y, ele->z, d);
                        ele = queues[d]->getNext();
                    }
                }

                /* Debug: extract minimal set */
                while (1)
                {
                    int num = 0;
                    queues[d]->reset();
                    ele = queues[d]->getNext();
                    while (ele != NULL)
                    {
                        int critical = 0;
                        setDataAt(ele->x, ele->y, ele->z, -1);

                        for (i = -1; i < 2; i++)
                        {
                            for (j = -1; j < 2; j++)
                            {
                                for (k = -1; k < 2; k++)
                                {
                                    if (i != 0 && j != 0 && k != 0)
                                    {
                                        continue;
                                    }
                                    int nx = ele->x + i;
                                    int ny = ele->y + j;
                                    int nz = ele->z + k;
                                    // if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteSheet( nx,ny,nz, facevol ) )
                                    if (getDataAt(nx, ny, nz) == d + 1 && !hasCompleteSheet(nx, ny, nz))
                                    {
                                        critical = 1;
                                        break;
                                    }
                                }
                                if (critical)
                                {
                                    break;
                                }
                            }
                            if (critical)
                            {
                                break;
                            }
                        }

                        if (critical)
                        {
                            setDataAt(ele->x, ele->y, ele->z, d);
                            ele = queues[d]->getNext();
                        }
                        else
                        {
                            setDataAt(ele->x, ele->y, ele->z, -d + 1);
                            if (d > 2)
                            {
                                queues[d - 1]->prepend(ele->x, ele->y, ele->z);
                            }
                            ele = queues[d]->remove();
                            num++;
                        }

                    }
#ifdef VERBOSE
                    printf("Non-minimal: %d\n", num);
#endif

                    if (num == 0)
                    {
                        break;
                    }
                }


                /* End of debugging */

                /*
                queue3->reset() ;
                ele = queue3->getNext() ;
                while ( ele != NULL )
                {
                setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
                ele = queue3->remove() ;
                }
                */
            }


            // Finally, threshold the volume
#ifdef VERBOSE
            //printf("Thresholding the volume to 0/1...\n") ;
#endif
            //threshold( -1, 1, 0, 0 ) ;
            threshold(0, 0, 1);

            delete facevol;
            delete fvol;
            delete queue2;
            delete queue3;
            for (d = -dis; d >= 2; d--) {
                delete queues[d];
            }
            delete[] queues;

            return -dis - 1;
        }

        void Volume::erodeSheetOld(int disthr)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(0.1f, -1, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2** queues = new GridQueue2*[64];

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (!hasCompleteSheet(i, j, k))
                    {
                        queue2->prepend(i, j, k);
                    }
                }
            }
            printf("Total %d nodes\n", queue2->getNumElements());

            // Start erosion
            gridQueueEle* ele;
            int dis = -1;
            while (queue2->getNumElements() > 0)
            {
                // First, set distance
                dis--;
                queues[-dis] = new GridQueue2();
                printf("Distance transform to %d...", dis);
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, dis);
                    queues[-dis]->prepend(ele->x, ele->y, ele->z);
                }
                printf("%d nodes\n", queues[-dis]->getNumElements());

                // Next, find next layer
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ele->x + neighbor6[m][0];
                        int ny = ele->y + neighbor6[m][1];
                        int nz = ele->z + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0 && !hasCompleteSheet(nx, ny, nz))
                        {
                            setDataAt(nx, ny, nz, 1);
                            queue3->prepend(nx, ny, nz);
                        }
                    }

                    ele = queue2->remove();
                }

                // Next, swap queues
                GridQueue2* temp = queue2;
                queue2 = queue3;
                queue3 = temp;
            }

            /* Deal with closed rings
            for ( i = 0 ; i < getSizeX() ; i ++ )
            for ( j = 0 ; j < getSizeY() ; j ++ )
            for ( k = 0 ; k < getSizeZ() ; k ++ )
            {
            if ( getDataAt( i, j, k ) == 0 )
            {
            setDataAt( i, j, k, dis - 1 ) ;
            }
            }
            */

            /* Find local minimum: to help determine erosion level */
            int cts[64];
            for (i = 0; i <= -dis; i++)
            {
                cts[i] = 0;
            }
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                double val = getDataAt(i, j, k);
                if (val < -1)
                {
                    int min = 1;
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = i + neighbor6[m][0];
                        int ny = j + neighbor6[m][1];
                        int nz = k + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) < val)
                        {
                            min = 0;
                            break;
                        }
                    }

                    if (min)
                    {
                        cts[(int)-val] ++;
                    }
                }
            }

            for (i = 2; i <= -dis; i++)
            {
                printf("Local minima: %d with %d nodes.\n", -i, cts[i]);
            }


            // return ;

            // Dilate back
            // Starting from nodes with distance - 2 - disthr

            for (int d = disthr + 1; d >= 2; d--)
            {
                delete queue3;
                queue3 = new GridQueue2();
                queues[d]->reset();
                while ((ele = queues[d]->getNext()) != NULL)
                {
                    int dilatable = 0;
                    // for ( int m = 0 ; m < 6 ; m ++ )
                    /*
                    for ( int mx = -1 ; mx < 2 ; mx ++ )
                    for ( int my = -1 ; my < 2 ; my ++ )
                    for ( int mz = -1 ; mz < 2 ; mz ++ )
                    {
                    if ( mx == 0 || my == 0 || mz == 0 )
                    {
                    int nx = ele->x + mx ; // neighbor6[m][0] ;
                    int ny = ele->y + my ; // neighbor6[m][1] ;
                    int nz = ele->z + mz ; // neighbor6[m][2] ;
                    if ( getDataAt( nx, ny, nz ) == - d - 1 )
                    {
                    dilatable = 1 ;
                    break ;
                    }
                    }
                    }
                    */
                    for (i = 0; i < 12; i++)
                    {
                        int flag = 1, flag2 = 0;
                        for (j = 0; j < 4; j++)
                        {
                            int nx = ele->x + sheetNeighbor[i][j][0];
                            int ny = ele->y + sheetNeighbor[i][j][1];
                            int nz = ele->z + sheetNeighbor[i][j][2];

                            double val = getDataAt(nx, ny, nz);

                            if (val > -d)
                            {
                                flag = 0;
                                break;
                            }
                            if (val == -d - 1)
                            {
                                flag2++;
                            }
                        }

                        if (flag && flag2)
                        {
                            dilatable = 1;
                            break;
                        }
                    }

                    if (!dilatable)
                    {
                        // setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
                        queue3->prepend(ele->x, ele->y, ele->z);
                        /*
                        setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
                        if ( d > 2 )
                        {
                        queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
                        }
                        */
                    }
                }

                queue3->reset();
                ele = queue3->getNext();
                while (ele != NULL)
                {
                    setDataAt(ele->x, ele->y, ele->z, -1);
                    ele = queue3->remove();
                }
            }

            // Finally, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(-1, 1, 0, 0);
        }



        void Volume::addNoise(float thr, float pos)
        {
            int i, j, k;
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);
            Volume* tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (tvol->getDataAt(i, j, k) >= 0 && isSimple(i, j, k))
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (tvol->getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            if (rand() < RAND_MAX * pos)
                            {
                                setDataAt(i, j, k, thr - 1);
                            }
                            break;
                        }
                    }
                }
            }

        }

        /************************************************************************/
        // A sequential thinning method for extracting 6-connected skeletons                
        // Properties: medial, topology preserving, shape preserving, noise resistance!
        // @param thr: threshold
        // @param type: 0 for curve preserving, 1 for surface preserving
        // @param noise: 0 for no-noise-reduction, n for level-n noise reduction
        /************************************************************************/
        void Volume::sequentialSkeleton(float thr, int type, int noise)
        {
            int i, j, k, m;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -1/0...\n");
#endif
            threshold(thr, -1, 0);

            // Values in the volume:
            // -1:	outside
            // 0:	unvisited
            // 1:	next layer
            // 2:	this layer - non complex nodes
            // 3:	complex nodes

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            gridQueueEle* ele;

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    for (m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            setDataAt(i, j, k, 1);
                            queue2->prepend(i, j, k);
                            break;
                        }
                    }
                }
            }
            // printf("Total %d nodes\n", queue2->getNumElements() ) ;


            // Perform erosion
            int dis = 0;
            int ox, oy, oz;
            int nx, ny, nz;
            while (queue2->getNumElements() > 0)
            {
                dis++;
                printf("At distance %d, there are %d nodes.\n", dis, queue2->getNumElements());

                // At the beginning of each iteration:
                // queue2 holds all nodes dis away from the background in city block metric
                // queue3 is empty

                // First, find out next layer and put into queue3
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    setDataAt(ox, oy, oz, 2);
                    for (m = 0; m < 6; m++)
                    {
                        nx = ox + neighbor6[m][0];
                        ny = oy + neighbor6[m][1];
                        nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            setDataAt(nx, ny, nz, 1);
                            queue3->prepend(nx, ny, nz);
                        }
                    }
                }

                // Next, find out complex nodes in this layer and remove from queue2
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if ((!isSimple(ox, oy, oz)) ||
                        (((type == 0 && isEdgeEnd(ox, oy, oz)) ||
                        (type == 1 && isFaceEnd(ox, oy, oz))) && !isNoise(ox, oy, oz, noise)))
                    {
                        setDataAt(ox, oy, oz, 3);
                        ele = queue2->remove();
                    }
                    else
                    {
                        ele = queue2->getNext();
                    }
                }

                // Now, the main loop: delete non-complex nodes until none
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    queue2->remove();
                    queue2->reset();

                    if (getDataAt(ox, oy, oz) != 2)
                    {
                        continue;
                    }

                    if ((((type == 0 && isEdgeEnd(ox, oy, oz)) ||
                        (type == 1 && isFaceEnd(ox, oy, oz))) && !isNoise(ox, oy, oz, noise)))
                    {
                        setDataAt(ox, oy, oz, 3);
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                    }

                    for (i = -1; i < 2; i++)
                    for (j = -1; j < 2; j++)
                    for (k = -1; k < 2; k++)
                    {
                        nx = ox + i;
                        ny = oy + j;
                        nz = oz + k;
                        int val = (int)(getDataAt(nx, ny, nz));
                        if (val > 1)
                        {
                            int complex = 0;
                            if ((!isSimple(nx, ny, nz)) ||
                                (((type == 0 && isEdgeEnd(nx, ny, nz)) ||
                                (type == 1 && isFaceEnd(nx, ny, nz))) && !isNoise(nx, ny, nz, noise)))
                            {
                                complex = 1;
                            }

                            if (val == 2 && complex)
                            {
                                // A non-complex node becomes complex
                                setDataAt(nx, ny, nz, 3);
                            }
                            else if (val == 3 && !complex)
                            {
                                // A complex node turns non-complex
                                setDataAt(nx, ny, nz, 2);
                                queue2->prepend(nx, ny, nz);
                            }

                        }
                    }

                    queue2->reset();
                }

                // Finally, swap queue3 and queue2
                GridQueue2* temp = queue3;
                queue3 = queue2;
                queue2 = temp;
            }


            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }

        void Volume::dumbsurfaceSkeleton(float thr)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            gridQueueEle* ele;


            while (1)
            {
                int n = 0;

                queue2->reset();
                for (i = 0; i < getSizeX(); i++)
                for (j = 0; j < getSizeY(); j++)
                for (k = 0; k < getSizeZ(); k++)
                {
                    if (getDataAt(i, j, k) == 0)
                    {
                        {
                            for (int m = 0; m < 6; m++)
                            {
                                if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                                {
                                    queue2->prepend(i, j, k);
                                    break;
                                }
                            }
                        }
                    }
                }

                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    int ox = ele->x;
                    int oy = ele->y;
                    int oz = ele->z;
                    if (isSimple(ox, oy, oz) && hasCompleteSheet(ox, oy, oz) == 1)
                    {
                        setDataAt(ox, oy, oz, -1);
                        n++;
                    }

                    ele = queue2->remove();
                }

                if (n == 0)
                {
                    break;
                }

                printf("%d simple nodes found.\n", n);
            }


            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }

        /* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
        *  Assuming the current volume has already been thresholded to 0/1
        */
        void Volume::surfaceSkeleton(Volume* grayvol, float lowthr, float highthr) {
            int i, j, k;
            threshold(0.5f, -MAX_ERODE, 0);

            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();

            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);
            int ct = 0;

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    float v = (float)grayvol->getDataAt(i, j, k);
                    if (v > highthr || v <= lowthr)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        ct++;

                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }


            // Perform erosion 
            int wid = MAX_ERODE;
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_SHEET
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {


                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif

                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_SHEET
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_SHEET <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                for (int cur = 1; cur < NOISE_DIS_SHEET; cur++)
                {
                    queue2->reset();
                    int count = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }


#endif

                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }


                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                int nowComplex = 0;

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;

                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }


#ifndef NOISE_DIS_SHEET
                    // if ( hasFeatureFace( ox, oy, oz ) )
                    if ((!isSimple(ox, oy, oz)) || isSheetEnd(ox, oy, oz))
                        // if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
#else
                    // if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
                    if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz, noisevol) || isHelixEnd(ox, oy, oz, noisevol))
                        // if ( isBertrandEndPoint( ox, oy, oz ) ) 
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;

                        nowComplex = 1;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                // setDataAt( nx, ny, nz, curwid + 1 ) ;
                                queue2->prepend(nx, ny, nz);
                            }
                        }

                        if (nowComplex)
                        {

                            // printf(" %d\n", score);
                        }
                    }

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

                if (numSimple == 0)
                {
                    if (queue2->getNumElements() > 0)
                    {
                        printf("*************************wierd**********************\n");
                    }
                    break;
                }
            }

            // Remove all internal voxels (contained in cells)

            queue4->reset();
            ele = queue4->getNext();
            while (ele != NULL)
            {
                ele = queue4->remove();
            }

            queue2->reset();
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0 && isInternal2(i, j, k) == 1)
                {
                    queue2->prepend(i, j, k);
                }
            }
            queue2->reset();
            ele = queue2->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;
                setDataAt(ox, oy, oz, -1);
                ele = queue2->remove();
            }



            // Finally, clean up
            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
            delete queue4;
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }

        void Volume::surfaceSkeleton(float thr)
        {
            int i, j, k;
            // First, threshold the volume
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
            printf("Initializing queue...\n");
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();

            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_SHEET
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 and queue4 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_SHEET
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_SHEET <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                for (int cur = 1; cur < NOISE_DIS_SHEET; cur++)
                {
                    queue2->reset();
                    int count = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }


#endif

                /* Commented for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_SHEET
                if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                int nowComplex = 0;

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    // printf("%d\n", score);
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_SHEET
                    // if ( hasFeatureFace( ox, oy, oz ) )
                    if ((!isSimple(ox, oy, oz)) || isSheetEnd(ox, oy, oz))
                        // if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
#else
                    // if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
                    if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz, noisevol) || isHelixEnd(ox, oy, oz, noisevol))
                        // if ( isBertrandEndPoint( ox, oy, oz ) ) 
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;

                        nowComplex = 1;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                        if (nowComplex)
                        {

                            // printf("Error: %d\n", score);
                        }
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_SHEET
                    ( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

                printf("%d complex, %d simple\n", numComplex, numSimple);

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
            printf("Thresholding the volume to 0/1...\n");
            threshold(0, 0, 1);

        }

        void Volume::surfaceSkeleton(float thr, Volume* svol)
        {
            int i, j, k;
            // First, threshold the volume
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
            printf("Initializing queue...\n");
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();

            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    if (svol->getDataAt(i, j, k) > 0)
                    {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());


            // Perform erosion 
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_SHEET
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 and queue4 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_SHEET
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_SHEET <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                for (int cur = 1; cur < NOISE_DIS_SHEET; cur++)
                {
                    queue2->reset();
                    int count = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }


#endif

                /* Commented for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_SHEET
                if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                int nowComplex = 0;

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    // printf("%d\n", score);
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_SHEET
                    // if ( hasFeatureFace( ox, oy, oz ) )
                    if ((!isSimple(ox, oy, oz)) || isSheetEnd(ox, oy, oz))
                        // if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
#else
                    // if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
                    if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz, noisevol) || isHelixEnd(ox, oy, oz, noisevol))
                        // if ( isBertrandEndPoint( ox, oy, oz ) ) 
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;

                        nowComplex = 1;
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                        if (nowComplex)
                        {

                            // printf("Error: %d\n", score);
                        }
                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_SHEET
                    ( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

                printf("%d complex, %d simple\n", numComplex, numSimple);

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
            printf("Thresholding the volume to 0/1...\n");
            threshold(0, 0, 1);

        }

        void Volume::surfaceSkeletonOld(float thr)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();

            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0)
                {
                    {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
#endif

            // Perform erosion 
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Start erosion to %d...\n", wid);
#endif
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_SHEET
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 and queue4 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif			

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_SHEET
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_SHEET <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                for (int cur = 1; cur < NOISE_DIS_SHEET; cur++)
                {
                    queue2->reset();
                    int count = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }


#endif

                /* Commented for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_SHEET
                if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();

                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    // printf("%d\n", score);
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_SHEET
                    // if ( hasFeatureFace( ox, oy, oz ) )
                    if ((!isSimple(ox, oy, oz)) || isSheetEnd(ox, oy, oz))
                        // if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
#else
                    // if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
                    if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz, noisevol) || isHelixEnd(ox, oy, oz, noisevol))
                        // if ( isBertrandEndPoint( ox, oy, oz ) ) 
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;

                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_SHEET
                    ( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }
#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif			

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
            delete queue;

        }

        void Volume::surfaceSkeletonPres(float thr, Volume * preserve)
        {
            int i, j, k;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Next, initialize the linked queue
#ifdef VERBOSE
            printf("Initializing queue...\n");
#endif
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();

            PriorityQueue <gridPoint, int> * queue = new PriorityQueue <gridPoint, int>(MAX_QUEUELEN);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= 0) {
                    if (preserve->getDataAt(i, j, k) > 0) {
                        setDataAt(i, j, k, MAX_ERODE);
                    }
                    else {
                        for (int m = 0; m < 6; m++)
                        {
                            if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                            {
                                // setDataAt( i, j, k, 1 ) ;
                                queue2->prepend(i, j, k);
                                break;
                            }
                        }
                    }
                }
            }
            int wid = MAX_ERODE;
#ifdef VERBOSE
            printf("Total %d nodes\n", queue2->getNumElements());
            printf("Start erosion to %d...\n", wid);
#endif


            // Perform erosion 
            gridQueueEle* ele;
            gridPoint* gp;
            int ox, oy, oz;
            int score;
            Volume* scrvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ());
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                scrvol->setDataAt(i, -1);
            }

#ifdef  NOISE_DIS_SHEET
            Volume* noisevol = new Volume(getSizeX(), getSizeY(), getSizeZ());
#endif

            for (int curwid = 1; curwid <= wid; curwid++)
            {
                // At the start of each iteration, 
                // queue2 and queue4 holds all the nodes for this layer
                // queue3 and queue are empty

                int numComplex = 0, numSimple = 0;
#ifdef VERBOSE
                printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid);
#endif			

                /*
                We first need to assign curwid + 1 to every node in this layer
                */
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    if (getDataAt(ox, oy, oz) == curwid)
                    {
                        ele = queue2->remove();
                    }
                    else
                    {
                        setDataAt(ox, oy, oz, curwid);
                        ele = queue2->getNext();
                    }
                }
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    queue2->prepend(ox, oy, oz);
                    ele = queue4->remove();
                }

                // Now queue2 holds all the nodes for this layer

#ifdef NOISE_DIS_SHEET
                /* Extra step: classify nodes in queue2 into noise and non-noise nodes */
                queue2->reset();

                // First run
                int flag = 0;
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    if (NOISE_DIS_SHEET <= 1)
                    {
                        noisevol->setDataAt(ox, oy, oz, 0);
                    }
                    else
                    {
                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) == 0)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                flag = 1;
                                break;
                            }
                        }
                        if (!flag)
                        {
                            noisevol->setDataAt(ox, oy, oz, 0);
                        }
                    }
                }

                for (int cur = 1; cur < NOISE_DIS_SHEET; cur++)
                {
                    queue2->reset();
                    int count = 0;

                    while ((ele = queue2->getNext()) != NULL)
                    {
                        ox = ele->x;
                        oy = ele->y;
                        oz = ele->z;

                        if (noisevol->getDataAt(ox, oy, oz) == 1)
                        {
                            continue;
                        }

                        flag = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (getDataAt(nx, ny, nz) > 0 && noisevol->getDataAt(nx, ny, nz) == 1)
                            {
                                noisevol->setDataAt(ox, oy, oz, 1);
                                count++;
                                break;
                            }
                        }
                    }

                    if (count == 0)
                    {
                        break;
                    }
                }


#endif

                /* Commented for debugging

                // First,
                // check for complex nodes in queue2
                // move them from queue2 to queue3
                queue2->reset() ;
                ele = queue2->getNext() ;
                while ( ele != NULL )
                {
                ox = ele->x ;
                oy = ele->y ;
                oz = ele->z ;

                // Check simple
                #ifndef NOISE_DIS_SHEET
                if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
                #else
                if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
                #endif
                {
                // Complex, set to next layer
                setDataAt( ox, oy, oz, curwid + 1 ) ;
                queue3->prepend( ox, oy, oz ) ;
                ele = queue2->remove() ;

                numComplex ++ ;
                }
                else
                {
                ele = queue2->getNext() ;
                }
                }
                */


                // Next,
                // Compute score for each node left in queue2
                // move them into priority queue
                queue2->reset();
                ele = queue2->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    // Compute score
                    score = getNumPotComplex(ox, oy, oz);
                    scrvol->setDataAt(ox, oy, oz, score);

                    // Push to queue
                    gp = new gridPoint;
                    gp->x = ox;
                    gp->y = oy;
                    gp->z = oz;
                    // queue->add( gp, -score ) ;
                    queue->add(gp, score);

                    ele = queue2->remove();
                }

                // Rename queue3 to be queue2, 
                // Clear queue3
                // From now on, queue2 holds nodes of next level
                delete queue2;
                queue2 = queue3;
                queue3 = new GridQueue2();


                // Next, start priority queue iteration
                while (!queue->isEmpty())
                {
                    // Retrieve the node with the highest score
                    queue->remove(gp, score);
                    ox = gp->x;
                    oy = gp->y;
                    oz = gp->z;
                    delete gp;
                    // printf("%d\n", score);
                    //				score = -score ;

                    // Ignore the node 
                    // if it has been processed before
                    // or it has an updated score
                    if (getDataAt(ox, oy, oz) != curwid || (int)scrvol->getDataAt(ox, oy, oz) != score)
                    {
                        continue;
                    }

                    /* Commented for debugging

                    // Remove this simple node
                    setDataAt( ox, oy, oz, -1 ) ;
                    numSimple ++ ;
                    // printf("Highest score: %d\n", score) ;
                    */

                    /* Added for debugging */
                    // Check simple 
#ifndef NOISE_DIS_SHEET
                    // if ( hasFeatureFace( ox, oy, oz ) )
                    if ((!isSimple(ox, oy, oz)) || isSheetEnd(ox, oy, oz))
                        // if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
#else
                    // if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
                    if (!isSimple(ox, oy, oz) || isSheetEnd(ox, oy, oz, noisevol) || isHelixEnd(ox, oy, oz, noisevol))
                        // if ( isBertrandEndPoint( ox, oy, oz ) ) 
#endif
                    {
                        // Complex, set to next layer
                        setDataAt(ox, oy, oz, curwid + 1);
                        queue4->prepend(ox, oy, oz);
                        numComplex++;

                    }
                    else
                    {
                        setDataAt(ox, oy, oz, -1);
                        numSimple++;

                    }
                    /* Adding ends */

                    // Move its neighboring unvisited node to queue2
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) == 0)
                        {
                            // setDataAt( nx, ny, nz, curwid + 1 ) ;
                            queue2->prepend(nx, ny, nz);
                        }
                    }

                    /* Commented for debugging

                    // Find complex nodes in its 3x3 neighborhood
                    // move them to queue2
                    for ( i = -1 ; i < 2 ; i ++ )
                    for ( j = -1 ; j < 2 ; j ++ )
                    for ( k = -1 ; k < 2 ; k ++ )
                    {
                    int nx = ox + i ;
                    int ny = oy + j ;
                    int nz = oz + k ;

                    // Check simple
                    if ( getDataAt( nx, ny, nz ) == curwid &&
                    // ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
                    #ifndef NOISE_DIS_SHEET
                    ( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
                    #else
                    ( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
                    #endif

                    {
                    // Complex, set to next layer
                    setDataAt( nx, ny, nz, curwid + 1 ) ;
                    queue2->prepend( nx, ny, nz ) ;
                    numComplex ++ ;
                    }
                    }
                    */

                    // Update scores for nodes in its 5x5 neighborhood
                    // insert them back into priority queue

                    for (i = -2; i < 3; i++)
                    for (j = -2; j < 3; j++)
                    for (k = -2; k < 3; k++)
                    {
                        int nx = ox + i;
                        int ny = oy + j;
                        int nz = oz + k;

                        if (getDataAt(nx, ny, nz) == curwid)
                        {
                            // Compute score
                            score = getNumPotComplex(nx, ny, nz);

                            if (score != (int)scrvol->getDataAt(nx, ny, nz))
                            {
                                // printf("Update\n") ;
                                scrvol->setDataAt(nx, ny, nz, score);
                                // Push to queue
                                gp = new gridPoint;
                                gp->x = nx;
                                gp->y = ny;
                                gp->z = nz;
                                // queue->add( gp, -score ) ;
                                queue->add(gp, score);
                            }
                        }
                    }


                }

#ifdef VERBOSE
                printf("%d complex, %d simple\n", numComplex, numSimple);
#endif			

                if (numSimple == 0)
                {
                    break;
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

            delete scrvol;
            delete queue;
            delete queue2;
            delete queue3;
            delete queue4;

        }

        /* Bertrand's parallel 6-connected surface thinning */
        void Volume::bertrandSurfaceSkeleton2(float thr)
        {
            int i, j, k;

            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);


            // Initialize queues
            printf("Initializing queues...\n");
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();
            Volume* fvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ(), 0);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 0)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) < 0)
                        {
                            fvol->setDataAt(i, j, k, 1);
                            queue2->prepend(i, j, k);
                            break;
                        }
                    }
                }
            }
            printf("Total %d nodes\n", queue2->getNumElements());

            // Start iteration
            int it = 0;
            gridQueueEle* ele;
            while (++it)
            {
                printf("Iteration %d... nodes in queue: %d...", it, queue2->getNumElements());

                // queue2 holds all boundary nodes

                int deleted = 0;
                for (i = 0; i < 6; i++)
                {
                    // At the beginning of each iteration, 
                    // queue2 holds all remaining boundary nodes
                    // queue3 is a deletable array, starting empty
                    // queue4 holds the candidates for next layer

                    queue2->reset();
                    ele = queue2->getNext();

                    // For each sub-iteration, go through queue2 first and find points to delete
                    while (ele != NULL)
                    {
                        int ox = ele->x;
                        int oy = ele->y;
                        int oz = ele->z;

                        if (isBertrandBorder(ox, oy, oz, i))
                        {
                            if (!isBertrandEndPoint(ox, oy, oz))
                            {
                                // Move this node to queue3
                                ele = queue2->remove();
                                queue3->prepend(ox, oy, oz);

                                // Move its neighboring unvisited node to queue4
                                for (int m = 0; m < 6; m++)
                                {
                                    int nx = ox + neighbor6[m][0];
                                    int ny = oy + neighbor6[m][1];
                                    int nz = oz + neighbor6[m][2];
                                    if (fvol->getDataAt(nx, ny, nz) == 0)
                                    {
                                        fvol->setDataAt(nx, ny, nz, 1);
                                        queue4->prepend(nx, ny, nz);
                                    }
                                }
                            }
                            else
                            {
                                ele = queue2->getNext();
                            }
                        }
                        else
                        {
                            ele = queue2->getNext();
                        }
                    }

                    // Now, queue2 holds all remaining nodes, 
                    // queue3 holds all deletable nodes, 
                    // and queue4 holds nodes to be added to the next layer

                    // Simultaneous deletion
                    if (queue3->getNumElements() == 0)
                    {
                        // break ;
                    }
                    else
                    {
                        queue3->reset();
                        ele = queue3->getNext();
                        while (ele != NULL)
                        {
                            setDataAt(ele->x, ele->y, ele->z, -1);
                            ele = queue3->remove();
                            deleted++;
                        }
                    }

                    // return ;
                }

                // After all sub-iterations
                // Move all queue4 nodes to queue2
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    queue2->prepend(ele->x, ele->y, ele->z);
                    ele = queue4->remove();
                }

                if (deleted == 0)
                {
                    printf("No more deletable nodes.\n");
                    break;
                }
                else
                {
                    printf("Deleted: %d\n", deleted);
                }
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);

        }



        void Volume::bertrandSurfaceSkeleton(float thr)
        {
            int i, j, k;

            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to -MAX_ERODE/0...\n");
#endif
            threshold(thr, -MAX_ERODE, 0);

            // Start iteration
            GridQueue2* queue2 = new GridQueue2();
            int dir = 0, ct = 0;
            while (++ct)
            {
                printf("Round %d, direction %d...", ct, dir);
                queue2->reset();

                for (i = 0; i < getSizeX(); i++)
                for (j = 0; j < getSizeY(); j++)
                for (k = 0; k < getSizeZ(); k++)
                {
                    if (getDataAt(i, j, k) >= 0)
                    {
                        if (isBertrandBorder(i, j, k, dir))
                        {
                            if (!isBertrandEndPoint(i, j, k))
                            {
                                queue2->prepend(i, j, k);
                            }
                        }
                    }

                }


                if (queue2->getNumElements() == 0)
                {
                    printf("Done.\n");
                    break;
                }
                else
                {
                    queue2->reset();
                    printf("%d nodes deleted.\n", queue2->getNumElements());
                    gridQueueEle* ele = queue2->getNext();
                    while (ele != NULL)
                    {
                        setDataAt(ele->x, ele->y, ele->z, -1);
                        ele = queue2->remove();
                    }
                }
                dir = (dir + 1) % 6;
            }

            // Finally, clean up
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(0, 0, 1);
        }


        /* Palagyi's parallel surface thinning */
        void Volume::palagyiSurfaceSkeleton(float thr)
        {
            int i, j, k;

            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(thr, 0, 1);

            // Next, initialize the templates
            printf("Initializing surface endpoints templates...\n");
            ThinningTemplate* US[6];

            int b0[] = { 12, 13 };
            int w0[] = { 2, 5, 8, 11, 14, 17, 20, 23, 26 };
            int ob0[] = { 10, 16 };
            int ob20[] = { 4, 22 };
            int nu[] = { 0 };
            US[0] = new ThinningTemplate(b0, 2, w0, 9, ob0, 2, ob20, 2, nu, 0, nu, 0);
            US[1] = new ThinningTemplate(US[0], 0, 1);
            /*
            int b01[] = {13,16} ;
            int w01[] = {0,1,2,9,10,11,18,19,20} ;
            int ob01[] = {12, 14} ;
            int ob201[] = {4, 22} ;
            US[1] = new ThinningTemplate( b01, 2, w01, 9, ob01, 2, ob201, 2, nu, 0, nu, 0 ) ;
            */

            int b1[] = { 12, 13, 16, 22 };
            int w1[] = { 2, 10, 11, 14 };
            int ow[] = { 1, 5 };
            US[2] = new ThinningTemplate(b1, 4, w1, 4, nu, 0, nu, 0, ow, 2, nu, 0);
            US[3] = new ThinningTemplate(US[2], 0);

            int b2[] = { 2, 12, 13, 16, 22 };
            int w2[] = { 10, 11, 14 };
            int op[] = { 1, 5 };
            US[4] = new ThinningTemplate(b2, 5, w2, 3, nu, 0, nu, 0, nu, 0, op, 2);
            US[5] = new ThinningTemplate(US[4], 0);

            ThinningTemplate * NE[6], *WD[6], *ES[6], *UW[6], *ND[6], *SW[6], *UN[6], *ED[6], *NW[6], *UE[6], *SD[6];

            for (i = 0; i < 6; i++)
            {
                SD[i] = new ThinningTemplate(US[i], 0, 1);
                ND[i] = new ThinningTemplate(SD[i], 0, 1);
                UN[i] = new ThinningTemplate(ND[i], 0, 1);

                ES[i] = new ThinningTemplate(US[i], 1, 1);
                NE[i] = new ThinningTemplate(ES[i], 2, 1);
                NW[i] = new ThinningTemplate(NE[i], 2, 1);
                SW[i] = new ThinningTemplate(NW[i], 2, 1);

                UE[i] = new ThinningTemplate(US[i], 2, 1);
                ED[i] = new ThinningTemplate(UE[i], 1, 1);
                WD[i] = new ThinningTemplate(ED[i], 1, 1);
                UW[i] = new ThinningTemplate(WD[i], 1, 1);
            }

            ThinningTemplate** alltemps[12] = { US, NE, WD, ES, UW, ND, SW, UN, ED, NW, UE, SD };

            // Initialize queues
            printf("Initializing queues...\n");
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            GridQueue2* queue4 = new GridQueue2();
            Volume* fvol = new Volume(this->getSizeX(), this->getSizeY(), this->getSizeZ(), 0);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 1)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        if (getDataAt(i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2]) == 0)
                        {
                            fvol->setDataAt(i, j, k, 1);
                            queue2->prepend(i, j, k);
                            break;
                        }
                    }
                }
            }
            printf("Total %d nodes\n", queue2->getNumElements());

            // Start iteration
            int it = 0;
            int vox[3][3][3];
            gridQueueEle* ele;
            while (queue2->getNumElements() > 0)
            {
                printf("Iteration %d... nodes in queue: %d...", it, queue2->getNumElements());

                // queue2 holds all boundary nodes

                int deleted = 0;
                for (i = 0; i < 12; i++)
                {
                    // At the beginning of each iteration, 
                    // queue2 holds all remaining boundary nodes
                    // queue3 is a deletable array, starting empty
                    // queue4 holds the candidates for next layer

                    queue2->reset();
                    ele = queue2->getNext();

                    // For each sub-iteration, go through queue2 first and find points to delete
                    while (ele != NULL)
                    {
                        int ox = ele->x;
                        int oy = ele->y;
                        int oz = ele->z;

                        // Check with templates
                        int match = 0;
                        for (int ci = -1; ci < 2; ci++)
                        for (int cj = -1; cj < 2; cj++)
                        for (int ck = -1; ck < 2; ck++)
                        {
                            vox[ci + 1][cj + 1][ck + 1] = (int)getDataAt(ox + ci, oy + cj, oz + ck);
                        }

                        for (j = 0; j < 6; j++)
                            // j = 1 ;
                        {
                            if (alltemps[i][j]->isMatch(vox))
                            {
                                /* Debug */
                                if (!isSimple2(vox))
                                {
                                    printf("Wrong! %d %d\n", i, j);
                                    for (int cci = 0; cci < 3; cci++)
                                    {
                                        for (int ccj = 0; ccj < 3; ccj++)
                                        {
                                            for (int cck = 0; cck < 3; cck++)
                                            {
                                                printf("%d ", vox[cci][ccj][cck]);
                                            }
                                            printf(" , ");
                                        }
                                        printf("\n");
                                    }
                                    exit(0);
                                }

                                // Move this node to queue3
                                ele = queue2->remove();
                                queue3->prepend(ox, oy, oz);

                                // Move its neighboring unvisited node to queue4
                                for (int m = 0; m < 6; m++)
                                {
                                    int nx = ox + neighbor6[m][0];
                                    int ny = oy + neighbor6[m][1];
                                    int nz = oz + neighbor6[m][2];
                                    if (fvol->getDataAt(nx, ny, nz) == 0)
                                    {
                                        fvol->setDataAt(nx, ny, nz, 1);
                                        queue4->prepend(nx, ny, nz);
                                    }
                                }

                                match = 1;
                                break;
                            }
                        }

                        if (match == 0)
                        {
                            ele = queue2->getNext();
                        }
                    }

                    // Now, queue2 holds all remaining nodes, 
                    // queue3 holds all deletable nodes, 
                    // and queue4 holds nodes to be added to the next layer

                    // Simultaneous deletion
                    queue3->reset();
                    ele = queue3->getNext();
                    while (ele != NULL)
                    {
                        setDataAt(ele->x, ele->y, ele->z, 0);
                        ele = queue3->remove();
                        deleted++;
                    }

                    // return ;
                }

                // After all sub-iterations
                // Move all queue4 nodes to queue2
                queue4->reset();
                ele = queue4->getNext();
                while (ele != NULL)
                {
                    queue2->prepend(ele->x, ele->y, ele->z);
                    ele = queue4->remove();
                }

                if (deleted == 0)
                {
                    printf("No more deletable nodes.\n");
                    break;
                }
                else
                {
                    printf("Deleted: %d\n", deleted);
                }
            }
        }

        /**
        * Normalize to a given range
        */
        void Volume::threshold(double thr)
        {
            threshold(thr, 0, 1, 0, true);
        }

        void Volume::threshold(double thr, int out, int in)
        {
            threshold(thr, out, in, out, true);
        }

        void Volume::threshold(double thr, int out, int in, int boundary)
        {
            threshold(thr, out, in, boundary, true);
        }

        void Volume::threshold(double thr, int out, int in, int boundary, bool markBoundary)
        {
            float val;
            for (int i = 0; i < getSizeX(); i++)
            for (int j = 0; j < getSizeY(); j++)
            for (int k = 0; k < getSizeZ(); k++)
            {
                val = (float)getDataAt(i, j, k);
                if (markBoundary) {
                    if (i > 1 && i < getSizeX() - 2 && j > 1 && j < getSizeY() - 2 && k > 1 && k < getSizeZ() - 2) {
                        if (val < thr) {
                            setDataAt(i, j, k, out);
                        }
                        else {
                            setDataAt(i, j, k, in);
                        }
                    }
                    else
                    {
                        setDataAt(i, j, k, boundary);
                    }
                }
                else {
                    if (val < thr) {
                        setDataAt(i, j, k, out);
                    }
                    else {
                        setDataAt(i, j, k, in);
                    }
                }
            }
        }

        void Volume::threshold2(double thr, int out, int in)
        {
            for (int i = 0; i < getSizeX(); i++)
            for (int j = 0; j < getSizeY(); j++)
            for (int k = 0; k < getSizeZ(); k++) {
                double val = getDataAt(i, j, k);
                if (val <= thr) {
                    setDataAt(i, j, k, out);
                }
                else {
                    setDataAt(i, j, k, in);
                }
            }
        }

        void Volume::smooth(float alpha)
        {
            VolumeData * smoothedData = new VolumeData(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, volData);

            for (int i = 1; i < getSizeX() - 1; i++)
            for (int j = 1; j < getSizeY() - 1; j++)
            for (int k = 1; k < getSizeZ() - 1; k++) {
                float v = (float)getDataAt(i - 1, j, k) +
                    (float)getDataAt(i + 1, j, k) +
                    (float)getDataAt(i, j - 1, k) +
                    (float)getDataAt(i, j + 1, k) +
                    (float)getDataAt(i, j, k - 1) +
                    (float)getDataAt(i, j, k + 1);
                smoothedData->SetDataAt(i, j, k, smoothedData->GetDataAt(i, j, k) * alpha + (1 - alpha) * v / 6);
            }
            delete volData;
            volData = smoothedData;
        }

        void Volume::normalize(double min, double max)
        {
            double imin = getMin();
            double imax = getMax();
            double irange = imax - imin;
            double range = max - min;

            int size = volData->GetMaxIndex();
            for (int i = 0; i < size; i++) {
                setDataAt(i, ((getDataAt(i) - (float)imin) / (float)irange) * (float)range + (float)min);
            }
        }

        void Volume::normalize(double min, double max, double thresh, double ithresh)
        {
            double imin = getMin();
            double imax = getMax();
            double irange1 = ithresh - imin;
            double irange2 = imax - ithresh;
            double range1 = thresh - min;
            double range2 = max - thresh;

            int size = volData->GetMaxIndex();
            for (int i = 0; i < size; i++) {
                if (getDataAt(i) < ithresh) {
                    setDataAt(i, ((getDataAt(i) - (float)imin) / (float)irange1) * (float)range1 + (float)min);
                }
                else
                {
                    setDataAt(i, (float)max - (((float)imax - getDataAt(i)) / (float)irange2) * (float)range2);
                }
            }
        }

        /* Set data at a pixel */

        Volume * Volume::getDataRange(int x, int y, int z, int radius) {
            Volume * range = new Volume(radius * 2 + 1, radius * 2 + 1, radius * 2 + 1);
            for (int xx = x - radius; xx <= x + radius; xx++) {
                for (int yy = y - radius; yy <= y + radius; yy++) {
                    for (int zz = z - radius; zz <= z + radius; zz++) {
                        range->setDataAt(xx - x + radius, yy - y + radius, zz - z + radius, getDataAt(xx, yy, zz));
                    }
                }
            }
            return range;
        }

        /* Get data at an interpolated voxel */
        double Volume::getInterpDataAt(double x, double y, double z)
        {
            /*
            double rad = getSizeX() / 4.0 ;
            double cent = ( getSizeX() - 1 ) / 2.0 ;

            double ox = x - cent ;
            double oy = y - cent ;
            double oz = z - cent ;

            double a = -0.3 ;
            double nx = ox ;
            double ny = cos( a ) * oy + sin( a ) * oz ;
            double nz = - sin( a ) * oy + cos( a ) * oz ;

            double b = 1.4 ;
            double nnx = cos( b ) * nx + sin( b ) * ny - 2;
            double nny = -sin( b ) * nx + cos ( b ) * ny - 1;
            double nnz = nz + 1;

            double dis = nnx * nnx + nny * nny ;
            return 10 - 10 * dis / ( rad * rad ) ;
            */

            double rvalue;
            int hx = (int)ceil(x);
            int lx = (int)floor(x);
            int hy = (int)ceil(y);
            int ly = (int)floor(y);
            int hz = (int)ceil(z);
            int lz = (int)floor(z);

            double x1 = x - lx, x2 = 1 - x1;
            double r1 = x2 * getDataAt(lx, ly, lz) + x1 * getDataAt(hx, ly, lz);
            double r2 = x2 * getDataAt(lx, ly, hz) + x1 * getDataAt(hx, ly, hz);
            double r3 = x2 * getDataAt(lx, hy, lz) + x1 * getDataAt(hx, hy, lz);
            double r4 = x2 * getDataAt(lx, hy, hz) + x1 * getDataAt(hx, hy, hz);

            double y1 = y - ly, y2 = 1 - y1;
            double r5 = y2 * r1 + y1 * r3;
            double r6 = y2 * r2 + y1 * r4;

            double z1 = z - lz, z2 = 1 - z1;
            rvalue = z2 * r5 + z1 * r6;

            return rvalue;
        }

        /* Rotation routine */
        void Volume::rotateX(double a)
        {
            int i;
            int sizeX = getSizeX(), sizeY = getSizeY(), sizeZ = getSizeZ();

            if (sizeX != sizeY || sizeX != sizeZ) {
                return;
            }

            VolumeData * newData = new VolumeData(sizeX, sizeY, sizeZ, 0, 0, 0, volData);

            double cent = (sizeX - 1) / 2.0;
            for (i = 0; i < sizeX; i++)
            for (int j = 0; j < sizeY; j++)
            for (int k = 0; k < sizeZ; k++)
            {
                double x = i - cent;
                double y = j - cent;
                double z = k - cent;

                double nx = x + cent;
                double ny = cos(a) * y + sin(a) * z + cent;
                double nz = -sin(a) * y + cos(a) * z + cent;

                if (nx < 0) {
                    nx = 0;
                }
                else if (nx > sizeX - 1) {
                    nx = sizeX - 1;
                }

                if (ny < 0) {
                    ny = 0;
                }
                else if (ny > sizeY - 1) {
                    ny = sizeY - 1;
                }

                if (nz < 0) {
                    nz = 0;
                }
                else if (nz > sizeZ - 1) {
                    nz = sizeZ - 1;
                }

                newData->SetDataAt(i, j, k, (float)getInterpDataAt(nx, ny, nz));
            }

            delete volData;
            volData = newData;
        }


        /* Write to file */
        void Volume::toMathematicaFile(char* fname)
        {
            FILE* fout = fopen(fname, "w");

            fprintf(fout, "{");
            for (int i = 0; i < getSizeX(); i++)
            {
                fprintf(fout, "{");
                for (int j = 0; j < getSizeY(); j++)
                {
                    fprintf(fout, "{");
                    for (int k = 0; k < getSizeZ(); k++)
                    {
                        fprintf(fout, "%.15f", getDataAt(i, j, k));
                        if (k < getSizeZ() - 1)
                        {
                            fprintf(fout, ",");
                        }
                    }
                    fprintf(fout, "}");
                    if (j < getSizeY() - 1)
                    {
                        fprintf(fout, ",\n");
                    }
                    else {
                        fprintf(fout, "\n");
                    }
                }
                fprintf(fout, "}");
                if (i < getSizeX() - 1)
                {
                    fprintf(fout, ",\n\n\n");
                }
                else {
                    fprintf(fout, "\n\n\n");
                }
            }
            fprintf(fout, "}");

            fclose(fout);

        }

        /* Write to file */
        void Volume::toMathematicaFile(char* fname, int lx, int hx, int ly, int hy, int lz, int hz)
        {
            FILE* fout = fopen(fname, "w");

            fprintf(fout, "{");
            for (int i = lx; i < hx; i++)
            {
                fprintf(fout, "{");
                for (int j = ly; j < hy; j++)
                {
                    fprintf(fout, "{");
                    for (int k = lz; k < hz; k++)
                    {
                        fprintf(fout, "%.15f", getDataAt(i, j, k));
                        if (k < hz - 1)
                        {
                            fprintf(fout, ",");
                        }
                    }
                    fprintf(fout, "}");
                    if (j < hy - 1)
                    {
                        fprintf(fout, ",");
                    }
                }
                fprintf(fout, "}");
                if (i < hx - 1)
                {
                    fprintf(fout, ",");
                }
            }
            fprintf(fout, "}");

            fclose(fout);

        }

        void Volume::toOFFCells(char* fname)
        {
            toOFFCells(fname, 0.0001f);
        }
        void Volume::toOFFCells2(char* fname)
        {
            toOFFCells2(fname, 0.0001f);
        }

        void Volume::toOFFCells2(char* fname, float thr)
        {
            int i, j, k;
            Volume* indvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), -1);

            // Get number of cells to write
            int numverts = 0, numfaces = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    indvol->setDataAt(i, j, k, numverts);
                    numverts++;

                    for (int mi = 0; mi < 3; mi++)
                    {
                        int find = mi * 4 + 3;
                        int isFace = 1;
                        for (int mj = 0; mj < 4; mj++)
                        {
                            int nx = i + sheetNeighbor[find][mj][0];
                            int ny = j + sheetNeighbor[find][mj][1];
                            int nz = k + sheetNeighbor[find][mj][2];

                            if (getDataAt(nx, ny, nz) < thr)
                            {
                                isFace = 0;
                                break;
                            }
                        }
                        if (isFace)
                        {
                            numfaces++;
                        }

                        int eind = mi * 2 + 1;
                        if (getDataAt(i + neighbor6[eind][0], j + neighbor6[eind][1], k + neighbor6[eind][2]) >= thr)
                        {
                            numfaces++;
                        }
                    }
                }
            }

            FILE* fin = fopen(fname, "w");
            fprintf(fin, "OFF\n");
            fprintf(fin, "%d %d 0\n", numverts, numfaces);

            // Write vertices
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    fprintf(fin, "%d %d %d\n", i, j, k);
                }
            }

            // Write faces
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    int thisvt = (int)(indvol->getDataAt(i, j, k));
                    for (int mi = 0; mi < 3; mi++)
                    {
                        int find = mi * 4 + 3;
                        int isFace = 1;
                        int vts[4];
                        for (int mj = 0; mj < 4; mj++)
                        {
                            int nx = i + sheetNeighbor[find][mj][0];
                            int ny = j + sheetNeighbor[find][mj][1];
                            int nz = k + sheetNeighbor[find][mj][2];

                            vts[mj] = (int)(indvol->getDataAt(nx, ny, nz));

                            if (getDataAt(nx, ny, nz) < thr)
                            {
                                isFace = 0;
                                break;
                            }
                        }
                        if (isFace)
                        {
                            fprintf(fin, "4 %d %d %d %d\n", vts[0], vts[1], vts[3], vts[2]);
                        }

                        int eind = mi * 2 + 1;
                        int mx = i + neighbor6[eind][0];
                        int my = j + neighbor6[eind][1];
                        int mz = k + neighbor6[eind][2];
                        int vt = (int)(indvol->getDataAt(mx, my, mz));
                        if (getDataAt(mx, my, mz) >= thr)
                        {
                            fprintf(fin, "4 %d %d %d %d\n", thisvt, thisvt, vt, vt);
                        }
                    }
                }
            }

            fclose(fin);
            delete indvol;
        }

        void Volume::toOFFCells(char* fname, float thr)
        {
            int i, j, k;
            // Volume* indvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), -1 ) ;

            // Get number of cells to write
            int numverts = 0, numfaces = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    numverts += 8;
                    for (int mi = 0; mi < 6; mi++)
                    {
                        if (getDataAt(i + neighbor6[mi][0], j + neighbor6[mi][1], k + neighbor6[mi][2]) < thr)
                        {
                            numfaces++;
                        }
                    }
                }
            }

            FILE* fin = fopen(fname, "w");
            fprintf(fin, "OFF\n");
            fprintf(fin, "%d %d 0\n", numverts, numfaces);

            // Write vertices
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    float ox = i - 0.5f;
                    float oy = j - 0.5f;
                    float oz = k - 0.5f;

                    for (int mi = 0; mi < 2; mi++)
                    for (int mj = 0; mj < 2; mj++)
                    for (int mk = 0; mk < 2; mk++)
                    {
                        fprintf(fin, "%f %f %f\n", ox + mi, oy + mj, oz + mk);
                    }
                }
            }

            // Write faces
            int ct = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) >= thr)
                {
                    for (int mi = 0; mi < 6; mi++)
                    {
                        if (getDataAt(i + neighbor6[mi][0], j + neighbor6[mi][1], k + neighbor6[mi][2]) < thr)
                        {
                            fprintf(fin, "4 %d %d %d %d\n", cubeFaces[mi][0] + ct, cubeFaces[mi][1] + ct, cubeFaces[mi][2] + ct, cubeFaces[mi][3] + ct);
                        }
                    }

                    ct += 8;
                }
            }

            fclose(fin);
            //delete indvol ;
        }

        void Volume::segment(float threshold, Volume* lowvol, Volume* highvol, char* mrcfile)
        {
            int i, j, k;
            Volume* segvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (lowvol->getDataAt(i, j, k) > 0)
                {
                    segvol->setDataAt(i, j, k, 1);
                }
                else if (highvol->getDataAt(i, j, k) > 0)
                {
                    segvol->setDataAt(i, j, k, 2);
                }
                else
                {
                    segvol->setDataAt(i, j, k, 0);
                }
            }

            writeSegmentation(threshold, segvol, NULL, mrcfile);
        }

        void Volume::segment(float threshold, Volume* vol, int maxDis, char* mrcfile)
        {
            int i, j;
            Volume* testvol = NULL;
            Volume* disvol = new Volume(getSizeX(), getSizeY(), getSizeZ());
            printf("Writing distance transform to %d levels.\n", maxDis);

            int totNodes = 0;
            int size = getSizeX() * getSizeY() * getSizeZ();
            for (i = maxDis; i >= 0; i--)
            {
                if (i == 1) continue;

                int nodes = 0;
                testvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, vol);
                testvol->erodeSheet(i);

                for (j = 0; j < size; j++)
                {
                    if (disvol->getDataAt(j) == 0 && testvol->getDataAt(j) > 0)
                    {
                        disvol->setDataAt(j, i + 1);
                        nodes++;
                    }
                }
                printf("Level %d has %d nodes.\n", i, nodes);
                totNodes += nodes;
                delete testvol;
            }
            printf("Totally %d nodes.\n", totNodes);

            writeSegmentation(threshold, disvol, NULL, mrcfile);
        }

        // Segment the volume using segvol
        // background voxels have values 0, others have values > 0
        void Volume::writeSegmentation(float threshold, Volume* segvol, char* txtfile, char* mrcfile)
        {
            printf("Start segmentation.\n");
            int i, j, k;
            Volume* vvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
            Volume* tvol1 = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
            Volume* tvol2 = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);

            // Initialize
            GridQueue2* queue = new GridQueue2();
            GridQueue2* queue2 = new GridQueue2();
            GridQueue2* queue3 = new GridQueue2();
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) < threshold || segvol->getDataAt(i, j, k) <= 0)
                {
                    continue;
                }

                vvol->setDataAt(i, j, k, 1);
                queue->prepend(i, j, k);
            }

            // Dilation
            printf("Dilation...");
            int ox, oy, oz;
            gridQueueEle* ele;
            while (queue->getNumElements() > 0)
            {
                // At the beginning
                // queue holds all nodes of previous layer
                // queue2 is empty, ready to fill with the current layer

                // First, fill current layer and assign values
                queue->reset();
                ele = queue->getNext();
                while (ele != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;
                    double seg = segvol->getDataAt(ox, oy, oz);
                    int isBorder = 0;

                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (getDataAt(nx, ny, nz) >= threshold && vvol->getDataAt(nx, ny, nz) == 0)
                        {
                            double ct = (int)tvol1->getDataAt(nx, ny, nz);
                            double val = tvol2->getDataAt(nx, ny, nz);

                            tvol1->setDataAt(nx, ny, nz, ct + 1);
                            tvol2->setDataAt(nx, ny, nz, val + seg);

                            if (ct == 0)
                            {
                                queue2->prepend(nx, ny, nz);
                            }
                        }
                        else if (getDataAt(nx, ny, nz) < threshold)
                        {
                            isBorder = 1;
                        }
                    }

                    if (isBorder)
                    {
                        queue3->prepend(ox, oy, oz);
                    }

                    ele = queue->remove();
                }

                // Normalize values of nodes in queue2
                queue2->reset();
                while ((ele = queue2->getNext()) != NULL)
                {
                    ox = ele->x;
                    oy = ele->y;
                    oz = ele->z;

                    double ct = (int)tvol1->getDataAt(ox, oy, oz);
                    double val = tvol2->getDataAt(ox, oy, oz);

                    if (ct == 0)
                    {
                        printf("Wrong! %f\n", ct);
                    }
                    segvol->setDataAt(ox, oy, oz, val / ct);
                    vvol->setDataAt(ox, oy, oz, 1);
                }

                // Finally, swap queues
                GridQueue2* temp = queue2;
                queue2 = queue;
                queue = temp;
            }
            printf("Done.\n");

            /* Debug
            for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
            {
            if ( getDataAt(i) >= threshold && segvol->getDataAt(i) == 0 )
            {
            segvol->setDataAt(i, 1) ;
            }
            else
            {
            segvol->setDataAt(i, 0) ;
            }
            }
            */

            // Write to MRC
            //printf("Writing to %s...", mrcfile) ;
            //segvol->toMRCFile( mrcfile ) ;


            // Growing out one layer for display
            printf("Growing out one layer...\n");
            // First, fill current layer and assign values
            queue2->reset();
            queue3->reset();
            ele = queue3->getNext();
            while (ele != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;
                double seg = segvol->getDataAt(ox, oy, oz);

                for (int mx = -1; mx < 2; mx++)
                for (int my = -1; my < 2; my++)
                for (int mz = -1; mz < 2; mz++)
                {
                    int nx = ox + mx; // neighbor6[m][0] ;
                    int ny = oy + my; // neighbor6[m][1] ;
                    int nz = oz + mz; // neighbor6[m][2] ;
                    if (vvol->getDataAt(nx, ny, nz) == 0)
                    {
                        double ct = (int)tvol1->getDataAt(nx, ny, nz);
                        double val = tvol2->getDataAt(nx, ny, nz);

                        tvol1->setDataAt(nx, ny, nz, ct + 1);
                        tvol2->setDataAt(nx, ny, nz, val + seg);

                        if (ct == 0)
                        {
                            queue2->prepend(nx, ny, nz);
                        }
                    }
                }
                ele = queue3->remove();
            }

            // Normalize values of nodes in queue2
            queue2->reset();
            while ((ele = queue2->getNext()) != NULL)
            {
                ox = ele->x;
                oy = ele->y;
                oz = ele->z;

                double ct = tvol1->getDataAt(ox, oy, oz);
                double val = tvol2->getDataAt(ox, oy, oz);

                if (ct == 0)
                {
                    printf("Wrong! %f\n", ct);
                }
                segvol->setDataAt(ox, oy, oz, val / ct);
            }

            printf("Writing...");
            segvol->toMRCFile(mrcfile);
            //			segvol->toMRCFile( "../colors.mrc" ) ;
            printf("Done.\n");

            printf("Segmentation...");
            for (i = 0; i < getSizeX() * getSizeY() * getSizeZ(); i++)
            {
                float segval = (float)segvol->getDataAt(i);
                // High values
                if (segval > 1.5f)
                {
                    tvol1->setDataAt(i, getDataAt(i));
                }
                else
                {
                    tvol1->setDataAt(i, -1);
                }

                // Low values
                if (segval < 1.5f && segval >= 1)
                {
                    tvol2->setDataAt(i, getDataAt(i));
                }
                else
                {
                    tvol2->setDataAt(i, -1);
                }
            }
            char nname[1024];
            sprintf(nname, "%s_sheet.mrc", mrcfile);
            tvol1->toMRCFile(nname);
            sprintf(nname, "%s_helix.mrc", mrcfile);
            tvol2->toMRCFile(nname);
            printf("Done.\n");
            return;

            /* Write to text
            if ( txtfile != NULL )
            {
            printf("Writing to %s...", txtfile) ;
            // Count border points
            queue->reset() ;
            for ( i = 0 ; i < getSizeX() ; i ++ )
            for ( j = 0 ; j < getSizeY() ; j ++ )
            for ( k = 0 ; k < getSizeZ() ; k ++ )
            {
            if ( getDataAt( i, j, k ) >= threshold )
            {
            for ( int m = 0 ; m < 6 ; m ++ )
            {
            if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < threshold )
            {
            queue->prepend( i, j, k ) ;
            break ;
            }
            }
            }
            }

            FILE* fout = fopen( txtfile, "w" ) ;
            fprintf( fout, "%d\n", queue->getNumElements() ) ;
            queue->reset() ;
            while( (ele=queue->getNext()) != NULL )
            {
            ox = ele->x ;
            oy = ele->y ;
            oz = ele->z ;
            fprintf( fout, "%d %d %d %f\n", ox, oy, oz, segvol->getDataAt(ox,oy,oz) ) ;
            }
            fclose( fout ) ;
            printf("Done.\n") ;
            }
            */
        }

        void Volume::floodFill(float thr)
        {
            int i;
            // First, threshold the volume
#ifdef VERBOSE
            printf("Thresholding the volume to 0/1...\n");
#endif
            threshold(thr, 0, 1);

            // Next, initialize queue
            Volume* tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
            GridQueue2* queue = new GridQueue2();
            gridQueueEle* ele;
            queue->prepend(0, 0, 0);
            tvol->setDataAt(0, 0, 0, -1);

            // Iteration
            printf("Flood filling...\n");
            queue->reset();
            ele = queue->getNext();
            int ct = 1;
            while (ele != NULL)
            {
                int ox = ele->x;
                int oy = ele->y;
                int oz = ele->z;
                queue->remove();

                for (int m = 0; m < 6; m++)
                {
                    int nx = ox + neighbor6[m][0];
                    int ny = oy + neighbor6[m][1];
                    int nz = oz + neighbor6[m][2];
                    if (nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ())
                    {
                        continue;
                    }
                    if (tvol->getDataAt(nx, ny, nz) == 0 && getDataAt(nx, ny, nz) == 0)
                    {
                        queue->prepend(nx, ny, nz);
                        tvol->setDataAt(nx, ny, nz, -1);
                        ct++;

                        if (ct % 100000 == 0)
                        {
                            printf("%d nodes processed.\n", ct);
                        }
                    }
                }

                queue->reset();
                ele = queue->getNext();
            }
            printf("Done.\n");

            // Done
            for (i = 0; i < getSizeX()*getSizeY()*getSizeZ(); i++)
            {
                this->setDataAt(i, tvol->getDataAt(i));
            }
        }

        void Volume::reduceComponent(int size)
        {
            int i, j, k;
            Volume* tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
            GridQueue2* queue = new GridQueue2();
            GridQueue2* queue2 = new GridQueue2();
            gridQueueEle* ele;
            int numC = 0, numBC = 0;

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 1 && tvol->getDataAt(i, j, k) == 0)
                {
                    numC++;
                    // Start searching for this component
                    queue->prepend(i, j, k);
                    queue2->prepend(i, j, k);
                    tvol->setDataAt(i, j, k, 1);
                    int ct = 1;

                    queue->reset();
                    ele = queue->getNext();
                    while (ele != NULL)
                    {
                        int ox = ele->x;
                        int oy = ele->y;
                        int oz = ele->z;
                        queue->remove();

                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ())
                            {
                                continue;
                            }
                            if (tvol->getDataAt(nx, ny, nz) == 0 && getDataAt(nx, ny, nz) == 1)
                            {
                                queue->prepend(nx, ny, nz);
                                queue2->prepend(nx, ny, nz);
                                tvol->setDataAt(nx, ny, nz, 1);
                                ct++;
                            }
                        }

                        queue->reset();
                        ele = queue->getNext();
                    }

                    if (ct < size)
                    {
                        // remove component

                        queue2->reset();
                        ele = queue2->getNext();
                        while (ele != NULL)
                        {
                            setDataAt(ele->x, ele->y, ele->z, 0);
                            ele = queue2->remove();
                        }
                        queue2->reset();

                    }
                    else
                    {
                        queue2 = new GridQueue2();
                        numBC++;
                    }
                }
            }

            printf("Total number of components: %d Remained: %d\n", numC, numBC);

        }

        void Volume::reduceComponent2(int num)
        {
            int i, j, k;
            Volume* tvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
            GridQueue2* queue = new GridQueue2();
            GridQueue2* queue2 = new GridQueue2();
            gridQueueEle* ele;
            int numC = 0, numBC = 0;
            int* tops = new int[num];
            int* topinds = new int[num];
            for (i = 0; i < num; i++)
            {
                tops[i] = 0;
                topinds[i] = -1;
            }

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 1 && tvol->getDataAt(i, j, k) == 0)
                {
                    numC++;
                    // Start searching for this component
                    queue->prepend(i, j, k);
                    queue2->prepend(i, j, k);
                    tvol->setDataAt(i, j, k, 1);
                    int ct = 1;

                    queue->reset();
                    ele = queue->getNext();
                    while (ele != NULL)
                    {
                        int ox = ele->x;
                        int oy = ele->y;
                        int oz = ele->z;
                        queue->remove();

                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ())
                            {
                                continue;
                            }
                            if (tvol->getDataAt(nx, ny, nz) == 0 && getDataAt(nx, ny, nz) == 1)
                            {
                                queue->prepend(nx, ny, nz);
                                queue2->prepend(nx, ny, nz);
                                tvol->setDataAt(nx, ny, nz, 1);
                                ct++;
                            }
                        }

                        queue->reset();
                        ele = queue->getNext();
                    }

                    queue2->reset();
                    ele = queue2->getNext();
                    while (ele != NULL)
                    {
                        ele = queue2->remove();
                    }
                    queue2->reset();

                    for (int ind = 0; ind < num; ind++)
                    {
                        if (ct > tops[ind])
                        {
                            for (int nind = num - 1; nind > ind; nind--)
                            {
                                tops[nind] = tops[nind - 1];
                                topinds[nind] = topinds[nind - 1];
                            }
                            tops[ind] = ct;
                            topinds[ind] = numC;
                            break;
                        }
                    }
                }
            }

            printf("%d components total, %d selected components:\n", numC, num);
            for (i = 0; i < num; i++)
            {
                printf("Id %d: size %d\n", topinds[i], tops[i]);
            }

            tvol->fill(0);

            numC = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (getDataAt(i, j, k) == 1 && tvol->getDataAt(i, j, k) == 0)
                {
                    numC++;
                    // Start searching for this component
                    queue->prepend(i, j, k);
                    queue2->prepend(i, j, k);
                    tvol->setDataAt(i, j, k, 1);
                    int ct = 1;

                    queue->reset();
                    ele = queue->getNext();
                    while (ele != NULL)
                    {
                        int ox = ele->x;
                        int oy = ele->y;
                        int oz = ele->z;
                        queue->remove();

                        for (int m = 0; m < 6; m++)
                        {
                            int nx = ox + neighbor6[m][0];
                            int ny = oy + neighbor6[m][1];
                            int nz = oz + neighbor6[m][2];
                            if (nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ())
                            {
                                continue;
                            }
                            if (tvol->getDataAt(nx, ny, nz) == 0 && getDataAt(nx, ny, nz) == 1)
                            {
                                queue->prepend(nx, ny, nz);
                                queue2->prepend(nx, ny, nz);
                                tvol->setDataAt(nx, ny, nz, 1);
                                ct++;
                            }
                        }

                        queue->reset();
                        ele = queue->getNext();
                    }

                    int removing = 1;
                    for (int ind = 0; ind < num; ind++)
                    {
                        if (topinds[ind] == numC)
                        {
                            removing = 0;
                            break;
                        }
                    }

                    if (removing)
                    {
                        // remove component

                        queue2->reset();
                        ele = queue2->getNext();
                        while (ele != NULL)
                        {
                            setDataAt(ele->x, ele->y, ele->z, 0);
                            ele = queue2->remove();
                        }
                        queue2->reset();

                    }
                    else
                    {
                        queue2 = new GridQueue2();
                        numBC++;
                    }
                }
            }

            printf("Total number of components: %d Remained: %d\n", numC, numBC);
            delete tvol;
        }

        void Volume::floodFillPQR(int offset)
        {
            int i;

            // Next, initialize queue
            GridQueue2* queue = new GridQueue2();
            gridQueueEle* ele;
            queue->prepend(0, 0, 0);

            // Iteration
            printf("Flood filling outside from  (0,0,0)...\n");
            int ct = 1;
            /*
            queue->reset() ;
            ele = queue->getNext() ;
            while( ele != NULL )
            {
            int ox = ele->x ;
            int oy = ele->y ;
            int oz = ele->z ;
            queue->remove() ;

            int boundary = 0 ;
            for ( int m = 0 ; m < 6 ; m ++ )
            {
            int nx = ox + neighbor6[m][0] ;
            int ny = oy + neighbor6[m][1] ;
            int nz = oz + neighbor6[m][2] ;
            if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
            {
            continue ;
            }
            if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 0 )
            {
            queue->prepend( nx, ny, nz ) ;
            tvol->setDataAt( nx, ny, nz, 1 ) ;
            ct ++ ;
            if ( ct % 100000 == 0 )
            {
            printf("%d nodes processed.\n", ct);
            }
            }
            else if ( getDataAt( nx, ny, nz ) == 1 )
            {
            boundary = 1 ;
            }
            }

            if ( boundary )
            {
            tvol->setDataAt( ox, oy, oz, 2 ) ;
            }

            queue->reset() ;
            ele = queue->getNext() ;
            }
            printf("Done.\n") ;
            for ( i = 0 ; i < getSizeX()*getSizeY()*getSizeZ() ; i ++ )
            {

            if ( tvol->getDataAt(i) == 2 )
            {
            setDataAt( i, 1 ) ;
            }
            }
            */


            // Find inside seed point
            int maxRounds = 1;
            for (int rounds = 0; rounds < maxRounds; rounds++)
            {

                int isSolid = 0;
                for (i = 0; i < getSizeX(); i++)
                {
                    if (getDataAt(i, getSizeY() / 2, getSizeZ() / 2) == 1)
                    {
                        isSolid = 1;
                    }
                    else if (isSolid == 1)
                    {
                        break;
                    }
                }

                Volume* invol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0);
                queue->prepend(i, getSizeY() / 2, getSizeZ() / 2);
                invol->setDataAt(i, getSizeY() / 2, getSizeZ() / 2, 1);
                printf("Flood filling inside from  (%d,%d,%d)...\n", i, getSizeY() / 2, getSizeZ() / 2);

                // Iteration
                queue->reset();
                ele = queue->getNext();
                ct = 1;
                while (ele != NULL)
                {
                    int ox = ele->x;
                    int oy = ele->y;
                    int oz = ele->z;
                    queue->remove();

                    int boundary = 0;
                    for (int m = 0; m < 6; m++)
                    {
                        int nx = ox + neighbor6[m][0];
                        int ny = oy + neighbor6[m][1];
                        int nz = oz + neighbor6[m][2];
                        if (nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ())
                        {
                            continue;
                        }
                        if (invol->getDataAt(nx, ny, nz) == 0 && getDataAt(nx, ny, nz) == 0)
                        {
                            queue->prepend(nx, ny, nz);
                            invol->setDataAt(nx, ny, nz, 1);
                            ct++;
                            if (ct % 100000 == 0)
                            {
                                printf("%d nodes processed.\n", ct);
                            }
                        }
                        else if (getDataAt(nx, ny, nz) == 1)
                        {
                            boundary = 1;
                        }
                    }

                    if (boundary)
                    {
                        invol->setDataAt(ox, oy, oz, 2);
                    }

                    queue->reset();
                    ele = queue->getNext();
                }
                printf("Done.\n");

                // Done
                for (i = 0; i < getSizeX()*getSizeY()*getSizeZ(); i++)
                {
                    /*
                    if ( tvol->getDataAt(i) == 2 )
                    {
                    setDataAt( i, 1 ) ;
                    }
                    */
                    if (invol->getDataAt(i) == 2)
                    {
                        setDataAt(i, 1);
                    }

                    /*
                    else if ( tvol->getDataAt(i) == 0 && invol->getDataAt(i) == 0 )
                    {
                    setDataAt( i, 1 ) ;
                    }
                    */
                }

                delete invol;

            }

            //		delete tvol ;

        }


        void Volume::writeDistances(char* fname, int maxDis)
        {
            int i, j, k;
            Volume* testvol = NULL;
            Volume* disvol = new Volume(getSizeX(), getSizeY(), getSizeZ());
            printf("Writing distance transform to %d levels.\n", maxDis);

            int totNodes = 0;
            int size = getSizeX() * getSizeY() * getSizeZ();
            float score = 10;
            for (i = maxDis; i >= 0; i--)
            {
                int nodes = 0;
                testvol = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);
                testvol->erodeSheet(i);

                for (j = 0; j < size; j++)
                {
                    if (disvol->getDataAt(j) == 0 && testvol->getDataAt(j) > 0)
                    {
                        disvol->setDataAt(j, i + 1);
                        nodes++;
                    }
                }
                printf("Level %d has %d nodes.\n", i, nodes);
                totNodes += nodes;
                score -= 0.5f;
                delete testvol;
            }
            printf("Totally %d nodes.\n", totNodes);

            //disvol->toMRCFile( "..\distance.mrc" ) ;
            FILE* fout = fopen(fname, "w");
            fprintf(fout, "%d\n", totNodes);
            int ct = 0;
            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                float val = (float)disvol->getDataAt(i, j, k);
                if (val > 0)
                {
                    fprintf(fout, "%d %d %d %f\n", i, j, k, val);
                    ct++;
                }
            }

            if (ct != totNodes)
            {
                printf("Counting wrong! %d %d\n", totNodes, ct);
            }
            fclose(fout);

        }

        void Volume::toPQRFile(char* fname, float spc, float minx, float miny, float minz, int padding)
        {
            FILE* fout = fopen(fname, "w");
            int i, j, k;

            for (i = 0; i < getSizeX(); i++)
            for (j = 0; j < getSizeY(); j++)
            for (k = 0; k < getSizeZ(); k++)
            {
                if (i < padding || i >= getSizeX() - padding || j < padding || j >= getSizeY() - padding || k < padding || k >= getSizeZ() - padding)
                {
                    continue;
                }
                float val = (float)this->getDataAt(i, j, k);
                if (val > 0)
                {
                    float x = (i - padding) * spc + minx;
                    float y = (j - padding) * spc + miny;
                    float z = (k - padding) * spc + minz;
                    fprintf(fout, "ATOM      1  X   DUM     1     %4.3f %4.3f %4.3f 0.000 0.000\n", x, y, z);
                }
            }

            fclose(fout);
        }

        void Volume::toMRCFile(char* fname)
        {
            FILE* fout = fopen(fname, "wb");

            // Write header
            int sizeX = getSizeX();
            int sizeY = getSizeY();
            int sizeZ = getSizeZ();
            fwrite(&sizeX, sizeof(int), 1, fout);
            fwrite(&sizeY, sizeof(int), 1, fout);
            fwrite(&sizeZ, sizeof(int), 1, fout);

            int mode = 2;
            fwrite(&mode, sizeof (int), 1, fout);

            int off[3] = { 0, 0, 0 };

            // Changed: MX == NX instead of MX = NX -1 (similar for Y and Z) because that is how EMAN2 and UCSF Chimera do it
            //int intv[3] = { getSizeX() - 1, getSizeY() - 1, getSizeZ() - 1 } ;
            int intv[3] = { getSizeX(), getSizeY(), getSizeZ() };

            fwrite(off, sizeof(int), 3, fout);
            fwrite(intv, sizeof(int), 3, fout);

            float cella[3] = { getSpacingX() * (float)(getSizeX()), getSpacingY() * (float)(getSizeY()), getSpacingZ() * (float)(getSizeZ()) };
            float cellb[3] = { 90, 90, 90 };
            fwrite(cella, sizeof(float), 3, fout);
            fwrite(cellb, sizeof(float), 3, fout);

            int cols[3] = { 1, 2, 3 };
            fwrite(cols, sizeof(int), 3, fout);

            double dmin = 100000, dmax = -100000;
            int i;
            int size = volData->GetMaxIndex();
            for (i = 0; i < size; i++) {
                float val = (float)getDataAt(i);
                if (val < dmin) {
                    dmin = val;
                }
                if (val > dmax) {
                    dmax = val;
                }
            }
            float ds[3] = { (float)dmin, (float)dmax, (float)0 };
            fwrite(ds, sizeof(float), 3, fout);

            int zero = 0;
            for (i = 23; i <= 49; i++)
            {
                fwrite(&zero, sizeof(int), 1, fout);
            }

            float origins[3];
            origins[0] = getOriginX() / getSpacingX() + 0.5f * (float)getSizeX();
            origins[1] = getOriginY() / getSpacingY() + 0.5f * (float)getSizeY();
            origins[2] = getOriginZ() / getSpacingZ() + 0.5f * (float)getSizeZ();

            fwrite(origins, sizeof(float), 3, fout);

            for (i = 53; i <= 256; i++)
            {
                fwrite(&zero, sizeof(int), 1, fout);
            }

            // Write contents
            for (int z = 0; z < getSizeZ(); z++)
            for (int y = 0; y < getSizeY(); y++)
            for (int x = 0; x < getSizeX(); x++)
            {
                float d = (float)getDataAt(x, y, z);
                fwrite(&d, sizeof(float), 1, fout);
            }

            fclose(fout);
        }

        // Returns the mean value of all the voxels
        float Volume::getMean()
        {
            int N = volData->GetMaxIndex();
            double mass = 0;
            for (int i = 0; i < N; i++)
                mass += getDataAt(i);
            float mean = (float)mass / N;
            return mean;
        }

        // Returns the mean value of all the surface voxels but no interior voxels
        float Volume::getEdgeMean()
        {
            int nx = getSizeX();
            int ny = getSizeY();
            int nz = getSizeZ();

            //Calculate the edge mean -- the average value of all the voxels on the surfaces (1 voxel) of the image
            double edge_sum = 0; //The sum of the values on the outer surfaces (1 voxel) of the image
            int num_voxels = 0;

            //sum the values of each voxel on the surfaces of the rectangular prism
            for (int i = 0; i < nx; i++)
            for (int j = 0; j<ny; j++)
            for (int k = 0; k<nz; k++)
            {
                if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || k == 0 || k == nz - 1)
                {
                    edge_sum += getDataAt(i, j, k);
                    num_voxels++;
                }
            }

            float edge_mean = (float)edge_sum / num_voxels;
            return edge_mean;
        }

        // Returns the population standard deviation of the values at all the voxels
        float Volume::getStdDev()
        {
            int N = volData->GetMaxIndex();

            //Calculate the standard deviation of all the voxels in the image
            double voxel_sum = 0;
            double voxel_squared_sum = 0;
            float val;

            for (int i = 0; i < N; i++)
            {
                val = (float)getDataAt(i);
                voxel_sum += val;
                voxel_squared_sum += val*val;
            }
            float std_dev = (float)sqrt((voxel_squared_sum - voxel_sum*voxel_sum / N) / N);
            return std_dev;
        }
        // Returns the center of mass of the image in pixels (not angstroms)
        Vector3DFloat Volume::getCenterOfMass()
        {
            int nx = getSizeX();
            int ny = getSizeY();
            int nz = getSizeZ();

            float mass = 0;
            float xmoment = 0;
            float ymoment = 0;
            float zmoment = 0;
            float val;

            for (int i = 0; i<nx; i++)
            for (int j = 0; j<ny; j++)
            for (int k = 0; k<nz; k++)
            {
                val = (float)getDataAt(i, j, k);
                mass += val;
                xmoment += i*val;
                ymoment += j*val;
                zmoment += k*val;
            }

            Vector3DFloat centerOfMass(xmoment / mass, ymoment / mass, zmoment / mass);
            return centerOfMass;
        }

        float* Volume::getArrayCopy(int padX, int padY, int padZ, float padValue) {
            return getVolumeData()->GetArrayCopy(padX, padY, padZ, padValue);
        }

        void Volume::buildHistogram(int binCount) {
            histogram.clear();
            for (int i = 0; i < binCount; i++) {
                histogram.push_back(0);
            }

            float minVal = getMin();
            float maxVal = getMax();
            float binSize = (maxVal - minVal) / (float)(binCount - 1);
            int binIx;
            for (unsigned int i = 0; i < getSizeX(); i++) {
                for (unsigned int j = 0; j < getSizeY(); j++) {
                    for (unsigned int k = 0; k < getSizeZ(); k++) {
                        binIx = (int)((getDataAt(i, j, k) - minVal) / binSize);
                        histogram[binIx]++;
                    }
                }
            }
        }

        int Volume::getHistogramBinValue(int binIx) {
            return histogram[binIx];

        }

    }
}

#endif
#endif
