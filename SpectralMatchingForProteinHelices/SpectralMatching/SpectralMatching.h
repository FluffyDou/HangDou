/*********************************************************************/
// SpectralMatching.h
// This is a class holding some helper functions, like line fitting.

// Hang 10/09/2014
/*********************************************************************/

#ifndef SPECTRALMATCHING_H
#define SPECTRALMATCHING_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#define EIGEN_DONT_PARALLELIZE

#include <random>
//#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <array>

#include "helper.h"
#include "Eigen\Dense"
#include "Eigen\Core"
//#include "Eigen\Eigenvalues"
#include "CPUTimer.h"
#include "PDBObject.h"
#include "HelixGeometry.h"
#include "PointCloud2D.h"
#include "ScoreMatrixDoubleLineDir.h"
#include "ScoreMatrixSingleLineDir.h"
#include "ScoreMatrix2DPointCloud.h"
#include "LineSegment.h"
#include "SMSolver2DP.h"
#include "SMSolverDLD.h"
#include "SMSolverSLD.h"
#include "Scripts.h"
#include "IPFP.h"

using Eigen::MatrixXd;

#endif