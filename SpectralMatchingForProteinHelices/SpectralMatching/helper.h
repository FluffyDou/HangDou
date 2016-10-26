/*********************************************************************/
// helper.h
// This is a class holding some helper functions, like line fitting.

// Hang 10/13/2014
/*********************************************************************/

#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <array>

#include "Eigen\Dense"

using Eigen::MatrixBase;
using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::Vector3d;
using Eigen::Matrix3f;
using Eigen::SelfAdjointEigenSolver;

typedef std::array<Vector3f,2> lineSeg;
typedef std::vector<Vector3f>  pointList;

#ifndef PI
#define PI 3.1415926585897f
#endif

const double indexEpsilon = 0.00001;

// Use least square to fit a 3D line from input points
lineSeg FitLine(pointList points);
// EuclideanDistance of two 3D points: EuDis
//float EuDis(Vector3f p1, Vector3f p2);
double EuDis(Vector3f p1, Vector3f p2);

// VectorAngle of two 3D vectors: VecAng
// Cast to double to increase precision
//float VecAng(Vector3f v1, Vector3f v2);
double VecAng(Vector3f v1, Vector3f v2);
double VecAng(Vector2f v1, Vector2f v2);

// Parse an input string by given character. Put the parsed result into a given vector
void ParseByCharater(std::string input, std::vector<std::string> & tokens, char flag);

std::string trimString(std::string const& str);

// Radian and Degree
double DToR(double degree);
double RToD(double Radian);

#endif