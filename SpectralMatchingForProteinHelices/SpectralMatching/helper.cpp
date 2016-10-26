/*********************************************************************/
// helper.pp
// This is a class holding some helper functions, like line fitting.

// Hang 10/13/2014
/*********************************************************************/

#include "helper.h"

// Degree to Radian
double DToR(double degree)
{
    return degree * PI / 180.0;
}

// Radian to Degree
double RToD(double Radian)
{
    return Radian * 180.0 / PI;
}

// Use least square to fit a 3D line from input points
lineSeg FitLine(pointList points)
{
    // Protect by checking if "points" has size 0. ToDo.

    // Compute the average location of the points
    Vector3f avg(0,0,0);
    for(auto it = points.begin(); it != points.end(); ++it)
    {
        avg += *it;
    }
    avg /= (float)points.size();

    Vector3f dp(0,0,0);
    float sX2 = 0, sY2 = 0, sZ2 = 0, sXY = 0, sXZ = 0, sYZ = 0;
    // Construct the covariance matrix
    for(auto it = points.begin(); it != points.end(); ++it)
    {
        dp = *it - avg;
        sX2 += (dp[0] * dp[0]);
        sY2 += (dp[1] * dp[1]);
        sZ2 += (dp[2] * dp[2]);
        sXY += (dp[0] * dp[1]);
        sXZ += (dp[0] * dp[2]);
        sYZ += (dp[1] * dp[2]);
    }

    // Matrix in Eigen format
    // Since this is a Hessian matrix, 
    // we actually only need to construct the lower triangle.
    Matrix3f covM;
    covM << sX2, sXY, sXZ,
            sXY, sY2, sYZ,
            sXZ, sYZ, sZ2;

    // Solve for principle eigenvector (line direction) through Eigen
    SelfAdjointEigenSolver<Matrix3f> es;
    es.compute(covM);
    Vector3f lineDir = es.eigenvectors().col(2);

    // Obtain the line segment --- project the first and last point on helix to the line direction
    Vector3f p1 = avg + lineDir*(lineDir.dot(*(points.begin()) - avg));
    Vector3f p2 = avg + lineDir*(lineDir.dot(*(points.end()-1) - avg));

    lineSeg ls;
    ls[0] = p1; ls[1] = p2;

    return ls;
}

// EuclideanDistance of two 3D points: EuDis
double EuDis(Vector3f p1, Vector3f p2)
{
    return (p1 - p2).cast<double>().norm();
}
//float EuDis(Vector3f p1, Vector3f p2)
//{
//    return (p1 - p2).norm();
//}

// VectorAngle of two 3D vectors: VecAng
// Cast to double to increase precision
double VecAng(Vector3f v1, Vector3f v2)
{    
    return acos( ((v1.cast<double>()).normalized()).dot((v2.cast<double>()).normalized()) );
}
// For 2D vectors
double VecAng(Vector2f v1, Vector2f v2)
{    
    return acos( ((v1.cast<double>()).normalized()).dot((v2.cast<double>()).normalized()) );
}

//float VecAng(Vector3f v1, Vector3f v2)
//{    
//    //return acos( ((v1.cast<double>()).normalized()).dot((v2.cast<double>()).normalized()) );
//    return acosf( (v1.normalized()).dot(v2.normalized()) );
//}

// Parse an input string by given character. Put the parsed result into a given vector
void ParseByCharater(std::string input, std::vector<std::string> & tokens, char flag)
{
    std::istringstream ss(input);
    std::string token;

    while(std::getline(ss, token, flag)) 
    {
        if(token.size() > 0)
        {
            tokens.push_back( trimString(token) );
        }        
    }
}

// Trim the space out of head and tail
std::string trimString(std::string const& str)
{
    std::size_t first = str.find_first_not_of(' ');    

    if( first == std::string::npos )
    {
        std::cout<<"Try to trim a space only string...";
        // Return an empty string
        return "";
    }

    std::size_t last  = str.find_last_not_of(' ');

    return str.substr(first, last-first+1);
}