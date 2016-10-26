/*********************************************************************/
// PointCloud2D.h
// This the class holding 2D point cloud, including the random point
// cloud generator with different distributions.

// Generate sample sets:

// Generate Module:
// Construct Module given point number

// Generate Data given Module: 
// Copy Module->Add noise->Add outliers->Rotate->Translate

// To generate Module:
// Add outliers in Module

// Hang 11/19/2014
/*********************************************************************/

#pragma once

//#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <array>
#include "Eigen\Dense"

using Eigen::Vector2f;

class PointCloud2D
{
public:

    PointCloud2D(void);

    PointCloud2D(std::vector<Vector2f> pointCloud2D);

    PointCloud2D(int ponitNum);

    ~PointCloud2D(void);

    /*** Setters ***/
    // This is for manually picked points --- assume no outliers
    void SetPointCloud2D(std::vector<Vector2f> pointCloud2D) 
    {
        m_pointCloud2D = pointCloud2D;
        m_cloudCenter  = Vector2f(0.0f, 0.0f);
        m_xyRange      = 0.0f;
        m_inliersNum   = pointCloud2D.size();
        m_outliersNum  = 0;

        // Compute the model center
        for(auto it = pointCloud2D.begin(); it != pointCloud2D.end(); ++it)
        {
            m_cloudCenter += *it;
            m_xyRange = (std::max)(m_xyRange, (std::max)(it->x(), it->y()));
        }
        m_cloudCenter /= (float)pointCloud2D.size();
    }

    //
    void SetPointCloud2D(PointCloud2D points) 
    {
        m_pointCloud2D = points.GetPointCloud(); 
        m_cloudCenter  = points.GetCloudCenter(); 
        m_xyRange      = points.GetXYRange();
        m_inliersNum   = points.GetInliersNum();
        m_outliersNum  = points.GetOutliersNum();
    }

    // Clear the object
    void Clear()
    {
        if( !m_pointCloud2D.empty() )
            m_pointCloud2D.clear();
        m_cloudCenter = Vector2f(0.0f, 0.0f);
        m_xyRange     = 0.0f;
        m_inliersNum  = 0;
        m_outliersNum = 0;
    }

    // Generate random 2d points uniformly
    void RandGenUniform(int pointNum, float gridUnit = 256.0f, float pointDense = 10.0f);

    // Add white noise to the point cloud
    void AddWhiteNoise(float miu = 0.0f, float sigma = 2.0f);

    // Rotate and translate the point cloud
    void TransForm(Vector2f translateXY, float rotateAngle);

    // Add uniformly distributed random outliers --- append to the point cloud list
    void AddOutLiersUniform(int pointNum = 0);

    /*** Getters ***/
    // Return the point
    const Vector2f & GetPoint(int i) const { return m_pointCloud2D[i]; }
    // Return the point cloud
    const std::vector<Vector2f> & GetPointCloud() const { return m_pointCloud2D; }
    //
    const Vector2f & GetCloudCenter() const { return m_cloudCenter; }
    //
    const float & GetXYRange() const { return m_xyRange; }
    //
    const int & GetInliersNum() const  { return m_inliersNum; }
    const int & GetOutliersNum() const { return m_outliersNum; }
    int GetPointNum()                  { return m_inliersNum + m_outliersNum; }

private:

    // Point list
    std::vector<Vector2f> m_pointCloud2D;

    // Point cloud center
    Vector2f m_cloudCenter;

    // x-y coordinate range of point cloud
    // We consider range on x is the same as y
    float m_xyRange;

    int m_inliersNum;
    int m_outliersNum;
};

