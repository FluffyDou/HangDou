/*********************************************************************/
// PointCloud2D.h
// This the class holding 2D point cloud, including the random point
// cloud generator with different distributions.

// Hang 11/19/2014
/*********************************************************************/

#include "PointCloud2D.h"


PointCloud2D::PointCloud2D(void)
{
    Clear();
}

PointCloud2D::PointCloud2D(int pointNum)
{
    RandGenUniform(pointNum);
}

PointCloud2D::PointCloud2D(std::vector<Vector2f> pointCloud2D)
{
    SetPointCloud2D(pointCloud2D);
}

PointCloud2D::~PointCloud2D(void)
{
}

// Generate random 2d points uniformly --- default block unit size: 256
void PointCloud2D::RandGenUniform(int pointNum, float gridUnit, float pointDense)
{
    // Allocate memory
    m_pointCloud2D = std::vector<Vector2f>(pointNum);
    m_cloudCenter  = Vector2f(0.0f, 0.0f);

    // Compute x-y coordinate range
    float xyRange = gridUnit * sqrtf( (float)pointNum / pointDense );
    m_xyRange =xyRange;

    // Set up uniformly distributed random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> uniformDis(0.0, xyRange);

    // Generate uniformly distributed random points
    for (int i = 0; i < pointNum; ++i) 
    {
        float x = uniformDis(gen);
        float y = uniformDis(gen);

        m_pointCloud2D[i] = Vector2f(x, y);
        // Update the point cloud center
        m_cloudCenter    += m_pointCloud2D[i];
    }

    m_cloudCenter /= (float)pointNum;

    // Record inlier number
    m_inliersNum  = pointNum;
    m_outliersNum = 0;
}


// Rotate and translate the point cloud
void PointCloud2D::TransForm(Vector2f translateXY, float rotateAngle)
{
    // New point cloud center
    Vector2f newCloudCenter(0.0f, 0.0f);

    // For each point, translate to cloud center first.
    // Rotate and than translate back. Finally translate again.
    for(auto it = m_pointCloud2D.begin(); it != m_pointCloud2D.end(); ++it)
    {
        *it = Eigen::Translation2f(translateXY) * Eigen::Translation2f(m_cloudCenter) * Eigen::Rotation2Df(rotateAngle) * Eigen::Translation2f(-m_cloudCenter) * (*it);
        // Update the new point cloud center
        newCloudCenter += *it;
    }

    m_cloudCenter = newCloudCenter/(float)m_pointCloud2D.size();
}


// Add white noise to the point cloud
void PointCloud2D::AddWhiteNoise(float miu, float sigma)
{
    if( 0 == sigma)
        return;

    // New point cloud center
    Vector2f newCloudCenter(0.0f, 0.0f);

    // Set up Gaussian distributed random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> normalDis(miu, sigma);

    // Loop to disturb each point
    for(auto it = m_pointCloud2D.begin(); it != m_pointCloud2D.end(); ++it)
    {
        float xShift = normalDis(gen);
        float yShift = normalDis(gen);

        *it += Vector2f(xShift, yShift);
        // Update the new point cloud center
        newCloudCenter += *it;
    }

    m_cloudCenter = newCloudCenter/(float)m_pointCloud2D.size();
}


// Add uniformly distributed random outliers --- append to the point cloud list
// ATTENTION! Add outliers before rotate and translate.
void PointCloud2D::AddOutLiersUniform(int pointNum)
{
    m_outliersNum = pointNum;

    if( 0 == pointNum )
        return;

    std::vector<Vector2f> outliers(pointNum);

    // Set up uniformly distributed random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> uniformDis(0.0, m_xyRange);

    for(auto it = outliers.begin(); it != outliers.end(); ++it)
    {
        float x = uniformDis(gen);
        float y = uniformDis(gen);

        *it = Vector2f(x,y);
    }

    m_pointCloud2D.insert(m_pointCloud2D.end(), outliers.begin(), outliers.end());

    // Update the model center
    m_cloudCenter = Vector2f(0.0f, 0.0f);
    for(auto it = m_pointCloud2D.begin(); it != m_pointCloud2D.end(); ++it)
    {
        m_cloudCenter += *it;
    }
    m_cloudCenter /= (float)m_pointCloud2D.size();
}