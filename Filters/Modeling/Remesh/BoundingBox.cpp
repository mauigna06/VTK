#include "BoundingBox.h"
#include "Face.h"
#include <algorithm>
#define EPSILON 1e-6

BoundingBox::BoundingBox():
min(Eigen::Vector3d::Zero()),
max(Eigen::Vector3d::Zero()),
extent(Eigen::Vector3d::Zero())
{
    
}

BoundingBox::BoundingBox(const Eigen::Vector3d& min0, const Eigen::Vector3d& max0):
min(min0),
max(max0)
{
    extent = max - min;
}

BoundingBox::BoundingBox(const Eigen::Vector3d& p):
min(p),
max(p)
{
    extent = max - min;
}

void BoundingBox::expandToInclude(const Eigen::Vector3d& p)
{
    if (min.x() > p.x()) min.x() = p.x();
    if (min.y() > p.y()) min.y() = p.y();
    if (min.z() > p.z()) min.z() = p.z();
    
    if (max.x() < p.x()) max.x() = p.x();
    if (max.y() < p.y()) max.y() = p.y();
    if (max.z() < p.z()) max.z() = p.z();
    
    extent = max - min;
}

void BoundingBox::expandToInclude(const BoundingBox& b)
{
    if (min.x() > b.min.x()) min.x() = b.min.x();
    if (min.y() > b.min.y()) min.y() = b.min.y();
    if (min.z() > b.min.z()) min.z() = b.min.z();
    
    if (max.x() < b.max.x()) max.x() = b.max.x();
    if (max.y() < b.max.y()) max.y() = b.max.y();
    if (max.z() < b.max.z()) max.z() = b.max.z();
    
    extent = max - min;
}

#include "BoundingBox.h"

#include <limits>

namespace vtkBotschKobbeltRemeshing
{

BoundingBox::BoundingBox()
{
    min.x() = std::numeric_limits<double>::max();
    min.y() = std::numeric_limits<double>::max();
    min.z() = std::numeric_limits<double>::max();

    max.x() = -std::numeric_limits<double>::max();
    max.y() = -std::numeric_limits<double>::max();
    max.z() = -std::numeric_limits<double>::max();

    extent = Eigen::Vector3d(0.0, 0.0, 0.0);
}

BoundingBox::BoundingBox(const Eigen::Vector3d& min0, const Eigen::Vector3d& max0)
{
    min = min0;
    max = max0;
    extent = max - min;
}

BoundingBox::BoundingBox(const Eigen::Vector3d& p)
{
    min = p;
    max = p;
    extent = Eigen::Vector3d(0.0, 0.0, 0.0);
}

void BoundingBox::expandToInclude(const Eigen::Vector3d& p)
{
    if (p.x() < min.x())
    {
        min.x() = p.x();
    }
    if (p.y() < min.y())
    {
        min.y() = p.y();
    }
    if (p.z() < min.z())
    {
        min.z() = p.z();
    }

    if (p.x() > max.x())
    {
        max.x() = p.x();
    }
    if (p.y() > max.y())
    {
        max.y() = p.y();
    }
    if (p.z() > max.z())
    {
        max.z() = p.z();
    }

    extent = max - min;
}

void BoundingBox::expandToInclude(const BoundingBox& b)
{
    if (b.min.x() < min.x())
    {
        min.x() = b.min.x();
    }
    if (b.min.y() < min.y())
    {
        min.y() = b.min.y();
    }
    if (b.min.z() < min.z())
    {
        min.z() = b.min.z();
    }

    if (b.max.x() > max.x())
    {
        max.x() = b.max.x();
    }
    if (b.max.y() > max.y())
    {
        max.y() = b.max.y();
    }
    if (b.max.z() > max.z())
    {
        max.z() = b.max.z();
    }

    extent = max - min;
}

int BoundingBox::maxDimension() const
{
    int k = 0;
    if (extent.y() > extent.x())
    {
        k = 1;
    }
    if (extent.z() > extent[k])
    {
        k = 2;
    }
    return k;
}

bool BoundingBox::intersect(const Eigen::Vector3d& p, const double dMin, double& dist) const
{
    const double& x = p.x();
    const double& y = p.y();
    const double& z = p.z();
    const double& minX = min.x();
    const double& minY = min.y();
    const double& minZ = min.z();
    const double& maxX = max.x();
    const double& maxY = max.y();
    const double& maxZ = max.z();

    dist = dMin;
    if (x < minX)
    {
        dist -= (minX - x) * (minX - x);
    }
    else if (x > maxX)
    {
        dist -= (x - maxX) * (x - maxX);
    }
    if (dist < 0)
    {
        return false;
    }

    if (y < minY)
    {
        dist -= (minY - y) * (minY - y);
    }
    else if (y > maxY)
    {
        dist -= (y - maxY) * (y - maxY);
    }
    if (dist < 0)
    {
        return false;
    }

    if (z < minZ)
    {
        dist -= (minZ - z) * (minZ - z);
    }
    else if (z > maxZ)
    {
        dist -= (z - maxZ) * (z - maxZ);
    }
    if (dist < 0)
    {
        return false;
    }

    return true;
}

} // namespace vtkBotschKobbeltRemeshing
