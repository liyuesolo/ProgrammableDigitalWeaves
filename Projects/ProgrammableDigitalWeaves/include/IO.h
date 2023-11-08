#ifndef IO_H
#define IO_H

#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "VecMatDef.h"

template<class T>
void extrudeMeshToObj(const std::vector<Vector<T, 2>>& points, 
    const std::vector<Vector<int, 2>>& edges,
    const std::string& filename, T z)
{
    using TV = Vector<T, 3>;

    std::ofstream out(filename);
    for (auto & point : points)
    {
        out << "v " << point[0] << " " << point[1] << " 0" << std::endl;
    }
    for (auto & point : points)
    {
        out << "v " << point[0] << " " << point[1] << " " << z << std::endl;
    }
    for (auto& edge: edges)
    {
        int vi = edge[0], vj = edge[1], vk = edge[0] + points.size(), vl = edge[1] + points.size();
        vi++; vj++; vk++; vl++;
        out << "f " << vi << " " << vj << " " << vl << std::endl;
        out << "f " << vl << " " << vk << " " << vi << std::endl;
    }
    out.close();
}

#endif