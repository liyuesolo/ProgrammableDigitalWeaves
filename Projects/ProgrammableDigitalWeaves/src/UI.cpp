#include "../include/UI.h"

void appendSphereMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double scale, Vector<double, 3> shift)
{
    Eigen::MatrixXd v_sphere;
    Eigen::MatrixXi f_sphere;

    igl::readOBJ("../../../Projects/ProgrammableDigitalWeaves/Data/sphere162.obj", v_sphere, f_sphere);

    v_sphere = v_sphere * scale;

    tbb::parallel_for(0, (int)v_sphere.rows(), [&](int row_idx){
        v_sphere.row(row_idx) += shift;
    });

    int n_vtx_prev = V.rows();
    int n_face_prev = F.rows();

    tbb::parallel_for(0, (int)f_sphere.rows(), [&](int row_idx){
        f_sphere.row(row_idx) += Eigen::Vector3i(n_vtx_prev, n_vtx_prev, n_vtx_prev);
    });

    V.conservativeResize(V.rows() + v_sphere.rows(), 3);
    F.conservativeResize(F.rows() + f_sphere.rows(), 3);

    V.block(n_vtx_prev, 0, v_sphere.rows(), 3) = v_sphere;
    F.block(n_face_prev, 0, f_sphere.rows(), 3) = f_sphere;
}

void appendSphereMeshWithColor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C,
        double scale, Vector<double, 3> shift)
{
    Eigen::MatrixXd v_sphere;
    Eigen::MatrixXi f_sphere;
    
    igl::readOBJ("../../../Projects/ProgrammableDigitalWeaves/Data/sphere162.obj", v_sphere, f_sphere);
    
    Eigen::MatrixXd c_sphere(f_sphere.rows(), f_sphere.cols());
    
    v_sphere = v_sphere * scale;

    tbb::parallel_for(0, (int)v_sphere.rows(), [&](int row_idx){
        v_sphere.row(row_idx) += shift;
    });

    int n_vtx_prev = V.rows();
    int n_face_prev = F.rows();

    tbb::parallel_for(0, (int)f_sphere.rows(), [&](int row_idx){
        f_sphere.row(row_idx) += Eigen::Vector3i(n_vtx_prev, n_vtx_prev, n_vtx_prev);
    });
    tbb::parallel_for(0, (int)f_sphere.rows(), [&](int row_idx){
        c_sphere.row(row_idx) = Eigen::Vector3d(1, 0, 0);
    });

    V.conservativeResize(V.rows() + v_sphere.rows(), 3);
    F.conservativeResize(F.rows() + f_sphere.rows(), 3);
    C.conservativeResize(C.rows() + f_sphere.rows(), 3);


    V.block(n_vtx_prev, 0, v_sphere.rows(), 3) = v_sphere;
    F.block(n_face_prev, 0, f_sphere.rows(), 3) = f_sphere;
    C.block(n_face_prev, 0, f_sphere.rows(), 3) = c_sphere;
    
}

void appendTorusMeshWithColor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C,
        double scale, Vector<double, 3> shift, Eigen::Matrix3d R)
{
    Eigen::MatrixXd v_sphere;
    Eigen::MatrixXi f_sphere;
    

    // igl::readOBJ("/home/yueli/Documents/ETH/WuKong/Projects/DigitalFabrics/Data/sphere.obj", v_sphere, f_sphere);
    igl::readOBJ("../../../Projects/ProgrammableDigitalWeaves/Data/torus.obj", v_sphere, f_sphere);
    
    Eigen::MatrixXd c_sphere(f_sphere.rows(), f_sphere.cols());
    
    v_sphere = v_sphere * scale;

    tbb::parallel_for(0, (int)v_sphere.rows(), [&](int row_idx){
        v_sphere.row(row_idx) = (R * v_sphere.row(row_idx).transpose()) + shift;
    });

    int n_vtx_prev = V.rows();
    int n_face_prev = F.rows();

    tbb::parallel_for(0, (int)f_sphere.rows(), [&](int row_idx){
        f_sphere.row(row_idx) += Eigen::Vector3i(n_vtx_prev, n_vtx_prev, n_vtx_prev);
    });

    tbb::parallel_for(0, (int)f_sphere.rows(), [&](int row_idx){
        c_sphere.row(row_idx) = Eigen::Vector3d(177, 156, 217) / 255.0;
    });

    V.conservativeResize(V.rows() + v_sphere.rows(), 3);
    F.conservativeResize(F.rows() + f_sphere.rows(), 3);
    C.conservativeResize(C.rows() + f_sphere.rows(), 3);


    V.block(n_vtx_prev, 0, v_sphere.rows(), 3) = v_sphere;
    F.block(n_face_prev, 0, f_sphere.rows(), 3) = f_sphere;
    C.block(n_face_prev, 0, f_sphere.rows(), 3) = c_sphere;
}

void removeSphereMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if(V.rows() > 0)
    {
        Eigen::MatrixXd v_sphere;
        Eigen::MatrixXi f_sphere;

        igl::readOBJ("../../../Projects/ProgrammableDigitalWeaves/Data/sphere162.obj", v_sphere, f_sphere);

        int n_vtx_prev = V.rows();
        int n_face_prev = F.rows();

        V.conservativeResize(V.rows() - v_sphere.rows(), 3);
        F.conservativeResize(F.rows() - f_sphere.rows(), 3);
    }
}

void removeSphereMeshWithColor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C)
{
    if(V.rows() > 0)
    {
        Eigen::MatrixXd v_sphere;
        Eigen::MatrixXi f_sphere;

        igl::readOBJ("../../../Projects/ProgrammableDigitalWeaves/Data/sphere162.obj", v_sphere, f_sphere);

        int n_vtx_prev = V.rows();
        int n_face_prev = F.rows();

        V.conservativeResize(V.rows() - v_sphere.rows(), 3);
        F.conservativeResize(F.rows() - f_sphere.rows(), 3);
        C.conservativeResize(C.rows() - f_sphere.rows(), 3);
    }
}
