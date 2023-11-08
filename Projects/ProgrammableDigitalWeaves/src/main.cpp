#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/unproject_on_plane.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
#include "igl/colormap.h"
#include <iostream>
#include "../include/UI.h"
#include "../include/EoLRodSim.h"
#include "../include/Homogenization.h"
#include "../include/HybridC2Curve.h"

#define T double
#define dim 3


bool USE_VIEWER = true;

EoLRodSim<T, dim> eol_sim;
HybridC2Curve<T, dim> hybrid_curve;
Homogenization<T, dim> homogenizer(eol_sim);

using TV = Vector<T, dim>;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

Eigen::MatrixXd nodes;


Eigen::MatrixXd V_drawing;
Eigen::MatrixXi F_drawing;

static bool tileUnit = false;
static bool showUnit = false;
static bool showStretching = false;
static bool show_index = false;
static bool show_original = false;
static bool per_yarn = true;
static bool slide = false;
static bool draw_unit = false;

static bool show_rest = false;

static bool drawing = true;
static bool editing = !drawing;
static bool show_data_points = true;

static float theta_pbc = 0;
static float strain = 1.0;
static int n_rod_per_yarn = 4;
static int modes = 9;

static bool show_target = true;
static bool show_tunnel = false;
static bool show_bc = true;
static bool target_drawn = false;
static bool perturb = true;
bool reset = false;
enum TestCase{
    DrawUnit, StaticSolve, BatchRendering, InverseDesign, StaticSolveIncremental
};

const char* test_case_names[] = {
    "DrawUnit", "StaticSolve", "BatchRendering", "InverseDesign", "StaticSolveIncremental"
};
using RGBMat = Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>;

int main(int argc, char *argv[])
{
    TestCase test_current = StaticSolve;

    double t = 0.0;

    int n_faces = 20;

    std::vector<HybridC2Curve<T, dim>> hybrid_curves;

    Eigen::VectorXd deformed_backup;

    Eigen::MatrixXd evectors;
    auto loadEigenVectors = [&]()
    {
        std::ifstream in("./eigen_vectors.txt");
        int row, col;
        in >> row >> col;
        evectors.resize(row, col);
        double entry;
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                in >> evectors(i, j);
        in.close();
    };

    igl::opengl::glfw::Viewer viewer;
    
    igl::opengl::glfw::imgui::ImGuiMenu menu;

    auto updateScreen = [&](igl::opengl::glfw::Viewer& viewer)
    {
        viewer.data().clear();
        if (draw_unit)
        {
            if(V_drawing.size())
            {
                viewer.data().set_mesh(V_drawing, F_drawing);
            }
        }
        else
        {
            if(tileUnit)
            {
                eol_sim.buildPeriodicNetwork(V, F, C, show_rest);
                viewer.data().set_mesh(V, F);     
                viewer.data().set_colors(C);
                return;
            }
            else
            {
                eol_sim.generateMeshForRendering(V, F, Vector<T, dim>::Zero(), show_rest);
                viewer.data().set_mesh(V, F); 
            }
            if(showUnit)
                viewer.data().set_colors(C);
            if (per_yarn)
            {
                
                int n_rods = 0;
                for (auto& rod : eol_sim.Rods)
                    n_rods += rod->numSeg();
                if (show_rest)
                    C.resize(n_rods * n_faces * 2, 3);
                else
                    C.resize(n_rods * n_faces, 3);
                tbb::parallel_for(0, n_rods, [&](int rod_idx){
                    for(int i = 0; i < n_faces; i++)
                    {
                        // if( int(std::floor(T(i) / 2)) % 2 == 0)
                        if( i % 2 == 0)
                            C.row(rod_idx * n_faces + i) = Eigen::Vector3d(0, 1, 0);
                        else
                            C.row(rod_idx * n_faces + i) = Eigen::Vector3d(0, 0, 1);
                        if (show_rest)
                            C.row(n_rods * n_faces + rod_idx * n_faces + i) = Eigen::Vector3d(1, 0, 0);
                    }
                    });
                if (test_current == InverseDesign)
                {
                    viewer.data().clear();
                    tbb::parallel_for(0, n_rods, [&](int rod_idx){
                    for(int i = 0; i < n_faces; i++)
                    {
                        // if( int(std::floor(T(i) / 2)) % 2 == 0)
                        C.row(rod_idx * n_faces + i) = Eigen::Vector3d(0, 1, 0);
                        if (show_rest)
                            C.row(n_rods * n_faces + rod_idx * n_faces + i) = Eigen::Vector3d(1, 0, 0);
                    }
                    });
                    if (show_target && !target_drawn)
                    {
                        for (auto target : eol_sim.targets)
                        {
                            appendSphereMeshWithColor(V, F, C, 0.03, target / eol_sim.unit);
                        }
                        
                        target_drawn = false;
                    }
                    else
                    {
                        target_drawn = false;
                    }
                }
                tbb::parallel_for(0, n_rods, [&](int rod_idx){
                    for(int i = 0; i < n_faces; i++)
                    {
                        // if( int(std::floor(T(i) / 2)) % 2 == 0)
                        C.row(rod_idx * n_faces + i) = Eigen::Vector3d(0, 0.3, 1);
                        if (show_rest)
                            C.row(n_rods * n_faces + rod_idx * n_faces + i) = Eigen::Vector3d(1, 0, 0);
                    }
                });
                if (show_tunnel)
                {
                    viewer.data().clear();
                    for (auto crossing : eol_sim.rod_crossings)
                    {
                        if (crossing->is_fixed)
                            continue;
                        int node_idx = crossing->node_idx;
                        TV pos;
                        auto rod = eol_sim.Rods[crossing->rods_involved.front()];
                        rod->x(crossing->node_idx, pos);
                        
                        for (int i = 0; i < crossing->rods_involved.size(); i++)
                        {
                            if (crossing->sliding_ranges[i].norm() < 1e-6)
                                continue;
                            Eigen::Matrix3d R;
                            if (i % 2 == 0) 
                                R = rotationMatrixFromEulerAngle(0.0, 0.0, M_PI/ 2.0);
                            if (i % 2 == 1) 
                                R = rotationMatrixFromEulerAngle(M_PI/ 2.0, 0.0, M_PI/ 2.0);
                            appendTorusMeshWithColor(V, F, C, 0.01, pos / eol_sim.unit, R);        
                        }
                    }
                }
                if (show_bc)
                {
                    viewer.data().clear();
                    for (auto data : eol_sim.boundary_spheres)
                    {
                        TV center = data.first;
                        T r = data.second;
                        appendSphereMeshWithColor(V, F, C, r / eol_sim.unit, center / eol_sim.unit);
                    }
                }
                viewer.data().set_mesh(V, F); 
                viewer.data().set_colors(C);
                if(tileUnit)
                {
                    eol_sim.getColorPerYarn(C, n_rod_per_yarn);
                    C.conservativeResize(F.rows(), 3);
                    tbb::parallel_for(0, eol_sim.n_rods * n_faces, [&](int i){
                        for(int j = 1; j < std::floor(F.rows()/eol_sim.n_rods/n_faces); j++)
                        {
                            C.row(j * eol_sim.n_rods * n_faces + i) = C.row(i);
                        }
                    });
                    viewer.data().set_colors(C);
                }
            }
            if(show_original && !tileUnit)
            {
                Eigen::MatrixXd X, x;
                eol_sim.getEulerianDisplacement(X, x);
                for (int i = 0; i < X.rows(); i++)
                {
                    // viewer.data().add_edges(X.row(i), x.row(i), Eigen::RowVector3d(1, 1, 1));
                }
                viewer.data().add_points(X, Eigen::RowVector3d(1,1,1));
                // viewer.data().add_points(x, Eigen::RowVector3d(0,0,0));  
            }
        }
        
    };

    viewer.callback_key_pressed = 
        [&](igl::opengl::glfw::Viewer & viewer,unsigned int key,int mods)->bool
    {
        if (key == ' ')
        {
            if (draw_unit)
            {
                std::vector<Vector<T, dim>> points_on_curve;
                hybrid_curve.getLinearSegments(points_on_curve);
                appendCylinderMesh<dim>(viewer, V_drawing, F_drawing, points_on_curve, true);
            }
            else if(test_current == StaticSolve)
            {
                eol_sim.advanceOneStep();   
                deformed_backup = eol_sim.deformed_states;
                
                // eol_sim.rest_states = eol_sim.deformed_states;
                // eol_sim.resetScene();
                
            }
            else if (test_current == StaticSolveIncremental)
            {
                for (int step = 0; step < eol_sim.incremental_steps; step++)
                {
                    eol_sim.staticSolveIncremental(step);
                    updateScreen(viewer);
                    igl::writeOBJ("output/" + std::to_string(step) + ".obj", V, F);
                    
                }    
            }
            updateScreen(viewer);
            return true;
        }
        else if (key == 'a')
        {
            viewer.core().is_animating = !viewer.core().is_animating;
            return true;
        }
        else if (key == '1')
        {
            modes -= 1;
            modes = (modes + evectors.cols()) % evectors.cols();
            std::cout << "modes " << evectors.cols() - modes << std::endl;
            return true;
        }
        else if (key == '2')
        {
            if (reset)
            {
                Eigen::VectorXd dq(eol_sim.W.cols());
                std::ifstream in("./eigen_vector.txt");
                double value;
                int cnt = 0;
                while(in >> value)
                {
                    dq[cnt++] = value;
                }
                in.close();
                std::cout << dq.norm() << std::endl;
                std::cout << (eol_sim.W * dq).norm() << std::endl;
                eol_sim.deformed_states = deformed_backup + eol_sim.W * dq;
            }
            else
            {
                eol_sim.deformed_states = deformed_backup;
            }
            if (tileUnit)
            {
                eol_sim.buildPeriodicNetwork(V, F, C, show_rest);
                viewer.data().set_mesh(V, F);     
                viewer.data().set_colors(C);
            }
            else
            {
                eol_sim.generateMeshForRendering(V, F, Vector<T, dim>::Zero(), show_rest);
                viewer.data().set_mesh(V, F); 
            }
            reset = !reset;
            return true;
        }
        return false;
    };

    
    

    int n_test_case = sizeof(test_case_names)/sizeof(const char*);
    
    int selected = -1;
    int rod_idx = -1;
    double u0 = 0.0, x0 = 0.0, y0 = 0.0;

    static TestCase test = BatchRendering;
    // TestCase test_current = StaticSolve; // set to be a different from above or change the above one to be a random one

    // test_current = InverseDesign;

    auto setupScene = [&](igl::opengl::glfw::Viewer& viewer)
    {   
        if (test_current == DrawUnit)
        {
            draw_unit = true;
            // viewer.core().camera_zoom = 0.1;
        }
        else if (test_current == StaticSolve || test_current == StaticSolveIncremental)
        {
            homogenizer.testOneSample();
            draw_unit = false;
        }
        else if (test_current == InverseDesign)
        {
            // eol_sim.inverse();
            homogenizer.initialize();
            // eol_sim.inverse();
            draw_unit = false;
        }
        updateScreen(viewer);
    };

    
    
    if (test_current == BatchRendering)
    {
        menu.callback_draw_viewer_menu = [&](){
            viewer.core().align_camera_center(viewer.data().V, viewer.data().F);
        };
    }
    else if (test_current == StaticSolve)
    {
        viewer.plugins.push_back(&menu);
        
        menu.callback_draw_viewer_menu = [&]()
        {   
            // menu.draw_viewer_menu();
            if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
            {   
                ImGui::Combo("TestCase", (int *)(&test_current), test_case_names, n_test_case);
                // if(test != test_current)
                // {
                //     test = test_current;
                //     setupScene(viewer);
                // }
            }
            if (ImGui::CollapsingHeader("PeriodicBC", ImGuiTreeNodeFlags_DefaultOpen))
            {   
                if (ImGui::DragFloat("Angle", &(eol_sim.theta), 0.f, 0.1f, M_PI * 2.f))
                {
                    eol_sim.resetScene();
                    Vector<T, dim> strain_dir, ortho_dir;
                    eol_sim.setUniaxialStrain(eol_sim.theta, 1.1, strain_dir, ortho_dir);
                    eol_sim.advanceOneStep();
                    updateScreen(viewer);
                }
                if (ImGui::DragFloat("Strain", &(strain), 1.f, 0.02f, 1.1f))
                {
                    eol_sim.resetScene();
                    Vector<T, dim> strain_dir, ortho_dir;
                    eol_sim.setUniaxialStrain(eol_sim.theta, strain, strain_dir, ortho_dir);
                    eol_sim.advanceOneStep();
                    updateScreen(viewer);
                }
                
                if (ImGui::Checkbox("Tunnel", &eol_sim.add_contact_penalty))
                {
                    if(eol_sim.add_contact_penalty)
                        eol_sim.releaseCrossing();
                    else
                        eol_sim.fixCrossing();
                    
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("ShowTunnel", &show_tunnel))
                {
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Shearing", &eol_sim.add_shearing))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Stretching", &eol_sim.add_stretching))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Bending", &eol_sim.add_bending))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Twisting", &eol_sim.add_twisting))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("RigidJoint", &eol_sim.add_rigid_joint))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("ShowRest", &show_rest))
                {
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("TileUnit", &tileUnit))
                {
                    updateScreen(viewer);
                }
                // if (ImGui::Checkbox("ShowIndex", &show_index))
                // {
                //     if(show_index)
                //     {
                //         for (int i = 0; i < eol_sim.n_nodes; i++)
                //             viewer.data().add_label(Eigen::Vector3d(eol_sim.q(0, i), eol_sim.q(1, i), 0), std::to_string(i));
                //         viewer.data().show_custom_labels = true;
                //     }
                // }
                // if (ImGui::Checkbox("ShowEulerianRest", &show_original))
                // {   
                //     updateScreen(viewer);
                // }
            }
            if (ImGui::CollapsingHeader("ColorScheme", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Checkbox("ShowStretching", &showStretching))
                {
                    viewer.data().clear();
                    viewer.data().set_mesh(V, F);
                    if (showStretching)
                    {
                        eol_sim.showStretching(C);
                        viewer.data().set_colors(C);
                    }   
                }
                if (ImGui::Checkbox("PerYarn", &per_yarn))
                {
                    updateScreen(viewer);
                }   
            }
            if (ImGui::Button("Solve", ImVec2(-1,0)))
            {
                eol_sim.advanceOneStep();
                
                updateScreen(viewer);
            }
            if (ImGui::Button("SaveMesh", ImVec2(-1,0)))
            {
                igl::writeOBJ("current_mesh.obj", V, F);
            }
            if (ImGui::Button("Reset", ImVec2(-1,0)))
            {
                eol_sim.resetScene();
                updateScreen(viewer);
            }
        };
    }
    else if (test_current == InverseDesign)
    {
        viewer.plugins.push_back(&menu);
        
        menu.callback_draw_viewer_menu = [&]()
        {   
            if (ImGui::CollapsingHeader("Inverse", ImGuiTreeNodeFlags_DefaultOpen))
            {   
                if (ImGui::Checkbox("ShowTarget", &show_target))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
            }
            if (ImGui::CollapsingHeader("Configurations", ImGuiTreeNodeFlags_DefaultOpen))
            {   
                if (ImGui::Checkbox("RegularizeEulerian", &eol_sim.add_eularian_reg))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Tunnel", &eol_sim.add_contact_penalty))
                {
                    if(eol_sim.add_contact_penalty)
                        eol_sim.releaseCrossing();
                    else
                        eol_sim.fixCrossing();
                    
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("ShowTunnel", &show_tunnel))
                {
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("ShowBCSphere", &show_bc))
                {
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Stretching", &eol_sim.add_stretching))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Bending", &eol_sim.add_bending))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Twisting", &eol_sim.add_twisting))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("RigidJoint", &eol_sim.add_rigid_joint))
                {
                    eol_sim.resetScene();
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("ShowRest", &show_rest))
                {
                    updateScreen(viewer);
                }
            }
            if (ImGui::Button("Solve", ImVec2(-1,0)))
            {
                eol_sim.advanceOneStep();
                
                updateScreen(viewer);
            }
            if (ImGui::Button("Reset", ImVec2(-1,0)))
            {
                eol_sim.resetScene();
                updateScreen(viewer);
            }
        };
    }
    else if (test_current == DrawUnit)
    {
        viewer.plugins.push_back(&menu);
        menu.callback_draw_viewer_menu = [&]()
        {
            if (ImGui::CollapsingHeader("CurveIO", ImGuiTreeNodeFlags_DefaultOpen))
            {
                float w = ImGui::GetContentRegionAvailWidth();
                float p = ImGui::GetStyle().FramePadding.x;
                if (ImGui::Button("Load##Curve##Data##Points", ImVec2((w-p)/2.f, 0)))
                {
                    std::string fname = igl::file_dialog_open();
                    if (fname.length() != 0)
                    {
                        hybrid_curve.data_points.clear();
                        std::ifstream in(fname);
                        double x, y;
                        while(in >> x >> y)
                        {
                            if constexpr (dim == 2)
                                hybrid_curve.data_points.push_back(Vector<T, dim>(x, y));
                            else
                                hybrid_curve.data_points.push_back(Vector<T, dim>(x, y, 0));
                        }
                        in.close();
                    }
                }
                ImGui::SameLine(0, p);
                if (ImGui::Button("Save##Curve##Data##Points", ImVec2((w-p)/2.f, 0)))
                {
                    std::string fname = igl::file_dialog_save();

                    if (fname.length() != 0)
                    {
                        std::ofstream out(fname);
                        for (auto pt : hybrid_curve.data_points)
                        {
                            out << pt.transpose() << std::endl;
                        }
                        out.close();
                    }
                }
            }

            if (ImGui::CollapsingHeader("Drawing Options", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::DragInt("SubDivision", &(hybrid_curve.sub_div), 1.f, 8, 64))
                {
                    std::vector<Vector<T, dim>> points_on_curve;
                    hybrid_curve.getLinearSegments(points_on_curve);
                    appendCylinderMesh<dim>(viewer, V_drawing, F_drawing, points_on_curve);
                    updateScreen(viewer);
                }
                if (ImGui::Checkbox("Drawing", &drawing))
                {
                    editing = !drawing;
                }
                if (ImGui::Checkbox("Editing", &editing))
                {
                    drawing = !editing;
                }
                if (ImGui::Checkbox("ShowDataPoints", &show_data_points))
                {

                }
                float w = ImGui::GetContentRegionAvailWidth();
                float p = ImGui::GetStyle().FramePadding.x;
                if (ImGui::Button("Add##Curve", ImVec2((w-p)/2.f, 0)))
                {
                    hybrid_curves.push_back(hybrid_curve);
                    hybrid_curve = HybridC2Curve<T, dim>();
                }
                ImGui::SameLine(0, p);
                if (ImGui::Button("Remove##Curve", ImVec2((w-p)/2.f, 0)))
                {
                    hybrid_curves.pop_back();
                }
            }
        };
    }
    
    auto draw_unit_func = [&](igl::opengl::glfw::Viewer& viewer, int button, int)->bool
    {
        if (!drawing)
            return false;
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        Eigen::Vector4f eye_n = (viewer.core().view).inverse().col(3);
        // Eigen::Vector3d point;
        // igl::unproject_on_plane(Eigen::Vector2d(x,y), viewer.core().proj*viewer.core().view, viewer.core().viewport, eye_n, point);
        
        if (button == 0) // left button
        {
            if constexpr (dim == 2)
                hybrid_curve.data_points.push_back(Vector<T, 2>(x, y));
            else
                hybrid_curve.data_points.push_back(Vector<T, 3>(x, y, 0));
            std::vector<Vector<T, dim>> points_on_curve;
            hybrid_curve.getLinearSegments(points_on_curve);
            appendCylinderMesh<dim>(viewer, V_drawing, F_drawing, points_on_curve);
            if (show_data_points)
            {
                for (auto pt : hybrid_curve.data_points)
                {
                    Eigen::Vector3d point;
                    igl::unproject_on_plane(pt, viewer.core().proj*viewer.core().view, viewer.core().viewport, eye_n, point);
                    appendSphereMesh(V_drawing, F_drawing, 0.05, point);
                }
            }
        }
        else if (button == 2) // right button
        {
            // removeSphereMesh(V_drawing, F_drawing);
            if(hybrid_curve.data_points.size())
            {
                hybrid_curve.data_points.pop_back();
                std::vector<Vector<T, dim>> points_on_curve;
                hybrid_curve.getLinearSegments(points_on_curve);
                appendCylinderMesh<dim>(viewer, V_drawing, F_drawing, points_on_curve, true);
                if (show_data_points)
                {
                    for (auto pt : hybrid_curve.data_points)
                    {
                        Eigen::Vector3d point;
                        igl::unproject_on_plane(pt, viewer.core().proj*viewer.core().view, viewer.core().viewport, eye_n, point);
                        appendSphereMesh(V_drawing, F_drawing, 0.05, point);
                    }
                }
            }
        }
        updateScreen(viewer);
        return false;
    };

    if (test_current == DrawUnit)
        viewer.callback_mouse_down = draw_unit_func;
    else
        viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
        {
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            if(eol_sim.new_frame_work)
            {
                for (auto& rod : eol_sim.Rods)
                {
                    for (int node_idx : rod->indices)
                    {
                        Vector<T, dim> xsim;
                        rod->x(node_idx, xsim);
                        Eigen::MatrixXd x3d(1, 3); x3d.setZero();
                        x3d.row(0).template segment<dim>(0) = xsim / eol_sim.unit;

                        Eigen::MatrixXd pxy(1, 3);

                        igl::project(x3d, viewer.core().view, viewer.core().proj, viewer.core().viewport, pxy);
                        
                        if(abs(pxy.row(0)[0]-x)<20 && abs(pxy.row(0)[1]-y)<20)
                        {
                            selected = node_idx;
                            rod_idx = rod->rod_id;
                            x0 = x;
                            y0 = y;
                            std::cout << "selected " << selected << std::endl;
                            return true;
                        }
                    }
                }
                return false;
            }
            else
            {
                Eigen::MatrixXd pxy = eol_sim.q.transpose().block(0, 0, eol_sim.n_nodes, 2) / eol_sim.unit;
                Eigen::MatrixXd rod_v(eol_sim.n_nodes, 3);
                rod_v.setZero();
                rod_v.block(0, 0, eol_sim.n_nodes, 2) = pxy;

                igl::project(rod_v, viewer.core().view, viewer.core().proj, viewer.core().viewport, pxy);

                for(int i=0; i<pxy.rows(); ++i)
                {
                    if(abs(pxy.row(i)[0]-x)<20 && abs(pxy.row(i)[1]-y)<20)
                    {
                        selected = i;
                        x0 = x;
                        y0 = y;
                        std::cout << "selected " << selected << std::endl;
                        return true;
                    }
                }
                return false;
            }
        };
    

    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	  {
		if(selected!=-1)
		{
			selected = -1;
            if (eol_sim.new_frame_work)
                eol_sim.rest_states = eol_sim.deformed_states;
            else
                eol_sim.q0 = eol_sim.q;
            // eol_sim.q0.transpose().block(0, 0, eol_sim.n_nodes, dim)  = eol_sim.q.transpose().block(0, 0, eol_sim.n_nodes, dim);
			return true;
		}
	    return false;
	  };

    if (test_current == DrawUnit)
    {
        
        viewer.callback_mouse_move =
        [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
        {
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;

            for (auto& rod : eol_sim.Rods)
            {
                
            }

            return false;
        };
    }
    else
        viewer.callback_mouse_move =
        [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
        {
            if(selected!=-1)
            {
                double x = viewer.current_mouse_x;
                double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            
                double delta_x = (x - x0) / viewer.core().viewport(2) * eol_sim.unit;
                double delta_y = (y - y0) / viewer.core().viewport(3) * eol_sim.unit;

                if (eol_sim.new_frame_work)
                {
                    // Vector<int, dim + 1> offset = eol_sim.Rods[rod_idx]->offset_map[selected];
                    // eol_sim.dirichlet_dof[eol_sim.Rods[rod_idx]->reduced_map[offset[0]]] = delta_x;
                    // eol_sim.dirichlet_dof[eol_sim.Rods[rod_idx]->reduced_map[offset[1]]] = delta_y;
                    // eol_sim.dirichlet_dof[eol_sim.Rods[rod_idx]->reduced_map[offset[2]]] = 0;


                    // eol_sim.Rods[0]->fixPointLagrangian(eol_sim.Rods[0]->indices.size() - 2, 
                    //         Vector<T, dim>(delta_x, delta_y, 0.0), 
                    //         eol_sim.dirichlet_dof);
                    // eol_sim.Rods[0]->fixPointLagrangian(eol_sim.Rods[0]->indices.size() - 1, 
                    //         Vector<T, dim>(delta_x, delta_y, 0.0), 
                    //         eol_sim.dirichlet_dof); 
                    
                    // 3d finger
                    // eol_sim.Rods[16]->fixPointLagrangian(eol_sim.Rods[16]->indices.size() - 1, 
                    //         Vector<T, dim>(0, delta_y, 0.0), 
                    //         eol_sim.dirichlet_dof);

                    //3d gripper
                    // eol_sim.Rods[1]->fixPointLagrangian(eol_sim.Rods[1]->indices.size() - 1, 
                    //         Vector<T, dim>(0, delta_y, 0.0), 
                    //         eol_sim.dirichlet_dof);

                    // eol_sim.Rods[7]->fixPointLagrangian(eol_sim.Rods[7]->indices.front(), 
                    //         Vector<T, dim>(0, delta_y, 0.0), 
                    //         eol_sim.dirichlet_dof);

                    T dz = perturb ? 1e-3 * eol_sim.unit : 0.0;
                    eol_sim.Rods[3]->fixPointLagrangian(eol_sim.Rods[3]->indices.back(), 
                            Vector<T, dim>(0, delta_y, perturb), 
                            eol_sim.dirichlet_dof);

                    eol_sim.advanceOneStep();
                    updateScreen(viewer);
                    perturb = false;
                    return true;
                }
                else
                {
                    // eol_sim.q(dim, 3) = u0;
                    Eigen::VectorXd delta_dof(4); delta_dof.setZero();
                    auto zero_delta = delta_dof;
                    Eigen::VectorXd mask_dof(4); mask_dof.setZero();
                    
                    // delta_dof(0) = delta_x * eol_sim.unit;
                    delta_dof(1) = delta_y * eol_sim.unit;
                    
                    // mask_dof(0) = 1;
                    mask_dof(1) = 1;
                    mask_dof(2) = 1;
                    mask_dof(3) = 1;
                    
                    eol_sim.dirichlet_data[selected] = std::make_pair(delta_dof, mask_dof);
                    eol_sim.advanceOneStep();
                    updateScreen(viewer);
                    std::cout << delta_x << " " << delta_y << std::endl;
                    return true;
                }
            }

            return false;
        };

    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
    {
        if(viewer.core().is_animating)
        {
            eol_sim.deformed_states = deformed_backup + eol_sim.W * evectors.col(modes) * std::sin(t);
            // Eigen::VectorXd dq = eol_sim.W.transpose() * (eol_sim.deformed_states - eol_sim.rest_states);
            // std::cout << eol_sim.computeTotalEnergy(dq) << std::endl;
            t += 0.1;
            if (tileUnit)
            {
                eol_sim.buildPeriodicNetwork(V, F, C, show_rest);
                viewer.data().set_mesh(V, F);     
                viewer.data().set_colors(C);
            }
            else
            {
                eol_sim.generateMeshForRendering(V, F, Vector<T, dim>::Zero(), show_rest);
                viewer.data().set_mesh(V, F); 
            }
            // std::cout << " per frame" << std::endl;
        }
        return false;
    };
    //================== Run GUI ==================
    
    
    // loadEigenVectors();

    int width = 800, height = 800;
    
    if (test_current == BatchRendering)
    {
        
    }
    else if (test_current == StaticSolve || test_current == InverseDesign || test_current == StaticSolveIncremental)
    {
        // if (test_current == InverseDesign)
            viewer.core().background_color.setOnes();
        viewer.data().set_face_based(true);
        viewer.data().shininess = 1.0;
        viewer.data().point_size = 25.0;
        setupScene(viewer);
        // viewer.data().show_lines = false;
        viewer.core().align_camera_center(V);
        viewer.core().animation_max_fps = 24.;

        viewer.launch();
    }
    else if (test_current == DrawUnit)
    {
        
        viewer.data().set_face_based(true);
        viewer.data().shininess = 1.0;
        // viewer.data().point_size = 25.0;
        setupScene(viewer);
        // viewer.launch();
    }

    //================== Run Diff Test ==================
    // eol_sim.buildPlanePeriodicBCScene3x3();
    // homogenizer.initialize();
    // eol_sim.buildPlanePeriodicBCScene3x3();
    // eol_sim.derivativeTest();
    // eol_sim.runDerivativeTest();

    // eol_sim.resetScene();
    // homogenizer.computeYoungsModulusPoissonRatioBatch();
    // homogenizer.fitComplianceFullTensor();
    // homogenizer.fitComplianceTensor();
    return 0;
}