#include "../include/EoLRodSim.h"
#include "igl/colormap.h"
#include "igl/readOBJ.h"

#include "../autodiff/EoLRodStretchingEnergy.h"

template<class T, int dim>
void EoLRodSim<T, dim>::generateMeshForRendering(Eigen::MatrixXd& V, Eigen::MatrixXi& F, TV shift, bool show_rest_shape)
{
    // std::cout << "generate mesh for rendering" << std::endl;
    int n_div = 10;
    
    T theta = 2.0 * EIGEN_PI / T(n_div);
    TV3Stack points = TV3Stack::Zero(3, n_div);

    // T visual_R = 0.01;
    
    // bottom face vertices
    for(int i = 0; i < n_div; i++)
        points.col(i) = TV3(visual_R * std::cos(theta * T(i)), 0.0, visual_R*std::sin(theta*T(i)));
    

    int n_rod_total = 0;
    for(auto& rod : Rods)
    {
        n_rod_total += rod->numSeg();
    }
    // std::cout << "n_rod_total: " << n_rod_total << std::endl;
    int rod_offset_v = n_div * 2;
    int rod_offset_f = n_div * 2;
    
    V.resize(n_rod_total * rod_offset_v, 3);
    V.setZero();
    F.resize(n_rod_total * rod_offset_f, 3);
    F.setZero();
    int rod_cnt = 0;
    
    for(auto& rod : Rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx){
            int rov = rod_cnt * rod_offset_v;
            int rof = rod_cnt * rod_offset_f;

            TV vtx_from_TV, vtx_to_TV;
            rod->x(node_i, vtx_from_TV);
            rod->x(node_j, vtx_to_TV);
            // rod->X(node_i, vtx_from_TV);
            // rod->X(node_j, vtx_to_TV);
            vtx_from_TV += shift;
            vtx_to_TV += shift;
            vtx_from_TV /= unit; vtx_to_TV /= unit;
            // std::cout << rod_cnt << " " <<  node_i << " " << vtx_from_TV.transpose() << " node j " << node_j << " " << vtx_to_TV.transpose() << std::endl;
            TV3 vtx_from = TV3::Zero();
            TV3 vtx_to = TV3::Zero();
            if constexpr (dim == 3)
            {
                vtx_from = vtx_from_TV;
                vtx_to = vtx_to_TV;
            }
            else
            {
                vtx_from = TV3(vtx_from_TV[0], vtx_from_TV[1], 0);
                vtx_to = TV3(vtx_to_TV[0], vtx_to_TV[1], 0);
            }

            TV3 normal_offset = TV3::Zero();

            vtx_from += normal_offset * R;
            vtx_to += normal_offset * R;
            
            TV3 axis_world = vtx_to - vtx_from;
            TV3 axis_local(0, axis_world.norm(), 0);

            
            TM3 R_local = Eigen::AngleAxis<T>(rod->reference_angles[rod_idx] + rod->reference_twist[rod_idx], axis_local.normalized()).toRotationMatrix();
            // TM3 R_local = Eigen::AngleAxis<T>(rod->reference_angles[rod_idx], axis_local.normalized()).toRotationMatrix();
        
            TM3 R = Eigen::Quaternion<T>().setFromTwoVectors(axis_world.normalized(), axis_local.normalized()).toRotationMatrix();
            
            // R =  R_local.transpose() * R;
                        
            for(int i = 0; i < n_div; i++)
            {
                for(int d = 0; d < 3; d++)
                {
                    V(rov + i, d) = points.col(i)[d];
                    V(rov + i+n_div, d) = points.col(i)[d];
                    if (d == 1)
                        V(rov + i+n_div, d) += axis_world.norm();
                }

                // central vertex of the top and bottom face
                if (rod_idx == 0)
                    V.row(rov + i) = ((V.row(rov + i) - Eigen::Vector3d(0, visual_R, 0).transpose()) * R).transpose() + vtx_from;
                else
                    V.row(rov + i) = (V.row(rov + i) * R).transpose() + vtx_from;
                if (rod_idx == rod->numSeg() - 1)
                    V.row(rov + i + n_div) = ((V.row(rov + i + n_div) + Eigen::Vector3d(0, visual_R, 0).transpose()) * R).transpose() + vtx_from;
                else
                    V.row(rov + i + n_div) = (V.row(rov + i + n_div) * R).transpose() + vtx_from;
                
                F.row(rof + i*2 ) = IV3(rov + i, rov + i+n_div, rov + (i+1)%(n_div));
                F.row(rof + i*2 + 1) = IV3(rov + (i+1)%(n_div), rov + i+n_div, rov + (i+1)%(n_div) + n_div);
            }
            rod_cnt++;
        });
        
    }
    if (!show_rest_shape)
        return;
    int v_off = n_rod_total * rod_offset_v;
    V.conservativeResize(V.rows() * 2, 3);
    F.conservativeResize(F.rows() * 2, 3);
    rod_cnt = 0;
    for(auto& rod : Rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx){
            int rov = v_off + rod_cnt * rod_offset_v;
            int rof = v_off + rod_cnt * rod_offset_f;

            TV vtx_from_TV, vtx_to_TV;
            rod->X(node_i, vtx_from_TV);
            rod->X(node_j, vtx_to_TV);
            vtx_from_TV /= unit; vtx_to_TV /= unit;
            // std::cout << rod_cnt << " " <<  node_i << " " << vtx_from_TV.transpose() << " node j " << node_j << " " << vtx_to_TV.transpose() << std::endl;
            TV3 vtx_from = TV3::Zero();
            TV3 vtx_to = TV3::Zero();
            if constexpr (dim == 3)
            {
                vtx_from = vtx_from_TV;
                vtx_to = vtx_to_TV;
            }
            else
            {
                vtx_from = TV3(vtx_from_TV[0], vtx_from_TV[1], 0);
                vtx_to = TV3(vtx_to_TV[0], vtx_to_TV[1], 0);
            }

            TV3 normal_offset = TV3::Zero();

            vtx_from += normal_offset * R;
            vtx_to += normal_offset * R;
            
            TV3 axis_world = vtx_to - vtx_from;
            TV3 axis_local(0, axis_world.norm(), 0);

            
            // TM3 R_local = Eigen::AngleAxis<T>(rod->reference_angles[rod_cnt], axis_local).toRotationMatrix();
        
            TM3 R = Eigen::Quaternion<T>().setFromTwoVectors(axis_world, axis_local).toRotationMatrix();
            
            // R = R * R_local;

            
            for(int i = 0; i < n_div; i++)
            {
                for(int d = 0; d < 3; d++)
                {
                    V(rov + i, d) = points.col(i)[d];
                    V(rov + i+n_div, d) = points.col(i)[d];
                    if (d == 1)
                        V(rov + i+n_div, d) += axis_world.norm();
                }

                // central vertex of the top and bottom face
                V.row(rov + i) = (V.row(rov + i) * R).transpose() + vtx_from;
                V.row(rov + i + n_div) = (V.row(rov + i + n_div) * R).transpose() + vtx_from;
                
            }
            rod_cnt++;
        });
    }
    F.block(n_rod_total * rod_offset_f, 0, n_rod_total * rod_offset_f, 3) = F.block(0, 0, n_rod_total * rod_offset_f, 3);
    for (int i = n_rod_total * rod_offset_f; i < F.rows(); i++)
    {
        F.row(i) += Eigen::Vector3i(n_rod_total * rod_offset_v, n_rod_total * rod_offset_v, n_rod_total * rod_offset_v);
    }
}

template<class T, int dim>
void EoLRodSim<T, dim>::buildMeshFromRodNetwork(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
    Eigen::Ref<const DOFStack> q_display, Eigen::Ref<const IV3Stack> rods_display, 
    Eigen::Ref<const TV3Stack> normal_tile)
{
    int n_div = 10;
    
    T theta = 2.0 * EIGEN_PI / T(n_div);
    TV3Stack points = TV3Stack::Zero(3, n_div);

    
    // bottom face vertices
    for(int i = 0; i < n_div; i++)
        points.col(i) = TV3(visual_R * std::cos(theta * T(i)), 0.0, visual_R*std::sin(theta*T(i)));
    
    int n_ros_draw = rods_display.cols();
    

    // int rod_offset_v = n_div * 2 + 2;
    // int rod_offset_f = n_div * 4;

    int rod_offset_v = n_div * 2;
    int rod_offset_f = n_div * 2;
    
    V.resize(n_ros_draw * rod_offset_v, 3);
    V.setZero();
    F.resize(n_ros_draw * rod_offset_f, 3);
    F.setZero();
    int rod_cnt = 0;
    
    tbb::parallel_for(0, n_ros_draw, [&](int rod_cnt){
        int rov = rod_cnt * rod_offset_v;
        int rof = rod_cnt * rod_offset_f;

        TV vtx_from_TV = q_display.col(rods_display.col(rod_cnt)[0]).template segment<dim>(0) / unit;
        TV vtx_to_TV = q_display.col(rods_display.col(rod_cnt)[1]).template segment<dim>(0) / unit;
        // TV vtx_from_TV = q_display.col(rods_display.col(rod_cnt)[0]).template segment<dim>(0) / 1;
        // TV vtx_to_TV = q_display.col(rods_display.col(rod_cnt)[1]).template segment<dim>(0) / 1;

        // int yarn_type = rods_display.col(rod_cnt)[2];
        // T u_from = q_display(dim + yarn_type, rods_display.col(rod_cnt)[0]);
        // T u_to = q_display(dim + yarn_type, rods_display.col(rod_cnt)[1]);

        // vtx_from_TV = vtx_from_TV + u_from * (vtx_to_TV - vtx_from_TV);
        // vtx_from_TV = vtx_from_TV + u_to * (vtx_to_TV - vtx_from_TV);

        TV3 vtx_from = TV3::Zero();
        TV3 vtx_to = TV3::Zero();
        if constexpr (dim == 3)
        {
            vtx_from = vtx_from_TV;
            vtx_to = vtx_to_TV;
        }
        else
        {
            // vtx_from = TV3(vtx_from_TV[0], 0, vtx_from_TV[1]);
            // vtx_to = TV3(vtx_to_TV[0], 0, vtx_to_TV[1]);
            vtx_from = TV3(vtx_from_TV[0], vtx_from_TV[1], 0);
            vtx_to = TV3(vtx_to_TV[0], vtx_to_TV[1], 0);
        }

        
        TV3 normal_offset = TV3::Zero();
        // if (rods_display.col(rod_cnt)[2] == WARP)
        //     normal_offset = normal_tile.col(rod_cnt);
        // else
        //     normal_offset = normal_tile.col(rod_cnt);

        vtx_from += normal_offset * R;
        vtx_to += normal_offset * R;
        
        TV3 axis_world = vtx_to - vtx_from;
        TV3 axis_local(0, axis_world.norm(), 0);

        
        // V(rov + n_div*2+1, 1) = axis_world.norm();
        
        // V.row(rov + n_div*2+1) = (V.row(rov + n_div*2+1) * R).transpose() + vtx_from;
        // V.row(rov + n_div*2) = vtx_from;
        
        for(int i = 0; i < n_div; i++)
        {
            for(int d = 0; d < 3; d++)
            {
                V(rov + i, d) = points.col(i)[d];
                V(rov + i+n_div, d) = points.col(i)[d];
                if (d == 1)
                    V(rov + i+n_div, d) += axis_world.norm();
            }

            // central vertex of the top and bottom face
            V.row(rov + i) = (V.row(rov + i) * R).transpose() + vtx_from;
            V.row(rov + i + n_div) = (V.row(rov + i + n_div) * R).transpose() + vtx_from;
            
            //top faces of the cylinder
            // F.row(rof + i) = IV3(rov + n_div*2, rov + i, rov + (i+1)%(n_div));
            //bottom faces of the cylinder
            // F.row(rof + i+n_div) = IV3(rov + n_div*2+1, rov + n_div + (i+1)%(n_div), rov + i + n_div);
            
            //side faces of the cylinder
            // F.row(rof + i*2 + 2 * n_div) = IV3(rov + i, rov + i+n_div, rov + (i+1)%(n_div));
            // F.row(rof + i*2 + 1 + 2 * n_div) = IV3(rov + (i+1)%(n_div), rov + i+n_div, rov + (i+1)%(n_div) + n_div);

            F.row(rof + i*2 ) = IV3(rov + i, rov + i+n_div, rov + (i+1)%(n_div));
            F.row(rof + i*2 + 1) = IV3(rov + (i+1)%(n_div), rov + i+n_div, rov + (i+1)%(n_div) + n_div);
        }

    });
}

template<class T, int dim>
void EoLRodSim<T, dim>::getEulerianDisplacement(Eigen::MatrixXd& X, Eigen::MatrixXd& x)
{
    int n_vtx = yarns.size() * (yarns[0].size() - 1);
    int cnt = 0;
    X.resize(n_vtx, 3); x.resize(n_vtx, 3);
    // tbb::parallel_for(0, (int)yarns.size(), [&](int i){
    for (int i = 0; i < (int)yarns.size(); i++){
        TV from = q.col(yarns[i][0]).template segment<dim>(0);
        TV to = q.col(yarns[i][yarns[i].size()-2]).template segment<dim>(0);

        int yarn_type = yarns[i][yarns[i].size()-1];

        for(int j = 0; j < yarns[i].size() - 1; j++)
        {
            int node = yarns[i][j];
            T u, u0;
            u = yarn_type == WARP ? q(dim + 0, node) : q(dim + 1, node);
            u0 = yarn_type == WARP ? q0(dim + 0, node) : q0(dim + 1, node);
            
            TV x_lag = from + u * (to - from);
            TV x_eul = from + u0 * (to - from);
            if constexpr (dim == 2)
            {
                X.row(cnt) = Eigen::RowVector3d(x_eul[0], x_eul[1], 0);
                x.row(cnt) = Eigen::RowVector3d(x_lag[0] + 1.0, x_lag[1], 0);
                cnt++;
            }
        }
    }
    // });
    assert(cnt == n_vtx);
}

template<class T, int dim>
void EoLRodSim<T, dim>::markSlidingRange(int idx, int dir, int depth, 
    std::vector<bool>& can_slide, int root)
{
    if (depth > slide_over_n_rods[dir] || idx == -1)
        return;

    // T rod_length = (q0.col(rods.col(0)(0)).template segment<dim>(0) - 
    //     q0.col(rods.col(0)(1)).template segment<dim>(0)).norm();
    
    // distance root/crossing node travels

    // std::vector<TV> X, X0, X_root, X0_root, dXdu, d2Xdu2;
    // getMaterialPositions(q0, {idx}, X, dir, dXdu, d2Xdu2, false, false );
    // getMaterialPositions(q, {idx}, X0, dir, dXdu, d2Xdu2, false, false );
    // getMaterialPositions(q0, {root}, X0_root, dir, dXdu, d2Xdu2, false, false );
    // getMaterialPositions(q, {root}, X_root, dir, dXdu, d2Xdu2, false, false );

    T root_sliding_dis = std::abs(q(dim + dir, root) - q0(dim + dir, root));
    
    T dis_to_root_rest_state = std::abs(q0(dim + dir, idx) - q0(dim + dir, root));

    T dis_to_root_current = std::abs(q(dim + dir, idx) - q(dim + dir, root));

    // T root_sliding_dis = (X_root[0] - X0_root[0]).norm();
    
    // T dis_to_root_rest_state = (X0[0] - X0_root[0]).norm();

    // T dis_to_root_current = (X[0] - X_root[0]).norm();

    T sliding_range = dir == 0 ? tunnel_u : tunnel_v;

    if (idx == root || root_sliding_dis < 1e-6)
        can_slide[idx * 2 + dir] = true;
    
    // if (root_sliding_dis > slide_over_n_rods[dir] * rod_length - 1e-6)
    if (root_sliding_dis > sliding_range - 1e-6)
    {
        if(dis_to_root_current > dis_to_root_rest_state)
            can_slide[idx * 2 + dir] = true;
    }
    else
    {
        can_slide[idx * 2 + dir] = true;
    }
    
    
    if(dir == 0)
    {
        markSlidingRange(connections(2, idx), dir, depth + 1, can_slide, root);
        markSlidingRange(connections(3, idx), dir, depth + 1, can_slide, root);
    }
    else
    {
        markSlidingRange(connections(0, idx), dir, depth + 1, can_slide, root);
        markSlidingRange(connections(1, idx), dir, depth + 1, can_slide, root);
    }   
}

template<class T, int dim>
void EoLRodSim<T, dim>::getColorPerYarn(Eigen::MatrixXd& C, int n_rod_per_yarn)
{
    int n_faces = 20;

    std::vector<bool> can_slide(n_nodes * 2, false);

    iterateSlidingNodes([&](int node_idx){
              
        if (dirichlet_data.find(node_idx) != dirichlet_data.end())
        {
            TVDOF mask = dirichlet_data[node_idx].first;
            if(!mask[dim])
                markSlidingRange(node_idx, 0, 0, can_slide, node_idx);
            if(!mask[dim+1])
                markSlidingRange(node_idx, 1, 0, can_slide, node_idx);
        }
        else
        {
            markSlidingRange(node_idx, 0, 0, can_slide, node_idx);
            markSlidingRange(node_idx, 1, 0, can_slide, node_idx);
            
        }
    });

    C.resize(n_rods * n_faces, 3);
    std::vector<Eigen::Vector3d> colors = {
        Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(0, 1, 0), 
        Eigen::Vector3d(0, 1, 1), Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(1, 0, 1)};
    int n_yarn = yarns.size();
    VectorXT delta_u(n_rods);
    delta_u.setZero();

    tbb::parallel_for(0, n_rods, [&](int rod_idx){
        int v0 = rods(0, rod_idx);
        int v1 = rods(1, rod_idx);
        T dU = (q0.col(v1).template segment<2>(dim) - q0.col(v0).template segment<2>(dim)).norm();
        T du = (q.col(v1).template segment<2>(dim) - q.col(v0).template segment<2>(dim)).norm();
        delta_u[rod_idx] = std::abs(du - dU);
    });
    
    // Eigen::MatrixXd sliding_color(n_rods, 3);
    // sliding_color.setZero();
    // std::cout << "color mapping" << std::endl;
    // igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_JET, delta_u, true, sliding_color);
    // std::cout << "color mapping done" << std::endl;
    tbb::parallel_for(0, n_rods, [&](int rod_idx){
        for(int i = 0; i < n_faces; i++)
        {

            if(can_slide[rods(0,rod_idx) * 2 + rods(2, rod_idx)] && can_slide[rods(1,rod_idx) * 2 + rods(2, rod_idx)])
                C.row(rod_idx * n_faces + i) = colors[0];
            else
                C.row(rod_idx * n_faces + i) = colors[2];

        }
            // C.row(rod_idx * 40 + i) = sliding_color.row(rod_idx);
            // C.row(rod_idx * 40 + i) = rod_idx % 2 == 0 ? Eigen::Vector3d::Ones() : Eigen::Vector3d::Zero();
            // C.row(rod_idx * 40 + i) = colors[yarn_map[rod_idx]];
    });
}

template<class T, int dim>
void EoLRodSim<T, dim>::showStretching(Eigen::MatrixXd& C)
{
    int n_faces = 20;
    int n_rod_total = 0;
    for(auto& rod : Rods)
    {
        n_rod_total += rod->numSeg();
    }
    C.resize(n_rod_total * n_faces, 3);

    VectorXT rod_energy(n_rod_total);
    rod_energy.setZero();

    int rod_cnt = 0;
    for (auto& rod : Rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx){
            // std::cout << "i " << node_i << " j " << node_j << std::endl;
            TV xi, xj, Xi, Xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;
            T V[1];
            if constexpr (dim == 2)
            {
                // #include "Maple/YarnStretchingV.mcg"
            }
            else if constexpr (dim == 3)
            {
                // #include "Maple/RodStretching3DV.mcg"
                V[0] = addStretchingEnergy();   
            }
            rod_energy[rod_cnt++] += V[0];
            // std::cout << "Rod " << node_i << "->" << node_j << ": " << V[0] << std::endl;
            std::cout << "Rod " << node_i << "->" << node_j << ": " << (xi - xj).norm() / (Xi - Xj).norm() << std::endl;
        });
    }

    if(rod_energy.maxCoeff() > 1e-4)
        rod_energy /= rod_energy.maxCoeff();
    else
        rod_energy.setZero();
    
    tbb::parallel_for(0, n_rod_total, [&](int rod_idx){
        for(int i = 0; i < n_faces; i++)
            C.row(rod_idx * n_faces + i) = Eigen::Vector3d(rod_energy[rod_idx], rod_energy[rod_idx], rod_energy[rod_idx]);
    });

}

template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;