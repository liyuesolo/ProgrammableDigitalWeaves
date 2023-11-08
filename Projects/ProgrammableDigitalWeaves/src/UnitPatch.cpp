#include <iostream>
#include <utility>
#include <fstream>
#include <unordered_map>
#include <chrono> 
#include "../include/UnitPatch.h"
#include "../include/HybridC2Curve.h"


#include "../include/IO.h"
#include "../include/GCodeGenerator.h"


static double ROD_A = 2.5e-4;
static double ROD_B = 2.5e-4;


#include <random>
#include <cmath>
std::random_device rd;
std::mt19937 gen( rd() );
std::uniform_real_distribution<> dis( 0.0, 1.0 );

static double zeta()
{
	return dis(gen);
}

template<class T, int dim>
void UnitPatch<T, dim>::addPoint(const TV& point, int& full_dof_cnt, int& node_cnt)
{
    deformed_states.conservativeResize(full_dof_cnt + dim);
    deformed_states.template segment<dim>(full_dof_cnt) = point;
    full_dof_cnt += dim;
    node_cnt++;
}

template<class T, int dim>
void UnitPatch<T, dim>::addCrossingPoint(std::vector<TV>& existing_nodes,
     const TV& point, int& full_dof_cnt, int& node_cnt)
{
    sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>())); 
    deformed_states.conservativeResize(full_dof_cnt + dim);
    deformed_states.template segment<dim>(full_dof_cnt) = point;
    existing_nodes.push_back(point);
    full_dof_cnt += dim;
    node_cnt++;
}

template<class T, int dim>
void UnitPatch<T, dim>::buildScene(int patch_type)
{
    if (patch_type == 0)
        build3DtestScene(4);
    else if (patch_type == 1)
        buildOneCrossScene(16);
    else if (patch_type == 2)
        buildGridScene(64);
    else if (patch_type == 3)
        buildOmegaScene(2);
    else if (patch_type == 4)
        buildStraightRodScene(16);
    else if (patch_type == 7)
        buildFingerScene(8);
    else if (patch_type == 8)
        buildShelterScene(32);
    else if (patch_type == 9)
        buildGripperScene(32);
    else if (patch_type == 10)
        buildGridLayoutGripper(32);
    else if (patch_type == 11)
        buildGridScene2(64);
    else if (patch_type == 12)
        buildSaddleScene(64);
    else if (patch_type == 13)
        buildTestJoint(32);
    else if (patch_type == 14)
        buildXJointsScene(32);
    else if (patch_type == 15)
        buildSquareCrossJointScene(32);
    else if (patch_type == 16)
        buildActiveTextileScene(32);
    else if (patch_type == 17)
        buildXJointsScene2(32);
    else if (patch_type == 18)
        buildInterlockingSquareScene(32);
    else if (patch_type == 19)
        buildDenseInterlockingSquareScene(8); // use 8 for printing
    else if (patch_type == 20)
        buildTestSceneJuan(16);
    else if (patch_type == 21)
        buildDenseInterlockingSquarePeriodicScene(32);
    else if (patch_type == 22)
        buildDomeScene(32);
    else if (patch_type == 23)
        buildShoeScene(8);
    else if (patch_type == 24)
        buildRandomPatchScene(4);
    else if (patch_type == 25)
        buildPeriodicCircleScene(32);
    else if (patch_type == 26)
        buildFullCircleScene(8);
    else if (patch_type == 27)
        buildShelterAcutationScene(32);
    else if (patch_type == 28)
        buildTestSceneJuan(32);
    else if (patch_type == 29)
        buildTennisBallWrapperScene(32);
    else if (patch_type == 30)
        buildActuationSingleStrandScene(64);
    else if (patch_type == 31)
        buildActuationSingleStrandSceneWithoutCrossing(32);
    else if (patch_type == 32)
        buildGridClosed(32);
}   

template<class T, int dim>
void UnitPatch<T, dim>::buildGridClosed(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        
        sim.visual_R = 0.005;

        sim.unit = 0.09;
        
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        int n_row = 20, n_col = 20;

        // push crossings first 
        T dy = 1.0 / n_row * sim.unit;
        T dx = 1.0 / n_col * sim.unit;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        

        std::vector<TV> nodal_positions;

        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                TV pt = TV(dx * col, dy * row, 0);
                addCrossingPoint(nodal_positions, pt, full_dof_cnt, node_cnt);
            }
        }

        for (int row = 0; row < n_row; row++)
        {
            std::vector<TV> passing_points;
            std::vector<int> passing_points_id;
            for (int col = 0; col < n_col; col++)
            {
                passing_points_id.push_back(row * n_col + col);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
        }

        for (int col = 0; col < n_col; col++)
        {
            std::vector<TV> passing_points;
            std::vector<int> passing_points_id;
            for (int row = 0; row < n_row; row++)
            {
                passing_points_id.push_back(row * n_col + col);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
        }

        for (auto rod : sim.Rods)
        {
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);
        }
        
        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        for (int row = 1; row < n_row - 1; row++)
        {
            for (int col = 1; col < n_col - 1; col++)
            {
                auto crossing = sim.rod_crossings[row * n_col + col];
                crossing->is_fixed = false;
                crossing->sliding_ranges[0] = Range::Ones();
            }
        }
        

        T r = 0.1 * sim.unit;
        TV center1, center2;
        sim.getCrossingPosition(0, center1);
        sim.getCrossingPosition(n_row * n_col - 1, center2);

        TV delta1 = TV(-0.3, -0.3, -1e-2) * sim.unit;

        auto circle1 = [r, center1, delta1](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = delta1;
            return (x - center1).norm() < r;
        };

        TV delta2 = TV(0.0, 0.0, 0) * sim.unit;
        auto circle2 = [r, center2, delta2](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = delta2;
            return (x - center2).norm() < r;

        };

        sim.fixRegionalDisplacement(circle1);
        sim.fixRegionalDisplacement(circle2);
        
        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }

        GCodeGenerator<T, dim>(sim, "closed_grid_sliding.gcode").generateGCodeClosedGrid(n_row, n_col, 0);
    }


}

template<class T, int dim>
void UnitPatch<T, dim>::buildActuationSingleStrandSceneWithoutCrossing(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.02;
        //sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

    

        // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
        std::vector<T> thetas;
        int theta_sub = 8;
        for (int i = 0; i< theta_sub; i++)
            thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

        std::cout << "thetas ";
        for(int i=0; i<thetas.size(); ++i)
            std::cout << " " << thetas[i];
        std::cout << std::endl;

        std::vector<bool> is_inner_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_center_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_outer_crossing = {true, true, true, true, true, true, true, true};


        std::vector<int> id_inner_crossing;
        std::vector<int> id_center_crossing;
        std::vector<int> id_outer_crossing;

        int cur_crossing_idx = 0;
        for(int i=0; i<is_inner_crossing.size(); ++i)
        {
            if(is_inner_crossing[i])
                id_inner_crossing.push_back(cur_crossing_idx++);
            else
                id_inner_crossing.push_back(-1);
        }
        for(int i=0; i<is_center_crossing.size(); ++i)
        {
            if(is_center_crossing[i])
                id_center_crossing.push_back(cur_crossing_idx++);
            else
                id_center_crossing.push_back(-1);
        }
        for(int i=0; i<is_outer_crossing.size(); ++i)
        {
            if(is_outer_crossing[i])
                id_outer_crossing.push_back(cur_crossing_idx++);
            else
                id_outer_crossing.push_back(-1);
        }

        //cur_crossing_idx++;

        std::vector<int> id_line_crossing_principal = {20, 12, 4, 0, 8, 16};
        std::vector<int> id_line_crossing_principal_horizontal = {22, 14, 6, 2, 10, 18};
        std::vector<std::vector<int>> id_line_crossings = {{1, 9, 17},{3, 11, 19}, {5, 13, 21}, {7, 15, 23}};

        std::vector<int> crossings_ids;
        for(int i=0; i<cur_crossing_idx; ++i)
            crossings_ids.push_back(i);

        std::vector<TV> cross_point_collection;

        for(int i=0; i<cur_crossing_idx; ++i)
                sim.rod_crossings.push_back(new RodCrossing<T, dim>(crossings_ids[i], std::vector<int>())); 

        int cur_idx = 0;

        std::vector<TV> points_inner;
        std::vector<TV> crossing_points_inner;
        std::vector<int> indices_inner;
        for(int i=0; i<thetas.size(); ++i)
            points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_inner.size(); ++i)
        {
            if(id_inner_crossing[i]!=-1)
            {
                addPoint(points_inner[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_inner[i]);
                crossing_points_inner.push_back(points_inner[i]);
                indices_inner.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_center;
        std::vector<TV> crossing_points_center;
        std::vector<int> indices_center;
        for(int i=0; i<thetas.size(); ++i)
            points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_center.size(); ++i)
        {
            if(id_center_crossing[i]!=-1)
            {
                addPoint(points_center[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_center[i]);
                crossing_points_center.push_back(points_center[i]);
                indices_center.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_outer;
        std::vector<TV> crossing_points_outer;
        std::vector<int> indices_outer;
        for(int i=0; i<thetas.size(); ++i)
            points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_outer.size(); ++i)
        {
            if(id_outer_crossing[i]!=-1)
            {
                addPoint(points_outer[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_outer[i]);
                crossing_points_outer.push_back(points_outer[i]);
                indices_outer.push_back(cur_idx++);
            }
        }
             
        //addPoint(TV(0,0,0), full_dof_cnt, node_cnt);

        std::vector<TV2> data_points;
        for(int i=0; i<points_inner.size(); ++i)
            data_points.push_back(points_inner[i].template head<2>());
        data_points.push_back(points_inner[0].template head<2>());

        addCurvedRod(data_points, crossing_points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
        int ii = 0;
        for (int i = 0; i < points_inner.size(); i++)
        {
            if(is_inner_crossing[i])
                addCrossingData(id_inner_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        data_points.clear();
        for(int i=0; i<points_center.size(); ++i)
            data_points.push_back(points_center[i].template head<2>());
        data_points.push_back(points_center[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_center.size(); i++)
        {
            if(is_center_crossing[i])
                addCrossingData(id_center_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        
        data_points.clear();
        for(int i=0; i<points_outer.size(); ++i)
            data_points.push_back(points_outer[i].template head<2>());
        data_points.push_back(points_outer[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_outer.size(); i++)
        {
            if(is_outer_crossing[i])
                addCrossingData(id_outer_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }
        
        TV to = TV(0, -2.5, 0.0) * sim.unit;
        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;

        passing_points_id = id_line_crossing_principal;
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        to = TV(2.0, 0.0, 0.0) * sim.unit;
        passing_points_id = id_line_crossing_principal_horizontal;
        passing_points.clear();
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        int idx =0;
        for (int i = 1; i < thetas.size(); i += 2)
        {
            passing_points = {points_inner[i], points_center[i], points_outer[i]};
            passing_points_id = {id_line_crossings[idx][0], id_line_crossings[idx][1], id_line_crossings[idx][2]};
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

            for (int j = 0; j < id_line_crossings[idx].size(); j++)
                addCrossingData(id_line_crossings[idx][j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

            ++idx;
        }

        for(int i=1; i<id_line_crossing_principal.size(); ++i)
        {
            sim.rod_crossings[id_line_crossing_principal[i]]->is_fixed = false;
            sim.rod_crossings[id_line_crossing_principal[i]]->sliding_ranges[1] = Range::Ones();
        }

        for(int i=1; i<id_line_crossing_principal_horizontal.size(); ++i)
        {
            sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->is_fixed = false;
            sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->sliding_ranges[1] = Range::Ones();
        }


        //sim.rod_crossings[24]->is_fixed = false;
        //sim.rod_crossings[24]->sliding_ranges[0] = Range::Ones();

        // for(int i=0; i<sim.rod_crossings.size(); ++i)
        //     std::cout << sim.rod_crossings[i]->sliding_ranges.size() << " shjdkasdk " << std::endl;

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        // Offset offset;
        // sim.Rods[3]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = 0.3 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        auto incrementalBC = [](EoLRodSim<T, dim>&_sim, int step)->void
        {
            Offset offset;
            _sim.Rods[3]->backOffsetReduced(offset);
            _sim.dirichlet_dof[offset[0]] = 0;
            _sim.dirichlet_dof[offset[1]] = 0.1 * _sim.unit;
            if (step == 0)
                _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
            else
                _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
        };

        sim.incremental_steps = 5;
        sim.incremental_bc = incrementalBC;

        // sim.Rods[4]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = - 0.25 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        // Offset offset;
        // sim.Rods[3]->getEntryReduced(sim.Rods[3]->indices[sim.Rods[3]->dof_node_location.back()], offset);
        // sim.dirichlet_dof[offset[3]] = 0.01 * sim.unit;

        sim.Rods[2]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset] = 0;
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-1] = 0;

        T r = 0.1 * sim.unit;
        TV center = TV(1.4*cos(thetas[0]), 1.4*sin(thetas[0]), 0) * sim.unit;
        TV center2 = TV(1.4*cos(thetas[2]), 1.4*sin(thetas[2]), 0) * sim.unit;

        auto circle = [r, center](const TV& x)->bool
        {
            return (x - center).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        sim.fixRegion(circle);
        //sim.fixRegion(circle2);

        sim.fixCrossing();

        sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));
        //sim.boundary_spheres.push_back(std::make_pair(center2, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildActuationSingleStrandScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = false;

        sim.ke = 1e-4;

        sim.unit = 0.02;
        //sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

    

        // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
        std::vector<T> thetas;
        int theta_sub = 8;
        for (int i = 0; i< theta_sub; i++)
            thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

        std::cout << "thetas ";
        for(int i=0; i<thetas.size(); ++i)
            std::cout << " " << thetas[i];
        std::cout << std::endl;

        std::vector<bool> is_inner_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_center_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_outer_crossing = {true, true, true, true, true, true, true, true};


        std::vector<int> id_inner_crossing;
        std::vector<int> id_center_crossing;
        std::vector<int> id_outer_crossing;

        int cur_crossing_idx = 0;
        for(int i=0; i<is_inner_crossing.size(); ++i)
        {
            if(is_inner_crossing[i])
                id_inner_crossing.push_back(cur_crossing_idx++);
            else
                id_inner_crossing.push_back(-1);
        }
        for(int i=0; i<is_center_crossing.size(); ++i)
        {
            if(is_center_crossing[i])
                id_center_crossing.push_back(cur_crossing_idx++);
            else
                id_center_crossing.push_back(-1);
        }
        for(int i=0; i<is_outer_crossing.size(); ++i)
        {
            if(is_outer_crossing[i])
                id_outer_crossing.push_back(cur_crossing_idx++);
            else
                id_outer_crossing.push_back(-1);
        }

        cur_crossing_idx++;

        std::vector<int> id_line_crossing_principal = {20, 12, 4, 24, 0, 8, 16};
        std::vector<int> id_line_crossing_principal_horizontal = {22, 14, 6, 24, 2, 10, 18};
        std::vector<std::vector<int>> id_line_crossings = {{1, 9, 17},{3, 11, 19}, {5, 13, 21}, {7, 15, 23}};

        std::vector<int> crossings_ids;
        for(int i=0; i<cur_crossing_idx; ++i)
            crossings_ids.push_back(i);

        std::vector<TV> cross_point_collection;

        for(int i=0; i<cur_crossing_idx; ++i)
                sim.rod_crossings.push_back(new RodCrossing<T, dim>(crossings_ids[i], std::vector<int>())); 

        int cur_idx = 0;

        std::vector<TV> points_inner;
        std::vector<TV> crossing_points_inner;
        std::vector<int> indices_inner;
        for(int i=0; i<thetas.size(); ++i)
            points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_inner.size(); ++i)
        {
            if(id_inner_crossing[i]!=-1)
            {
                addPoint(points_inner[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_inner[i]);
                crossing_points_inner.push_back(points_inner[i]);
                indices_inner.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_center;
        std::vector<TV> crossing_points_center;
        std::vector<int> indices_center;
        for(int i=0; i<thetas.size(); ++i)
            points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_center.size(); ++i)
        {
            if(id_center_crossing[i]!=-1)
            {
                addPoint(points_center[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_center[i]);
                crossing_points_center.push_back(points_center[i]);
                indices_center.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_outer;
        std::vector<TV> crossing_points_outer;
        std::vector<int> indices_outer;
        for(int i=0; i<thetas.size(); ++i)
            points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_outer.size(); ++i)
        {
            if(id_outer_crossing[i]!=-1)
            {
                addPoint(points_outer[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_outer[i]);
                crossing_points_outer.push_back(points_outer[i]);
                indices_outer.push_back(cur_idx++);
            }
        }
             
        addPoint(TV(0,0,0), full_dof_cnt, node_cnt);

        std::vector<TV2> data_points;
        for(int i=0; i<points_inner.size(); ++i)
            data_points.push_back(points_inner[i].template head<2>());
        data_points.push_back(points_inner[0].template head<2>());

        addCurvedRod(data_points, crossing_points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
        int ii = 0;
        for (int i = 0; i < points_inner.size(); i++)
        {
            if(is_inner_crossing[i])
                addCrossingData(id_inner_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        data_points.clear();
        for(int i=0; i<points_center.size(); ++i)
            data_points.push_back(points_center[i].template head<2>());
        data_points.push_back(points_center[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_center.size(); i++)
        {
            if(is_center_crossing[i])
                addCrossingData(id_center_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        
        data_points.clear();
        for(int i=0; i<points_outer.size(); ++i)
            data_points.push_back(points_outer[i].template head<2>());
        data_points.push_back(points_outer[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_outer.size(); i++)
        {
            if(is_outer_crossing[i])
                addCrossingData(id_outer_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }
        
        TV to = TV(0, -2.2, 0.0) * sim.unit;
        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;

        passing_points_id = id_line_crossing_principal;
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        to = TV(2.0, 0.0, 0.0) * sim.unit;
        passing_points_id = id_line_crossing_principal_horizontal;
        passing_points.clear();
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        int idx =0;
        for (int i = 1; i < thetas.size(); i += 2)
        {
            passing_points = {points_inner[i], points_center[i], points_outer[i]};
            passing_points_id = {id_line_crossings[idx][0], id_line_crossings[idx][1], id_line_crossings[idx][2]};
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

            for (int j = 0; j < id_line_crossings[idx].size(); j++)
                addCrossingData(id_line_crossings[idx][j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

            ++idx;
        }

        for (auto rod : sim.Rods)
        {
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);
        }
        
        for(int i=1; i<id_line_crossing_principal.size(); ++i)
        {
            sim.rod_crossings[id_line_crossing_principal[i]]->is_fixed = false;
            sim.rod_crossings[id_line_crossing_principal[i]]->sliding_ranges[1] = Range::Ones();

            auto crossing = sim.rod_crossings[id_line_crossing_principal[i]];
            auto rod = sim.Rods[crossing->rods_involved[0]];
            int cnt = 0;
            for (int idx : rod->dof_node_location)
            {
                if(rod->indices[idx] == crossing->node_idx)
                    rod->fixed_by_crossing[cnt] = false;     
                cnt ++;
            }
        }

        

        // for(int i=1; i<id_line_crossing_principal_horizontal.size(); ++i)
        // {
        //     sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->is_fixed = false;
        //     sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->sliding_ranges[1] = Range::Ones();
        // }


        sim.rod_crossings[24]->is_fixed = false;
        sim.rod_crossings[24]->sliding_ranges[0] = Range::Ones();

        // for(int i=0; i<sim.rod_crossings.size(); ++i)
        //     std::cout << sim.rod_crossings[i]->sliding_ranges.size() << " shjdkasdk " << std::endl;

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        // Offset offset;
        // sim.Rods[3]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = 0.3 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        auto incrementalBC = [](EoLRodSim<T, dim>&_sim, int step)->void
        {
            Offset offset;
            _sim.Rods[3]->backOffsetReduced(offset);
            _sim.dirichlet_dof[offset[0]] = 0;
            _sim.dirichlet_dof[offset[1]] = 0.025 * _sim.unit;
            if (step == 0)
                _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
            else
                _sim.dirichlet_dof[offset[2]] = 0.0 * _sim.unit;
        };

        sim.incremental_steps = 25;
        sim.incremental_bc = incrementalBC;

        // sim.Rods[4]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = - 0.25 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        // Offset offset;
        // sim.Rods[3]->getEntryReduced(sim.Rods[3]->indices[sim.Rods[3]->dof_node_location.back()], offset);
        // sim.dirichlet_dof[offset[3]] = 0.01 * sim.unit;
        for (int n = 0; n < 5; n++)
        {
            sim.Rods[2]->fixPointLagrangian(0 + n, TV::Zero(), sim.dirichlet_dof);
            sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 1 - n, TV::Zero(), sim.dirichlet_dof);    
            sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + n] = 0;
            sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg() - 1 - n] = 0;
        }
        // sim.Rods[2]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-1] = 0;

        T r = 0.1 * sim.unit;
        TV center = TV(1.4*cos(thetas[0]), 1.4*sin(thetas[0]), 0) * sim.unit;
        TV center2 = TV(1.4*cos(thetas[2]), 1.4*sin(thetas[2]), 0) * sim.unit;

        auto circle = [r, center](const TV& x)->bool
        {
            return (x - center).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        // sim.fixRegion(circle);
        //sim.fixRegion(circle2);

        sim.fixCrossing();

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));
        //sim.boundary_spheres.push_back(std::make_pair(center2, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r ) * sim.unit;
            
        }
        // GCodeGenerator<T, dim>(sim, "single_strand.gcode").generateGCodeSingleStrand();
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildTennisBallWrapperScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.04;
        //sim.unit = 1.0;
        sim.visual_R = 0.01;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        auto horizontal = [](double r, double L, double t){
            TV cg;
            cg[0] = r * cos(t);
            cg[1] = r * sin(t) - sqrt(-L * L + r * r);
            cg[2] = 0;
            return cg;
        };

        auto vertical = [](double r, double L, double t){
            TV cg;
            cg[0] = r * cos(t) - sqrt(-L * L + r * r);
            cg[1] = r * sin(t);
            cg[2] = 0;
            return cg;
        };

        auto crossings_f = []()
        {
            std::vector<std::vector<double>> cg1(3, std::vector<double>(3));

            cg1[0][0] = 0.1254203428e1;
            cg1[0][1] = 0.1322764991e1;
            cg1[0][2] = 0.1379390878e1;
            cg1[1][0] = 1.280360672;
            cg1[1][1] = 0.1340519331e1;
            cg1[1][2] = 0.1391579417e1;
            cg1[2][0] = 1.320641700;
            cg1[2][1] = 1.370605522;
            cg1[2][2] = 0.1413913841e1;


            return cg1;
        };

        std::vector<TV> points;
        points.push_back(TV(-1,0,0));
        points.push_back(TV(1,0,0));
        points.push_back(TV(0,1,0));
        points.push_back(TV(0,-1,0));

        std::vector<double> radius = {1.3, 1.5, 1.8};
        std::vector<std::vector<double>> crossings_t = crossings_f();

        for(int i=0; i<3; ++i)
        {
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][0]));
            points.back()[0] *= -1.0;
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][1]));
            points.back()[0] *= -1.0;
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][2]));
            points.back()[0] *= -1.0;
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][2]));
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][1]));
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][0]));
        }

        for(int i=2; i>=0; --i)
        {
            points.push_back(-horizontal(radius[i], 1.0, crossings_t[i][0]));
            points.push_back(-horizontal(radius[i], 1.0, crossings_t[i][1]));
            points.push_back(-horizontal(radius[i], 1.0, crossings_t[i][2]));
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][2]));
            points.back()[1] *= -1.0;
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][1]));
            points.back()[1] *= -1.0;
            points.push_back(horizontal(radius[i], 1.0, crossings_t[i][0]));
            points.back()[1] *= -1.0;
        }

        for(int i=0; i<points.size(); ++i)
            std::cout << points[i].transpose() << std::endl;

        for(int i=0; i<points.size(); ++i)
        {
            points[i] *= sim.unit;
            addPoint(points[i], full_dof_cnt, node_cnt);
        }

        for(int i=0; i<points.size(); ++i)
             sim.rod_crossings.push_back(new RodCrossing<T, dim>(i, std::vector<int>())); 

        for(int i=0; i<6; ++i)
        {
            std::vector<int> passing_idx = {0};
            for(int j=0; j<6; ++j)
                passing_idx.push_back(4 + i*6 + j);
            passing_idx.push_back(1);

            std::vector<TV2> data_points;
            std::vector<TV> crossing_points;
            for(int j=0; j<passing_idx.size(); ++j)
            {
                data_points.push_back(points[passing_idx[j]].template head<2>());
                crossing_points.push_back(points[passing_idx[j]]);
            }

            addCurvedRod(data_points, crossing_points, passing_idx, sub_div, full_dof_cnt, node_cnt, rod_cnt, false);

            for (int j = 0; j < passing_idx.size(); j++)
                addCrossingData(passing_idx[j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);
        }

        for(int i=0; i<6; ++i)
        {
            std::vector<int> passing_idx = {2};
            for(int j=0; j<6; ++j)
                passing_idx.push_back(4 + j*6 + i);
            passing_idx.push_back(3);

            std::vector<TV2> data_points;
            std::vector<TV> crossing_points;
            for(int j=0; j<passing_idx.size(); ++j)
            {
                data_points.push_back(points[passing_idx[j]].template head<2>());
                crossing_points.push_back(points[passing_idx[j]]);
            }

            addCurvedRod(data_points, crossing_points, passing_idx, sub_div, full_dof_cnt, node_cnt, rod_cnt, false);

            for (int j = 0; j < passing_idx.size(); j++)
                addCrossingData(passing_idx[j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);
        }

        for (auto rod : sim.Rods)
        {
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);
        }


        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
        for(int i=4; i<points.size(); ++i)
        {
            sim.rod_crossings[i]->is_fixed = false;
            sim.rod_crossings[i]->sliding_ranges[0] = Range::Ones();
            auto crossing = sim.rod_crossings[i];
            auto rod = sim.Rods[crossing->rods_involved[1]];
            int cnt = 0;
            for (int idx : rod->dof_node_location)
            {
                if(rod->indices[idx] == crossing->node_idx)
                    rod->fixed_by_crossing[cnt] = false;     
                cnt ++;
            }
        }

    
        Offset offset;
        sim.Rods[3]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] =  0.2 * sim.unit;
        sim.dirichlet_dof[offset[1]] = 0;
        sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        // sim.Rods[3]->getEntryReduced(0, offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = -0.1 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        // sim.Rods[2]->getEntryReduced(6, offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = 0.1 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        //sim.Rods[3]->getEntryReduced(6, offset);
        sim.Rods[3]->backOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = -0.2 * sim.unit;
        sim.dirichlet_dof[offset[1]] = 0;
        sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        sim.Rods[6]->backOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = 0.2 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        sim.Rods[6]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        // sim.Rods[4]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0.5 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        // // Offset offset;
        // // sim.Rods[3]->getEntryReduced(sim.Rods[3]->indices[sim.Rods[3]->dof_node_location.back()], offset);
        // // sim.dirichlet_dof[offset[3]] = 0.01 * sim.unit;

        // sim.Rods[2]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-1] = 0;

        T r = 0.1 * sim.unit;
        TV center = TV(-1.0, 0.0, 0) * sim.unit;
        TV center2 = TV(1.0, 0.0, 0) * sim.unit;

        // auto circle = [r, center](const TV& x)->bool
        // {
        //     return (x - center).norm() < r;
        // };

        // auto circle2 = [r, center2](const TV& x)->bool
        // {
        //     return (x - center2).norm() < r;
        // };

        // sim.fixRegion(circle);
        // sim.fixRegion(circle2);

        sim.fixCrossing();

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));
        // sim.boundary_spheres.push_back(std::make_pair(center2, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
        GCodeGenerator<T, dim>(sim, "tennis.gcode").generateGCodeShelter();
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildShelterAcutationScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = false;

        sim.ke = 1e-4;

        sim.unit = 0.02;
        //sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

    

        // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
        std::vector<T> thetas;
        int theta_sub = 8;
        for (int i = 0; i< theta_sub; i++)
            thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

        std::cout << "thetas ";
        for(int i=0; i<thetas.size(); ++i)
            std::cout << " " << thetas[i];
        std::cout << std::endl;

        std::vector<bool> is_inner_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_center_crossing = {true, true, true, true, true, true, true, true};
        std::vector<bool> is_outer_crossing = {true, true, true, true, true, true, true, true};


        std::vector<int> id_inner_crossing;
        std::vector<int> id_center_crossing;
        std::vector<int> id_outer_crossing;

        int cur_crossing_idx = 0;
        for(int i=0; i<is_inner_crossing.size(); ++i)
        {
            if(is_inner_crossing[i])
                id_inner_crossing.push_back(cur_crossing_idx++);
            else
                id_inner_crossing.push_back(-1);
        }
        for(int i=0; i<is_center_crossing.size(); ++i)
        {
            if(is_center_crossing[i])
                id_center_crossing.push_back(cur_crossing_idx++);
            else
                id_center_crossing.push_back(-1);
        }
        for(int i=0; i<is_outer_crossing.size(); ++i)
        {
            if(is_outer_crossing[i])
                id_outer_crossing.push_back(cur_crossing_idx++);
            else
                id_outer_crossing.push_back(-1);
        }

        std::vector<int> id_line_crossing_principal = {20, 12, 4, 0, 8, 16};
        std::vector<int> id_line_crossing_principal_horizontal = {22, 14, 6, 2, 10, 18};
        std::vector<std::vector<int>> id_line_crossings = {{1, 9, 17},{3, 11, 19}, {5, 13, 21}, {7, 15, 23}};

        std::vector<int> crossings_ids;
        for(int i=0; i<cur_crossing_idx; ++i)
            crossings_ids.push_back(i);

        std::vector<TV> cross_point_collection;

        for(int i=0; i<cur_crossing_idx; ++i)
                sim.rod_crossings.push_back(new RodCrossing<T, dim>(crossings_ids[i], std::vector<int>())); 

        int cur_idx = 0;

        std::vector<TV> points_inner;
        std::vector<TV> crossing_points_inner;
        std::vector<int> indices_inner;
        for(int i=0; i<thetas.size(); ++i)
            points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_inner.size(); ++i)
        {
            if(id_inner_crossing[i]!=-1)
            {
                addPoint(points_inner[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_inner[i]);
                crossing_points_inner.push_back(points_inner[i]);
                indices_inner.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_center;
        std::vector<TV> crossing_points_center;
        std::vector<int> indices_center;
        for(int i=0; i<thetas.size(); ++i)
            points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_center.size(); ++i)
        {
            if(id_center_crossing[i]!=-1)
            {
                addPoint(points_center[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_center[i]);
                crossing_points_center.push_back(points_center[i]);
                indices_center.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_outer;
        std::vector<TV> crossing_points_outer;
        std::vector<int> indices_outer;
        for(int i=0; i<thetas.size(); ++i)
            points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_outer.size(); ++i)
        {
            if(id_outer_crossing[i]!=-1)
            {
                addPoint(points_outer[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_outer[i]);
                crossing_points_outer.push_back(points_outer[i]);
                indices_outer.push_back(cur_idx++);
            }
        }
                
        std::vector<TV2> data_points;
        for(int i=0; i<points_inner.size(); ++i)
            data_points.push_back(points_inner[i].template head<2>());
        data_points.push_back(points_inner[0].template head<2>());

        addCurvedRod(data_points, crossing_points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
        int ii = 0;
        for (int i = 0; i < points_inner.size(); i++)
        {
            if(is_inner_crossing[i])
                addCrossingData(id_inner_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        data_points.clear();
        for(int i=0; i<points_center.size(); ++i)
            data_points.push_back(points_center[i].template head<2>());
        data_points.push_back(points_center[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_center.size(); i++)
        {
            if(is_center_crossing[i])
                addCrossingData(id_center_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        
        data_points.clear();
        for(int i=0; i<points_outer.size(); ++i)
            data_points.push_back(points_outer[i].template head<2>());
        data_points.push_back(points_outer[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_outer.size(); i++)
        {
            if(is_outer_crossing[i])
                addCrossingData(id_outer_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }
        
        TV to = TV(0, -2.2, 0.0) * sim.unit;
        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;

        passing_points_id = id_line_crossing_principal;
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        to = TV(2.2, 0.0, 0.0) * sim.unit;
        passing_points_id = id_line_crossing_principal_horizontal;
        passing_points.clear();
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        int idx =0;
        for (int i = 1; i < thetas.size(); i += 2)
        {
            passing_points = {points_inner[i], points_center[i], points_outer[i]};
            passing_points_id = {id_line_crossings[idx][0], id_line_crossings[idx][1], id_line_crossings[idx][2]};
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

            for (int j = 0; j < id_line_crossings[idx].size(); j++)
                addCrossingData(id_line_crossings[idx][j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

            ++idx;
        }

        for (auto rod : sim.Rods)
        {
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);
        }

        for(int i=1; i<id_line_crossing_principal.size(); ++i)
        {
            sim.rod_crossings[id_line_crossing_principal[i]]->is_fixed = false;
            sim.rod_crossings[id_line_crossing_principal[i]]->sliding_ranges[1] = Range::Ones();

            auto crossing = sim.rod_crossings[id_line_crossing_principal[i]];
            
            auto rod = sim.Rods[crossing->rods_involved[1]];
            int cnt = 0;
            for (int idx : rod->dof_node_location)
            {
                if(rod->indices[idx] == crossing->node_idx)
                    rod->fixed_by_crossing[cnt] = false;     
                cnt ++;
            }
            
        }

        for(int i=1; i<id_line_crossing_principal_horizontal.size(); ++i)
        {
            sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->is_fixed = false;
            sim.rod_crossings[id_line_crossing_principal_horizontal[i]]->sliding_ranges[1] = Range::Ones();

            auto crossing = sim.rod_crossings[id_line_crossing_principal_horizontal[i]];
            // std::cout << crossing->node_idx << std::endl;
            auto rod = sim.Rods[crossing->rods_involved[1]];
            int cnt = 0;
            for (int idx : rod->dof_node_location)
            {
                if(rod->indices[idx] == crossing->node_idx)
                    rod->fixed_by_crossing[cnt] = false;     
                cnt ++;
            }
        }

        // for(int i=0; i<sim.rod_crossings.size(); ++i)
        //     std::cout << sim.rod_crossings[i]->sliding_ranges.size() << " shjdkasdk " << std::endl;

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        // Offset offset;
        // sim.Rods[3]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = 0.5 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        // sim.Rods[4]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0.5 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        auto incrementalBC = [](EoLRodSim<T, dim>&_sim, int step)->void
        {
            Offset offset;
            if (step < 8)
            {
                _sim.Rods[3]->backOffsetReduced(offset);
                _sim.dirichlet_dof[offset[0]] = 0;
                _sim.dirichlet_dof[offset[1]] = -0.1 * _sim.unit;
                if (step == 0)
                    _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
                else
                    _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
            }
            else
            {
                _sim.Rods[4]->backOffsetReduced(offset);
                _sim.dirichlet_dof[offset[0]] = 0.05 * _sim.unit;
                _sim.dirichlet_dof[offset[1]] = 0.0;
                if (step == 0)
                    _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
                else
                    _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
                for (int n = 0; n < 5; n++)
                {
                    _sim.Rods[2]->fixPointLagrangian(0 + n, TV::Zero(), _sim.dirichlet_dof);
                    _sim.Rods[2]->fixPointLagrangian(_sim.Rods[2]->indices.size() - 1 - n, TV::Zero(), _sim.dirichlet_dof);    
                    _sim.dirichlet_dof[_sim.Rods[2]->theta_reduced_dof_start_offset + n] = 0;
                    _sim.dirichlet_dof[_sim.Rods[2]->theta_reduced_dof_start_offset + _sim.Rods[2]->numSeg() - 1 - n] = 0;

                    int loc = _sim.Rods[2]->dof_node_location[2];
                    _sim.Rods[2]->fixPointLagrangian(loc - n, TV::Zero(), _sim.dirichlet_dof);
                    _sim.Rods[2]->fixPointLagrangian(loc + n, TV::Zero(), _sim.dirichlet_dof);
                    _sim.dirichlet_dof[_sim.Rods[2]->theta_reduced_dof_start_offset + loc + n] = 0;
                    _sim.dirichlet_dof[_sim.Rods[2]->theta_reduced_dof_start_offset + loc - 1 - n] = 0;
                }
            }
            
        };

        sim.incremental_steps = 11;
        sim.incremental_bc = incrementalBC;

        // Offset offset;
        // sim.Rods[3]->getEntryReduced(sim.Rods[3]->indices[sim.Rods[3]->dof_node_location.back()], offset);
        // sim.dirichlet_dof[offset[3]] = 0.01 * sim.unit;

        // sim.Rods[2]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-1] = 0;

        for (int n = 0; n < 5; n++)
        {
            sim.Rods[2]->fixPointLagrangian(0 + n, TV::Zero(), sim.dirichlet_dof);
            sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 1 - n, TV::Zero(), sim.dirichlet_dof);    
            sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + n] = 0;
            sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg() - 1 - n] = 0;

            // int loc = sim.Rods[2]->dof_node_location[2];
            // sim.Rods[2]->fixPointLagrangian(loc - n, TV::Zero(), sim.dirichlet_dof);
            // sim.Rods[2]->fixPointLagrangian(loc + n, TV::Zero(), sim.dirichlet_dof);
            // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + loc + n] = 0;
            // sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + loc - 1 - n] = 0;

        }

        
        T r = 0.3 * sim.unit;
        TV center = TV(1.4*cos(thetas[0]), 1.4*sin(thetas[0]), 0) * sim.unit;
        TV center2 = TV(1.4*cos(thetas[2]), 1.4*sin(thetas[2]), 0) * sim.unit;

        auto circle = [r, center](const TV& x)->bool
        {
            return (x - center).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        // sim.fixRegionAvoidRod(circle, 3);
        // sim.fixRegionAvoidRod(circle2, 4);

        sim.fixCrossing();

        

        // sim.boundary_spheres.push_back(std::make_pair(center, r ));
        // sim.boundary_spheres.push_back(std::make_pair(center2, r ));

        
        

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            if (crossing->is_fixed)
                continue;
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }

        GCodeGenerator<T, dim>(sim, "shelter.gcode").generateShelterScene();
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildFullCircleScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = true;
        sim.add_pbc_twisting = true;
        sim.add_pbc = true;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.01;
        // sim.unit = 1.0;
        sim.visual_R = 0.012;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };        

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        int n_row = 5, n_col = 5;

        T r = 0.02; 
        T d = 0.9 * 2.0 * r;
        TV bottom_left_center = TV(0, 0, 0);

        // add top and right crossings
        for (int row = 0; row < n_row - 1; row++)
        {
            
            for (int col = 0; col < n_col - 1; col++)
            {
                TV center = bottom_left_center + TV(d * row, d * col, 0);
                TV next_right = bottom_left_center + TV(d * row, d * (col + 1), 0);
                TV next_top = bottom_left_center + TV(d * (row + 1), d * col, 0);

                TV ixn0_right, ixn1_right;
                TV ixn0_top, ixn1_top;

                // circleCircleIntersection(center, r, next_right, r, ixn0_right, ixn1_right);
                // circleCircleIntersection(center, r, next_right, r, ixn0_right, ixn1_right);

            }    
        }
        
        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        Offset offset;
        sim.Rods[0]->frontOffsetReduced(offset);
        for (int d = 0; d < dim; d++) sim.dirichlet_dof[offset[d]] = 0;

        sim.fixCrossing();

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildPeriodicCircleScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = true;
        sim.add_pbc_twisting = true;
        sim.add_pbc = true;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.01;
        // sim.unit = 1.0;
        sim.visual_R = 0.012;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };        

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        T r = 0.55 * sim.unit;

        TV center0 = TV(0, 0, 0) * sim.unit;
        TV center1 = TV(1, 0, 0) * sim.unit;
        TV center2 = TV(1, 1, 0) * sim.unit;
        TV center3 = TV(0, 1, 0) * sim.unit;
        
        std::vector<TV> centers = {center0, center1, center2, center3};

        // add 8 boundary points counterclock wise
        for (int i = 0; i < centers.size(); i++)
        {
            TV current = centers[i], next = centers[(i + 1) % centers.size()];
            TV v0 = next - r  / sim.unit * (next - current);
            TV v1 = current + r / sim.unit * (next - current);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
        }
        

        TV vtx8, vtx9, vtx10, vtx11, dummy;
        circleCircleIntersection(center0, r, center1, r, vtx8, dummy);
        circleCircleIntersection(center1, r, center2, r, vtx9, dummy);
        circleCircleIntersection(center2, r, center3, r, vtx10, dummy);
        circleCircleIntersection(center3, r, center0, r, vtx11, dummy);

        addCrossingPoint(nodal_positions, vtx8, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, vtx9, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, vtx10, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, vtx11, full_dof_cnt, node_cnt);


        auto addCurvedRodFromIDs = [&](const std::vector<int>& ids)
        {
            std::vector<TV> passing_points;
            std::vector<TV2> data_points;
            for (int id : ids)
            {
                passing_points.push_back(nodal_positions[id]);
                data_points.push_back(passing_points.back().template head<2>());
            }
            addCurvedRod(data_points, passing_points, ids, sub_div, full_dof_cnt, node_cnt, rod_cnt, false);
            for (int i = 0; i < ids.size(); i++)
            {
                addCrossingData(ids[i], rod_cnt - 1, sim.Rods[rod_cnt-1]->dof_node_location[i]);
            }
        };

        addCurvedRodFromIDs({0, 8, 9, 3});
        addCurvedRodFromIDs({2, 9, 10, 5});
        addCurvedRodFromIDs({4, 10, 11, 7});
        addCurvedRodFromIDs({6, 11, 8, 1});
        

        auto setPBCData = [&](int rod_idx0, int rod_idx1, int direction, bool unique, bool reverse)
        {
            auto rod0 = sim.Rods[rod_idx0];
            auto rod1 = sim.Rods[rod_idx1];
            Offset end0, end1;
            if (reverse)
            {
                rod0->backOffset(end0); rod1->frontOffset(end1);
            }
            else
            {
                rod0->frontOffset(end0); rod1->backOffset(end1);
            }
            if (unique)
                sim.pbc_pairs_reference[direction] = std::make_pair(std::make_pair(end0, end1), 
                        std::make_pair(rod0->rod_id, rod1->rod_id));
            sim.pbc_pairs.push_back(std::make_pair(direction, std::make_pair(end0, end1)));

            Offset a, b;
            if (reverse)
            {
                rod0->getEntryByLocation(rod1->indices.size() - 2, a);
                rod1->getEntryByLocation(1, b); 
                sim.pbc_bending_pairs.push_back({end0, a, b, end1});
                sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});
            }
            else
            {
                rod0->getEntryByLocation(1, a); 
                rod1->getEntryByLocation(rod1->indices.size() - 2, b);
                sim.pbc_bending_pairs.push_back({end0, a, b, end1});
                sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});
            }
            
        };

        // now we set the periodic data
        setPBCData(3, 0, 0, true, false);
        setPBCData(2, 1, 0, false, true);

        setPBCData(0, 1, 1, true, false);
        setPBCData(3, 2, 1, false, true);

        for (auto crossing : sim.rod_crossings)
        {
            if (crossing->node_idx < 8)
                continue;
            // crossing->is_fixed = false;
            // crossing->sliding_ranges[0] = Range(0.9, 0.9);
        }
        
        // auto crossing = sim.rod_crossings[8];
        // crossing->is_fixed = false;
        // crossing->sliding_ranges[0] = Range::Ones();
        // crossing = sim.rod_crossings[11];
        // crossing->is_fixed = false;
        // crossing->sliding_ranges[0] = Range::Ones();


        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        Offset offset;
        sim.Rods[0]->frontOffsetReduced(offset);
        for (int d = 0; d < dim; d++) sim.dirichlet_dof[offset[d]] = 0;

        sim.fixCrossing();

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
        // GCodeGenerator<T, dim>(this->sim, "circle_patch_every_other.gcode").circlePatchGCode(3, 5, 3, false);
        // GCodeGenerator<T, dim>(this->sim, "circle_patch.gcode").circlePatchGCode(5, 6, 1, true);
        // GCodeGenerator<T, dim>(this->sim, "circle_patch_test.gcode").circlePatchGCode(3, 4, 5, true);
        // GCodeGenerator<T, dim>(this->sim, "circle_patch_fused.gcode").circlePatchGCode(5, 6, 2, true);
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildRandomPatchScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        // sim.unit = 0.05;
        sim.unit = 1.0;
        sim.visual_R = 0.012;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };



        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        TV bottom_left = TV(0, 0, 0) * sim.unit;
        TV top_right = TV(1, 1, 0) * sim.unit;

        int n_sample = 2;
        TV bottom_right = TV(1, 0, 0) * sim.unit;
        TV top_left = TV(0, 1, 0) * sim.unit;

        T dx = (bottom_right - bottom_left).norm() / (1 + n_sample);
        T dy = (top_left - bottom_left).norm() / (1 + n_sample);

        addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);

        T epsilon = 0;//0.05 * sim.unit;
        for (int i = 0; i < n_sample + 1; i++)
        {
            TV node = TV(bottom_left[0] + epsilon+ dx * i + dx * zeta(), bottom_left[1], 0);
            addCrossingPoint(nodal_positions, node, full_dof_cnt, node_cnt);
            node = TV(top_right[0], bottom_left[1] + epsilon + dy * i + dy * zeta(), 0);
            addCrossingPoint(nodal_positions, node, full_dof_cnt, node_cnt);
            node = TV(top_right[0] - (epsilon + dx * i + dx * zeta()), top_right[1], 0);
            addCrossingPoint(nodal_positions, node, full_dof_cnt, node_cnt);
            node = TV(bottom_left[0], top_right[1] - (epsilon + dy * i + dy * zeta()), 0);
            addCrossingPoint(nodal_positions, node, full_dof_cnt, node_cnt);
        }

        std::vector<std::vector<int>> random_locations(4, std::vector<int>(n_sample + 1));
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < n_sample + 1; j++)
            {
                random_locations[i][j] = j;    
            }
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(random_locations[i].begin(), random_locations[i].end(), std::default_random_engine(seed));
        }

        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;
        
        int n_new_points_edge = n_sample + 1;

        // compute crossings and rod 

        using Edge = std::pair<TV2, TV2>;

        std::vector<Edge> existing_edges;

        struct TemporaryRod
        {
            std::vector<TV> points;
            std::vector<int> ids;
        };

        std::vector<TemporaryRod> temporary_rods;
        // four edges
        for (int j = 0; j < n_new_points_edge; j++)
        {
            for (int edge = 0; edge < 4; edge++)
            {
                int next_edge = (edge + 1) % 4;
                // int next_edge = (edge + 2);
                
                // int node_from = 4 + random_locations[edge][j] * 4 + edge;
                // int node_to = 4 + random_locations[next_edge][j] * 4 + next_edge;

                int node_from = 4 + j * 4 + edge;
                int node_to = 4 + j * 4 + next_edge;

                // std::cout << "sample " << j << " at edge " << edge << std::endl;
                std::cout << "from " << node_from << " to " << node_to << std::endl;
                
                TV from = nodal_positions[node_from];
                TV to = nodal_positions[node_to];
                
                TV2 from2d = from.template head<2>();
                TV2 to2d = to.template head<2>();
                std::vector<TV> intersections;
                int crossing_cnt_temp = sim.rod_crossings.size();
                for (Edge edge : existing_edges)
                {
                    TV2 intersection;
                    if (lineSegementsIntersect2D(from2d, to2d, edge.first, edge.second, intersection))
                    {
                        if ((intersection - from2d).norm() < 1e-6 || (intersection - to2d).norm() < 1e-6)
                            continue;
                        std::cout << "intersect at " << intersection.transpose() << std::endl;
                        intersections.push_back(TV(intersection[0], intersection[1], 0.0));
                        addCrossingPoint(nodal_positions, intersections.back(), full_dof_cnt, node_cnt);
                    }
                }
                std::cout << "intersections: " << intersections.size() << std::endl;
                TemporaryRod temp_rod;
                temp_rod.points.push_back(from);
                temp_rod.ids.push_back(node_from);
                for (TV intersection : intersections)
                {
                    temp_rod.points.push_back(intersection);
                    temp_rod.ids.push_back(crossing_cnt_temp++);
                }      
                temp_rod.points.push_back(to);
                temp_rod.ids.push_back(node_to);
                std::sort(temp_rod.points.begin(), temp_rod.points.end(), 
                    [&](const TV & a, const TV & b) -> bool
                    { 
                    return (a - from).norm() > (b - from).norm(); 
                });
                std::sort(temp_rod.ids.begin(), temp_rod.ids.end(), 
                    [&](const int a, const int b) -> bool
                    { 
                    return (nodal_positions[a] - from).norm() > (nodal_positions[b] - from).norm(); 
                });
                existing_edges.push_back(std::make_pair(from2d, to2d));
                temporary_rods.push_back(temp_rod);
            }
        }
        std::cout << "Intersection done" << std::endl;
        
        for (TemporaryRod temp_rod : temporary_rods)
        {
            

            addAStraightRod(temp_rod.points.front(), temp_rod.points.back(), temp_rod.points, temp_rod.ids, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < temp_rod.ids.size(); i++)
                addCrossingData(temp_rod.ids[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
        }

        {
            passing_points_id.push_back(0);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            for (int i = 0; i < n_sample + 1; i++)
            {
                passing_points_id.push_back(4 + i * 4);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            passing_points_id.push_back(1);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
            
            passing_points.clear();
            passing_points_id.clear();

            passing_points_id.push_back(1);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            for (int i = 0; i < n_sample + 1; i++)
            {
                passing_points_id.push_back(4 + i * 4 + 1);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            passing_points_id.push_back(2);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
            passing_points.clear();
            passing_points_id.clear();

            passing_points_id.push_back(2);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            for (int i = 0; i < n_sample + 1; i++)
            {
                passing_points_id.push_back(4 + i * 4 + 2);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            passing_points_id.push_back(3);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
            passing_points.clear();
            passing_points_id.clear();

            passing_points_id.push_back(3);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            for (int i = 0; i < n_sample + 1; i++)
            {
                passing_points_id.push_back(4 + i * 4 + 3);
                passing_points.push_back(nodal_positions[passing_points_id.back()]);
            }
            passing_points_id.push_back(0);
            passing_points.push_back(nodal_positions[passing_points_id.back()]);
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_id.size(); i++)
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
            passing_points.clear();
            passing_points_id.clear();
        }

        
        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        sim.fixCrossing();

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        // for (auto& crossing : sim.rod_crossings)
        // {
        //     Offset off;
        //     sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
        //     T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
        //     int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
        //     sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        // }
    }
}


template<class T, int dim>
void UnitPatch<T, dim>::buildShoeScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        // sim.unit = 0.05;
        sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };



        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;
        T theta = M_PI / 4.0;

        TV v0 = TV(0.5, -0.5, 0.0) * sim.unit;
        TV v1 = TV(0.5, 0.0, 0.0) * sim.unit;
        addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);

        TV v2 = TV(0.5 - 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta), 0.0) * sim.unit;
        TV v3 = TV(0.5 + 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta), 0.0) * sim.unit;
        addCrossingPoint(nodal_positions, v2, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v3, full_dof_cnt, node_cnt);
        // TV v4 = 0.5 * (v0 + v1);
        // addCrossingPoint(nodal_positions, v4, full_dof_cnt, node_cnt);
        
        std::vector<TV2> data_points;
        data_points.push_back(TV2(0.5, -0.5) * sim.unit);
        data_points.push_back(TV2(0, 0) * sim.unit);
        data_points.push_back(TV2(0.5 - 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        data_points.push_back(TV2(0.5, 0.0) * sim.unit);
        data_points.push_back(TV2(0.5 + 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        data_points.push_back(TV2(1.0, 0.0) * sim.unit);
        data_points.push_back(TV2(0.5, -0.5) * sim.unit);        

        std::vector<TV> passing_points = {v0, v2, v1, v3};
        std::vector<int> passing_points_id = {0, 2, 1, 3};

        addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        for (int i = 0; i < passing_points_id.size(); i++)
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
        
        // data_points.clear();
        // data_points.push_back(TV2(0.5 - 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        // data_points.push_back(v4.template head<2>());
        // data_points.push_back(TV2(0.5 + 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        // passing_points = {v2, v4, v3};
        // passing_points_id = {2, 4, 3};
        // addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, false);
        // for (int i = 0; i < passing_points_id.size(); i++)
        //     addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        // passing_points = {v1, v4, v0};
        // passing_points_id = {1, 4, 0};
        // TV to = TV(0.5, -1.5, 0.0) * sim.unit;
        // addAStraightRod(v1, to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        // for (int i = 0; i < passing_points_id.size(); i++)
        //     addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
        
        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        // auto crossing = sim.rod_crossings[0];
        // crossing->is_fixed = false;
        // crossing->sliding_ranges[1] = Range::Ones();
        
    
        // Offset offset;
        // sim.Rods[1]->backOffsetReduced(offset);

        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = 0.6 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0.01 * sim.unit;

        // T r = 0.1 * sim.unit;
        // TV center;
        // sim.getCrossingPosition(0, center);

        // auto circle = [r, center](const TV& x)->bool
        // {
        //     return (x - center).norm() < r;
        // };

        

        // sim.fixRegion(circle);

        sim.fixCrossing();

        // sim.Rods[0]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        // sim.Rods[0]->fixPointLagrangian(sim.Rods[0]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        
        // sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset + sim.Rods[0]->numSeg()-1] = 0;

        // sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildDomeScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.05;
        //sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;
        
        std::vector<T> thetas;
        int theta_sub = 8;
        for (int i = 0; i< theta_sub; i++)
            thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

        int quater = thetas.size()/4;

        std::vector<TV> points_on_circle;
        for(int i=0; i<thetas.size(); ++i)
            points_on_circle.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);

        for(int i=0; i<thetas.size(); ++i)
            points_on_circle.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
        
        for(int i=0; i<thetas.size(); ++i)
            points_on_circle.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);

        std::vector<TV> nodal_positions;
        addCrossingPoint(nodal_positions, points_on_circle[0], full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, points_on_circle[thetas.size()/2], full_dof_cnt, node_cnt);

        for (int i = 0; i < 3; i++)
        {
            
            addCrossingPoint(nodal_positions, points_on_circle[i * thetas.size() + quater], full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, points_on_circle[i * thetas.size() + quater * 3], full_dof_cnt, node_cnt);
        }

        std::vector<TV> passsing_points;
        std::vector<int> passsing_points_id;
        for (int i = 0; i < 3; i++)
        {
            
            if (i == 0)
            {
                passsing_points = {nodal_positions[0], nodal_positions[2 + i*2], nodal_positions[1], nodal_positions[2 + i*2 + 1]};
                passsing_points_id = { 0, 2 + i * 2, 1, 2 + i * 2 + 1 };
            }
            else
            {
                passsing_points = {nodal_positions[2 + i*2], nodal_positions[2 + i*2 + 1]};
                passsing_points_id = { 2 + i * 2, 2 + i * 2 + 1 };
            }
            
            std::vector<TV2> data_points;
            for (int j = 0; j < thetas.size(); j++)    
                data_points.push_back(points_on_circle[i*thetas.size() + j].template head<2>());
            data_points.push_back(points_on_circle[i*thetas.size()].template head<2>());

            addCurvedRod(data_points, passsing_points, passsing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
            for (int j = 0; j < passsing_points_id.size(); j++)
                addCrossingData(passsing_points_id[j], rod_cnt - 1, sim.Rods[rod_cnt-1]->dof_node_location[j]);
        }
        
        passsing_points = {nodal_positions[0], nodal_positions[1]};
        passsing_points_id = {0, 1};
        TV drag = TV(0, -3.0, 0) * sim.unit;
        addAStraightRod(drag, passsing_points.back(), passsing_points, passsing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int j = 0; j < passsing_points_id.size(); j++)
            addCrossingData(passsing_points_id[j], rod_cnt - 1, sim.Rods[rod_cnt-1]->dof_node_location[j]);
        
        passsing_points.resize(6);
        passsing_points_id.resize(6);
        for (int i = 0; i < 3; i++)
        {
            passsing_points[i] = nodal_positions[2 + i * 2 + 0];
            passsing_points[5-i] = nodal_positions[2 + i * 2 + 1];
            passsing_points_id[i] = 2 + i * 2;
            passsing_points_id[5 - i] = 2 + i * 2 + 1;
        }

        addAStraightRod(passsing_points.front(), passsing_points.back(), passsing_points, passsing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int j = 0; j < passsing_points_id.size(); j++)
            addCrossingData(passsing_points_id[j], rod_cnt - 1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        auto crossing = sim.rod_crossings[0];
        crossing->is_fixed = false;
        crossing->sliding_ranges[1] = Range(1, 1);

    
        Offset offset;
        sim.Rods[3]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -0.5 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0;

        T r = 0.1 * sim.unit;
        TV center= TV(0, -1.4, 0) * sim.unit;

        auto circle = [r, center](const TV& x)->bool
        {
            return (x - center).norm() < r;
        };

        sim.fixRegion(circle);

        sim.fixCrossing();

        sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildTestSceneJuan(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.01;
        //sim.unit = 1.0;
        sim.visual_R = 0.02;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

    

        // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
        std::vector<T> thetas;
        int theta_sub = 8;
        for (int i = 0; i< theta_sub; i++)
            thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

        std::cout << "thetas ";
        for(int i=0; i<thetas.size(); ++i)
            std::cout << " " << thetas[i];
        std::cout << std::endl;

        std::vector<bool> is_inner_crossing = {true, true, false, true, true, true, false, true};
        std::vector<bool> is_center_crossing = {true, true, false, true, true, true, false, true};
        std::vector<bool> is_outer_crossing = {true, true, false, true, true, true, false, true};

        std::vector<int> id_inner_crossing = {0, 1, -1, 2, 3, 4, -1, 5};
        std::vector<int> id_center_crossing = {6, 7, -1, 8, 9, 10, -1, 11};
        std::vector<int> id_outer_crossing = {12, 13, -1, 14, 15, 16, -1, 17};

        std::vector<int> id_line_crossing_principal = {15, 9, 3, 0, 6, 12};
        std::vector<std::vector<int>> id_line_crossings = {{1, 7, 13}, {2, 8, 14}, {4, 10, 16}, {5, 11, 17}};

        std::vector<int> crossings_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

        std::vector<TV> cross_point_collection;

        for(int i=0; i<18; ++i)
                sim.rod_crossings.push_back(new RodCrossing<T, dim>(crossings_ids[i], std::vector<int>())); 

        int cur_idx = 0;

        std::vector<TV> points_inner;
        std::vector<TV> crossing_points_inner;
        std::vector<int> indices_inner;
        for(int i=0; i<thetas.size(); ++i)
            points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_inner.size(); ++i)
        {
            if(id_inner_crossing[i]!=-1)
            {
                addPoint(points_inner[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_inner[i]);
                crossing_points_inner.push_back(points_inner[i]);
                indices_inner.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_center;
        std::vector<TV> crossing_points_center;
        std::vector<int> indices_center;
        for(int i=0; i<thetas.size(); ++i)
            points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_center.size(); ++i)
        {
            if(id_center_crossing[i]!=-1)
            {
                addPoint(points_center[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_center[i]);
                crossing_points_center.push_back(points_center[i]);
                indices_center.push_back(cur_idx++);
            }
        }

        std::vector<TV> points_outer;
        std::vector<TV> crossing_points_outer;
        std::vector<int> indices_outer;
        for(int i=0; i<thetas.size(); ++i)
            points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
        for(int i=0; i<points_outer.size(); ++i)
        {
            if(id_outer_crossing[i]!=-1)
            {
                addPoint(points_outer[i], full_dof_cnt, node_cnt);
                cross_point_collection.push_back(points_outer[i]);
                crossing_points_outer.push_back(points_outer[i]);
                indices_outer.push_back(cur_idx++);
            }
        }
                
        TV to = TV(0, -2.1, 0.0) * sim.unit;

        std::vector<TV2> data_points;
        for(int i=0; i<points_inner.size(); ++i)
            data_points.push_back(points_inner[i].template head<2>());
        data_points.push_back(points_inner[0].template head<2>());

        addCurvedRod(data_points, crossing_points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
        int ii = 0;
        for (int i = 0; i < points_inner.size(); i++)
        {
            if(is_inner_crossing[i])
                addCrossingData(id_inner_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        data_points.clear();
        for(int i=0; i<points_center.size(); ++i)
            data_points.push_back(points_center[i].template head<2>());
        data_points.push_back(points_center[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_center.size(); i++)
        {
            if(is_center_crossing[i])
                addCrossingData(id_center_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }

        
        data_points.clear();
        for(int i=0; i<points_outer.size(); ++i)
            data_points.push_back(points_outer[i].template head<2>());
        data_points.push_back(points_outer[0].template head<2>());
        
        addCurvedRod(data_points, crossing_points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

        ii = 0;
        for (int i = 0; i < points_outer.size(); i++)
        {
            if(is_outer_crossing[i])
                addCrossingData(id_outer_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
        }
        
        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;

        passing_points_id = id_line_crossing_principal;
        for(int i=0; i<passing_points_id.size(); ++i)
            passing_points.push_back(cross_point_collection[passing_points_id[i]]);

        addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        // std::cout << "rod I want " << rod_cnt - 1 << std::endl;

        for (int i = 0; i < id_line_crossing_principal.size(); i++)
            addCrossingData(id_line_crossing_principal[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

        int idx =0;
        for (int i = 1; i < thetas.size(); i += 2)
        {
            passing_points = {points_inner[i], points_center[i], points_outer[i]};
            passing_points_id = {id_line_crossings[idx][0], id_line_crossings[idx][1], id_line_crossings[idx][2]};
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

            for (int j = 0; j < id_line_crossings[idx].size(); j++)
                addCrossingData(id_line_crossings[idx][j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

            ++idx;
        }

        sim.rod_crossings[3]->is_fixed = false;
        sim.rod_crossings[3]->sliding_ranges[1] = Range::Ones();

        sim.rod_crossings[9]->is_fixed = false;
        sim.rod_crossings[9]->sliding_ranges[1] = Range::Ones();

        // sim.rod_crossings[15]->is_fixed = false;
        // sim.rod_crossings[15]->sliding_ranges[1] = Range::Ones();

        sim.rod_crossings[12]->is_fixed = false;
        sim.rod_crossings[12]->sliding_ranges[1] = Range::Ones();

        sim.rod_crossings[6]->is_fixed = false;
        sim.rod_crossings[6]->sliding_ranges[1] = Range::Ones();

        sim.rod_crossings[0]->is_fixed = false;
        sim.rod_crossings[0]->sliding_ranges[1] = Range::Ones();
        // for(int i=0; i<sim.rod_crossings.size(); ++i)
        //     std::cout << sim.rod_crossings[i]->sliding_ranges.size() << " shjdkasdk " << std::endl;

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        // Offset offset;
        // sim.Rods[3]->backOffsetReduced(offset);
        // // sim.dirichlet_dof[offset[0]] = 0;
        // // sim.dirichlet_dof[offset[1]] = 0.5 * sim.unit;
        // // sim.dirichlet_dof[offset[2]] = 0.001 * sim.unit;

        auto incrementalBC = [](EoLRodSim<T, dim>&_sim, int step)->void
        {
            Offset offset;
            _sim.Rods[3]->backOffsetReduced(offset);
            _sim.dirichlet_dof[offset[0]] = 0;
            _sim.dirichlet_dof[offset[1]] = 0.1 * _sim.unit;
            if (step == 0)
                _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
            else
                _sim.dirichlet_dof[offset[2]] = 0.001 * _sim.unit;
        };

        sim.incremental_steps = 5;
        sim.incremental_bc = incrementalBC;

        // Offset offset;
        // sim.Rods[3]->getEntryReduced(sim.Rods[3]->dof_node_location.back(), offset);
        // sim.dirichlet_dof[offset[3]] = 0.2 * sim.unit;

        sim.Rods[2]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(2, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[2]->fixPointLagrangian(sim.Rods[2]->indices.size() - 3, TV::Zero(), sim.dirichlet_dof);
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset] = 0;
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + 1] = 0;
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-1] = 0;
        sim.dirichlet_dof[sim.Rods[2]->theta_reduced_dof_start_offset + sim.Rods[2]->numSeg()-2] = 0;

        T r = 0.1 * sim.unit;
        TV center= TV(1.4*cos(thetas[0]), 1.4*sin(thetas[0]), 0) * sim.unit;

        auto circle = [r, center](const TV& x)->bool
        {
            return (x - center).norm() < r;
        };

        sim.fixRegion(circle);

        sim.fixCrossing();

        sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
    }
    // if constexpr (dim == 3)
    // {
    //     auto unit_yarn_map = sim.yarn_map;
    //     sim.yarn_map.clear();
        
    //     clearSimData();

    //     sim.add_rotation_penalty = false;
    //     sim.add_pbc_bending = false;
    //     sim.add_pbc_twisting = false;
    //     sim.add_pbc = false;

    //     sim.add_contact_penalty=true;
    //     sim.new_frame_work = true;
    //     sim.add_eularian_reg = true;

    //     sim.ke = 1e-4;

    //     sim.unit = 0.05;
    //     //sim.unit = 1.0;
    //     sim.visual_R = 0.02;

    //     auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
    //     {
    //         sim.rod_crossings[crossing_idx]->is_fixed = true;
    //         sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
    //         sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
    //         sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
    //     };

    //     std::vector<Eigen::Triplet<T>> w_entry;
    //     int full_dof_cnt = 0;
    //     int node_cnt = 0;
    //     int rod_cnt = 0;

    

    //     // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
    //     std::vector<T> thetas;
    //     int theta_sub = 8;
    //     for (int i = 0; i< theta_sub; i++)
    //         thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

    //     std::cout << "thetas ";
    //     for(int i=0; i<thetas.size(); ++i)
    //         std::cout << " " << thetas[i];
    //     std::cout << std::endl;

    //     std::vector<bool> is_inner_crossing = {true, true, false, true, true, true, false, true};
    //     std::vector<bool> is_center_crossing = {true, true, false, true, false, true, false, true};
    //     std::vector<bool> is_outer_crossing = {true, true, false, true, false, true, false, true};

    //     std::vector<int> id_inner_crossing = {0, 1, -1, 2, 3, 4, -1, 5};
    //     std::vector<int> id_center_crossing = {6, 7, -1, 8, -1, 9, -1, 10};
    //     std::vector<int> id_outer_crossing = {11, 12, -1, 13, -1, 14, -1, 15};

    //     std::vector<int> id_line_crossing_principal = {3, 0, 6, 11};
    //     std::vector<std::vector<int>> id_line_crossings = {{1, 7, 12}, {2, 8, 13}, {4, 9, 14}, {5, 10, 15}};

    //     std::vector<int> crossings_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    //     for(int i=0; i<16; ++i)
    //         sim.rod_crossings.push_back(new RodCrossing<T, dim>(crossings_ids[i], std::vector<int>())); 

    //     int cur_idx = 0;

    //     std::vector<TV> points_inner;
    //     std::vector<TV> crossing_points_inner;
    //     std::vector<int> indices_inner;
    //     for(int i=0; i<thetas.size(); ++i)
    //         points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
    //     for(int i=0; i<points_inner.size(); ++i)
    //     {
    //         if(id_inner_crossing[i]!=-1)
    //         {
    //             addPoint(points_inner[i], full_dof_cnt, node_cnt);
    //             crossing_points_inner.push_back(points_inner[i]);
    //             indices_inner.push_back(cur_idx++);
    //         }
    //     }

    //     std::vector<TV> points_center;
    //     std::vector<TV> crossing_points_center;
    //     std::vector<int> indices_center;
    //     for(int i=0; i<thetas.size(); ++i)
    //         points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
    //     for(int i=0; i<points_center.size(); ++i)
    //     {
    //         if(id_center_crossing[i]!=-1)
    //         {
    //             addPoint(points_center[i], full_dof_cnt, node_cnt);
    //             crossing_points_center.push_back(points_center[i]);
    //             indices_center.push_back(cur_idx++);
    //         }
    //     }

    //     std::vector<TV> points_outer;
    //     std::vector<TV> crossing_points_outer;
    //     std::vector<int> indices_outer;
    //     for(int i=0; i<thetas.size(); ++i)
    //         points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
    //     for(int i=0; i<points_outer.size(); ++i)
    //     {
    //         if(id_outer_crossing[i]!=-1)
    //         {
    //             addPoint(points_outer[i], full_dof_cnt, node_cnt);
    //             crossing_points_outer.push_back(points_outer[i]);
    //             indices_outer.push_back(cur_idx++);
    //         }
    //     }
                
    //     TV to = TV(0, -3.0, 0.0) * sim.unit;

    //     std::vector<TV2> data_points;
    //     for(int i=0; i<points_inner.size(); ++i)
    //         data_points.push_back(points_inner[i].template head<2>());
    //     data_points.push_back(points_inner[0].template head<2>());

    //     addCurvedRod(data_points, crossing_points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
    //     int ii = 0;
    //     for (int i = 0; i < points_inner.size(); i++)
    //     {
    //         if(is_inner_crossing[i])
    //             addCrossingData(id_inner_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
    //     }

    //     data_points.clear();
    //     for(int i=0; i<points_center.size(); ++i)
    //         data_points.push_back(points_center[i].template head<2>());
    //     data_points.push_back(points_center[0].template head<2>());
        
    //     addCurvedRod(data_points, crossing_points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

    //     ii = 0;
    //     for (int i = 0; i < points_center.size(); i++)
    //     {
    //         if(is_center_crossing[i])
    //             addCrossingData(id_center_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
    //     }

        
    //     data_points.clear();
    //     for(int i=0; i<points_outer.size(); ++i)
    //         data_points.push_back(points_outer[i].template head<2>());
    //     data_points.push_back(points_outer[0].template head<2>());
        
    //     addCurvedRod(data_points, crossing_points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

    //     ii = 0;
    //     for (int i = 0; i < points_outer.size(); i++)
    //     {
    //         if(is_outer_crossing[i])
    //             addCrossingData(id_outer_crossing[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[ii++]);
    //     }
        
    //     std::vector<TV> passing_points;
    //     std::vector<int> passing_points_id;
    //     passing_points.push_back(crossing_points_inner[3]);
    //     passing_points.push_back(crossing_points_inner[0]); passing_points.push_back(crossing_points_center[0]); passing_points.push_back(crossing_points_outer[0]);

    //     passing_points_id.push_back(3);
    //     passing_points_id.push_back(0); passing_points_id.push_back(6); passing_points_id.push_back(11);

    //     addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

    //     // std::cout << "rod I want " << rod_cnt - 1 << std::endl;

    //     for (int i = 0; i < id_line_crossing_principal.size(); i++)
    //         addCrossingData(id_line_crossing_principal[i], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[i]);

    //     int idx =0;
    //     for (int i = 1; i < thetas.size(); i += 2)
    //     {
    //         passing_points = {points_inner[i], points_center[i], points_outer[i]};
    //         passing_points_id = {id_line_crossings[idx][0], id_line_crossings[idx][1], id_line_crossings[idx][2]};
    //         addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

    //         for (int j = 0; j < id_line_crossings[idx].size(); j++)
    //             addCrossingData(id_line_crossings[idx][j], rod_cnt-1, sim.Rods[rod_cnt-1]->dof_node_location[j]);

    //         ++idx;
    //     }


    //     // for(int i=0; i<sim.rod_crossings.size(); ++i)
    //     //     std::cout << sim.rod_crossings[i]->sliding_ranges.size() << " shjdkasdk " << std::endl;

    //     // sim.rod_crossings[11]->is_fixed = false;
    //     // sim.rod_crossings[11]->sliding_ranges[1] = Range::Ones();

    //     // sim.rod_crossings[6]->is_fixed = false;
    //     // sim.rod_crossings[6]->sliding_ranges[1] = Range::Ones();

    //     // sim.rod_crossings[0]->is_fixed = false;
    //     // sim.rod_crossings[0]->sliding_ranges[1] = Range::Ones();

    //     int dof_cnt = 0;
    //     markCrossingDoF(w_entry, dof_cnt);
        
    //     for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
    //     appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
    //     sim.rest_states = deformed_states;
        
    //     sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
    //     sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
    //     for (auto& rod : sim.Rods)
    //     {
    //         rod->fixEndPointEulerian(sim.dirichlet_dof);
    //         rod->setupBishopFrame();
    //     }
        
    
    //     Offset offset;
    //     sim.Rods[3]->backOffsetReduced(offset);
    //     sim.dirichlet_dof[offset[0]] = 0;
    //     sim.dirichlet_dof[offset[1]] = -0.05 * sim.unit;
    //     sim.dirichlet_dof[offset[2]] = 0;

    //     T r = 0.1 * sim.unit;
    //     TV center= TV(0, -1.4, 0) * sim.unit;

    //     auto circle = [r, center](const TV& x)->bool
    //     {
    //         return (x - center).norm() < r;
    //     };

    //     sim.fixRegion(circle);

    //     sim.fixCrossing();

    //     sim.boundary_spheres.push_back(std::make_pair(center, r * 0.5));

    //     sim.perturb = VectorXT::Zero(sim.W.cols());

    //     // for (auto& crossing : sim.rod_crossings)
    //     // {
    //     //     Offset off;
    //     //     sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
    //     //     T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
    //     //     int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
    //     //     sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
    //     // }
    // }


}

// template<class T, int dim>
// void UnitPatch<T, dim>::buildTestSceneJuan(int sub_div)
// {
//     if constexpr (dim == 3)
//     {
//         auto unit_yarn_map = sim.yarn_map;
//         sim.yarn_map.clear();
        
//         clearSimData();

//         sim.add_rotation_penalty = false;
//         sim.add_pbc_bending = false;
//         sim.add_pbc_twisting = false;
//         sim.add_pbc = false;

//         sim.add_contact_penalty=true;
//         sim.new_frame_work = true;
//         sim.add_eularian_reg = true;

//         sim.ke = 1e-4;

//         sim.unit = 0.05;
//         sim.visual_R = 0.02;

//         auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
//         {
//             sim.rod_crossings[crossing_idx]->is_fixed = true;
//             sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
//             sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
//             sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
//         };

//         std::vector<Eigen::Triplet<T>> w_entry;
//         int full_dof_cnt = 0;
//         int node_cnt = 0;
//         int rod_cnt = 0;

    

//         // std::vector<T> thetas = {0, M_PI/3.0, M_PI/2.0, 2.0*M_PI/3.0, M_PI, 4.0*M_PI/3.0, 3.0*M_PI/2.0, 5.0*M_PI/3.0};
//         std::vector<T> thetas;
//         int theta_sub = 8;
//         for (int i = 0; i< theta_sub; i++)
//             thetas.push_back(2.0 * M_PI / theta_sub * ((i + int(0.75 * theta_sub)) % theta_sub));

//         std::vector<TV> points_inner;
//         for(int i=0; i<thetas.size(); ++i)
//             points_inner.push_back(TV(0.6*cos(thetas[i]), 0.6*sin(thetas[i]), 0) * sim.unit);
//         for(int i=0; i<points_inner.size(); ++i)
//             addPoint(points_inner[i], full_dof_cnt, node_cnt);
//         std::vector<int> indices_inner = {0, 1, 2, 3, 4, 5, 6, 7};

//         std::vector<TV> points_center;
//         for(int i=0; i<thetas.size(); ++i)
//             points_center.push_back(TV(1.0*cos(thetas[i]), 1.0*sin(thetas[i]), 0) * sim.unit);
//         for(int i=0; i<points_center.size(); ++i)
//             addPoint(points_center[i], full_dof_cnt, node_cnt);
//         std::vector<int> indices_center = {8, 9, 10, 11, 12, 13, 14, 15};

//         std::vector<TV> points_outer;
//         for(int i=0; i<thetas.size(); ++i)
//             points_outer.push_back(TV(1.4*cos(thetas[i]), 1.4*sin(thetas[i]), 0) * sim.unit);
//         for(int i=0; i<points_outer.size(); ++i)
//             addPoint(points_outer[i], full_dof_cnt, node_cnt);
//         std::vector<int> indices_outer = {16, 17, 18, 19, 20, 21, 22, 23};
                

//         std::vector<TV2> data_points;
//         for(int i=0; i<points_inner.size(); ++i)
//             data_points.push_back(points_inner[i].template head<2>());
//         data_points.push_back(points_inner[0].template head<2>());


//         addCurvedRod(data_points, points_inner, indices_inner, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
//         data_points.clear();
//         for(int i=0; i<points_center.size(); ++i)
//             data_points.push_back(points_center[i].template head<2>());
//         data_points.push_back(points_center[0].template head<2>());
        
//         addCurvedRod(data_points, points_center, indices_center, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

//         data_points.clear();
//         for(int i=0; i<points_outer.size(); ++i)
//             data_points.push_back(points_outer[i].template head<2>());
//         data_points.push_back(points_outer[0].template head<2>());
        
//         addCurvedRod(data_points, points_outer, indices_outer, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
//         std::vector<TV> passing_points;
//         std::vector<int> passing_points_id;
//         passing_points.push_back(points_inner[theta_sub/2]);
//         passing_points.push_back(points_inner[0]); passing_points.push_back(points_center[0]); passing_points.push_back(points_outer[0]);

//         passing_points_id.push_back(indices_inner[theta_sub/2]);
//         passing_points_id.push_back(indices_inner[0]); passing_points_id.push_back(indices_center[0]); passing_points_id.push_back(indices_outer[0]);
//         TV to = TV(0, -3.0, 0.0) * sim.unit;
//         addAStraightRod(passing_points.front(), to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

//         int left = 0.75 * theta_sub;
//         int right = 0.25 * theta_sub;

//         passing_points = {points_outer[left], points_center[left], points_inner[left], points_inner[right], points_center[right], points_outer[right]};
//         passing_points_id = {indices_outer[left], indices_center[left], indices_inner[left], indices_inner[right], indices_center[right], indices_outer[right]};
//         addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

//         for (int i = 1; i < thetas.size(); i += 2)
//         {
//             passing_points = {points_outer[i], points_center[i], points_inner[i]};
//             passing_points_id = {indices_outer[i], indices_center[i], indices_inner[i]};
//             addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
//         }

        

//         // for (int i = 0; i < indices_center.size(); i++)
//         // {
//         //     sim.rod_crossings.push_back(new RodCrossing<T, dim>(i, {}));
//         // }
        
        
//         // data_points.clear();
//         // for(int i=0; i<points_center.size(); ++i)
//         //     data_points.push_back(points_outer[i].template head<2>());
//         // data_points.push_back(points_outer[0].template head<2>());
//         // std::vector<TV> nodal_positions;

//         // TV v0 = TV(0.5, 0, 0.0) * sim.unit;
//         // TV v1 = TV(0.0, 0.5, 0.0) * sim.unit;
//         // TV v2 = TV(0.5, 1, 0.0) * sim.unit;
//         // TV v3 = TV(1, 0.5, 0.0) * sim.unit;
//         // addPoint(v0, full_dof_cnt, node_cnt);
//         // addPoint(v1, full_dof_cnt, node_cnt);
//         // addPoint(v2, full_dof_cnt, node_cnt);
//         // addPoint(v3, full_dof_cnt, node_cnt);

//         // TV v4 = TV(0.5, 0.2, 0.0) * sim.unit;
//         // TV v5 = TV(0.2, 0.5, 0.0) * sim.unit;
//         // TV v6 = TV(0.5, 0.8, 0.0) * sim.unit;
//         // TV v7 = TV(0.8, 0.5, 0.0) * sim.unit;
//         // addPoint(v4, full_dof_cnt, node_cnt);
//         // addPoint(v5, full_dof_cnt, node_cnt);
//         // addPoint(v6, full_dof_cnt, node_cnt);
//         // addPoint(v7, full_dof_cnt, node_cnt);
        
//         // T theta = M_PI / 4.0;
//         // std::vector<TV2> data_points;

//         // data_points.push_back(TV2(0.5, 0) * sim.unit);
//         // data_points.push_back(TV2(0, 0.5) * sim.unit);
//         // data_points.push_back(TV2(0.5, 1) * sim.unit);
//         // data_points.push_back(TV2(1, 0.5) * sim.unit);
//         // data_points.push_back(TV2(0.5, 0) * sim.unit);

//         // std::vector<TV> passing_points = {v0, v1, v2, v3};
//         // std::vector<int> passing_points_id = {0, 1, 2, 3};

//         // addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);

//         // passing_points = {v4, v5, v6, v7};
//         // passing_points_id = {4, 5, 6, 7};
//         // data_points.clear();
//         // data_points.push_back(TV2(0.5, 0.2) * sim.unit);
//         // data_points.push_back(TV2(0.2, 0.5) * sim.unit);
//         // data_points.push_back(TV2(0.5, 0.8) * sim.unit);
//         // data_points.push_back(TV2(0.8, 0.5) * sim.unit);
//         // data_points.push_back(TV2(0.5, 0.2) * sim.unit);
//         // addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);
        
//         // passing_points = {v1, v3};
//         // passing_points_id = {1, 3};

//         // addAStraightRod(v1, v3, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

//         // passing_points = {v2, v0};
//         // passing_points_id = {2, 0};
//         // TV to = TV(0.5, -0.5, 0) * sim.unit;
//         // addAStraightRod(v2, to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

//         // RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(0, {0, 2});
//         // crossing->on_rod_idx[0] = sim.Rods[0]->dof_node_location[0];
//         // crossing->on_rod_idx[2] = sim.Rods[2]->dof_node_location[1];
//         // crossing->is_fixed = false;
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // crossing->sliding_ranges.push_back(Range(1, 1));
//         // sim.rod_crossings.push_back(crossing);

//         // crossing = new RodCrossing<T, dim>(1, {0, 1});
//         // crossing->on_rod_idx[0] = sim.Rods[0]->dof_node_location[1];
//         // crossing->on_rod_idx[1] = sim.Rods[1]->dof_node_location[0];
//         // crossing->is_fixed = true;
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // sim.rod_crossings.push_back(crossing);

//         // crossing = new RodCrossing<T, dim>(2, {0, 2});
//         // crossing->on_rod_idx[0] = sim.Rods[0]->dof_node_location[2];
//         // crossing->on_rod_idx[2] = sim.Rods[2]->dof_node_location[0];
//         // crossing->is_fixed = true;
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // sim.rod_crossings.push_back(crossing);

//         // crossing = new RodCrossing<T, dim>(3, {0, 1});
//         // crossing->on_rod_idx[0] = sim.Rods[0]->dof_node_location[3];
//         // crossing->on_rod_idx[1] = sim.Rods[1]->dof_node_location[1];
//         // crossing->is_fixed = true;
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // crossing->sliding_ranges.push_back(Range(0, 0));
//         // sim.rod_crossings.push_back(crossing);

//         for (auto& rod : sim.Rods)
//             rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

//         int dof_cnt = 0;
//         markCrossingDoF(w_entry, dof_cnt);
        
//         for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
//         appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
//         sim.rest_states = deformed_states;
        
//         sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
//         sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
//         for (auto& rod : sim.Rods)
//         {
//             rod->fixEndPointEulerian(sim.dirichlet_dof);
//             rod->setupBishopFrame();
//         }
        
    
//         // Offset offset;
//         // sim.Rods[2]->backOffsetReduced(offset);
//         // sim.dirichlet_dof[offset[0]] = 0;
//         // sim.dirichlet_dof[offset[1]] = -0.5 * sim.unit;
//         // sim.dirichlet_dof[offset[2]] = 0;

//         // T r = 0.1 * sim.unit;
//         // TV center1, center2;
//         // sim.getCrossingPosition(0, center1);

//         // auto circle1 = [r, center1](const TV& x)->bool
//         // {
//         //     return (x - center1).norm() < r;
//         // };

//         // sim.fixRegion(circle1);
//         // sim.boundary_spheres.push_back(std::make_pair(center1, r * 0.5));
//         // sim.fixCrossing();
//         // sim.perturb = VectorXT::Zero(sim.W.cols());
//         // for (auto& crossing : sim.rod_crossings)
//         // {
//         //     Offset off;
//         //     sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
//         //     T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
//         //     int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
//         //     sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
//         // }
//     }
// }

template<class T, int dim>
void UnitPatch<T, dim>::buildDenseInterlockingSquarePeriodicScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = true;
        sim.add_pbc_twisting = true;
        sim.add_pbc = true;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;
        sim.visual_R = 0.005;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;


        T square_width = 0.012; 
        T overlap = square_width * 0.25;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        auto addInnerSquare = [&](const TV& bottom_left)
        {
            TV v0, v1;

            // add bottom left corner
            addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
            v0 = bottom_left + TV(overlap, 0, 0);
            v1 = bottom_left + TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            // add bottom right corner
            TV bottom_right = bottom_left + TV(square_width, 0, 0);
            addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
            v0 = bottom_right - TV(overlap, 0, 0);
            v1 = bottom_right + TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            TV top_right = bottom_left + TV(square_width, square_width, 0);
            addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
            v0 = top_right - TV(overlap, 0, 0);
            v1 = top_right - TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            TV top_left = bottom_left + TV(0, square_width, 0);
            addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);
            v0 = top_left + TV(overlap, 0, 0);
            v1 = top_left - TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
        };

        TV x0 = TV(overlap, overlap, 0);
        TV x1 = TV(square_width - overlap, overlap, 0);
        TV x2 = TV(square_width - overlap, square_width - overlap, 0);
        TV x3 = TV(overlap, square_width - overlap, 0);


        addCrossingPoint(nodal_positions, x0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, x1, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, x2, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, x3, full_dof_cnt, node_cnt);
        
        TV bottom_left = TV::Zero();
        addInnerSquare(bottom_left);

        auto addRodFromIds = [&](std::vector<int> passing_points_ids)
        {
            std::vector<TV> passing_points;
            for (int idx : passing_points_ids)
                passing_points.push_back(nodal_positions[idx]);
            addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_ids, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = 0; i < passing_points_ids.size(); i++)
                addCrossingData(passing_points_ids[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);
            
        };
        
        auto addRodFrom2Ids = [&](int idx0, int idx1, TV extra, bool extra_node_first)
        {
            if (extra_node_first)
            {
                addAStraightRod(extra, nodal_positions[idx0], {nodal_positions[idx1], nodal_positions[idx0]}, {idx1, idx0}, sub_div, full_dof_cnt, node_cnt, rod_cnt);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[0]);
                addCrossingData(idx0, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            }
            else
            {
                addAStraightRod(nodal_positions[idx0], extra, {nodal_positions[idx0], nodal_positions[idx1]}, {idx0, idx1}, sub_div, full_dof_cnt, node_cnt, rod_cnt);
                addCrossingData(idx0, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[0]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            }
            
        };

        T distance = square_width * 0.5 - overlap;
        for (int corner = 0; corner < 4; corner++)
        {
            int x_shift = 4 + corner * 3 + 1;
            int y_shift = 4 + corner * 3 + 2;
            TV extra;
            if (corner == 0)
            {
                extra = nodal_positions[corner] + (nodal_positions[x_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, x_shift, extra, true);
                extra = nodal_positions[corner] + (nodal_positions[y_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, y_shift, extra, true);
            }
            else if (corner == 1)
            {
                extra = nodal_positions[corner] + (nodal_positions[x_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, x_shift, extra, true);
                extra = nodal_positions[corner] + (nodal_positions[y_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, y_shift, extra, false);
            }
            else if (corner == 2)
            {
                extra = nodal_positions[corner] + (nodal_positions[x_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, x_shift, extra, false);
                extra = nodal_positions[corner] + (nodal_positions[y_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, y_shift, extra, false);
            }
            else if (corner == 3)
            {
                extra = nodal_positions[corner] + (nodal_positions[x_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, x_shift, extra, false);
                extra = nodal_positions[corner] + (nodal_positions[y_shift] - nodal_positions[corner]).normalized() * square_width * 0.5;
                addRodFrom2Ids(corner, y_shift, extra, true);
            }
            
        }
        
        addRodFromIds({4, 5, 8, 7});
        addRodFromIds({7, 9, 12, 10});
        addRodFromIds({10, 11, 14, 13});
        addRodFromIds({13, 15, 6, 4});
        

        auto rod0 = sim.Rods[1];
        auto rod1 = sim.Rods[3];
        Offset end0, end1;
        rod0->frontOffset(end0); rod1->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(rod0->rod_id, rod1->rod_id));
        sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));

        Offset a, b;
        rod0->getEntryByLocation(1, a); rod1->getEntryByLocation(rod1->indices.size() - 2, b);
        sim.pbc_bending_pairs.push_back({end0, a, b, end1});
        sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});
        
        
        rod0 = sim.Rods[7];
        rod1 = sim.Rods[5];
        rod0->frontOffset(end0); rod1->backOffset(end1);
        sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        rod0->getEntryByLocation(1, a); rod1->getEntryByLocation(rod1->indices.size() - 2, b);
        sim.pbc_bending_pairs.push_back({end0, a, b, end1});
        sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});


        rod0 = sim.Rods[0];
        rod1 = sim.Rods[6];
        rod0->frontOffset(end0); rod1->backOffset(end1);
        sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(rod0->rod_id, rod1->rod_id));
        sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
        rod0->getEntryByLocation(1, a); rod1->getEntryByLocation(rod1->indices.size() - 2, b);
        sim.pbc_bending_pairs.push_back({end0, a, b, end1});
        sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});

        rod0 = sim.Rods[2];
        rod1 = sim.Rods[4];
        rod0->frontOffset(end0); rod1->backOffset(end1);
        sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
        rod0->getEntryByLocation(1, a); rod1->getEntryByLocation(rod1->indices.size() - 2, b);
        sim.pbc_bending_pairs.push_back({end0, a, b, end1});
        sim.pbc_bending_pairs_rod_id.push_back({rod0->rod_id, rod0->rod_id, rod1->rod_id, rod1->rod_id});

        enum FixingType
        {
            ALL, AllowX, OppositeDirection, Allow1XFix1X, FixOneCorner
        };
        T dx = (square_width * 0.5 - overlap) / square_width * 1.4;
        T dy = overlap / square_width * 0.92;

        FixingType type = AllowX;

        for (int corner = 0; corner < 4; corner++)
        {
            if (type == ALL)
                break;
            if (type == FixOneCorner && corner == 0)
            {
                continue;
            }
            int x_shift = 4 + corner * 3 + 1;
            int y_shift = 4 + corner * 3 + 2;

            auto crossing = sim.rod_crossings[x_shift];
            if (type == AllowX || type == OppositeDirection || type == Allow1XFix1X)
            {
                crossing->is_fixed = false;
                crossing->sliding_ranges[1] = Range(dx, dx);
            }
            crossing = sim.rod_crossings[y_shift];
            if (type == AllowX)
            {
                crossing->is_fixed = false;
                crossing->sliding_ranges[0] = Range(dx, dx);
            }
            else if (type == OppositeDirection || FixOneCorner)
            {
                crossing->is_fixed = false;
                crossing->sliding_ranges[1] = Range(dy, dy);
            }
        }
        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        // std::cout << "mark dof" << std::endl;
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
        Offset offset;
        sim.Rods[0]->frontOffsetReduced(offset);
        for (int d = 0; d < dim; d++) sim.dirichlet_dof[offset[d]] = 0;
        
        sim.fixCrossing();
        // std::cout << "fix crossing" << std::endl;
        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
        // std::cout << sim.rod_crossings[38]->is_fixed << std::endl;
  
        // GCodeGenerator<T, dim>(sim, "sliding_blocks.gcode").slidingBlocksGCode(n_row, n_col, 0);
    }
    // std::cout << "done" << std::endl;
}

template<class T, int dim>
void UnitPatch<T, dim>::buildDenseInterlockingSquareScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;
        sim.visual_R = 0.005;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;


        T square_width = 0.012; 
        T overlap = square_width * 0.3;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        auto addSquare = [&](const TV& bottom_left)
        {
            addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
            TV bottom_right = bottom_left + TV(square_width, 0, 0);
            addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
            TV top_right = bottom_left + TV(square_width, square_width, 0);
            addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
            TV top_left = bottom_left + TV(0, square_width, 0);
            addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);
        };


        auto addInnerSquare = [&](const TV& bottom_left)
        {
            TV v0, v1;

            // add bottom left corner
            addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
            v0 = bottom_left + TV(overlap, 0, 0);
            v1 = bottom_left + TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            // add bottom right corner
            TV bottom_right = bottom_left + TV(square_width, 0, 0);
            addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
            v0 = bottom_right - TV(overlap, 0, 0);
            v1 = bottom_right + TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            TV top_right = bottom_left + TV(square_width, square_width, 0);
            addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
            v0 = top_right - TV(overlap, 0, 0);
            v1 = top_right - TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
            TV top_left = bottom_left + TV(0, square_width, 0);
            addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);
            v0 = top_left + TV(overlap, 0, 0);
            v1 = top_left - TV(0, overlap, 0);
            addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
            addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
            
        };

        auto addRodsForASquare = [&](int bottom_left_node_idx)
        {
            TV bottom_left = nodal_positions[bottom_left_node_idx];
            TV bottom_right = nodal_positions[bottom_left_node_idx + 1];
            TV top_right = nodal_positions[bottom_left_node_idx + 2];
            TV top_left = nodal_positions[bottom_left_node_idx + 3];

            addAStraightRod(bottom_left, bottom_right, 
                {bottom_left, bottom_right}, {bottom_left_node_idx, bottom_left_node_idx + 1},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            
            addCrossingData(bottom_left_node_idx, rod_cnt - 1, 0);
            addCrossingData(bottom_left_node_idx + 1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            addAStraightRod(bottom_right, top_right, 
                {bottom_right, top_right}, {bottom_left_node_idx + 1, bottom_left_node_idx + 2},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(bottom_left_node_idx + 1, rod_cnt - 1, 0);
            addCrossingData(bottom_left_node_idx + 2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            
            addAStraightRod(top_right, top_left, 
                {top_right, top_left}, {bottom_left_node_idx + 2, bottom_left_node_idx + 3},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(bottom_left_node_idx + 2, rod_cnt - 1, 0);
            addCrossingData(bottom_left_node_idx + 3, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            addAStraightRod(top_left, bottom_left, 
                {top_left, bottom_left}, {bottom_left_node_idx + 3, bottom_left_node_idx},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(bottom_left_node_idx + 3, rod_cnt - 1, 0);
            addCrossingData(bottom_left_node_idx, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        };

        int n_row = 5, n_col = 5;
        
        T length_x = square_width * n_col + (square_width - 2.0 * overlap) * (n_col - 1);
        T length_y = length_x;


        // "longer" rows
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                T left = col * 2.0 * (square_width - overlap);
                T bottom = row * 2.0 * (square_width - overlap);
                TV bottom_left = TV(left, bottom, 0.0);
                addSquare(bottom_left);
            }    
        }

        for (int row = 0; row < n_row - 1; row++)
        {
            for (int col = 0; col < n_col - 1; col++)
            {
                T left = square_width - overlap + col * 2.0 * (square_width - overlap);
                T bottom = square_width - overlap + row * 2.0 * (square_width - overlap);
                TV bottom_left = TV(left, bottom, 0.0);
                // std::cout << bottom_left.transpose() << std::endl;
                addInnerSquare(bottom_left);
            }    
        }

        // std::cout << nodal_positions.size() << std::endl;

        // add Boundary rods first
        // these are the top row and bottom row

        for (int col = 0; col < n_col; col++)
        {
            int idx0 = (0 * n_col + col) * 4; // 4 is four nodes per square
            int idx1 = ((n_row - 1) * n_col + col) * 4; // 4 is four nodes per square
            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1 + 3];
            TV v0_next = nodal_positions[idx0 + 1], v1_next = nodal_positions[idx1 + 2]; 
            addAStraightRod(v0, v0_next, 
                {v0, v0_next}, {idx0, idx0 + 1},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0, rod_cnt - 1, 0);
            addCrossingData(idx0 + 1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        }

        for (int col = n_col - 1; col > -1; col--)
        {
            
            int idx1 = ((n_row - 1) * n_col + col) * 4; // 4 is four nodes per square
            TV v1 = nodal_positions[idx1 + 2];
            TV v1_next = nodal_positions[idx1 + 3]; 
            
            addAStraightRod(v1, v1_next, 
                {v1, v1_next}, {idx1 + 2, idx1 + 3},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx1 + 2, rod_cnt - 1, 0);
            addCrossingData(idx1 + 3, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        }
        
        for (int col = 0; col < n_col; col++)
        {
            auto rod0 = sim.Rods[col];
            auto rod1 = sim.Rods[n_col + n_col - col - 1];
            Offset end0, end1;
            // std::cout << rod0->indices.front() << " " << rod1->indices.back() << std::endl;
            rod0->frontOffset(end0); rod1->backOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
            if (col == 0)
            {
                sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(rod0->rod_id, rod1->rod_id));
            }
            rod0->backOffset(end0); rod1->frontOffset(end1);
            // std::cout << rod0->indices.back() << " " << rod1->indices.front() << std::endl;
            sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
        }
        

        // these are the left and right most rows
        for (int row = n_row - 1; row > -1; row--)
        {
            int idx0 = (row * n_col + 0) * 4; // 4 is four nodes per square
            int idx1 = (row * n_col + n_col - 1) * 4; // 4 is four nodes per square
            
            TV v0 = nodal_positions[idx0 + 3];
            TV v0_next = nodal_positions[idx0 + 0];

            addAStraightRod(v0, v0_next, 
                {v0, v0_next}, {idx0 + 3, idx0},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0 + 3, rod_cnt - 1, 0);
            addCrossingData(idx0, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        }

        
        for (int row = 0; row < n_row; row++)
        {
            
            int idx1 = (row * n_col + n_col - 1) * 4; // 4 is four nodes per square
            
            TV v1 = nodal_positions[idx1 + 1];
            TV v1_next = nodal_positions[idx1 + 2]; 

            addAStraightRod(v1, v1_next, 
                {v1, v1_next}, {idx1 + 1, idx1 + 2},
                sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx1 + 1, rod_cnt - 1, 0);
            addCrossingData(idx1 + 2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        }

        for (int row = 0; row < n_row; row++)
        {
            auto rod0 = sim.Rods[row + 2 * n_col];
            auto rod1 = sim.Rods[2 * n_col + 2 * n_row - row - 1];

            Offset end0, end1;
            rod0->frontOffset(end0); rod1->backOffset(end1);
            // std::cout << rod0->indices.front() << " " << rod1->indices.back() << std::endl;
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
            if (row == 0)
            {
                sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(rod0->rod_id, rod1->rod_id));
            }
            rod0->backOffset(end0); rod1->frontOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
            // std::cout << rod0->indices.back() << " " << rod1->indices.front() << std::endl;
        }
        
        for (int col = 0; col < n_col - 1; col++)
        {
            // vertical ones
            for (int row = 0; row < n_row - 1; row++)    
            {
                int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
                int idx1 = (row * n_col + col) * 4 + 2; 

                int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
                int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

                TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                TV v_middle1 = nodal_positions[idx_middle1];
                TV v_middle2;
                if (row == 0)
                {
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (row == n_row - 2)
                {
                    idx0 = ((row + 1) * n_col + col) * 4 + 1; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col) * 4 + 2; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

                }
                if (row != n_row - 2)
                {
                    idx0 = ((row + 1) * n_col + col) * 4 + 1; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col) * 4 + 2; 
                    
                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10;
                    idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                    
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
            }
        }

        // add inner rods
        for (int col = 0; col < n_col - 1; col++)
        {
            // vertical ones
            for (int row = 0; row < n_row - 1; row++)
            {
                //these are right column of each square
                int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
                int idx1 = (row * n_col + col) * 4 + 2; 

                int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
                int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

                TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                TV v_middle1 = nodal_positions[idx_middle1];
                TV v_middle2;
                if (row == 0)
                {
                    

                    idx0 = (row * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                    idx1 = (row * n_col + col + 1) * 4 + 3; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 4;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (row == n_row - 2)
                {
                    

                    idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col + 1) * 4 + 3; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

                }
                if (row != n_row - 2)
                {
                    

                    idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col + 1) * 4 + 3; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7;
                    idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 4;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
            }
        }

        
        

        // add inner rods
        for (int row = 0; row < n_row - 1; row++)
        {
            // vertical ones
            for (int col = 0; col < n_col - 1; col++)
            {
                //these are right column of each square
                int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
                int idx1 = (row * n_col + col) * 4 + 2; 

                int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
                int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

                TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                TV v_middle1 = nodal_positions[idx_middle1];
                TV v_middle2;
                

                if (col == 0)
                {
                    idx0 = (row * n_col + col) * 4 + 3; // 4 is four nodes per square
                    idx1 = (row * n_col + col) * 4 + 2; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 2;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (col == n_col - 2)
                {
                    idx0 = (row * n_col + col + 1) * 4 + 3; // 4 is four nodes per square
                    idx1 = (row * n_col + col + 1) * 4 + 2; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (n_col - 2 != col)
                {
                    idx0 = (row * n_col + col + 1) * 4 + 3; // 4 is four nodes per square
                    idx1 = (row * n_col + col + 1) * 4 + 2; 
                    
                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5;
                    idx_middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col + 1) * 12 + 2;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                    
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
            }
        }




        // add inner rods
        for (int row = 0; row < n_row - 1; row++)
        {
            // vertical ones
            for (int col = 0; col < n_col - 1; col++)
            {
                //these are right column of each square
                int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
                int idx1 = (row * n_col + col) * 4 + 2; 

                int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
                int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

                TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                TV v_middle1 = nodal_positions[idx_middle1];
                TV v_middle2;
                
                if (col == 0)
                {
                    

                    idx0 = ((row + 1) * n_col + col) * 4 + 0; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col) * 4 + 1; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 11;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (col == n_col - 2)
                {
                    
                    idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col + 1) * 4 + 1; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
                if (n_col - 2 != col)
                {
                    

                    idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                    idx1 = ((row + 1) * n_col + col + 1) * 4 + 1; 

                    idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8;
                    idx_middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col + 1) * 12 + 11;
                    v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                    v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                    addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                    addCrossingData(idx0, rod_cnt - 1, 0);
                    addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                    addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                    addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                }
            }
        }



        // std::cout << "inner" << std::endl;
        for (int row = 0; row < n_row - 1; row++)
        {
            for (int col = 0; col < n_col - 1; col++)
            {
                int idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
                int idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 3; 
                int middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1; 
                int middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 4; 

                
                TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                TV v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
                addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

                idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 3;
                idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 6; 
                middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5; 
                middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8; 

                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
                addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

                idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 6;
                idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 9; 
                middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7; 
                middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10; 

                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
                addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
                
                idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 9;
                idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 0; 
                middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 11; 
                middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 2; 

                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
                addAStraightRod(v0, v1, 
                        {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                        sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
        }

        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

        // for (int row = 0; row < n_row - 1; row++)
        // {
        //     for (int col = 0; col < n_col - 1; col++)
        //     {
        //         int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
        //         for (int corner = 0; corner < 4; corner++)
        //         {
                    
        //             auto crossing = sim.rod_crossings[base + corner * 3 + 1];
        //             crossing->is_fixed = false;
        //             crossing->sliding_ranges[1] = Range::Ones();
        //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[1] = false;
        //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[2] = false;

        //             sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[1] = false;

        //             crossing = sim.rod_crossings[base + corner * 3 + 2];
        //             crossing->is_fixed = false;
        //             crossing->sliding_ranges[0] = Range::Ones();
        //             sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[1] = false;

        //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[1] = false;
        //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[2] = false;

        //         }
        //     }
        // }

        

        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        // std::cout << "mark dof" << std::endl;
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
        if (sim.add_pbc)
            sim.Rods[0]->fixPointLagrangianByID(0, TV::Zero(), Mask::Ones(), sim.dirichlet_dof);

        TV bottom_left, top_right;
        sim.computeBoundingBox(bottom_left, top_right);

        TV shear_x_right = TV(0.1, 0.0, 0.0) * length_x;
        TV shear_x_left = TV(0.0, 0.0, 0) * sim.unit;


        T rec_width = 0.015 * sim.unit;

        auto rec1 = [bottom_left, top_right, shear_x_left, rec_width](
            const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = shear_x_left;
            if (x[0] < bottom_left[0] + rec_width 
                )
                return true;
            return false;
        };

        auto rec2 = [bottom_left, top_right, shear_x_right, rec_width](
            const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = shear_x_right;

            if (x[0] > top_right[0] - rec_width)
                return true;
            return false;
        };

        if (!sim.add_pbc)
        {
            sim.fixRegionalDisplacement(rec2);
            sim.fixRegionalDisplacement(rec1);
        }

        sim.fixCrossing();
        // std::cout << "fix crossing" << std::endl;
        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            
        }
        // std::cout << sim.rod_crossings[38]->is_fixed << std::endl;
  
        // GCodeGenerator<T, dim>(sim, "sliding_blocks_"+ std::to_string(n_row) + "x" + std::to_string(n_col) +".gcode").slidingBlocksGCode(n_row, n_col, 0);
        // GCodeGenerator<T, dim>(sim, "sliding_blocks_1X1Y_fuse_corner.gcode").slidingBlocksGCode(n_row, n_col, 3, true);
        // GCodeGenerator<T, dim>(sim, "sliding_blocks_1X1Y_bar.gcode").slidingBlocksGCode(n_row, n_col, 1, true);
        // GCodeGenerator<T, dim>(sim, "sliding_blocks_fused.gcode").slidingBlocksGCode(n_row, n_col, 2, true);
        // GCodeGenerator<T, dim>(sim, "sliding_blocks_allX_bar.gcode").slidingBlocksGCode(n_row, n_col, 0, true);
    }
    // std::cout << "done" << std::endl;
}

template<class T, int dim>
void UnitPatch<T, dim>::buildInterlockingSquareScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        auto addSquare = [&](const TV& bottom_left, T width)
        {
            
            sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
            addPoint(bottom_left, full_dof_cnt, node_cnt);
            nodal_positions.push_back(bottom_left);

            sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
            TV bottom_right = bottom_left + TV(width, 0, 0);
            addPoint(bottom_right, full_dof_cnt, node_cnt);
            nodal_positions.push_back(bottom_right);

            sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
            TV top_right = bottom_left + TV(width, width, 0);
            addPoint(top_right, full_dof_cnt, node_cnt);
            nodal_positions.push_back(top_right);

            sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
            TV top_left = bottom_left + TV(0, width, 0);
            addPoint(top_left, full_dof_cnt, node_cnt);
            nodal_positions.push_back(top_left);
        };

        TV s0 = TV(0, 0, 0) * sim.unit;
        addSquare(s0, 1.0 * sim.unit);

        TV s1 = TV(0.75, 0.75, 0) * sim.unit;
        addSquare(s1, 1.0 * sim.unit);

        TV crossing0 = TV(1.0, 0.75, 0.0) * sim.unit;
        sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
        addPoint(crossing0, full_dof_cnt, node_cnt);
        nodal_positions.push_back(crossing0);

        sim.rod_crossings.push_back(new RodCrossing<T, dim>(node_cnt, std::vector<int>()));
        TV crossing1 = TV(0.75, 1.0, 0.0) * sim.unit;
        addPoint(crossing1, full_dof_cnt, node_cnt);
        nodal_positions.push_back(crossing1);

        
        std::vector<TV> passing_points = {nodal_positions[0], nodal_positions[1]};
        std::vector<int> passing_points_id = {0, 1};

        addAStraightRod(passing_points.front(), passing_points.back(), 
            passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        
        passing_points = {nodal_positions[1], nodal_positions[8], nodal_positions[2]};
        passing_points_id = {1, 8, 2};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[2], nodal_positions[9], nodal_positions[3]};
        passing_points_id = {2, 9, 3};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[3], nodal_positions[0]};
        passing_points_id = {3, 0};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[4], nodal_positions[8], nodal_positions[5]};
        passing_points_id = {4, 8, 5};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[5], nodal_positions[6]};
        passing_points_id = {5, 6};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[6], nodal_positions[7]};
        passing_points_id = {6, 7};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        passing_points = {nodal_positions[7], nodal_positions[9], nodal_positions[4]};
        passing_points_id = {7, 9, 4};
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        for (int i = 0; i < passing_points_id.size(); i++) 
            addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

        // for (auto& crossing : sim.rod_crossings)
        // {
        //     std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
        //     crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
        // }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
        TV center0;
        T r = 0.05 * sim.unit;
        sim.Rods[0]->x(sim.Rods[0]->indices.front(), center0);
        auto circle0 = [r, center0](const TV& x)->bool
        {
            return (x - center0).norm() < r;
        };

        sim.fixRegion(circle0);

        Offset offset;
        sim.Rods[5]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0.3 * sim.unit;
        sim.dirichlet_dof[offset[1]] = 0;
        sim.dirichlet_dof[offset[2]] = 0;

        sim.Rods[5]->backOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0.3 * sim.unit;
        sim.dirichlet_dof[offset[1]] = 0;
        sim.dirichlet_dof[offset[2]] = 0;

        auto crossing = sim.rod_crossings[8];
        crossing->is_fixed = false;
        crossing->sliding_ranges[1] = Range::Ones();

        crossing = sim.rod_crossings[9];
        crossing->is_fixed = false;
        crossing->sliding_ranges[0] = Range::Ones();

        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            // sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
        sim.boundary_spheres.push_back(std::make_pair(center0, r));
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildActiveTextileScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        enum ActuationType
        {
            OneHorizontal, HorizontalFree, LongRec, BothFree, BothExtend
        };

        ActuationType type = BothFree;
        bool center = false;

        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        int n_row = 11, n_col = 11;

        // push crossings first 
        T dy = 1.0 / n_row * sim.unit;
        T dx = 1.0 / n_col * sim.unit;

        if (type == LongRec)
        {
            dy = 4.0 / n_row * sim.unit;
            n_row = 11; n_col = 11;
        }

        int half_row = (n_row - 1) / 2;
        int half_col = (n_col - 1) / 2;
        
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                TV x = TV(row * dx, col * dy, 0);
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                addPoint(x, full_dof_cnt, node_cnt);
                nodal_positions.push_back(x);
            }
        }

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        for (int row = 0; row < n_row; row++)
        {
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            int start = row * (n_col);
            int end = row * n_col + n_col;

            for (int i = start; i < end; i++)
            {
                passing_points.push_back(nodal_positions[i]);
                passing_points_id.push_back(i);
            }

            TV from, to;
            int sub_div_extend;

            TV extend_to = passing_points.back() + 
                (passing_points.back() - passing_points.front()) * 0.3;
            TV extend = passing_points.front() - 
                    (passing_points.back() - passing_points.front()) * 0.3;
            if (type == OneHorizontal || type == HorizontalFree || type == LongRec)
            {
                from = passing_points.front();
                sub_div_extend = sub_div;
                to = passing_points.back();
            }
            else if (type == BothExtend)
            {
                from = row == 0 || row == n_row - 1 ? passing_points.front() : extend;
                to = row == 0 || row == n_row - 1 ? passing_points.back() : extend_to;
                sub_div_extend = row == 0 || row == n_row - 1 ? sub_div : sub_div * 2;
            }
            else
            {
                
                from = row == 0 || row == n_row - 1 ? passing_points.front() : extend;
                if (center)
                    from = row == half_row  ? extend : passing_points.front();
                sub_div_extend = row == 0 || row == n_row - 1 ? sub_div : sub_div * 2;
                to = passing_points.back();
            }

            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
            for (int i = start; i < end; i++)
            {
                addCrossingData(i, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i - start]);
            }
        }
        
        for (int col = 0; col < n_col; col++)
        {
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            for (int i = 0; i < n_row; i++)
            {
                int middle = i * n_col + col;
                passing_points.push_back(nodal_positions[middle]);
                passing_points_id.push_back(middle);
            }

            TV extend = passing_points.back() + 
                (passing_points.back() - passing_points.front()) * 0.3;

            TV extend_from = passing_points.front() - 
                (passing_points.back() - passing_points.front()) * 0.3;

            TV to, from;
            int sub_div_extend;

            if (type == OneHorizontal || type == HorizontalFree || type == LongRec)
            {
                to = passing_points.back();
                from = passing_points.front();
                sub_div_extend = sub_div;
            }
            else if (type == BothFree)
            {
                to = col == 0 || col == n_col - 1 ? passing_points.back() : extend;
                if (center)
                    to = col == half_col ? extend : passing_points.back();
                from = passing_points.front();
                sub_div_extend = sub_div * 2;
            }
            else if (type == BothExtend)
            {
                to = col == 0 || col == n_col - 1 ? passing_points.back() : extend;
                from = col == 0 || col == n_col - 1 ? passing_points.front() : extend_from;
                sub_div_extend = col == 0 || col == n_col - 1 ? sub_div : sub_div * 2;
            }
            else
            {
                
            }

            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div_extend, full_dof_cnt, node_cnt, rod_cnt);
            
            int cnt = 0;
            for (int idx : passing_points_id)
            {
                addCrossingData(idx, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[cnt++]);
            }
            
        }
        
        

        T dv = 1.0 / (n_row - 1);
        T du = 1.0 / (n_col - 1);
        
        for (auto& crossing : sim.rod_crossings)
        {
            std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
            crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
        }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
        
        Offset offset;
        if (type == OneHorizontal)
        {
            for (int col = 1; col < n_col; col++)
            {    
                int node_idx = col * n_col + half_row;
                auto crossing = sim.rod_crossings[node_idx];
                crossing->is_fixed = false;
                // crossing->sliding_ranges[0] = Range(0.4, 0.4);

                // crossing->sliding_ranges[1] = Range(0.24 * du, 0.24 * du);
                crossing->sliding_ranges[1] = Range::Ones() * 1e6;
            }
        }
        else if (type == HorizontalFree)
        {
            for (int col = 1; col < n_col; col++)
            {
                for (int row = 1; row < n_row - 1; row ++)
                {
                    
                    int node_idx = col * n_col + row;
                    // std::cout << node_idx << std::endl;
                    auto crossing = sim.rod_crossings[node_idx];
                    crossing->is_fixed = false;
                    crossing->sliding_ranges[1] = Range::Ones() * 1e6;
                    
                    // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
                }    
            }
            
            // std::cout << node_idx << std::endl;
            auto crossing = sim.rod_crossings[1];
            crossing->is_fixed = false;
            crossing->sliding_ranges[1] = Range::Ones() * 1e6;
            crossing = sim.rod_crossings[3];
            crossing->is_fixed = false;
            crossing->sliding_ranges[1] = Range::Ones() * 1e6;

            sim.rod_crossings[21]->is_fixed = true;
            sim.rod_crossings[23]->is_fixed = true;
        }
        else if (type == LongRec)
        {
            for (int col = 1; col < n_col; col++)
            {
                for (int row = 1; row < n_row - 1; row ++)
                {
                    
                    int node_idx = col * n_col + row;
                    // std::cout << node_idx << std::endl;
                    auto crossing = sim.rod_crossings[node_idx];
                    crossing->is_fixed = false;
                    crossing->sliding_ranges[1] = Range::Ones() * 1e6;
                    // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
                }    
            }
        }
        else if (type == BothFree)
        {
            for (int col = 1; col < n_col; col++)
            {
                for (int row = 1; row < n_row - 1; row ++)
                {
                    
                    int node_idx = col * n_col + row;
                    // std::cout << node_idx << std::endl;
                    auto crossing = sim.rod_crossings[node_idx];
                    // if (row == half_row)
                        crossing->is_fixed = false;
                    crossing->sliding_ranges[1] = Range::Ones() * 1e6;
                    // std::cout << crossing->rods_involved[1] << std::endl;
                    // if (row == half_row)
                    sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[col] = false;
                    // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
                }    
            }
            // auto crossing = sim.rod_crossings[5];
            // crossing->is_fixed = false;
            // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
            // sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[0] = false;
            // crossing = sim.rod_crossings[10];
            // crossing->is_fixed = false;
            // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
            // sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[0] = false;
            // crossing = sim.rod_crossings[15];
            // crossing->is_fixed = false;
            // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
            // sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[0] = false;

            for (int col = 1; col < n_col - 1; col++)
            {
                auto crossing = sim.rod_crossings[col * n_row];

                crossing->is_fixed = false;
                crossing->sliding_ranges[0] = Range::Ones() * 1e6;
                sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[0] = false;
            }

            // sim.rod_crossings[6]->is_fixed = false;
            // sim.rod_crossings[7]->is_fixed = false;
            // sim.rod_crossings[8]->is_fixed = false;

            // sim.rod_crossings[11]->is_fixed = false;
            // sim.rod_crossings[12]->is_fixed = false;
            // sim.rod_crossings[13]->is_fixed = false;


        }
        else if (type == BothExtend)
        {
            for (int col = 1; col < n_col; col++)
            {
                for (int row = 1; row < n_row - 1; row ++)
                {
                    
                    int node_idx = col * n_col + row;
                    std::cout << node_idx << std::endl;
                    auto crossing = sim.rod_crossings[node_idx];
                    crossing->is_fixed = false;
                    crossing->sliding_ranges[1] = Range::Ones() * 1e6;
                    // crossing->sliding_ranges[0] = Range::Ones() * 1e6;
                }    
            }
            auto crossing = sim.rod_crossings[9];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range::Ones() * 1e6;
            crossing = sim.rod_crossings[10];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range::Ones() * 1e6;
            crossing = sim.rod_crossings[19];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range::Ones() * 1e6;
        }

        

        T r = 0.05 * sim.unit;
        TV center1, center2, center3, center4;
        sim.Rods[n_row-1]->x(sim.Rods[n_row-1]->indices.front(), center1);
        sim.Rods[n_row-1]->x(sim.Rods[n_row-1]->indices.back(), center2);

        sim.Rods[0]->x(sim.Rods[0]->indices.front(), center3);
        sim.Rods[0]->x(sim.Rods[0]->indices.back(), center4);

        // sim.Rods[n_row-1]->fixEndPointLagrangian(sim.dirichlet_dof);
        // sim.Rods[0]->fixPointLagrangianByID(0, TV::Zero(), Mask::Ones(), sim.dirichlet_dof);

        auto circle1 = [r, center1](const TV& x)->bool
        {
            return (x - center1).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        auto circle3 = [r, center3](const TV& x)->bool
        {
            return (x - center3).norm() < r;
        };

        sim.fixRegion(circle1);
        sim.fixRegion(circle2);
        sim.fixRegion(circle3);

        sim.boundary_spheres.push_back(std::make_pair(center1, r));
        sim.boundary_spheres.push_back(std::make_pair(center2, r));
        sim.boundary_spheres.push_back(std::make_pair(center3, r));

        if (type == BothFree)
        {
            int sub_type = 4;
            if(sub_type == 0)
            {
                sim.Rods[half_col]->frontOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0;
                sim.dirichlet_dof[offset[1]] = -0.1 * sim.unit;
                sim.dirichlet_dof[offset[2]] = 0;    

                sim.Rods[n_row + half_col]->backOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0.1 * sim.unit;
                sim.dirichlet_dof[offset[1]] = 0;
                sim.dirichlet_dof[offset[2]] = 0;
            }
            else if (sub_type == 1)
            {
                sim.Rods[n_col - 2]->frontOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0;
                sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
                sim.dirichlet_dof[offset[2]] = 0;    

                sim.Rods[n_row + 1]->backOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0.2 * sim.unit;
                sim.dirichlet_dof[offset[1]] = 0;
                sim.dirichlet_dof[offset[2]] = 0;
            }
            else if (sub_type == 2)
            {
                sim.Rods[1]->frontOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0;
                sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
                sim.dirichlet_dof[offset[2]] = 0;    

                sim.Rods[n_row + n_col - 2]->backOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0.2 * sim.unit;
                sim.dirichlet_dof[offset[1]] = 0;
                sim.dirichlet_dof[offset[2]] = 0;
            }
        }
        else if (type == BothExtend)
        {
            int sub_type = 0;
            if(sub_type == 0)
            {
                sim.Rods[half_col]->frontOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0;
                sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
                sim.dirichlet_dof[offset[2]] = 0;    

                sim.Rods[n_row + half_col]->backOffsetReduced(offset);
                sim.dirichlet_dof[offset[0]] = 0.2 * sim.unit;
                sim.dirichlet_dof[offset[1]] = 0;
                sim.dirichlet_dof[offset[2]] = 0;
            }
        }
        

        // sim.Rods[n_row + 1]->frontOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = -0.2 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0;

        // sim.Rods[n_row + n_col - 2]->frontOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = -0.2 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0;
        // sim.dirichlet_dof[offset[2]] = 0;
        

        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            // sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
        // VectorXT dq(sim.W.cols());
        // dq.setZero();
        // sim.checkHessianPD(dq);
        // GCodeGenerator<T, dim>(sim, "acitive_textile2_all_fused.gcode").activeTexticleGCode2(true);
        // TV target = TV(0.0133826, 0.0601781, 0.0324977);

        TV traj_center;
        sim.getCrossingPosition(half_col * n_row + half_row, traj_center);
        TV tip;
        sim.getCrossingPosition(n_row - 1, tip);
        T r_traj = (tip - traj_center).norm() * 0.8;
        T d_theta = M_PI / 2.0 / 5.0;
        TV up = traj_center + TV(0, 0, r_traj);
        for (T theta = d_theta; theta < M_PI / 2.0 + 1e-4; theta+=d_theta)
        {
            TV target;
            target[2] = r_traj * std::sin(theta);
            
            target.template segment<2>(0) = (traj_center + std::cos(theta) * (tip - traj_center)).template segment<2>(0);
            sim.targets.push_back(target);
        }
        
        // sim.targets.push_back(target);
        
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildSquareCrossJointScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<TV> nodal_positions;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                TV x = TV(0.5 * i, 0.5 * j, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }
        }

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                TV x = TV(0.25 + 0.5 * i, 0.25 + 0.5 * j, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }   
        }

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };
        
        std::vector<TV> passing_points;
        std::vector<int> passing_points_id;

        // passing_points = { nodal_positions[0], nodal_positions[1], nodal_positions[2] };
        // passing_points_id = { 0, 1, 2 };
        // addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        // addCrossingData(0, rod_cnt - 1, 0);
        // addCrossingData(1, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        // addCrossingData(2, rod_cnt - 1, sim.Rods.back()->numSeg());
        
        passing_points = { nodal_positions[3], nodal_positions[9], nodal_positions[1] };
        passing_points_id = { 3, 9, 1 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(3, rod_cnt - 1, 0);
        addCrossingData(9, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(1, rod_cnt - 1, sim.Rods.back()->numSeg());

        passing_points = { nodal_positions[1], nodal_positions[10], nodal_positions[5] };
        passing_points_id = { 1, 10, 5 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(1, rod_cnt - 1, 0);
        addCrossingData(10, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(5, rod_cnt - 1, sim.Rods.back()->numSeg());

        passing_points = { nodal_positions[3], nodal_positions[11], nodal_positions[7] };
        passing_points_id = { 3, 11, 7 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(3, rod_cnt - 1, 0);
        addCrossingData(11, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(7, rod_cnt - 1, sim.Rods.back()->numSeg());

        passing_points = { nodal_positions[7], nodal_positions[12], nodal_positions[5] };
        passing_points_id = { 7, 12, 5 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(7, rod_cnt - 1, 0);
        addCrossingData(12, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(5, rod_cnt - 1, sim.Rods.back()->numSeg());
        


        passing_points = { nodal_positions[0], nodal_positions[9], nodal_positions[4], nodal_positions[12], nodal_positions[8] };
        passing_points_id = { 0, 9, 4, 12, 8 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(0, rod_cnt - 1, 0);
        addCrossingData(9, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(4, rod_cnt - 1, sim.Rods.back()->dof_node_location[2]);
        addCrossingData(12, rod_cnt - 1, sim.Rods.back()->dof_node_location[3]);
        addCrossingData(8, rod_cnt - 1, sim.Rods.back()->numSeg());

        passing_points = { nodal_positions[6], nodal_positions[11], nodal_positions[4], nodal_positions[10], nodal_positions[2] };
        passing_points_id = { 6, 11, 4, 10, 2 };
        addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        addCrossingData(6, rod_cnt - 1, 0);
        addCrossingData(11, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        addCrossingData(4, rod_cnt - 1, sim.Rods.back()->dof_node_location[2]);
        addCrossingData(10, rod_cnt - 1, sim.Rods.back()->dof_node_location[3]);
        addCrossingData(2, rod_cnt - 1, sim.Rods.back()->numSeg());

        // passing_points = { nodal_positions[6], nodal_positions[7], nodal_positions[8] };
        // passing_points_id = { 6, 7, 8 };
        // addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        // addCrossingData(6, rod_cnt - 1, 0);
        // addCrossingData(7, rod_cnt - 1, sim.Rods.back()->dof_node_location[1]);
        // addCrossingData(8, rod_cnt - 1, sim.Rods.back()->numSeg());

        for (auto& crossing : sim.rod_crossings)
        {
            std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
            crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
        }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        // sim.Rods[0]->fixEndPointLagrangian(sim.dirichlet_dof);
        
        Offset offset;

        // sim.Rods[2]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0.3 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        // sim.Rods[0]->backOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;

        T r = 0.05 * sim.unit;

        TV delta0 = TV(0, 0, 0) * sim.unit;
        TV center0;
        sim.Rods[0]->x(sim.Rods[0]->indices.back(), center0);    
        auto circle0 = [r, center0, delta0](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta0;
            mask.setConstant(true);
            return (x - center0).norm() < r;
        };

        TV delta1 = TV(0.3, 0, 0) * sim.unit;
        TV center1;
        sim.Rods[2]->x(sim.Rods[2]->indices.back(), center1);    
        auto circle1 = [r, center1, delta1](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta1;
            mask.setConstant(true);
            return (x - center1).norm() < r;
        };


        sim.fixRegionalDisplacement(circle0);
        sim.fixRegionalDisplacement(circle1);

        
        // sim.Rods[5]->backOffsetReduced(offset);
        // // for (int d = 0; d < dim; d++) sim.dirichlet_dof[offset[d]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[0]] = -0.1 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;
        
        // sim.Rods[6]->frontOffsetReduced(offset);
        // // for (int d = 0; d < dim; d++) sim.dirichlet_dof[offset[d]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[0]] = -0.1 * sim.unit;
        // sim.dirichlet_dof[offset[1]] = 0 * sim.unit;
        // sim.dirichlet_dof[offset[2]] = 0 * sim.unit;
        

        for (int crossing_idx = 9; crossing_idx < 13; crossing_idx++)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = false;
            // sim.rod_crossings[crossing_idx]->sliding_ranges[1] = Range(0.2, 0.2);
            sim.rod_crossings[crossing_idx]->sliding_ranges[0] = Range(0.2, 0.2);
        }

        sim.rod_crossings[4]->is_fixed = false;
        sim.rod_crossings[4]->sliding_ranges[0] = Range(0.2, 0.2);
        // sim.rod_crossings[4]->sliding_ranges[1] = Range(0.2, 0.2);
        
        
        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildXJointsScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = true;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        int n_row = 6, n_col = 6;

        T dx = 1.0 / n_col;
        T dy = 1.0 / n_row;

        std::vector<int> middle_nodes;

        std::vector<TV> nodal_positions;
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                TV x = TV(col * dx, row * dy, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }
            if (row == n_row - 1)
                    continue;
            for (int col = 0; col < n_col - 1; col++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                middle_nodes.push_back(node_cnt);
                TV x = TV(col * dx + 0.5 * dx, row * dy + 0.5 * dy, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }
        }

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        for (int row = 0; row < n_row; row++)
        {
            int start = row * (2 * n_col - 1);
            int end = row * (2 * n_col - 1) + n_col;
            std::vector<TV> passing_points = std::vector<TV>(
                nodal_positions.begin() + row * (2 * n_col - 1),
                nodal_positions.begin() + row * (2 * n_col - 1) + n_col
            );
            
            std::vector<int> passing_points_id;
            for (int i = start; i < end; i++) passing_points_id.push_back(i);
            addAStraightRod(nodal_positions[start], nodal_positions[end - 1],
                passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
            );
            int rod_idx = rod_cnt - 1;
            
            for (int i = start; i < end; i++) 
                addCrossingData(i, rod_idx, sim.Rods[rod_idx]->dof_node_location[i - start]);
            
            // addCrossingData(start, rod_idx, 0);
            // addCrossingData(end-1, rod_idx, sim.Rods[rod_idx]->numSeg());
            
            
        }
        
        for (int row = 0; row < n_row - 1; row++)
        {
            for (int col = 0; col < n_col - 1; col++)
            {
                int start = row * (2 * n_col - 1) + col;
                int end = (row + 1) * (2 * n_col - 1) + col + 1;
                int middle = start + n_col;

                std::vector<TV> passing_points = { nodal_positions[start], 
                                                    nodal_positions[middle], 
                                                    nodal_positions[end]
                                                };                
                
                std::vector<int> passing_points_id = { start, middle, end };

                addAStraightRod(passing_points.front(), passing_points.back(),
                    passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
                );

                int rod_idx = rod_cnt - 1;
                // std::cout << "dia" << std::endl;
                addCrossingData(start, rod_idx, 0);
                addCrossingData(end, rod_idx, sim.Rods[rod_idx]->numSeg());
                addCrossingData(middle, rod_idx, sim.Rods[rod_idx]->dof_node_location[1]);
                // std::cout << "push" << std::endl;

                start = row * (2 * n_col - 1) + col + 1;
                end = (row + 1) * (2 * n_col - 1) + col;
                middle = start + n_col - 1;

                passing_points = { nodal_positions[start], 
                                    nodal_positions[middle], 
                                    nodal_positions[end]
                                };                
                
                passing_points_id = { start, middle, end };

                addAStraightRod(passing_points.front(), passing_points.back(),
                    passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
                );

                rod_idx = rod_cnt - 1;
                addCrossingData(start, rod_idx, 0);
                addCrossingData(end, rod_idx, sim.Rods[rod_idx]->numSeg());
                addCrossingData(middle, rod_idx, sim.Rods[rod_idx]->dof_node_location[1]);
            }
        }

        for (auto& crossing : sim.rod_crossings)
        {
            std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
            crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
        }

        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        // for (int crossing_idx : {3, 4, 8, 9})
        // {
        //     auto crossing = sim.rod_crossings[crossing_idx];
        //     crossing->is_fixed = false;
        //     crossing->sliding_ranges[0] = Range(0.4, 0.4);
        // }

        // for (int crossing_idx : {6, 7, 8, 9, 10})
        // {
        //     auto crossing = sim.rod_crossings[crossing_idx];
        //     crossing->is_fixed = false;
        //     crossing->sliding_ranges[0] = Range(0.4, 0.4);
        // }

        // for (int crossing_idx : {17, 18, 19, 20, 21})
        // {
        //     auto crossing = sim.rod_crossings[crossing_idx];
        //     crossing->is_fixed = false;
        //     crossing->sliding_ranges[0] = Range(0.4, 0.4);
        // }

        // for (int crossing_idx : {12, 13, 14, 15})
        // {
        //     auto crossing = sim.rod_crossings[crossing_idx];
        //     crossing->is_fixed = false;
        //     crossing->sliding_ranges[0] = Range(0.4, 0.4);
        // }
        for (int row = 0; row < n_row; row++)
        {
            auto rod = sim.Rods[row];
            Offset end0, end1;
            rod->frontOffset(end0); rod->backOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        }

        for (int col = 0; col < n_col; col++)
        {
            auto rod1 = sim.Rods[n_row - 1];
            auto rod0 = sim.Rods[0];
            Offset end0, end1;
            rod0->getEntryByLocation(col, end0);
            rod1->getEntryByLocation(col, end1);
            sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
            if (col == 0)
            {
                sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(0, n_row-1));
            }
        }


        for (int idx : middle_nodes)
        {
            auto crossing = sim.rod_crossings[idx];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range(1, 1);
        }

        
        sim.fixCrossing();

        // sim.Rods[0]->fixEndPointLagrangian(sim.dirichlet_dof);
        // sim.Rods[1]->fixEndPointLagrangian(sim.dirichlet_dof);
        
        Offset offset;

        

        Offset end0, end1;
        sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0, 0));

        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[0]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[1]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[2]]] = 0;

        // for (int d = 0; d < dim; d++)
        // {
        //     sim.dirichlet_dof[sim.rod_crossings[0]->reduced_dof_offset + d] = 0;
        // }
        

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildXJointsScene2(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = true;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        int n_row = 3, n_col = 3;

        T dx = 1.0 / n_col;
        T dy = 1.0 / n_row;

        std::vector<int> middle_nodes;

        std::vector<TV> nodal_positions;
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                TV x = TV(col * dx, row * dy, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }
            if (row == n_row - 1)
                    continue;
            for (int col = 0; col < n_col - 1; col++)
            {
                RodCrossing<T, dim>* crossing = new RodCrossing<T, dim>(node_cnt, std::vector<int>());
                sim.rod_crossings.push_back(crossing);
                middle_nodes.push_back(node_cnt);
                TV x = TV(col * dx + 0.5 * dx, row * dy + 0.5 * dy, 0.0) * sim.unit;
                nodal_positions.push_back(x);
                addPoint(x, full_dof_cnt, node_cnt);
            }
        }

        auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
        {
            sim.rod_crossings[crossing_idx]->is_fixed = true;
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
            sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
            sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
        };

        for (int row = 0; row < n_row; row++)
        {
            int start = row * (2 * n_col - 1);
            int end = row * (2 * n_col - 1) + n_col;
            std::vector<TV> passing_points = std::vector<TV>(
                nodal_positions.begin() + row * (2 * n_col - 1),
                nodal_positions.begin() + row * (2 * n_col - 1) + n_col
            );
            
            std::vector<int> passing_points_id;
            for (int i = start; i < end; i++) passing_points_id.push_back(i);
            addAStraightRod(nodal_positions[start], nodal_positions[end - 1],
                passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
            );
            int rod_idx = rod_cnt - 1;
            
            for (int i = start; i < end; i++) 
                addCrossingData(i, rod_idx, sim.Rods[rod_idx]->dof_node_location[i - start]);
            
        }
        
        // for (int row = 0; row < n_row - 1; row++)
        // {
        //     for (int col = 0; col < n_col - 1; col++)
        //     {
        //         int start = row * (2 * n_col - 1) + col;
        //         int end = (row + 1) * (2 * n_col - 1) + col + 1;
        //         int middle = start + n_col;

        //         std::vector<TV> passing_points = { nodal_positions[start], 
        //                                             nodal_positions[middle], 
        //                                             nodal_positions[end]
        //                                         };                
                
        //         std::vector<int> passing_points_id = { start, middle, end };

        //         addAStraightRod(passing_points.front(), passing_points.back(),
        //             passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
        //         );

        //         int rod_idx = rod_cnt - 1;
        //         // std::cout << "dia" << std::endl;
        //         addCrossingData(start, rod_idx, 0);
        //         addCrossingData(end, rod_idx, sim.Rods[rod_idx]->numSeg());
        //         addCrossingData(middle, rod_idx, sim.Rods[rod_idx]->dof_node_location[1]);
        //         // std::cout << "push" << std::endl;

        //         start = row * (2 * n_col - 1) + col + 1;
        //         end = (row + 1) * (2 * n_col - 1) + col;
        //         middle = start + n_col - 1;

        //         passing_points = { nodal_positions[start], 
        //                             nodal_positions[middle], 
        //                             nodal_positions[end]
        //                         };                
                
        //         passing_points_id = { start, middle, end };

        //         addAStraightRod(passing_points.front(), passing_points.back(),
        //             passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
        //         );

        //         rod_idx = rod_cnt - 1;
        //         addCrossingData(start, rod_idx, 0);
        //         addCrossingData(end, rod_idx, sim.Rods[rod_idx]->numSeg());
        //         addCrossingData(middle, rod_idx, sim.Rods[rod_idx]->dof_node_location[1]);
        //     }
        // }

        for (int col = 0; col < n_col - 2; col++)
        {
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            int start = col;
            passing_points_id.push_back(start);
            passing_points.push_back(nodal_positions[start]);
            for (int row = 0; row < n_row; row++)
            {
                int x0 = row * (2 * n_col - 1) + col;
                int middle1 = x0 + n_col;
                if (row != 0 && row != n_row - 1)
                {
                    passing_points_id.push_back(x0);
                    passing_points.push_back(nodal_positions[x0]);    
                }
                passing_points_id.push_back(middle1);
                passing_points.push_back(nodal_positions[middle1]);
            }
            int end = (n_row-1) * (2 * n_col - 1) + col + 1;
            passing_points_id.push_back(end);
            passing_points.push_back(nodal_positions[end]);
            addAStraightRod(passing_points.front(), passing_points.back(),
                    passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt
            );
            for (int i = 0; i < passing_points_id.size(); i++)
            {
                addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt-1]->dof_node_location[i]);
            }
            
        }
        

        for (auto& crossing : sim.rod_crossings)
        {
            std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
            crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
        }

        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        for (int row = 0; row < n_row; row++)
        {
            auto rod = sim.Rods[row];
            Offset end0, end1;
            rod->frontOffset(end0); rod->backOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        }

        for (int col = 0; col < n_col; col++)
        {
            auto rod1 = sim.Rods[n_row - 1];
            auto rod0 = sim.Rods[0];
            Offset end0, end1;
            rod0->getEntryByLocation(col, end0);
            rod1->getEntryByLocation(col, end1);
            sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
            if (col == 0)
            {
                sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), 
                    std::make_pair(0, n_row-1));
            }
        }


        for (int idx : middle_nodes)
        {
            auto crossing = sim.rod_crossings[idx];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range(1, 1);
        }

        
        sim.fixCrossing();

        
        Offset offset;

        

        Offset end0, end1;
        sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0, 0));

        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[0]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[1]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[2]]] = 0;

        // for (int d = 0; d < dim; d++)
        // {
        //     sim.dirichlet_dof[sim.rod_crossings[0]->reduced_dof_offset + d] = 0;
        // }
        

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildTestJoint(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-6;

        sim.unit = 0.09;

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        TV x0 = TV(0, 0, 0) * sim.unit;
        TV x1 = TV(1, 0, 0) * sim.unit;
        TV x2 = TV(0, 1, 0) * sim.unit;
        TV x3 = TV(1, 1, 0) * sim.unit;
        TV x4 = TV(0.5, 0.5, 0) * sim.unit;

        addPoint(x0, full_dof_cnt, node_cnt);
        addPoint(x1, full_dof_cnt, node_cnt);
        addPoint(x2, full_dof_cnt, node_cnt);
        addPoint(x3, full_dof_cnt, node_cnt);
        addPoint(x4, full_dof_cnt, node_cnt);

        std::vector<TV> passing_points = {x0, x1};
        std::vector<int> passing_points_id = {0, 1};

        addAStraightRod(x0, x1, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        
        passing_points = {x1, x4, x2};
        passing_points_id = {1, 4, 2};
        addAStraightRod(x1, x2, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        passing_points = {x2, x3};
        passing_points_id = {2, 3};
        addAStraightRod(x2, x3, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        passing_points = {x0, x4, x3};
        passing_points_id = {0, 4, 3};
        addAStraightRod(x0, x3, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);


        sim.rod_crossings.push_back(new RodCrossing<T, dim>(0, {0, 3}));
        sim.rod_crossings.back()->is_fixed = true;
        sim.rod_crossings.back()->on_rod_idx[0] = 0;
        sim.rod_crossings.back()->on_rod_idx[3] = 0;
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));


        sim.rod_crossings.push_back(new RodCrossing<T, dim>(1, {0, 1}));
        sim.rod_crossings.back()->is_fixed = true;
        sim.rod_crossings.back()->on_rod_idx[0] = sim.Rods[0]->numSeg();
        sim.rod_crossings.back()->on_rod_idx[1] = 0;
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));


        sim.rod_crossings.push_back(new RodCrossing<T, dim>(2, {1, 2}));
        sim.rod_crossings.back()->is_fixed = true;
        sim.rod_crossings.back()->on_rod_idx[1] = sim.Rods[0]->numSeg();
        sim.rod_crossings.back()->on_rod_idx[2] = 0;
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));

        sim.rod_crossings.push_back(new RodCrossing<T, dim>(3, {2, 3}));
        sim.rod_crossings.back()->is_fixed = true;
        sim.rod_crossings.back()->on_rod_idx[2] = sim.Rods[2]->numSeg();
        sim.rod_crossings.back()->on_rod_idx[3] = sim.Rods[3]->numSeg();
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));


        sim.rod_crossings.push_back(new RodCrossing<T, dim>(4, {1, 3}));
        sim.rod_crossings.back()->is_fixed = true;
        sim.rod_crossings.back()->on_rod_idx[1] = sim.Rods[1]->dof_node_location[1];
        sim.rod_crossings.back()->on_rod_idx[3] = sim.Rods[3]->dof_node_location[1];
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0.4, 0.4));
        sim.rod_crossings.back()->sliding_ranges.push_back(Range(0, 0));

        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        
        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        sim.fixCrossing();
        Offset offset;
        sim.Rods[0]->fixEndPointLagrangian(sim.dirichlet_dof);
        // sim.Rods[2]->fixEndPointLagrangian(sim.dirichlet_dof);

        sim.Rods[2]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = -0.1 * sim.unit;
        sim.dirichlet_dof[offset[1]] = -0.1 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0.0 * sim.unit;

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildGridLayoutGripper(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-2;

        sim.unit = 0.1;

        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<Entry> w_entry;

        // int n_row = 3, n_col = 4;
        int n_row = 5, n_col = 7;

        //n_row +=2 
        //n_col = 4, 7, 13, 25

        int half_row = std::floor(n_row / 2);
        int first_quator_col = n_col % 2 == 1 ? std::ceil(n_col / 4)  + 1 : std::ceil(n_col / 4);
        int second_quator_col = first_quator_col * 2;

        T dx = 1.0 / T(n_col);
        T dy = 1.0 / T(n_row);

        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                if (row > half_row && col > first_quator_col && col < second_quator_col)
                    continue;
                TV pt = TV(col * dx, row * dy, 0.0) * sim.unit;
                addPoint(pt, full_dof_cnt, node_cnt);
            }
        }

        int temp_cnt  = node_cnt;
        
        // for (int row = 0; row <= half_row; row++)
        // {
        //     TV pt = TV(0.5, row * dy, 0.0) * sim.unit;
        //     addPoint(pt, full_dof_cnt, node_cnt);
        // }

        auto node = [&](int idx)
        {
            return deformed_states.template segment<dim>(idx * dim);
        };

        std::unordered_map<int, int> offset_map;

        // add horizontal rods
        int offset_cnt = 0;
        
        for (int row = 0; row < n_row; row++)
        {
            std::vector<TV> passing_points;
            std::vector<int> passing_points_id;
            int row_offset_cnt = 0;
            for (int col = 0; col < n_col; col++)
            {
                if (row > half_row && col > first_quator_col && col < second_quator_col)
                {
                    offset_cnt ++;
                    offset_map[row * n_col + col] = offset_cnt;
                    row_offset_cnt++;
                    continue;
                }
                offset_map[row * n_col + col] = offset_cnt;
                passing_points.push_back(node(row * n_col + col - offset_cnt));
                passing_points_id.push_back(row * n_col + col - offset_cnt);
            }
            if (row <= half_row)
            {   
                addAStraightRod(passing_points.front(), passing_points.back(), 
                    passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
                
            }
            else
            {
                
                addAStraightRod(passing_points.front(), passing_points[first_quator_col], 
                    std::vector<TV>(passing_points.begin(), passing_points.begin() + first_quator_col + 1), 
                    std::vector<int>(passing_points_id.begin(), passing_points_id.begin() + first_quator_col + 1), 
                    (int)std::floor(sub_div / 3), full_dof_cnt, node_cnt, rod_cnt);

                
                int start = second_quator_col - row_offset_cnt;
                addAStraightRod(passing_points[start], passing_points.back(), 
                    std::vector<TV>(passing_points.begin() + start, passing_points.end()), 
                    std::vector<int>(passing_points_id.begin() + start, passing_points_id.end()), 
                    (int)std::floor(sub_div / 3), full_dof_cnt, node_cnt, rod_cnt);
            }
        }
        int n_h_rod = rod_cnt;
        
        
        // add vertical rods
        for (int col = 0; col < n_col; col++)
        {
            std::vector<TV> passing_points;
            std::vector<int> passing_points_id;

            for (int row = 0; row < n_row; row++)
            {
                offset_cnt = offset_map[row * n_col + col];
                if (row > half_row && col > first_quator_col && col < second_quator_col)
                {
                    continue;
                }
                passing_points.push_back(node(row * n_col + col - offset_cnt));
                passing_points_id.push_back(row * n_col + col - offset_cnt);
            }
            
            if ((col <= first_quator_col || col >= second_quator_col))
            {
                addAStraightRod(passing_points.front(), passing_points.back(), 
                    passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
            }
            else
            {
                TV from = passing_points.front() - TV(0.0, 0.5, 0) * sim.unit;
                // addAStraightRod(passing_points.front(), passing_points.back(), 
                //     passing_points, passing_points_id, sub_div / 2, full_dof_cnt, node_cnt, rod_cnt);
                addAStraightRod(from, passing_points.back(), 
                    passing_points, passing_points_id, sub_div / 2, full_dof_cnt, node_cnt, rod_cnt);
            }

        }
        
        offset_cnt = 0;

        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), true);

        // std::vector<int> hard_coded_sliding_node = {3, 9, 10, 11, 8, 15, 12, 19};
        
        int row_offset_cnt = 0;
        std::unordered_map<int, int> crossing_id_map;
    
        for (int row = 0; row < n_row; row++)
        {
            int col_offset_cnt = 0;
            if (row > half_row)
                row_offset_cnt++;
            for (int col = 0; col < n_col; col++)
            {
                
                if (col > first_quator_col && col < second_quator_col && row > half_row)
                    col_offset_cnt++;
                if (row > half_row && col > first_quator_col && col < second_quator_col)
                {
                    offset_cnt++;
                    continue;
                }
                
                int node_id = row * n_col + col - offset_cnt;
                int h_rod;
                if (row <= half_row)
                    h_rod = row;
                else if (col <= first_quator_col && row > half_row)
                    h_rod = half_row + row_offset_cnt;
                else if (col >= second_quator_col && row > half_row)
                {
                    if (col == second_quator_col)
                        row_offset_cnt++;
                    h_rod = half_row + row_offset_cnt;
                }

                int v_rod = col;

                std::vector<int> rods_involved = {h_rod, n_h_rod + v_rod};
                // std::cout << "node " << node_id << " hrod " << h_rod << " v_rod " << n_h_rod + v_rod << " col_offset_cnt " << col_offset_cnt << std::endl;
                // std::cout << "node " << node_id << " hrod " << h_rod << " v_rod " << v_rod << " col_offset_cnt " << col_offset_cnt << std::endl;
                RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(node_id, rods_involved);

                // if (row < half_row && col > first_quator_col && col < second_quator_col)
                // if (std::find(hard_coded_sliding_node.begin(), hard_coded_sliding_node.end(), node_id) != hard_coded_sliding_node.end())
                // {
                //     crossing->is_fixed = false;
                //     crossing->sliding_ranges.push_back(Range(0.1, 0.1));
                //     crossing->sliding_ranges.push_back(Range(0.2, 0.2));
                // }
                // else

                {
                    crossing->is_fixed = true;
                    crossing->sliding_ranges.push_back(Range(0, 0));
                    crossing->sliding_ranges.push_back(Range(0, 0));
                }

                if (row <= half_row)
                    crossing->on_rod_idx[h_rod] = sim.Rods[h_rod]->dof_node_location[v_rod];
                else if (row > half_row && col <= first_quator_col)
                    crossing->on_rod_idx[h_rod] = sim.Rods[h_rod]->dof_node_location[v_rod];
                else if (row > half_row && col >= second_quator_col)
                    crossing->on_rod_idx[h_rod] = sim.Rods[h_rod]->dof_node_location[v_rod - second_quator_col];

                if (sim.Rods[n_h_rod + v_rod]->dof_node_location.size() > row)
                    crossing->on_rod_idx[n_h_rod + v_rod] = sim.Rods[n_h_rod + v_rod]->dof_node_location[row];
                
                // std::cout << h_rod << " "  << v_rod << " " << v_rod - second_quator_col << std::endl;
                // std::cout << n_h_rod + v_rod << " "  << row << std::endl;
                sim.rod_crossings.push_back(crossing);
                crossing_id_map[node_id] = sim.rod_crossings.size() - 1;
                
            }
        }


        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        

        for (auto& crossing : sim.rod_crossings)
            if(!crossing->is_fixed)
                std::cout << "free" << std::endl;

        Offset offset;

        sim.Rods[10]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0;

        // for (int i = 3; i < 4; i++)
        // {
        //     sim.Rods[0]->getEntry(i, offset);
        //     for (int d = 0; d < dim; d++)
        //         sim.dirichlet_dof[sim.Rods[0]->reduced_map[offset[d]]] = 0;
        //     int loc = sim.Rods[0]->dof_node_location[3];
        //     sim.Rods[0]->getEntryByLocation(loc + 1, offset);
        //     for (int d = 0; d < dim; d++)
        //         sim.dirichlet_dof[sim.Rods[0]->reduced_map[offset[d]]] = 0;
        //     sim.Rods[0]->getEntryByLocation(loc - 1, offset);
        //     for (int d = 0; d < dim; d++)
        //         sim.dirichlet_dof[sim.Rods[0]->reduced_map[offset[d]]] = 0;

        //     sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset+loc] = 0;
        //     sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset+loc-1] = 0;

            
        // }   

        sim.Rods[0]->getEntry(3, offset);
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[offset[2]]] = 0;

        sim.Rods[0]->fixEndPointLagrangian(sim.dirichlet_dof);
        

        RodCrossing<T, dim>* crossing = sim.rod_crossings[crossing_id_map[3]];
        crossing->is_fixed = false;
        // crossing->sliding_ranges[0] = Range(0.2, 0.2);
        crossing->sliding_ranges[0] = Range(0, 0);
        crossing->sliding_ranges[1] = Range(0.25, 0.25);
        // crossing->sliding_ranges[1] = Range(1.0, 1.0) * 1e-4;

        crossing = sim.rod_crossings[crossing_id_map[10]];
        crossing->is_fixed = false;
        // crossing->sliding_ranges[0] = Range(0.2, 0.2);
        crossing->sliding_ranges[0] = Range(0, 0);
        crossing->sliding_ranges[1] = Range(0.15, 0.15);
        // crossing->sliding_ranges[1] = Range(1.0, 1.0) * 1e-4;

        for (int idx : {15, 19, 22, 25, 8, 12, 11, 9})
        {
            crossing = sim.rod_crossings[crossing_id_map[idx]];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range(0.0, 0.0);
            crossing->sliding_ranges[1] = Range(0.0, 0.1);

            // crossing->sliding_ranges[1] = Range(1.0, 1.0) * 1e-4;
        }

        for (int idx : {16, 18})
        {
            crossing = sim.rod_crossings[crossing_id_map[idx]];
            crossing->is_fixed = false;
            crossing->sliding_ranges[0] = Range(0.0, 0.0);
            // crossing->sliding_ranges[1] = Range(0.0, 0.2); // this is the correct value
            crossing->sliding_ranges[1] = Range(0.0, 0.1);
            // crossing->sliding_ranges[1] = Range(1.0, 1.0) * 1e-4;
        }

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            // dq[z_off] += 0.001 * r * unit;
            
        }
        
        sim.fixCrossing();

        T r = 0.1 * sim.unit;
        TV center1, center2;
        sim.Rods.front()->x(sim.Rods[0]->indices.front(), center1);
        sim.Rods.front()->x(sim.Rods[0]->indices.back(), center2);


        auto circle1 = [r, center1](const TV& x)->bool
        {
            return (x - center1).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        sim.fixRegion(circle1);
        sim.fixRegion(circle2);

        // GCodeGenerator<T, dim>(this->sim, "gripper_straight.gcode").generateGCodeFromRodsGridGripperHardCoded();
        // GCodeGenerator<T, dim>(this->sim, "gripper_straight_fused.gcode").generateGCodeFromRodsFixedGridGripperHardCoded();

        // GCodeGenerator<T, dim>(sim, "test_multi_layer.gcode").crossingTest();
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildGripperScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;
        sim.unit = 0.05;

        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<Entry> w_entry;

        TV v0 = TV(0.5, 0.0, 0.0) * sim.unit;
        TV v1 = TV(0.5, -0.5, 0.0) * sim.unit;
        addPoint(v0, full_dof_cnt, node_cnt);
        addPoint(v1, full_dof_cnt, node_cnt);
        
        T theta = M_PI / 4.0;
        std::vector<TV2> data_points;
        data_points.push_back(TV2(0.5, -0.5) * sim.unit);
        data_points.push_back(TV2(0, 0) * sim.unit);
        data_points.push_back(TV2(0.5 - 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        data_points.push_back(TV2(0.5, 0.0) * sim.unit);
        data_points.push_back(TV2(0.5 + 0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta)) * sim.unit);
        data_points.push_back(TV2(1.0, 0.0) * sim.unit);
        data_points.push_back(TV2(0.5, -0.5) * sim.unit);        

        std::vector<TV> passing_points = {v0, v1};
        std::vector<int> passing_points_id = {0, 1};

        addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, true);


        TV to = TV(0.5, -1, 0.0) * sim.unit;
        addAStraightRod(v0, to, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);

        for (int i = 0; i < 2; i++)
        {
            RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(i, {0, 1});
            crossing->on_rod_idx[0] = sim.Rods[0]->dof_node_location[1 - i];
            crossing->on_rod_idx[1] = sim.Rods[1]->dof_node_location[i];
            if (i == 0)
            // if(true)
            {
                crossing->is_fixed = true;
                crossing->sliding_ranges.push_back(Range(0, 0));
                crossing->sliding_ranges.push_back(Range(0, 0));
            }
            else
            {
                crossing->is_fixed = false;
                crossing->sliding_ranges.push_back(Range(0, 0));
                crossing->sliding_ranges.push_back(Range(0.4, 0.4));
            }
            sim.rod_crossings.push_back(crossing);
        }

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        Offset offset;
        sim.Rods[0]->getEntry(1, offset);
        for (int d = 0; d < dim; d++)
        {        
            sim.dirichlet_dof[sim.Rods[0]->reduced_map[offset[d]]] = 0;
        }

        sim.fixCrossing();
        sim.Rods[0]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[0]->fixPointLagrangian(sim.Rods[0]->indices.size() - 2, TV::Zero(), sim.dirichlet_dof);
        
        sim.Rods[1]->backOffsetReduced(offset);

        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -0.2 * sim.unit;
        sim.dirichlet_dof[offset[2]] = 0;

        sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset] = 0;
        sim.dirichlet_dof[sim.Rods[0]->theta_reduced_dof_start_offset + sim.Rods[0]->numSeg()-1] = 0;

        // GCodeGenerator<T, dim>(this->sim, "gripper.gcode").generateGCodeFromRodsCurveGripperHardCoded();
        

        sim.perturb = VectorXT::Zero(sim.W.cols());

        for (auto& crossing : sim.rod_crossings)
        {
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
        }
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildFingerScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;
        
        sim.ke = 1e-4;

        clearSimData();
        
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;

        std::vector<Entry> w_entry;

        std::vector<std::vector<TV>> passing_points(6, std::vector<TV>());
        std::vector<std::vector<int>> passing_points_id(6, std::vector<int>());

        for (int k = 0; k < 2; k++)
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    TV v0 = TV(0.25 * (j+1), 1.0 - 0.25 * i, -0.1 * T(k)) * sim.unit;
                    addPoint(v0, full_dof_cnt, node_cnt);
                    passing_points_id[3*k + j].push_back( 12 * k + i * 3 + j);
                    passing_points[3*k + j].push_back(v0);
                }
            }
        }
        
        
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 4; i++)
            {
                TV from = deformed_states.template segment<dim>((12 * j + i * 3) * dim);
                TV middle = deformed_states.template segment<dim>((12 * j + i * 3 + 1) * dim);
                TV to = deformed_states.template segment<dim>((12 * j + i * 3 + 2) * dim);
                std::vector<TV> points = { from, middle, to };
                std::vector<int> ids = { 12 * j + i * 3, 12 * j + i * 3 + 1, 12 * j + i * 3 + 2 };
                addAStraightRod(from, to, points, ids, 
                    sub_div, full_dof_cnt, node_cnt, rod_cnt);            
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                TV from = deformed_states.template segment<dim>((12 + i * 3 + j) * dim);
                TV to = deformed_states.template segment<dim>((i * 3 + j) * dim);
                std::vector<TV> points = { from, to };
                std::vector<int> ids = { 12 + i * 3 + j , i * 3 + j};
                if (j == 1)
                    continue;
                addAStraightRod(from, to, points, ids, 
                    sub_div, full_dof_cnt, node_cnt, rod_cnt);            
            }
        }

        // addAStraightRod(passing_points[0].front(), passing_points[0].back(), passing_points[0], passing_points_id[0], 
        //     sub_div * 2, full_dof_cnt, node_cnt, rod_cnt, false);

        TV from = passing_points[1].front();
        TV to = TV(0.5, 0, 0) * sim.unit;
        addAStraightRod(from, to, passing_points[1], passing_points_id[1], 
            sub_div * 2, full_dof_cnt, node_cnt, rod_cnt);
        
        for (int i = 3; i < 6; i++)
        {
            addAStraightRod(passing_points[i].front(), passing_points[i].back(), passing_points[i], passing_points_id[i], 
                sub_div * 2, full_dof_cnt, node_cnt, rod_cnt);
        }

        for (int k = 0; k < 2; k++)
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int rod0 = 4 * k + i; // x direction
                    int rod1 = j == 0 ? 8 + i * 2 : 8 + i * 2 + 1; // z direction
                    int rod2 = k == 0 ? 16 : 17 + j; // y direction
                    
                    std::vector<int> rods_involved;
                    if (k == 0 && j != 1)
                        rods_involved = {rod0, rod1};
                    else if (j == 1)
                        rods_involved = {rod0, rod2};
                    else
                        rods_involved = {rod0, rod1, rod2};

                    RodCrossing<T, dim>* crossing = 
                            new RodCrossing<T, dim>(k * 12 + i * 3 + j, rods_involved);
                                
                    crossing->undeformed_twist.push_back(Vector<T, 2>(0, 0)); 
                    if (j != 1)
                        crossing->undeformed_twist.push_back(Vector<T, 2>(0, 0)); 
                    if (k==1 || (k ==0 && j == 1))
                        crossing->undeformed_twist.push_back(Vector<T, 2>(0, 0)); 
                    
                    crossing->on_rod_idx[rod0] = sim.Rods[rod0]->dof_node_location[j];
                    
                    if (j != 1)
                        crossing->on_rod_idx[rod1] = sim.Rods[rod1]->dof_node_location[1] - sim.Rods[rod1]->dof_node_location[k];
                    
                    if (k==1 || (k ==0 && j == 1))
                        crossing->on_rod_idx[rod2] = sim.Rods[rod2]->dof_node_location[i];
                    
                    
                    if (i > 0 && j==1 && k == 0)
                    // if (i == 3 && j==1 && k == 0)
                    // if(false)
                    {
                        crossing->is_fixed = false;
                        crossing->sliding_ranges.push_back(Range(0, 0));
                        crossing->sliding_ranges.push_back(Range(0.2, 0.2));
                        // std::cout << "free" << std::endl;
                    }
                    else
                    {
                        crossing->is_fixed = true;
                        crossing->sliding_ranges.push_back(Range(0, 0));
                        if (j != 1)
                            crossing->sliding_ranges.push_back(Range(0, 0));
                        if (k==1 || (k ==0 && j == 1))
                            crossing->sliding_ranges.push_back(Range(0, 0));
                    }
                    
                    sim.rod_crossings.push_back(crossing);
                }
                
            }
        }

    
        
        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }

        Offset offset;
        sim.Rods[16]->getEntry(10, offset);
        for (int d = 0; d < dim; d++)
        {        
            sim.dirichlet_dof[sim.Rods[16]->reduced_map[offset[d]]] = 0;
        }

    
        sim.Rods[16]->backOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -sim.unit * 0.1;
        sim.dirichlet_dof[offset[2]] = 0;

        sim.Rods[3]->fixEndPointLagrangian(sim.dirichlet_dof);
        
        // for (int i = 17; i < 20; i++)
        // {
        //     sim.Rods[i]->backOffsetReduced(offset);
            
        //     sim.dirichlet_dof[offset[0]] = 0;
        //     sim.dirichlet_dof[offset[1]] = 0;
        //     sim.dirichlet_dof[offset[2]] = 0;
        // }
        // for (int i = 0; i < sim.Rods[16]->numSeg(); i++)
        // {
        //     sim.dirichlet_dof[sim.Rods[16]->theta_reduced_dof_start_offset + i] = 0;
        // }

        
        sim.fixCrossing();
    }
}


// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'



template<class T, int dim>
void UnitPatch<T, dim>::cropTranslationalUnitByparallelogram(const std::vector<std::vector<TV2>>& input_points,
    std::vector<TV2>& output_points, const TV2& top_left, const TV2& top_right,
    const TV2& bottom_right, const TV2& bottom_left, std::vector<Vector<int, 2>>& edge_pairs,
    std::unordered_map<int, std::vector<int>>& crossing_tracker,
    std::vector<std::vector<Vector<int, 2>>>& boundary_pairs,
    std::vector<std::vector<int>>& boundary_pair_rod_idx)
{
    if constexpr (dim == 3)
    {
        crossing_tracker.clear();

        std::vector<TV2> parallogram = {top_left, top_right, bottom_right, bottom_left};

        using Edge = std::pair<TV2, TV2>;

        std::vector<Edge> edges;

        boundary_pairs.resize(4, std::vector<Vector<int, 2>>());
        boundary_pair_rod_idx.resize(4, std::vector<int>());
        
        for (auto one_tile : input_points)
        {
            // -2 because the tile vertices loop back to the first one
            for (int i = 0; i < one_tile.size() - 2; i++)
            {
                const TV2 xi = one_tile[i];
                const TV2 xj = one_tile[i+ 1];
                
                bool xi_inside = insidePolygon(parallogram, xi);
                bool xj_inside = insidePolygon(parallogram, xj);

                // both points are inside the parallelogram
                if (xi_inside && xj_inside)
                {
                    
                    Edge xij = std::make_pair(xi, xj);

                    auto find_edge_iter = std::find_if(edges.begin(), edges.end(), [&xij](Edge e)
                        {   
                            return (((e.first - xij.first).norm() < 1e-6) && ((e.second - xij.second).norm() < 1e-6)) || 
                                (((e.first - xij.second).norm() < 1e-6) && ((e.second - xij.first).norm() < 1e-6));
                        }
                    );

                    bool new_edge = find_edge_iter == edges.end();
                    
                    if(new_edge)
                    {
                        edges.push_back(std::make_pair(xi, xj));
                        // std::cout << xi.transpose() << " " << xj.transpose() << std::endl;
                        auto find_xi_iter = std::find_if(output_points.begin(), output_points.end(), 
                            [&xi](const TV2 x)->bool
                               { return (x - xi).norm() < 1e-6; }
                             );
                        int xi_idx = -1, xj_idx = -1;
                    
                        if (find_xi_iter == output_points.end())
                        {
                            // xi is a new vtx
                            output_points.push_back(xi);
                            xi_idx = int(output_points.size()) - 1;
                            //pre push this edge
                            crossing_tracker[xi_idx] = {int(edge_pairs.size())};
                        }
                        else
                        {
                            int index = std::distance(output_points.begin(), find_xi_iter);
                            if (crossing_tracker.find(index) == crossing_tracker.end())
                            {
                                crossing_tracker[index] = {int(edge_pairs.size())};
                            }
                            else
                            {
                                crossing_tracker[index].push_back(int(edge_pairs.size()));
                            }
                            xi_idx = index;
                        }

                        auto find_xj_iter = std::find_if(output_points.begin(), output_points.end(), 
                                [&xj](const TV2 x)->bool
                               { return (x - xj).norm() < 1e-6; }
                        );
                        if (find_xj_iter == output_points.end())
                        {
                            output_points.push_back(xj);
                            xj_idx = int(output_points.size()) - 1;
                            crossing_tracker[xj_idx] = {int(edge_pairs.size())};
                        }
                        else
                        {
                            int index = std::distance(output_points.begin(), find_xj_iter);
                            if (crossing_tracker.find(index) == crossing_tracker.end())
                            {
                                crossing_tracker[index] = {int(edge_pairs.size())};
                            }
                            else
                            {
                                crossing_tracker[index].push_back(int(edge_pairs.size()));
                            }
                            xj_idx = index;
                        }
                        
                        edge_pairs.push_back(Vector<int,2>(xi_idx, xj_idx));
                    }
                }
                else if(!xi_inside && xj_inside)
                {
                    
                    // std::cout << "One is inside" << std::endl;
                    Edge xij = std::make_pair(xi, xj);
                    auto find_edge_iter = std::find_if(edges.begin(), edges.end(), [&xij](Edge e)
                        {   
                            return (((e.first - xij.first).norm() < 1e-6) && ((e.second - xij.second).norm() < 1e-6)) || 
                                (((e.first - xij.second).norm() < 1e-6) && ((e.second - xij.first).norm() < 1e-6));
                        }
                    );

                    bool new_edge = find_edge_iter == edges.end();

                    if(new_edge)
                    {
                        edges.push_back(std::make_pair(xi, xj));
                        TV2 intersection;
                        int xj_idx = -1;
                        bool intersected = false;
                        int intersecting_edge = -1;
                        if (lineSegementsIntersect2D(xi, xj, parallogram[0], parallogram[1], intersection))
                        {
                            intersected = true;
                            intersecting_edge = 0;
                        }
                        else if(lineSegementsIntersect2D(xi, xj, parallogram[1], parallogram[2], intersection))
                        {
                            intersected = true;
                            intersecting_edge = 1;
                        }
                        else if(lineSegementsIntersect2D(xi, xj, parallogram[2], parallogram[3], intersection))
                        {
                            intersected = true;
                            intersecting_edge = 2;
                        }
                        else if (lineSegementsIntersect2D(xi, xj, parallogram[3], parallogram[0], intersection))
                        {
                            intersected = true;
                            intersecting_edge = 3;
                        }
                        if (intersected)
                        {
                            
                            output_points.push_back(intersection);
                            int xi_idx = output_points.size() - 1;
                            crossing_tracker[xi_idx] = {xi_idx};

                            auto find_xj_iter = std::find_if(output_points.begin(), output_points.end(), [&xj](const TV2 x)->bool
                               { return (x - xj).norm() < 1e-6; }
                             );
                            if (find_xj_iter == output_points.end())
                            {
                                output_points.push_back(xj);
                                xj_idx = int(output_points.size()) - 1;
                                crossing_tracker[xj_idx] = {int(edge_pairs.size())};
                            }
                            else
                            {
                                int index = std::distance(output_points.begin(), find_xj_iter);
                                if (crossing_tracker.find(index) == crossing_tracker.end())
                                {
                                    crossing_tracker[index] = {int(edge_pairs.size())};
                                }
                                else
                                {
                                    crossing_tracker[index].push_back(int(edge_pairs.size()));
                                }
                                xj_idx = index;
                            }

                            boundary_pairs[intersecting_edge].push_back(Vector<int,2>(xj_idx, xi_idx));
                            boundary_pair_rod_idx[intersecting_edge].push_back(edge_pairs.size());
                            edge_pairs.push_back(Vector<int,2>(xj_idx, xi_idx));
                        }
                    }
                    
                }
                else if(!xj_inside && xi_inside)
                {
                    // ignored due to duplicated edges
                }
                else
                {
                    // continue;
                }
            }
        }
    
        
    }
    
}

template<class T, int dim>
void UnitPatch<T, dim>::buildStraightRodScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.new_frame_work = true;

        clearSimData();

        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;
        
        TV from = TV(0, 0.5, 0) * sim.unit;
        TV to = TV(2, 0.5001, 0.001) * sim.unit;

        std::vector<int> passing_points_id;
        std::vector<TV> passing_points;

        addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
        
        int dof_cnt;
        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
    
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());


        int cnt = 0;
        for (auto& rod : sim.Rods)
        {
            std::cout << rod->kt << std::endl;
            // rod->kt  =0 ;
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            sim.dirichlet_dof[rod->theta_reduced_dof_start_offset] = 0;
            // sim.dirichlet_dof[rod->theta_reduced_dof_start_offset + rod->indices.size()-1] = 0;
            rod->setupBishopFrame();
            Offset end0, end1;
            rod->frontOffset(end0); rod->backOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        }

        Offset end0, end1;
        sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0, 0));

        sim.Rods[0]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
        sim.Rods[0]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);

        // sim.Rods[0]->fixPointLagrangian(sim.Rods[0]->indices.size() - 2, 
        //         TV(-0.3, 0.01, 0.01) * sim.unit, 
        //         sim.dirichlet_dof);
        // sim.Rods[0]->fixPointLagrangian(sim.Rods[0]->indices.size() - 1, 
        //         TV(-0.3, 0.01, 0.01) * sim.unit, 
        //         sim.dirichlet_dof); 


    }
}

template<class T, int dim>
void UnitPatch<T, dim>::addCurvedRod(const std::vector<TV2>& data_points,
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div, int& full_dof_cnt, int& node_cnt, int& rod_cnt, bool closed)
{
    if (passing_points.size() != passing_points_id.size())
        std::cout << " passing_points.size() != passing_points_id.size() " << std::endl;

    int first_node_idx = node_cnt;
    int sub_div_2 = sub_div / 2;
    HybridC2Curve<T, 2>* curve = new HybridC2Curve<T, 2>(sub_div);
    for (const auto& pt : data_points)
        curve->data_points.push_back(pt);
    
    std::vector<TV2> points_on_curve;
    curve->sampleCurves(points_on_curve);
    // for (auto data_pt : data_points)
    //     std::cout << data_pt.transpose() << " | ";
    // std::cout << std::endl;
    // for (auto data_pt : points_on_curve)
    //     std::cout << data_pt.transpose() << " | ";
    // std::cout << std::endl;
    // // std::cout << points_on_curve.size() << std::endl;
    // std::getchar();

    std::unordered_map<int, int> dof_node_location;
    if (closed)
        deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size() - 1 - passing_points_id.size()) * (dim + 1));
    else
        deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size() - passing_points_id.size()) * (dim + 1));
    

    
    Rod<T, dim>* rod = new Rod<T, dim>(deformed_states, sim.rest_states, rod_cnt, closed, ROD_A, ROD_B);
    std::unordered_map<int, Offset> offset_map;
    std::vector<int> node_index_list;
    std::vector<T> data_points_discrete_arc_length;
    int full_dof_before = full_dof_cnt;
    int not_found_cnt = 0;
    for (int i = 0; i < points_on_curve.size(); i++)
    {
        
        TV2 pt = points_on_curve[i];
        TV pt_search;
        pt_search.template segment<2>(0) = pt;
        // std::cout << "pt on curve " << pt.transpose() << std::endl;
        //if points already added as crossings
        auto find_node_iter = std::find_if(passing_points.begin(), passing_points.end(), 
            [&pt_search](TV pt_in_vec)
            { return (pt_search - pt_in_vec).norm() < 1e-8; }
        );
        
        if (find_node_iter == passing_points.end())
        {
            if (closed && i == points_on_curve.size() - 1)
            {
                node_index_list.push_back(first_node_idx);
                break;
            }
            offset_map[node_cnt] = Offset::Zero();
            node_index_list.push_back(node_cnt);
            //push Lagrangian DoF
            for (int d = 0; d < dim; d++) 
            {
                deformed_states[full_dof_cnt] = pt_search[d];
                offset_map[node_cnt][d] = full_dof_cnt++;  
            }
            // push Eulerian DoF
            deformed_states[full_dof_cnt] = T(i) * (curve->data_points.size() - 1) / (points_on_curve.size() - 1);
            offset_map[node_cnt][dim] = full_dof_cnt++;
            node_cnt++;
        }
        else
        {
            not_found_cnt++;
            int dof_loc = std::distance(passing_points.begin(), find_node_iter);
            
            // std::cout << passing_points[dof_loc].transpose() << std::endl;

            node_index_list.push_back(passing_points_id[dof_loc]);
            dof_node_location[dof_loc] = i;
            if (closed && i == points_on_curve.size() - 1)
                dof_node_location[dof_loc] = 0;
        }
           
    }
    // checking scene
    // if (not_found_cnt != passing_points.size() + int(closed))
        // std::cout << "not_found_cnt " << not_found_cnt << " should be: passing_points.size() " << passing_points.size() + int(closed) << std::endl;
    int dof_added = full_dof_cnt - full_dof_before;

    int dof_added_should_be = (points_on_curve.size() - passing_points.size() - int(closed)) * (dim + 1);
    if (dof_added != dof_added_should_be)
        std::cout << "after Lagrangian dof_added " << dof_added << " dof_added should be " << dof_added_should_be << std::endl;
    
    // now we add the Eulerian Dof of the passing points
    deformed_states.conservativeResize(full_dof_cnt + passing_points.size());

    for (int i = 0; i < passing_points.size(); i++)
    {
        offset_map[passing_points_id[i]] = Offset::Zero();
        for (int d = 0; d < dim; d++)
            offset_map[passing_points_id[i]][d] = passing_points_id[i] * dim + d;
        
        deformed_states[full_dof_cnt] = T(dof_node_location[i]) * (curve->data_points.size() - 1) / (points_on_curve.size() - 1);
        
        offset_map[passing_points_id[i]][dim] = full_dof_cnt++; 
    }

    dof_added = full_dof_cnt - full_dof_before;
    dof_added_should_be += passing_points.size();
    if (dof_added != dof_added_should_be)
        std::cout << "after eulerian:  dof_added " << dof_added << " dof_added should be " << dof_added_should_be << std::endl;
    
    rod->offset_map = offset_map;
    rod->indices = node_index_list;
    
    for(int i = 0; i < curve->data_points.size(); i++)
    {
        int node_idx = rod->indices[i * sub_div_2];        
        data_points_discrete_arc_length.push_back(deformed_states[offset_map[node_idx][dim]]);
    }

    Vector<T, dim + 1> q0, q1;
    rod->frontDoF(q0); rod->backDoF(q1);

    DiscreteHybridCurvature<T, dim>* rest_state_rod0 = new DiscreteHybridCurvature<T, dim>(q0, q1);
    sim.curvature_functions.push_back(rest_state_rod0);
    rest_state_rod0->setData(curve, data_points_discrete_arc_length);
    

    rod->rest_state = rest_state_rod0;
    std::vector<int> ordered_location;
    for (auto item : dof_node_location)
        ordered_location.push_back(item.second);
    std::sort(ordered_location.begin(), ordered_location.end());
    
    
    rod->dof_node_location = ordered_location;

    sim.Rods.push_back(rod);

    rod_cnt++;
    
}

template<class T, int dim>
void UnitPatch<T, dim>::buildOmegaScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.new_frame_work = true;

        clearSimData();

        std::vector<Entry> w_entry;
    
        int full_dof_cnt = 0;
        int node_cnt = 0;
        int rod_cnt = 0;
        std::vector<TV2> data_points;

        data_points.push_back(TV2(0.0, 0.0) * sim.unit);
        data_points.push_back(TV2(0.4, 0.0) * sim.unit);
        data_points.push_back(TV2(0.2, 0.3) * sim.unit);
        data_points.push_back(TV2(0.8, 0.3) * sim.unit);
        data_points.push_back(TV2(0.6, 0.0) * sim.unit);
        data_points.push_back(TV2(1.0, 0.0) * sim.unit);
        
        std::vector<TV> passing_points = {};
        std::vector<int> passing_points_id = {};

        addCurvedRod(data_points, passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt, false);

        int dof_cnt = 0;
        // markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;
    
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;
        
        int cnt = 0;
        for (auto& rod : sim.Rods)
        {            
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            // sim.dirichlet_dof[rod->theta_reduced_dof_start_offset] = M_PI/8.0;
            // sim.dirichlet_dof[rod->theta_reduced_dof_start_offset] = 0.0;
            rod->setupBishopFrame();
            Offset end0, end1;
            rod->frontOffset(end0); rod->backOffset(end1);
            sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        }

        sim.rest_states = sim.deformed_states;
        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());

        Offset end0, end1;
        sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0,0));

        Offset ob, of;
        sim.Rods[0]->backOffsetReduced(ob);
        sim.Rods[0]->frontOffsetReduced(of);


        sim.dirichlet_dof[of[0]] = 0.0 * sim.unit;
        sim.dirichlet_dof[of[1]] = 0.0 * sim.unit;
        sim.dirichlet_dof[of[2]] = 0.0 * sim.unit;

        sim.dirichlet_dof[ob[0]] = 0.2 * sim.unit;
        sim.dirichlet_dof[ob[1]] = 0.0 * sim.unit;
        sim.dirichlet_dof[ob[2]] = 0.1 * sim.unit;

        
    }
}

// this assumes passing points to be pushed before everything
template<class T, int dim>
void UnitPatch<T, dim>::addAStraightRod(const TV& from, const TV& to, 
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div, int& full_dof_cnt, int& node_cnt, int& rod_cnt)
{
    
    std::unordered_map<int, Offset> offset_map;

    std::vector<TV> points_on_curve;
    std::vector<int> rod_indices;
    std::vector<int> key_points_location_rod;
    addStraightYarnCrossNPoints(from, to, passing_points, passing_points_id,
                                sub_div, points_on_curve, rod_indices,
                                key_points_location_rod, node_cnt);
                   

    deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size()) * (dim + 1));

    Rod<T, dim>* rod = new Rod<T, dim>(deformed_states, sim.rest_states, rod_cnt, false, ROD_A, ROD_B);

    for (int i = 0; i < points_on_curve.size(); i++)
    {
        offset_map[node_cnt] = Offset::Zero();
        //push Lagrangian DoF    
        deformed_states.template segment<dim>(full_dof_cnt) = points_on_curve[i];
        for (int d = 0; d < dim; d++)
        {
            offset_map[node_cnt][d] = full_dof_cnt++;  
        }
        // push Eulerian DoF
        deformed_states[full_dof_cnt] = (points_on_curve[i] - from).norm() / (to - from).norm();
        offset_map[node_cnt][dim] = full_dof_cnt++;
        node_cnt++;
    }
    
    deformed_states.conservativeResize(full_dof_cnt + passing_points.size());

    for (int i = 0; i < passing_points.size(); i++)
    {
        deformed_states[full_dof_cnt] = (passing_points[i] - from).norm() / (to - from).norm();
        offset_map[passing_points_id[i]] = Offset::Zero();
        offset_map[passing_points_id[i]][dim] = full_dof_cnt++; 
        Vector<int, dim> offset_dof_lag;
        for (int d = 0; d < dim; d++)
        {
            offset_dof_lag[d] = passing_points_id[i] * dim + d;
        }
        offset_map[passing_points_id[i]].template segment<dim>(0) = offset_dof_lag;
    }
    
    rod->offset_map = offset_map;
    rod->indices = rod_indices;
    Vector<T, dim + 1> q0, q1;
    rod->frontDoF(q0); rod->backDoF(q1);

    rod->rest_state = new LineCurvature<T, dim>(q0, q1);
    
    rod->dof_node_location = key_points_location_rod;
    
    sim.Rods.push_back(rod);
    rod_cnt++;
}

template<class T, int dim>
void UnitPatch<T, dim>::appendThetaAndJointDoF(std::vector<Entry>& w_entry, 
    int& full_dof_cnt, int& dof_cnt)
{
    // for (auto& rod : sim.Rods)
    // {
    //     rod->theta_dof_start_offset = full_dof_cnt;
    //     rod->theta_reduced_dof_start_offset = dof_cnt;
    //     deformed_states.conservativeResize(full_dof_cnt + rod->indices.size() - 1);
    //     for (int i = 0; i < rod->indices.size() - 1; i++)
    //         w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
    //     deformed_states.template segment(rod->theta_dof_start_offset, 
    //         rod->indices.size() - 1).setZero();
    // }   

    deformed_states.conservativeResize(full_dof_cnt + sim.rod_crossings.size() * dim);
    deformed_states.template segment(full_dof_cnt, sim.rod_crossings.size() * dim).setZero();

    for (auto& crossing : sim.rod_crossings)
    {
        crossing->dof_offset = full_dof_cnt;
        crossing->reduced_dof_offset = dof_cnt;
        for (int d = 0; d < dim; d++)
        {
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        }
    }

    for (auto& rod : sim.Rods)
    {
        rod->theta_dof_start_offset = full_dof_cnt;
        rod->theta_reduced_dof_start_offset = dof_cnt;
        deformed_states.conservativeResize(full_dof_cnt + rod->indices.size() - 1);
        for (int i = 0; i < rod->indices.size() - 1; i++)
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        deformed_states.template segment(rod->theta_dof_start_offset, 
            rod->indices.size() - 1).setZero();
    }
}


template<class T, int dim>
void UnitPatch<T, dim>::buildGridScene2(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = false;

        sim.ke = 1e-4;

        sim.unit = 0.09;
        sim.visual_R = 0.0035;
        
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;

        int n_row = 20, n_col = 20;

        // push crossings first 
        T dy = 1.0 / n_row * sim.unit;
        T dx = 1.0 / n_col * sim.unit;
        
        //num of crossing
        deformed_states.resize(n_col * n_row * dim);
        
        std::unordered_map<int, Offset> crossing_offset_copy;

        auto getXY = [=](int row, int col, T& x, T& y)
        {
            if (row == 0) y = 0.5 * dy;
            else if (row == n_row) y = n_row * dy;
            else y = 0.5 * dy + (row ) * dy;
            if (col == 0) x = 0.5 * dx;
            else if (col == n_col) x = n_col * dx;
            else x = 0.5 * dx + (col ) * dx;
        };


        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                T x, y;
                getXY(row, col, x, y);
                deformed_states.template segment<dim>(node_cnt * dim) = TV(x, y, 0);
                
                full_dof_cnt += dim;
                node_cnt ++;       
            }
        }

        int rod_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            T x0 = 0.0, x1 = 1.0 * sim.unit;
            T x, y;
            
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }

            getXY(row, 0, x, y);

            TV from = TV(x0, y, 0);
            TV to = TV(x1, y, 0);
        
            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (int col = 0; col < n_col; col++)
        {
            T y0 = 0.0, y1 = 1.0 * sim.unit;
            T x, y;
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            getXY(0, col, x, y);
            for (int row = 0; row < n_row; row++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }
            
            TV from = TV(x, y0, 0);
            TV to = TV(x, y1, 0);

            addAStraightRod(from, to, passing_points, passing_points_id, sub_div, 
                            full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        T dv = 1.0 / n_row;
        T du = 1.0 / n_col;

        int odd_even_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(node_idx, {row, n_row + col});

                crossing->on_rod_idx[row] = sim.Rods[row]->dof_node_location[col];
                crossing->on_rod_idx[n_row + col] = sim.Rods[n_row + col]->dof_node_location[row];
                
                // if (odd_even_cnt % 2 == 0)
                //     crossing->is_fixed = true;

                // if (row % 2 == 0)
                //     crossing->is_fixed = true;

                // if (row ==  col)
                //     crossing->is_fixed = true;

                // if (row != col)
                //     crossing->is_fixed = true;

                // crossing->is_fixed = true;

                // if (row == 0 || row == n_row - 1 || col == 0 || col == n_col - 1)                    
                    crossing->is_fixed = true;

                sim.Rods[row]->fixed_by_crossing[col] = false;
                sim.Rods[n_row + col]->fixed_by_crossing[row] = false;
                // if (col % 2 == 0)
                {
                    // crossing->sliding_ranges.push_back(Range(0, 0));    
                    // crossing->sliding_ranges.push_back(Range(1.0/20.0 - 1e3, 1.0/20.0 - 1e3));
                    // crossing->sliding_ranges.push_back(Range(1, 1));    
                    crossing->sliding_ranges.push_back(Range(1, 1));    
                    crossing->sliding_ranges.push_back(Range(0, 0));    
                }
                // else
                // {
                //     crossing->sliding_ranges.push_back(Range(0.02, 0.02));    
                //     crossing->sliding_ranges.push_back(Range(0, 0));
                // }
                
                sim.rod_crossings.push_back(crossing);
                odd_even_cnt++;
            }
        }    

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;

        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        

        T r = 0.05 * sim.unit;
        TV center1, center2;
        sim.getCrossingPosition(0, center1);
        sim.getCrossingPosition(n_row * n_col - 1, center2);

        TV delta1 = TV(-0.05, -0.05, -1e-2) * sim.unit;

        auto circle1 = [r, center1, delta1](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = delta1;
            return (x - center1).norm() < r;
        };

        TV delta2 = TV(0.0, 0.0, 0) * sim.unit;
        auto circle2 = [r, center2, delta2](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = delta2;
            return (x - center2).norm() < r;

        };
        
        TV bottom_left, top_right;
        sim.computeBoundingBox(bottom_left, top_right);
        TV shear_y_left = TV(0.0, 0.5, 0.1) * sim.unit;
        TV shear_y_right = TV(0.0, 0.0, 0) * sim.unit;
        
        T rec_width = 0.1 * sim.unit;

        auto rec1 = [bottom_left, top_right, shear_y_left, rec_width](
            const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = shear_y_left;
            T one_third = (top_right[1] - bottom_left[1]) / 3.0;
            if (x[0] < bottom_left[0] + rec_width 
                // &&
                // (x[1] > bottom_left[1] + one_third && x[1] < bottom_left[1] + one_third * 2)
                )
                return true;
            return false;
        };

        auto rec2 = [bottom_left, top_right, shear_y_right, rec_width](
            const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            mask = Vector<bool, dim>(true, true, true);
            delta = shear_y_right;
            T one_third = (top_right[1] - bottom_left[1]) / 3.0;

            if (x[0] > top_right[0] - rec_width 
                // && 
                // (x[1] > bottom_left[1] + one_third && x[1] < bottom_left[1] + one_third * 2)
                )
                return true;
            return false;
        };

        sim.fixRegionalDisplacement(circle1);
        sim.fixRegionalDisplacement(circle2);

        // sim.fixRegionalDisplacement(rec1);
        // sim.fixRegionalDisplacement(rec2);


        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        {
            if (crossing->is_fixed)
                return;
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // sim.perturb[z_off] += 0.001 * r * sim.unit;
        }
        // VectorXT dq = VectorXT::Zero(sim.W.cols());
        // sim.checkHessianPD(dq);
        // GCodeGenerator<T, dim>(this->sim, "fabric_patch.gcode").generateGCodeFromRodsGridHardCoded();
    }
}


template<class T, int dim>
void UnitPatch<T, dim>::buildSaddleScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;

        sim.ke = 1e-4;

        sim.unit = 0.09;
        
        
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;

        int n_row = 9, n_col = 9;

        // push crossings first 
        T dy = 1.0 / n_row * sim.unit;
        T dx = 1.0 / n_col * sim.unit;
        
        //num of crossing
        deformed_states.resize(n_col * n_row * dim);
        
        std::unordered_map<int, Offset> crossing_offset_copy;

        auto getXY = [=](int row, int col, T& x, T& y)
        {
            if (row == 0) y = 0.5 * dy;
            else if (row == n_row) y = n_row * dy;
            else y = 0.5 * dy + (row ) * dy;
            if (col == 0) x = 0.5 * dx;
            else if (col == n_col) x = n_col * dx;
            else x = 0.5 * dx + (col ) * dx;
        };


        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                T x, y;
                getXY(row, col, x, y);
                deformed_states.template segment<dim>(node_cnt * dim) = TV(x, y, 0);
                
                full_dof_cnt += dim;
                node_cnt ++;       
            }
        }

        int rod_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            T x0 = 0.0, x1 = 1.0 * sim.unit;
            T x, y;
            
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }

            getXY(row, 0, x, y);

            TV from = TV(x0, y, 0);
            TV to = TV(x1, y, 0);
        
            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (int col = 0; col < n_col; col++)
        {
            T y0 = 0.0, y1 = 1.0 * sim.unit;
            T x, y;
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            getXY(0, col, x, y);
            for (int row = 0; row < n_row; row++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }
            
            TV from = TV(x, y0, 0);
            TV to = TV(x, y1, 0);

            addAStraightRod(from, to, passing_points, passing_points_id, sub_div, 
                            full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        T dv = 1.0 / n_row;
        T du = 1.0 / n_col;

        int half_row = (n_row - 1) / 2;
        int half_col = (n_col - 1) / 2;
        
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(node_idx, {row, n_row + col});

                // crossing->sliding_ranges = { Range(0.4 * du, 0.4 * du), Range(0.4 * dv, 0.4 * dv)};
                

                crossing->on_rod_idx[row] = sim.Rods[row]->dof_node_location[col];
                crossing->on_rod_idx[n_row + col] = sim.Rods[n_row + col]->dof_node_location[row];
                
                // if (row < n_row - 1 && (col == n_col - 1 || col == 0))
                // if (row < n_row - 1 && col == (n_col - 1) / 2)
                // if (row < n_row - 1 && col > 0 && col < n_col - 1)
                // if (row < n_row - 1 && col == (n_col - 1) / 2 || node_idx == 7)
                // if (false)
                // if (true)
                // if (row > 0 && row < n_row - 1 && col > 0 && col < n_col - 1)
                // if (row > 0 && row < n_row - 1)
                // {
                //     crossing->is_fixed = false;
                //     sim.Rods[row]->fixed_by_crossing[col] = false;
                //     sim.Rods[n_row + col]->fixed_by_crossing[row] = false;
                //     crossing->sliding_ranges.push_back(Range(0, 0));
                //     crossing->sliding_ranges.push_back(Range(0.4 * dv, 0.4 * dv));
                    
                // }
                // else
                // {
                //     crossing->is_fixed = true;
                    
                //     sim.Rods[row]->fixed_by_crossing[col] = true;
                //     sim.Rods[n_row + col]->fixed_by_crossing[row] = true;

                //     crossing->sliding_ranges.push_back(Range(0, 0));
                //     crossing->sliding_ranges.push_back(Range(0, 0));
                // }

                // if ((row == 0 || row == n_row - 1)  && (col == 0 || col == n_col - 1 || col == half_col))
                if ((row == 0 || row == n_row - 1)  && (col == 0 || col == n_col - 1))
                // if(true)
                {
                    crossing->is_fixed = true;
                    
                    sim.Rods[row]->fixed_by_crossing[col] = true;
                    sim.Rods[n_row + col]->fixed_by_crossing[row] = true;

                    crossing->sliding_ranges.push_back(Range(0, 0));
                    crossing->sliding_ranges.push_back(Range(0, 0));
                    
                }
                else
                {
                    crossing->is_fixed = false;
                    sim.Rods[row]->fixed_by_crossing[col] = false;
                    sim.Rods[n_row + col]->fixed_by_crossing[row] = false;
                    crossing->sliding_ranges.push_back(Range(0, 0));
                    crossing->sliding_ranges.push_back(Range(0.4 * dv, 0.4 * dv));
                }
                sim.rod_crossings.push_back(crossing);
            }
        }    

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;

        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        //drag
        Offset offset;
        

        T r = 0.1 * sim.unit;
        TV center1, center2, center3, center4, center5, center6;
        sim.Rods.front()->x(sim.Rods[0]->indices.front(), center1);
        sim.Rods.front()->x(sim.Rods[0]->indices.back(), center2);

        sim.Rods[n_row - 1]->x(sim.Rods[n_row - 1]->indices.front(), center3);
        sim.Rods[n_row - 1]->x(sim.Rods[n_row - 1]->indices.back(), center4);

        
        sim.Rods[half_row]->x(sim.Rods[half_row]->indices.front(), center5);
        sim.Rods[half_row]->x(sim.Rods[half_row]->indices.back(), center6);

        std::vector<std::function<bool(const TV&)>> circles;
        
        T move_x = 0.05;
        T move_y = 0.01;

        TV delta1 = TV(move_x, move_y, 0.2) * sim.unit;

        auto circle1 = [r, center1, delta1](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta1;
            mask.setConstant(true);
            return (x - center1).norm() < r;
        };

        TV delta2  = TV(-move_x, move_y, 0.2) * sim.unit;

        auto circle2 = [r, center2, delta2](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta2;
            mask.setConstant(true);
            return (x - center2).norm() < r;
        };

        TV delta3 = TV(move_x, -move_y, 0.2) * sim.unit;
        auto circle3 = [r, center3, delta3](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta3;
            mask.setConstant(true);
            return (x - center3).norm() < r;
        };

        TV delta4 = TV(-move_x, -move_y, 0.2) * sim.unit;
        auto circle4 = [r, center4, delta4](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta4;
            mask.setConstant(true);
            return (x - center4).norm() < r;
        };

        TV delta5 = TV(0.0, 0.0, -0.1) * sim.unit;

        auto circle5 = [r, center5, delta5](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta5;
            mask.setConstant(true);
            return (x - center5).norm() < r;
        };

        auto circle6 = [r, center6, delta5](const TV& x, TV& delta, Vector<bool, dim>& mask)->bool
        {
            delta = delta5;
            mask.setConstant(true);
            return (x - center6).norm () < r;
        };

        sim.fixCrossing();

        // sim.fixRegionalDisplacement(circle1);
        // sim.fixRegionalDisplacement(circle2);
        // sim.fixRegionalDisplacement(circle3);
        // sim.fixRegionalDisplacement(circle4);

        // sim.fixRegionalDisplacement(circle5);
        // sim.fixRegionalDisplacement(circle6);

        // sim.Rods.back()->frontOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = -sim.unit * 0.0;
        // sim.dirichlet_dof[offset[2]] = 0;

        

        sim.Rods[half_row]->fixPointLagrangianByID(36, TV(0.05, 0, -0.1) * sim.unit, Mask(true, false, true), sim.dirichlet_dof);
        sim.Rods[half_row]->fixPointLagrangianByID(44, TV(-0.05, 0, -0.1) * sim.unit, Mask(true, false, true), sim.dirichlet_dof);

        // sim.Rods[half_row]->fixPointLagrangianByID(40, TV(0, 0, 0) * sim.unit, Mask(true, false, true), sim.dirichlet_dof);

        sim.Rods[n_row + half_col]->fixPointLagrangianByID(4, TV(0.0, -0.05, 0.1) * sim.unit, Mask(false, true, true), sim.dirichlet_dof);
        sim.Rods[n_row + half_col]->fixPointLagrangianByID(76, TV(0.0, 0.05, 0.1) * sim.unit, Mask(false, true, true), sim.dirichlet_dof);

        // sim.Rods[n_row + (int)((n_col - 1)/2)]->frontOffsetReduced(offset);
        // sim.dirichlet_dof[offset[0]] = 0;
        // sim.dirichlet_dof[offset[1]] = -sim.unit * 0.4;
        // sim.dirichlet_dof[offset[2]] = 0;    

        Mask maskY(false, true, true);

        move_y = 0.0;
        T move_z = 0.0;

        // sim.dirichlet_dof[sim.rod_crossings[0]->reduced_dof_offset + 1] = -M_PI / 4.0;
        // sim.dirichlet_dof[sim.rod_crossings[0]->reduced_dof_offset] = M_PI / 8.0;
        sim.fixCrossingLagrangian(0, TV(0, move_y, 0)* sim.unit, maskY);

        // sim.dirichlet_dof[sim.rod_crossings[n_row - 1]->reduced_dof_offset + 1] = M_PI / 4.0;
        // sim.dirichlet_dof[sim.rod_crossings[n_row - 1]->reduced_dof_offset] = M_PI / 8.0;
        sim.fixCrossingLagrangian(n_row - 1, TV(0, move_y, 0)* sim.unit, maskY);


        // sim.dirichlet_dof[sim.rod_crossings[sim.rod_crossings.size() - n_col]->reduced_dof_offset + 1] = -M_PI / 4.0;
        sim.fixCrossingLagrangian(sim.rod_crossings.size() - n_col, TV(0, -move_y, move_z)* sim.unit, maskY);

        // sim.dirichlet_dof[sim.rod_crossings.back()->reduced_dof_offset + 1] = M_PI / 4.0;
        sim.fixCrossingLagrangian(sim.rod_crossings.size() - 1, TV(0, -move_y, move_z) * sim.unit, maskY);

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }
        
        
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildShelterScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        clearSimData();

        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.add_pbc_twisting = false;
        sim.add_pbc = false;

        sim.add_contact_penalty=true;
        sim.new_frame_work = true;
        sim.add_eularian_reg = true;


        sim.unit = 0.03;
        
        
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;

        int n_row = 4, n_col = 4;

        // push crossings first 
        T dy = 1.0 / n_row * sim.unit;
        T dx = 1.0 / n_col * sim.unit;
        
        //num of crossing
        deformed_states.resize(n_col * n_row * dim);
        
        std::unordered_map<int, Offset> crossing_offset_copy;

        auto getXY = [=](int row, int col, T& x, T& y)
        {
            if (row == 0) y = 0.5 * dy;
            else if (row == n_row) y = n_row * dy;
            else y = 0.5 * dy + (row ) * dy;
            if (col == 0) x = 0.5 * dx;
            else if (col == n_col) x = n_col * dx;
            else x = 0.5 * dx + (col ) * dx;
        };


        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                T x, y;
                getXY(row, col, x, y);
                deformed_states.template segment<dim>(node_cnt * dim) = TV(x, y, 0);
                
                full_dof_cnt += dim;
                node_cnt ++;       
            }
        }

        int rod_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            T x0 = 0.0, x1 = 1.0 * sim.unit;
            T x, y;
            
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }

            getXY(row, 0, x, y);

            TV from = TV(x0, y, 0);
            TV to = TV(x1, y, 0);
        
            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (int col = 0; col < n_col; col++)
        {
            T y0 = 0.0, y1 = 1.0 * sim.unit;
            T x, y;
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            getXY(0, col, x, y);
            for (int row = 0; row < n_row; row++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }
            
            TV from = TV(x, y0, 0);
            TV to = TV(x, y1, 0);

            addAStraightRod(from, to, passing_points, passing_points_id, sub_div, 
                            full_dof_cnt, node_cnt, rod_cnt);
            
        }
        
        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

        T dv = 1.0 / n_row;
        T du = 1.0 / n_col;

        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(node_idx, {row, n_row + col});

                // crossing->sliding_ranges = { Range(0.4 * du, 0.4 * du), Range(0.4 * dv, 0.4 * dv)};
                

                crossing->on_rod_idx[row] = sim.Rods[row]->dof_node_location[col];
                crossing->on_rod_idx[n_row + col] = sim.Rods[n_row + col]->dof_node_location[row];
                
                // if (row < n_row - 1 && (col == n_col - 1 || col == 0))
                if (row < n_row - 1 && (col == n_col - 1 || col == 0))
                // if (false)
                {
                    crossing->is_fixed = false;
                    sim.Rods[row]->fixed_by_crossing[col] = false;
                    sim.Rods[n_row + col]->fixed_by_crossing[row] = false;
                    crossing->sliding_ranges.push_back(Range(0, 0));
                    crossing->sliding_ranges.push_back(Range(0.2, 0.01));
                }
                else
                {
                    crossing->is_fixed = true;
                    
                    sim.Rods[row]->fixed_by_crossing[col] = true;
                    sim.Rods[n_row + col]->fixed_by_crossing[row] = true;

                    crossing->sliding_ranges.push_back(Range(0, 0));
                    crossing->sliding_ranges.push_back(Range(0, 0));
                }
                sim.rod_crossings.push_back(crossing);
            }
        }    

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        // std::cout << sim.W << std::endl;

        
        for (auto& rod : sim.Rods)
        {
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            rod->setupBishopFrame();
        }
        
    
        //drag
        Offset offset;
        sim.Rods.front()->fixEndPointLagrangian(sim.dirichlet_dof);

        // for (int i = 0; i < n_col; i ++)
        //     sim.Rods.front()->fixPointLagrangianByID(i, TV::Zero(), sim.dirichlet_dof);

        // sim.dirichlet_dof[sim.Rods.front()->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[sim.Rods.front()->theta_reduced_dof_start_offset + sim.Rods.front()->numSeg() - 1] = 0;

        T r = 0.05 * sim.unit;
        TV center1, center2;
        sim.Rods.front()->x(sim.Rods[0]->indices.front(), center1);
        sim.Rods.front()->x(sim.Rods[0]->indices.back(), center2);


        auto circle1 = [r, center1](const TV& x)->bool
        {
            return (x - center1).norm() < r;
        };

        auto circle2 = [r, center2](const TV& x)->bool
        {
            return (x - center2).norm() < r;
        };

        sim.fixRegion(circle1);
        sim.fixRegion(circle2);

        sim.Rods.back()->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -sim.unit * 0.2;
        sim.dirichlet_dof[offset[2]] = 0;

        sim.Rods[n_row]->frontOffsetReduced(offset);
        sim.dirichlet_dof[offset[0]] = 0;
        sim.dirichlet_dof[offset[1]] = -sim.unit * 0.2;
        sim.dirichlet_dof[offset[2]] = 0;

        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            // sim.perturb[z_off] += 0.001 * r * sim.unit;
            // sim.perturb[z_off] += 0.001 * sim.unit;
            
        }

        GCodeGenerator<T, dim>(this->sim, "shelter.gcode").generateGCodeFromRodsShelterHardCoded();
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildGridScene(int sub_div)
{
    
    if constexpr (dim == 3)
    {
        
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        
        sim.add_rotation_penalty = true;
        sim.add_pbc_bending = true;
        sim.add_pbc = true;
        sim.add_pbc_twisting = true;
        sim.add_contact_penalty=true;
        sim.new_frame_work = true;

        sim.ke = 1e-2;
        
        sim.unit = 0.08;
        
        sim.visual_R = 0.002;
        // sim.unit = 4;
        // std::cout << sim.unit << std::endl;
        
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;

        int n_row = 20, n_col = 20;

        T length = 1.0 / 20.0 * sim.unit;

        // push crossings first 
        T dy = length;
        T dx = length;
        
        //num of crossing
        deformed_states.resize(n_col * n_row * dim);
        
        std::unordered_map<int, Offset> crossing_offset_copy;

        auto getXY = [=](int row, int col, T& x, T& y)
        {
            if (row == 0) y = 0.5 * dy;
            else if (row == n_row) y = n_row * dy;
            else y = 0.5 * dy + (row ) * dy;
            if (col == 0) x = 0.5 * dx;
            else if (col == n_col) x = n_col * dx;
            else x = 0.5 * dx + (col ) * dx;
        };

        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                T x, y;
                getXY(row, col, x, y);
                deformed_states.template segment<dim>(node_cnt * dim) = TV(x, y, 0);
                
                full_dof_cnt += dim;
                node_cnt ++;       
            }
        }

        int rod_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            // T x0 = 0.0, x1 = 1.0 * sim.unit;
            T x0 = 0.0, x1 = dx * n_row;
            T x, y;
            
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }

            getXY(row, 0, x, y);

            TV from = TV(x0, y, 0);
            TV to = TV(x1, y, 0);
        
            addAStraightRod(from, to, passing_points, passing_points_id, 
                sub_div, full_dof_cnt, node_cnt, rod_cnt);
        }
        
        for (int col = 0; col < n_col; col++)
        {
            // T y0 = 0.0, y1 = 1.0 * sim.unit;
            T y0 = 0.0, y1 = dy * n_col;
            T x, y;
            std::vector<int> passing_points_id;
            std::vector<TV> passing_points;
            getXY(0, col, x, y);
            for (int row = 0; row < n_row; row++)
            {
                int node_idx = row * n_col + col;
                passing_points_id.push_back(node_idx);
                passing_points.push_back(deformed_states.template segment<dim>(node_idx * dim));
            }
            
            TV from = TV(x, y0, 0);
            TV to = TV(x, y1, 0);

            addAStraightRod(from, to, passing_points, passing_points_id, sub_div, 
                            full_dof_cnt, node_cnt, rod_cnt);
        }

        for (auto& rod : sim.Rods)
            rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);        

        T dv = 1.0 / n_row;
        T du = 1.0 / n_col;
        int odd_even_cnt = 0;
        for (int row = 0; row < n_row; row++)
        {
            for (int col = 0; col < n_col; col++)
            {
                int node_idx = row * n_col + col;
                RodCrossing<T, dim>* crossing = 
                    new RodCrossing<T, dim>(node_idx, {row, n_row + col});

                // crossing->sliding_ranges = { 
                //                             // Range(0.4 * du, 0.4 * du), 
                //                             Range(0.0, 0.0), 
                //                             // Range(0.4 * dv, 0.4 * dv)
                //                             Range(1, 1)
                //                             };
                // if (col == 0)
                {
                    crossing->sliding_ranges = { 
                                            Range(0.0, 0.0), 
                                            Range(1, 1)
                                            };
                }
                // else
                {
                    // crossing->sliding_ranges = { 
                    //                         Range(1, 1),
                    //                         Range(0.0, 0.0)
                    //                         };
                }
                crossing->on_rod_idx[row] = sim.Rods[row]->dof_node_location[col];
                crossing->on_rod_idx[n_row + col] = sim.Rods[n_row + col]->dof_node_location[row];
                
                
                // if (odd_even_cnt % 2 == 0)
                // if (row % 2 == 0)
                    // crossing->is_fixed = true;
                // else
                    // crossing->is_fixed = false;
                // if (row != 0 && col != 0)
                //     crossing->is_fixed = true;
                // if (row == col)
                    // crossing->is_fixed = true;
                // if (row == 0 && col == 0)
                    // crossing->is_fixed = true;
                sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[row] = false;
                sim.rod_crossings.push_back(crossing);
                odd_even_cnt++;
            }
        }    
        for (int row = 0; row < n_row; row++)
            for (int i = 0; i < sim.Rods[row]->dof_node_location.size(); i++)
                sim.Rods[row]->fixed_by_crossing[i] = false;

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);

        for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
        
        appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
        
        sim.rest_states = deformed_states;

        
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
        
        

        int cnt = 0;
        for (auto& rod : sim.Rods)
        {
            
            rod->fixEndPointEulerian(sim.dirichlet_dof);
            
            rod->setupBishopFrame();
            cnt++;
            Offset end0, end1;
            rod->frontOffset(end0); rod->backOffset(end1);
            if (rod->rod_id < n_row)
                sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
            else
                sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
            
            Offset a, b;
            rod->getEntryByLocation(1, a); rod->getEntryByLocation(rod->indices.size() - 2, b);
            sim.pbc_bending_pairs.push_back({end0, a, b, end1});
            sim.pbc_bending_pairs_rod_id.push_back({rod->rod_id, rod->rod_id, rod->rod_id, rod->rod_id});
        }

        Offset end0, end1;
        sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
        sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0, 0));
        
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[0]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[1]]] = 0;
        sim.dirichlet_dof[sim.Rods[0]->reduced_map[end0[2]]] = 0;

        sim.Rods[n_row]->frontOffset(end0); sim.Rods[n_row]->backOffset(end1);
        sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), std::make_pair(n_row, n_row));


        sim.fixCrossing();

        sim.perturb = VectorXT::Zero(sim.W.cols());
        for (auto& crossing : sim.rod_crossings)
        // for (int i = 0; i < 10; i++)
        {
            // auto crossing = rod_crossings[i];
            Offset off;
            sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[dim-1]];
            // sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
            // break;
            sim.perturb[z_off] += 0.001 * r * sim.unit;
            
        }

        GCodeGenerator<T, dim>(this->sim, "grid_new.gcode").buildGridScene(n_row, n_col, 0);

        // GCodeGenerator<T, dim>(this->sim, "dense_patch_free.gcode").generateGCodeFromRodsGridHardCoded(n_row, n_col, 2);
        // GCodeGenerator<T, dim>(this->sim, "dense_patch_free_change_boundary.gcode").generateGCodeFromRodsGridHardCoded(n_row, n_col, 3);
        // GCodeGenerator<T, dim>(this->sim, "dense_patch_free_change_boundary.gcode").generateGCodeFromRodsGridHardCoded(n_row, n_col, 1);
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::buildOneCrossScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        int sub_div_2 = sub_div / 2;
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.new_frame_work = true;

        clearSimData();
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        std::unordered_map<int, Offset> offset_map;
        
        TV from(0.0, 0.5, 0.0);
        TV to(1.0, 0.5, 0.0);
        from *= sim.unit; to *= sim.unit;

        TV center = TV(0.5, 0.5, 0.0) * sim.unit;
        int center_id = 0;
        deformed_states.resize(dim);
        deformed_states.template segment<dim>(full_dof_cnt) = center;
        offset_map[node_cnt] = Offset::Zero();
        for (int d = 0; d < dim; d++) offset_map[node_cnt][d] = full_dof_cnt++;
        node_cnt++;
        auto center_offset = offset_map[center_id];

        std::vector<TV> points_on_curve;
        std::vector<int> rod0;
        std::vector<int> key_points_location_rod0;

        addStraightYarnCrossNPoints(from, to, {center}, {0}, sub_div, points_on_curve, rod0, key_points_location_rod0, node_cnt);

        deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size()) * (dim + 1));

        Rod<T, dim>* r0 = new Rod<T, dim>(deformed_states, sim.rest_states, 0, false, ROD_A, ROD_B);

        for (int i = 0; i < points_on_curve.size(); i++)
        {
            offset_map[node_cnt] = Offset::Zero();
            
            //push Lagrangian DoF
            
            deformed_states.template segment<dim>(full_dof_cnt) = points_on_curve[i];
            
            for (int d = 0; d < dim; d++)
            {
                offset_map[node_cnt][d] = full_dof_cnt++;    
            }
            // push Eulerian DoF
            deformed_states[full_dof_cnt] = (points_on_curve[i] - from).norm() / (to - from).norm();
            offset_map[node_cnt][dim] = full_dof_cnt++;
            node_cnt++;
        }
        deformed_states.conservativeResize(full_dof_cnt + 1);
        deformed_states[full_dof_cnt] = (center - from).norm() / (to - from).norm();
        offset_map[center_id][dim] = full_dof_cnt++;

        r0->offset_map = offset_map;
        r0->indices = rod0;

        Vector<T, dim + 1> q0, q1;
        r0->frontDoF(q0); r0->backDoF(q1);
        r0->rest_state = new LineCurvature<T, dim>(q0, q1);
        
        r0->dof_node_location = key_points_location_rod0;
        sim.Rods.push_back(r0);

        offset_map.clear();
        
        TV rod1_from(0.5, 0.0, 0.0);
        TV rod1_to(0.5, 1.0, 0.0);
        rod1_from *= sim.unit; rod1_to *= sim.unit;

        points_on_curve.clear();
        points_on_curve.resize(0);
        std::vector<int> rod1;
        std::vector<int> key_points_location_rod1;

        addStraightYarnCrossNPoints(rod1_from, rod1_to, {center}, {0}, sub_div, points_on_curve, rod1, key_points_location_rod1, node_cnt);

        deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size()) * (dim + 1));

        Rod<T, dim>* r1 = new Rod<T, dim>(deformed_states, sim.rest_states, 1, false, ROD_A, ROD_B);
        for (int i = 0; i < points_on_curve.size(); i++)
        {
            offset_map[node_cnt] = Offset::Zero();
            //push Lagrangian DoF
            deformed_states.template segment<dim>(full_dof_cnt) = points_on_curve[i];
            // std::cout << points_on_curve[i].transpose() << std::endl;
            for (int d = 0; d < dim; d++)
            {
                offset_map[node_cnt][d] = full_dof_cnt++;    
            }
            // push Eulerian DoF
            deformed_states[full_dof_cnt] = (points_on_curve[i] - rod1_from).norm() / (rod1_to - rod1_from).norm();
            offset_map[node_cnt][dim] = full_dof_cnt++;
            node_cnt++;
        }

        deformed_states.conservativeResize(full_dof_cnt + 1);

        deformed_states[full_dof_cnt] = (center - rod1_from).norm() / (rod1_to - rod1_from).norm();
        offset_map[center_id] = Offset::Zero();
        offset_map[center_id].template segment<dim>(0) = center_offset.template segment<dim>(0);
        offset_map[center_id][dim] = full_dof_cnt++;

        r1->offset_map = offset_map;
        r1->indices = rod1;

        r1->frontDoF(q0); r1->backDoF(q1);
        r1->rest_state = new LineCurvature<T, dim>(q0, q1);
        
        r1->dof_node_location = key_points_location_rod1;
        sim.Rods.push_back(r1);

        RodCrossing<T, dim>* rc0 = new RodCrossing<T, dim>(0, {0, 1});
        rc0->sliding_ranges = { Range(0.2, 0.2), Range(0.2, 0.2)};
        rc0->on_rod_idx[0] = key_points_location_rod0[0];
        rc0->on_rod_idx[1] =  key_points_location_rod1[0];
        sim.rod_crossings.push_back(rc0);

        

        int dof_cnt = 0;
        markCrossingDoF(w_entry, dof_cnt);
        r0->markDoF(w_entry, dof_cnt);
        r1->markDoF(w_entry, dof_cnt);

        r0->theta_dof_start_offset = full_dof_cnt;
        r0->theta_reduced_dof_start_offset = dof_cnt;        
        int theta_reduced_dof_offset0 = dof_cnt;
        deformed_states.conservativeResize(full_dof_cnt + r0->indices.size() - 1);
        for (int i = 0; i < r0->indices.size() - 1; i++)
        {
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        }   
        deformed_states.template segment(r0->theta_dof_start_offset, 
            r0->indices.size() - 1).setZero();

        r1->theta_dof_start_offset = full_dof_cnt;
        
        int theta_reduced_dof_offset1 = dof_cnt;
        r1->theta_reduced_dof_start_offset = dof_cnt;
        deformed_states.conservativeResize(full_dof_cnt + r1->indices.size() - 1);
        for (int i = 0; i < r1->indices.size() - 1; i++)
        {
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        }   
        deformed_states.template segment(r1->theta_dof_start_offset, 
            r1->indices.size() - 1).setZero();

        deformed_states.conservativeResize(full_dof_cnt + sim.rod_crossings.size() * dim);
        deformed_states.template segment(full_dof_cnt, sim.rod_crossings.size() * dim).setZero();

        for (auto& crossing : sim.rod_crossings)
        {
            crossing->dof_offset = full_dof_cnt;
            crossing->reduced_dof_offset = dof_cnt;
            for (int d = 0; d < dim; d++)
            {
                w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
            }
        }
        
        sim.rest_states = sim.deformed_states;
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());

        // std::cout << "r0->theta_dof_start_offset " << r0->theta_dof_start_offset << " sim.W.cols() " << sim.W.cols() << std::endl;
        

        Offset ob, of;
        r0->backOffsetReduced(ob);
        r0->frontOffsetReduced(of);

        // std::cout << ob.transpose() << " " << of.transpose() << std::endl;
        r0->fixEndPointEulerian(sim.dirichlet_dof);
        r1->fixEndPointEulerian(sim.dirichlet_dof);

        // r1->fixEndPointLagrangian(sim.dirichlet_dof);

        // sim.fixCrossing();

        sim.dirichlet_dof[ob[0]] = -0.3 * sim.unit;
        sim.dirichlet_dof[ob[1]] = 0.3 * sim.unit;
        // sim.dirichlet_dof[ob[2]] = 0;
        sim.dirichlet_dof[ob[2]] = 0.3 * sim.unit;


        sim.dirichlet_dof[theta_reduced_dof_offset0] = 0;
        sim.dirichlet_dof[theta_reduced_dof_offset1] = 0;

        Offset ob1, of1;
        r1->backOffsetReduced(ob1);
        r1->frontOffsetReduced(of1);


        sim.dirichlet_dof[ob1[0]] = 0.15 * sim.unit;
        sim.dirichlet_dof[ob1[1]] = -0.2 * sim.unit;
        sim.dirichlet_dof[ob1[2]] = -0.1 * sim.unit;

        for (int d = 0; d < dim; d++)
        {
            sim.dirichlet_dof[of[d]] = 0;
            sim.dirichlet_dof[ob1[d]] = 0;
            sim.dirichlet_dof[of1[d]] = 0;
        }

        sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][0]]] = 0.0 * sim.unit;
        sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][1]]] = 0.0 * sim.unit;
        sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][2]]] = 0.0 * sim.unit;
        
        for (auto& rod : sim.Rods)
        {
            rod->setupBishopFrame();
        }
        
    }
}
template<class T, int dim>
void UnitPatch<T, dim>::build3DtestScene(int sub_div)
{
    if constexpr (dim == 3)
    {
        int sub_div_2 = sub_div / 2;
        auto unit_yarn_map = sim.yarn_map;
        sim.yarn_map.clear();
        sim.add_rotation_penalty = false;
        sim.add_pbc_bending = false;
        sim.new_frame_work = true;

        clearSimData();
        std::vector<Eigen::Triplet<T>> w_entry;
        int full_dof_cnt = 0;
        int node_cnt = 0;
        
        TV from(0.0, 0.5, 0.0);
        TV to(1.0, 0.5, 0.0);
        from *= sim.unit; to *= sim.unit;

        std::vector<TV> points_on_curve;
        std::vector<int> rod0;
        std::vector<int> dummy;

        addStraightYarnCrossNPoints(from, to, {}, {}, sub_div, points_on_curve, rod0, dummy, 0);

        // std::cout << points_on_curve.size() << " " << rod0.size() << std::endl;

        deformed_states.resize((points_on_curve.size()) * (dim + 1));

        Rod<T, dim>* r0 = new Rod<T, dim>(deformed_states, sim.rest_states, 0, false, ROD_A, ROD_B);

        std::unordered_map<int, Offset> offset_map;
        std::vector<int> node_index_list;

        std::vector<T> data_points_discrete_arc_length;
        
        for (int i = 0; i < points_on_curve.size(); i++)
        {
            offset_map[i] = Offset::Zero();
            node_cnt++;
            node_index_list.push_back(i);
            //push Lagrangian DoF
            
            deformed_states.template segment<dim>(full_dof_cnt) = points_on_curve[i];
            // std::cout << points_on_curve[i].transpose() << std::endl;
            for (int d = 0; d < dim; d++)
            {
                offset_map[i][d] = full_dof_cnt++;    
            }
            // push Eulerian DoF
            deformed_states[full_dof_cnt] = (points_on_curve[i] - from).norm() / (to - from).norm();
            
            offset_map[i][dim] = full_dof_cnt++;
            
        }

        r0->offset_map = offset_map;
        r0->indices = node_index_list;
        // for (int idx : node_index_list)
        //     std::cout << idx << " ";
        // std::cout << std::endl;
        Vector<T, dim + 1> q0, q1;
        r0->frontDoF(q0); r0->backDoF(q1);
        

        r0->rest_state = new LineCurvature<T, dim>(q0, q1);
        
        r0->dof_node_location = {};
        sim.Rods.push_back(r0);

        int dof_cnt = 0;
        
        r0->markDoF(w_entry, dof_cnt);
        r0->theta_dof_start_offset = full_dof_cnt;
        
        int theta_reduced_dof_offset = dof_cnt;
        deformed_states.conservativeResize(full_dof_cnt + r0->indices.size() - 1);
        for (int i = 0; i < r0->indices.size() - 1; i++)
        {
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        }   
        deformed_states.template segment(r0->theta_dof_start_offset, 
            r0->indices.size() - 1).setZero();
        
        sim.rest_states = sim.deformed_states;
        sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
        sim.W.setFromTriplets(w_entry.begin(), w_entry.end());

        // std::cout << "r0->theta_dof_start_offset " << r0->theta_dof_start_offset << " sim.W.cols() " << sim.W.cols() << std::endl;
        

        Offset ob, of;
        r0->backOffsetReduced(ob);
        r0->frontOffsetReduced(of);

        // std::cout << ob.transpose() << " " << of.transpose() << std::endl;
        r0->fixEndPointEulerian(sim.dirichlet_dof);
        

        sim.dirichlet_dof[ob[0]] = -0.3 * sim.unit;
        sim.dirichlet_dof[ob[1]] = 0.3 * sim.unit;
        // sim.dirichlet_dof[ob[2]] = 0;


        for (int i = theta_reduced_dof_offset; i < dof_cnt; i++)
        {
            sim.dirichlet_dof[i] = 0;
            break;
            // sim.dirichlet_dof[i] = T(i) * M_PI / 4;
        }

        // sim.dirichlet_dof[dof_cnt-1] = M_PI / 2.0;
        // sim.dirichlet_dof[dof_cnt-1] = M_PI;
        // sim.dirichlet_dof[dof_cnt-1] = 0;

        
        // sim.dirichlet_dof[ob[0]] = -0.3 * sim.unit;
        // sim.dirichlet_dof[ob[1]] = 0.1 * sim.unit;
        // sim.dirichlet_dof[ob[2]] = 0.0 * sim.unit;

        for (int d = 0; d < dim; d++)
        {
            sim.dirichlet_dof[of[d]] = 0;
        }

        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][0]]] = 0.0 * sim.unit;
        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][1]]] = 0.0 * sim.unit;
        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[1][2]]] = 0.0 * sim.unit;


        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[node_cnt-2][0]]] = -0.19 * sim.unit;
        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[node_cnt-2][1]]] = 0.18 * sim.unit;
        // sim.dirichlet_dof[r0->reduced_map[r0->offset_map[node_cnt-2][2]]] = 0.18 * sim.unit;
        
        for (auto& rod : sim.Rods)
        {
            rod->setupBishopFrame();
        }
        
        // std::ifstream in("./testdq.txt");
        // for(int i =0; i < deformed_states.rows(); i++)
        //     in>> deformed_states[i];
        
        // in.close();

        
        // VectorXT dq = sim.W.transpose() * deformed_states;
        // sim.testGradient(dq);
        // sim.testHessian(dq);

        // for (auto& it : sim.dirichlet_dof)
        //     std::cout << it.first << " " << it.second << std::endl;
        // std::cout << deformed_states << std::endl;
        // std::cout << sim.W << std::endl;
        // std::cout << sim.W.rows() << " " << sim.W.cols() << " " << deformed_states.rows() << " " << r0->theta_dof_start_offset<< std::endl;
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::markCrossingDoF(std::vector<Eigen::Triplet<T>>& w_entry,
        int& dof_cnt)
{
    for (auto& crossing : sim.rod_crossings)
    {
        int node_idx = crossing->node_idx;
        // std::cout << "node " << node_idx << std::endl;
        std::vector<int> rods_involved = crossing->rods_involved;

        Offset entry_rod0; 
        sim.Rods[rods_involved.front()]->getEntry(node_idx, entry_rod0);

        // push Lagrangian dof first
        for (int d = 0; d < dim; d++)
        {
            
            for (int rod_idx : rods_involved)
            {
                // std::cout << "rods involved " << rod_idx << std::endl;
                // if (node_idx == 21)
                //     std::cout << "rods involved " << rod_idx << std::endl;
                sim.Rods[rod_idx]->reduced_map[entry_rod0[d]] = dof_cnt;
            }    
            w_entry.push_back(Entry(entry_rod0[d], dof_cnt++, 1.0));
        }
        
        // push Eulerian dof for all rods
        for (int rod_idx : rods_involved)
        {
            // std::cout << "dim on rod " <<  rod_idx << std::endl;
            sim.Rods[rod_idx]->getEntry(node_idx, entry_rod0);
            // std::cout << "dim dof on rod " <<  entry_rod0[dim] << std::endl;
            sim.Rods[rod_idx]->reduced_map[entry_rod0[dim]] = dof_cnt;
            w_entry.push_back(Entry(entry_rod0[dim], dof_cnt++, 1.0));
        }
        // std::getchar();
        
    }
}

template<class T, int dim>
void UnitPatch<T, dim>::clearSimData()
{
    sim.kc = 1e8;
    sim.add_pbc = true;

    if(sim.disable_sliding)
    {
        sim.add_shearing = true;
        sim.add_eularian_reg = false;
        sim.k_pbc = 1e8;
        sim.k_strain = 1e8;
    }
    else
    {
        sim.add_shearing = false;
        sim.add_eularian_reg = true;
        sim.ke = 1e-4;    
        sim.k_yc = 1e8;
    }
    sim.k_pbc = 1e4;
    sim.k_strain = 1e7;
    sim.kr = 1e3;
    
    // sim.pbc_ref_unique.clear();
    // sim.dirichlet_data.clear();
    // sim.pbc_ref.clear();
    // sim.pbc_bending_pairs.clear();
    sim.yarns.clear();
}

// assuming passing points sorted long from to to direction
template<class T, int dim>
void UnitPatch<T, dim>::addStraightYarnCrossNPoints(const TV& from, const TV& to,
    const std::vector<TV>& passing_points, 
    const std::vector<int>& passing_points_id, int sub_div,
    std::vector<TV>& sub_points, std::vector<int>& node_idx, 
    std::vector<int>& key_points_location, 
    int start, bool pbc)
{
    
    int cnt = 1;
    if(passing_points.size())
    {
        if ((from - passing_points[0]).norm() < 1e-6 )
        {
            node_idx.push_back(passing_points_id[0]);
            cnt = 0;
        }
        else
        {
            node_idx.push_back(start);
            sub_points.push_back(from);
        }
    }
    else
    {
        node_idx.push_back(start);
        sub_points.push_back(from);
    }
    
    T length_yarn = (to - from).norm();
    TV length_vec = (to - from).normalized();
    
    TV loop_point = from;
    TV loop_left = from;
    for (int i = 0; i < passing_points.size(); i++)
    {
        if ((from - passing_points[i]).norm() < 1e-6 )
        {
            key_points_location.push_back(0);
            continue;
        }
        T fraction = (passing_points[i] - loop_point).norm() / length_yarn;
        int n_sub_nodes = std::ceil(fraction * sub_div);
        T length_sub = (passing_points[i] - loop_point).norm() / T(n_sub_nodes);
        for (int j = 0; j < n_sub_nodes - 1; j++)
        {
            sub_points.push_back(loop_left + length_sub * length_vec);
            loop_left = sub_points.back();
            node_idx.push_back(start + cnt);
            cnt++;
        }
        node_idx.push_back(passing_points_id[i]);
        key_points_location.push_back(cnt + i);
        loop_point = passing_points[i];
        loop_left = passing_points[i];
    }
    if (passing_points.size())
    {
        if ((passing_points.back() - to).norm() < 1e-6)
        {
            
            return;
        }
    }
    T fraction;
    int n_sub_nodes;
    T length_sub;
    if( passing_points.size() )
    {
        fraction = (to - passing_points.back()).norm() / length_yarn;
        n_sub_nodes = std::ceil(fraction * sub_div);
        length_sub = (to - passing_points.back()).norm() / T(n_sub_nodes);
    }
    else
    {
        n_sub_nodes = sub_div + 1;
        length_sub = (to - from).norm() / T(sub_div);
    }
    for (int j = 0; j < n_sub_nodes - 1; j++)
    {
        if (j == 0)
        {
            if(passing_points.size())
            {
                sub_points.push_back(passing_points.back() + length_sub * length_vec);
                loop_left = sub_points.back();
            }
        }
        else
        {
            sub_points.push_back(loop_left + length_sub * length_vec);
            loop_left = sub_points.back();
        }
        if(passing_points.size() == 0 && j == 0)
            continue;
        node_idx.push_back(start + cnt);
        cnt++;
    }
    node_idx.push_back(start + cnt);
    sub_points.push_back(to);
}


template class UnitPatch<double, 3>;
template class UnitPatch<double, 2>;   