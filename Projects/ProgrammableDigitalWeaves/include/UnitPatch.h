#ifndef UNIT_PATCH_H
#define UNIT_PATCH_H

#include "EoLRodSim.h"

template<class T, int dim>
class EoLRodSim;

template<class T, int dim>
class UnitPatch
{
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;

    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;

    
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;

    using Offset = Vector<int, dim + 1>;
    using Range = Vector<T, 2>;
    using Mask = Vector<bool, dim>;
    using Mask2 = Vector<bool, 2>;


private:
    EoLRodSim<T, dim>& sim;

    VectorXT& deformed_states = sim.deformed_states;
    

public:
    UnitPatch(EoLRodSim<T, dim>& eol_sim) : sim(eol_sim) {}
    ~UnitPatch() {}

    void buildScene(int patch_type);
    void buildGridClosed(int sub_div);

    //Mechanisim
    void buildFingerScene(int sub_div);
    void buildTennisBallWrapperScene(int sub_div);
    void buildShoeScene(int sub_div);

    void buildActuationSingleStrandScene(int sub_div);
    void buildActuationSingleStrandSceneWithoutCrossing(int sub_div);

    void buildPeriodicCircleScene(int sub_div);
    void buildFullCircleScene(int sub_div);

    void buildShelterAcutationScene(int sub_div);

    void buildRandomPatchScene(int sub_div);

    void buildSaddleScene(int sub_div);
    void buildGripperScene(int sub_div);
    void buildGridLayoutGripper(int sub_div);
    void buildSquareCrossJointScene(int sub_div);

    void buildInterlockingSquareScene(int sub_div);
    void buildDenseInterlockingSquarePeriodicScene(int sub_div);

    void buildDenseInterlockingSquareScene(int sub_div);
    
    void buildActiveTextileScene(int sub_div);

    void buildTestSceneJuan(int sub_div);

    void buildDomeScene(int sub_div);
    
    
    void buildTestJoint(int sub_div);
    void buildXJointsScene(int sub_div);
    void buildXJointsScene2(int sub_div);

    void buildShelterScene(int sub_div);
    void buildGridScene2(int sub_div);

    void build3DtestScene(int sub_div);
    void buildOneCrossScene(int sub_div);
    void buildGridScene(int sub_div);

    void buildOmegaScene(int sub_div);
    void buildStraightRodScene(int sub_div);

private:
    
    void cropTranslationalUnitByparallelogram(const std::vector<std::vector<TV2>>& input_points,
    std::vector<TV2>& output_points, const TV2& top_left, const TV2& top_right,
    const TV2& bottom_right, const TV2& bottom_left, std::vector<Vector<int, 2>>& edge_pairs,
    std::unordered_map<int, std::vector<int>>& crossing_tracker,
    std::vector<std::vector<Vector<int, 2>>>& boundary_pairs,
    std::vector<std::vector<int>>& boundary_pair_rod_idx);

    void clearSimData();
    
    void appendThetaAndJointDoF(std::vector<Entry>& w_entry, 
        int& full_dof_cnt, int& dof_cnt);
    
    void addAStraightRod(const TV& from, const TV& to, 
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div,
        int& full_dof_cnt, int& node_cnt, int& rod_cnt);
    
    void addCurvedRod(const std::vector<TV2>& data_points,
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div, int& full_dof_cnt, int& node_cnt, int& rod_cnt, bool closed);

    void addStraightYarnCrossNPoints(const TV& from, const TV& to,
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div,
        std::vector<TV>& sub_points, std::vector<int>& node_idx,
        std::vector<int>& key_points_location,
        int start, bool pbc = false);


    void markCrossingDoF(
        std::vector<Eigen::Triplet<T>>& w_entry,
        int& dof_cnt);
    
    void addPoint(const TV& point, int& full_dof_cnt, int& node_cnt);

    void addCrossingPoint(std::vector<TV>& existing_nodes, 
        const TV& point, int& full_dof_cnt, int& node_cnt);
};

#endif