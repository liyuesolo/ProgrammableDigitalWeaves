#include "../include/EoLRodSim.h"
#include "../include/UnitPatch.h"

template<class T, int dim>
void EoLRodSim<T, dim>::buildPeriodicNetwork(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
    Eigen::MatrixXd& C, bool show_rest)
{
    V.resize(0, 0); F.resize(0, 0);
    int n_tile_x = 2, n_tile_y = 2;
    int n_faces = 20;

    if constexpr (dim == 3)
    {
        TV xi = deformed_states.template segment<dim>(pbc_pairs_reference[0].first.first[0]);
        TV xj = deformed_states.template segment<dim>(pbc_pairs_reference[0].first.second[0]);
        TV xk = deformed_states.template segment<dim>(pbc_pairs_reference[1].first.first[0]);
        TV xl = deformed_states.template segment<dim>(pbc_pairs_reference[1].first.second[0]);

        TV ref0_shift = xj - xi, ref1_shift = xl - xk;
        // for (int i = -n_tile_x + 1; i < n_tile_x; i++)
        for (int i = 0; i < n_tile_x; i++)
        {
            for (int j = 0; j < n_tile_y; j++)
            {
                Eigen::MatrixXd V_tile;
                Eigen::MatrixXi F_tile;
                
                generateMeshForRendering(V_tile, F_tile, (T(i) * ref0_shift + T(j) * ref1_shift), show_rest);
                for (int row = 0; row < F_tile.rows(); row++)
                {
                    for (int d = 0; d < dim; d++)
                    {
                        F_tile(row, d) += V.rows();
                    }
                }
                // std::cout << "V " << V_tile.row(0) << " F " << F_tile.row(0) << std::endl;
                
                int v_row_current = V.rows();
                // std::cout << "V rows " << V.rows() << " V_tile rows " << V_tile.rows() << std::endl;
                V.conservativeResize(v_row_current + V_tile.rows(), 3);
                V.block(v_row_current, 0, V_tile.rows(), 3) = V_tile;
                
                int f_row_current = F.rows();
                // std::cout << "F rows " << F.rows() << " F_tile rows " << F_tile.rows() << std::endl;
                F.conservativeResize(f_row_current+ F_tile.rows(), 3);
                F.block(f_row_current, 0, F_tile.rows(), 3) = F_tile;
            }
        }
        
    }
    
    C.resize(F.rows(), 3);
    tbb::parallel_for(0, int(F.rows()), [&](int i){
        // if( i % 2 == 0)
        //     C.row(i) = Eigen::Vector3d(0, 1, 0);
        // else
        //     C.row(i) = Eigen::Vector3d(0, 0, 1);
        C.row(i) = Eigen::Vector3d(0, 0.3, 1);
    });
    
    // std::cout << V << std::endl;
    // std::cout << F << std::endl;
    // std::cout << C << std::endl;
}

template<class T, int dim>
void EoLRodSim<T, dim>::buildSceneFromUnitPatch(int patch_id)
{
    UnitPatch<T, dim> unit_patch(*this);
    unit_patch.buildScene(patch_id);
}

template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;

//not working for float yet
// template class EoLRodSim<float>;