#include "../include/EoLRodSim.h"

#include "../autodiff/EoLRodBendingAndTwistingPBC.h"


template<class T, int dim>
void EoLRodSim<T, dim>::buildPBCData(std::vector<TV>& data, 
    const std::vector<Offset>& offsets, 
    const std::vector<int>& rod_id,
    std::vector<TV>& dXdu, std::vector<TV>& d2Xdu2,
    bool g, bool f)
{
    data.resize(8); dXdu.resize(4); d2Xdu2.resize(4);

    for (int i = 0; i < 4; i++)
    {
        data[i] = deformed_states.template segment<dim>(offsets[i][0]);
        TV pos, dpos, ddpos;
        T u = deformed_states[offsets[i][dim]];
        Rods[rod_id[i]]->rest_state->getMaterialPos(u, pos, dpos, ddpos, g, f);
        data[4 + i] = pos;
        dXdu[i] = dpos;
        d2Xdu2[i] = ddpos;
    }
    
}

template<class T, int dim>
T EoLRodSim<T, dim>::add3DPBCBendingAndTwistingEnergy(bool bending, bool twisting)
{
    T energy = 0.0;
    if constexpr (dim == 3)
    {
        iteratePBCBendingPairs([&](const std::vector<Offset>& offsets, 
            const std::vector<int>& rod_id)
        {
            std::vector<TV> data, dXdu, d2Xdu2;
            buildPBCData(data, offsets, rod_id, dXdu, d2Xdu2, false, false);
            
            TV referenceNormal1 = Rods[rod_id.front()]->reference_frame_us.front();
            TV referenceNormal2 = Rods[rod_id.back()]->reference_frame_us.back();

            TV referenceTangent1 = Rods[rod_id.front()]->prev_tangents.front();
            TV referenceTangent2 = Rods[rod_id.back()]->prev_tangents.back();

            T theta0 = Rods[rod_id.front()]->reference_angles[0];
            T theta1 = Rods[rod_id.back()]->reference_angles[int(Rods[rod_id.back()]->indices.size()) - 2];
            
            T reference_twist = Rods[rod_id.front()]->reference_twist[0];
            
            Matrix<T, 2, 2> bending_coeffs = bending ? Rods[rod_id.front()]->bending_coeffs : Matrix<T, 2, 2>::Zero();
            T kt = twisting ? Rods[rod_id.front()]->kt : 0.0;
            energy += computeRodBendingAndTwistPBCEnergy(
                bending_coeffs, 
                kt, 0.0, 
                referenceNormal1, referenceTangent1, 
                referenceNormal2, referenceTangent2,
                reference_twist, data[0], data[1], data[2], data[3],
                data[4], data[5], data[6], data[7], theta0, theta1
            );
        });
    }
    return energy;
}

template<class T, int dim>
void EoLRodSim<T, dim>::add3DPBCBendingAndTwistingForce(Eigen::Ref<VectorXT> residual, bool bending, bool twisting)
{
    VectorXT residual_cp = residual;
    if constexpr (dim == 3)
    {
        iteratePBCBendingPairs([&](const std::vector<Offset>& offsets, 
            const std::vector<int>& rod_id)
        {
            std::vector<TV> data, dXdu, d2Xdu2;
            buildPBCData(data, offsets, rod_id, dXdu, d2Xdu2, true, false);
            
            TV referenceNormal1 = Rods[rod_id.front()]->reference_frame_us.front();
            TV referenceNormal2 = Rods[rod_id.back()]->reference_frame_us.back();

            TV referenceTangent1 = Rods[rod_id.front()]->prev_tangents.front();
            TV referenceTangent2 = Rods[rod_id.back()]->prev_tangents.back();

        

            T theta0 = Rods[rod_id.front()]->reference_angles[0];
            T theta1 = Rods[rod_id.back()]->reference_angles[int(Rods[rod_id.back()]->indices.size()) - 2];
            
            T reference_twist = Rods[rod_id.front()]->reference_twist[0];

            Matrix<T, 2, 2> bending_coeffs = bending ? Rods[rod_id.front()]->bending_coeffs : Matrix<T, 2, 2>::Zero();
            T kt = twisting ? Rods[rod_id.front()]->kt : 0.0;
            
            Vector<T, 26> F;
            F.setZero();
            computeRodBendingAndTwistPBCEnergyGradient(
                bending_coeffs, 
                kt, 0.0, 
                referenceNormal1, referenceTangent1, 
                referenceNormal2, referenceTangent2,
                reference_twist, data[0], data[1], data[2], data[3],
                data[4], data[5], data[6], data[7], theta0, theta1, F
            );
            // std::cout << "F" << std::endl;
            // std::cout << F << std::endl;
            // std::cout << "x" << std::endl;
            // // for(TV& x : data)
            // //     std::cout << x.transpose() << std::endl;
            // std::cout << (data[1] - data[0]).normalized().transpose() << std::endl;
            // std::cout << (data[3] - data[2]).normalized().transpose() << std::endl;
            // std::cout << referenceTangent1.transpose() << " " << referenceTangent2.transpose() << std::endl;
            // std::getchar();

            for (int i = 0; i < 4; i++)
            {
                residual.template segment<dim>(offsets[i][0]) += -F.template segment<dim>(i * dim);
                residual[offsets[i][dim]] += -F.template segment<dim>(4*dim + i * dim).dot(dXdu[i]);
            }
            residual[Rods[rod_id.front()]->theta_dof_start_offset] += -F[24];
            residual[Rods[rod_id.back()]->theta_dof_start_offset + 
                int(Rods[rod_id.back()]->indices.size()) - 2] += -F[25];
            
        });
    }
    if(print_force_mag)
        std::cout << "pbc bending and twist norm: " << (residual - residual_cp).norm() << std::endl;
}

template<class T, int dim>
void EoLRodSim<T, dim>::add3DPBCBendingAndTwistingK(std::vector<Entry>& entry_K, bool bending, bool twisting)
{
    if constexpr (dim == 3)
    {
        iteratePBCBendingPairs([&](const std::vector<Offset>& offsets, 
            const std::vector<int>& rod_id)
        {
            std::vector<TV> data, dXdu, d2Xdu2;
            buildPBCData(data, offsets, rod_id, dXdu, d2Xdu2, true, false);
            
            TV referenceNormal1 = Rods[rod_id.front()]->reference_frame_us.front();
            TV referenceNormal2 = Rods[rod_id.back()]->reference_frame_us.back();

            TV referenceTangent1 = Rods[rod_id.front()]->prev_tangents.front();
            TV referenceTangent2 = Rods[rod_id.back()]->prev_tangents.back();

            T theta0 = Rods[rod_id.front()]->reference_angles[0];
            T theta1 = Rods[rod_id.back()]->reference_angles[int(Rods[rod_id.back()]->indices.size()) - 2];
            
            T reference_twist = Rods[rod_id.front()]->reference_twist[0];

            Matrix<T, 2, 2> bending_coeffs = bending ? Rods[rod_id.front()]->bending_coeffs : Matrix<T, 2, 2>::Zero();
            T kt = twisting ? Rods[rod_id.front()]->kt : 0.0;
            
            Matrix<T, 26, 26> J;
            J.setZero();
            computeRodBendingAndTwistPBCEnergyHessian(
                bending_coeffs, 
                kt, 0.0, 
                referenceNormal1, referenceTangent1, 
                referenceNormal2, referenceTangent2,
                reference_twist, data[0], data[1], data[2], data[3],
                data[4], data[5], data[6], data[7], theta0, theta1, J
            );

            std::vector<int> theta_dofs = {
                Rods[rod_id.front()]->theta_dof_start_offset, 
                Rods[rod_id.back()]->theta_dof_start_offset + int(Rods[rod_id.back()]->indices.size()) - 2};

            for(int k = 0; k < offsets.size(); k++)
                for(int l = 0; l < offsets.size(); l++)
                    for(int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                        {
                            entry_K.push_back(Entry(offsets[k][i], offsets[l][j], J(k*dim + i, l * dim + j)));
                                
                                entry_K.push_back(Entry(offsets[k][i], offsets[l][dim], J(k*dim + i, 4 * dim + l * dim + j) * dXdu[l][j]));

                                entry_K.push_back(Entry(offsets[k][dim], offsets[l][j] , J(4 * dim + k * dim + i, l * dim + j) * dXdu[k][i]));

                                // dX/du ^T d2e/dX2 dXdu                    
                                entry_K.push_back(Entry(offsets[k][dim], 
                                                        offsets[l][dim], 
                                                        J(4 * dim + k*dim + i, 4 * dim + l * dim + j) * dXdu[l][j] * dXdu[k][i]));
                        }
            for(int k = 0; k < offsets.size(); k++)
               for (int j = 0; j < 2; j++)
                {
                    for(int i = 0; i < dim; i++)
                    {
                        entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], theta_dofs[j], J(k*dim + i, 24 + j)));
                        entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][i], J(24 + j, k * dim + i)));

                        entry_K.push_back(Eigen::Triplet<T>(offsets[k][dim], theta_dofs[j], J(4 * dim + k * dim + i, 24 + j) * dXdu[k][i]));
                        entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][dim], J(24 + j, 4 * dim + k * dim + i) * dXdu[k][i]));
                    }
                }
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        entry_K.push_back(Eigen::Triplet<T>(theta_dofs[i], theta_dofs[j], J(24 + i, 24 + j)));

        });
    }
}

template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;