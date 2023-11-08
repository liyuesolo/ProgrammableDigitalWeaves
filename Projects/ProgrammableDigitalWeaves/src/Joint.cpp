#include "../include/EoLRodSim.h"
// #include "JointBendingAndTwisting.h"
// #include "EoLRodQuaternionBendingAndTwisting.h"
#include "../autodiff/EoLRodEulerAngleBendingAndTwisting.h"
#include "../autodiff/EoLRodBendingAndTwisting.h"

template<class T, int dim>
void EoLRodSim<T, dim>::addJointBendingAndTwistingK(std::vector<Entry>& entry_K)
{
    if constexpr (dim == 3)
    {
        for (auto& crossing : rod_crossings)
        {
            int node_i = crossing->node_idx;
            
            std::vector<int> rods_involved = crossing->rods_involved;
            std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

            std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
            Vector<T, 3> omega = crossing->omega;
            Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;

            int cnt = 0;
            for (int rod_idx : rods_involved)
            {
                Offset offset_i, offset_j, offset_k;
                TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
                
                auto rod = Rods[rod_idx];
                rod->getEntry(node_i, offset_i);
                rod->x(node_i, xi); rod->XdX(node_i, Xi, dXi);

                int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
                int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
                int edge1 = 0, edge0 = 0;
                if (node_j != -1)
                {
                    rod->getEntry(node_j, offset_j);
                    rod->x(node_j, xj);
                    rod->XdX(node_j, Xj, dXj);
                    edge1 = on_rod_idx[rod_idx];
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge1 = 0;
                }
                if (node_k != -1)
                {
                    rod->getEntry(node_k, offset_k);
                    rod->x(node_k, xk);
                    rod->XdX(node_k, Xk, dXk);
                    edge0 = on_rod_idx[rod_idx]-1;
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge0 = rod->numSeg() - 1;
                }

                T theta_edge0 = rod->reference_angles[edge0];
                T theta_edge1 = rod->reference_angles[edge1];

                Matrix<T, 2, 2> B = rod->bending_coeffs ;
                T kt = rod->kt;
                T reference_twist_edge0 = rod->reference_twist[edge0];
                T reference_twist_edge1 = rod->reference_twist[edge1];

                TV referenceNormal1 = rod->reference_frame_us[edge0];
                TV referenceNormal2 = rod->reference_frame_us[edge1];

                TV referenceTangent1 = rod->prev_tangents[edge0];
                TV referenceTangent2 = rod->prev_tangents[edge1];            
                
                TV rest_tangent1 = rod->rest_tangents[edge0];
                TV rest_tangent2 = rod->rest_tangents[edge1];
                TV rest_normal1 = rod->rest_normals[edge0];
                TV rest_normal2 = rod->rest_normals[edge1];

                if (crossing->is_fixed)
                {
                    Matrix<T, 16, 16> J0, J1;
                    J0.setZero(); J1.setZero();
                    // xi is the rigid body
                    
                    int dof_theta0 = rod->theta_dof_start_offset + edge0;
                    int dof_theta1 = rod->theta_dof_start_offset + edge1;

                    Offset omega_dof;
                    omega_dof[0] = crossing->dof_offset;
                    omega_dof[1] = crossing->dof_offset + 1;
                    omega_dof[2] = crossing->dof_offset + 2;

                    std::vector<Offset> offsets1 = { offset_i, offset_j };
                    std::vector<Offset> offsets2 = { offset_i, offset_k };
                    std::vector<TV> dXdu1 = { dXi, dXj };
                    std::vector<TV> dXdu2 = { dXi, dXk };
                    
                    int entry_cnt = 0;
                    if (node_j != -1)
                    {
                        computeRodEulerAngleBendingAndTwistEnergyRBFirstHessian(B, kt, 0.0, 
                            rest_tangent1, rest_normal1, 
                            referenceTangent2, referenceNormal2, 
                            reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1, J0);
                        
                        J0 *= 0.5;
                        for(int k = 0; k < offsets1.size(); k++)
                        {
                            //dx/dx
                            //dx/dX dX/dx
                            for(int l = 0; l < offsets1.size(); l++)
                                for (int i = 0; i < dim; i++)
                                    for (int j = 0; j < dim; j++)
                                    {
                                        entry_K.push_back(Eigen::Triplet<T>(offsets1[k][i], offsets1[l][j], J0(k*dim + i, l * dim + j)));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets1[k][i], offsets1[l][dim], J0(k*dim + i, 2 * dim + l * dim + j) * dXdu1[l][j]));
                                        // entry_K.push_back(Eigen::Triplet<T>(offsets1[l][dim], offsets1[k][i], J0(2 * dim + l * dim + j, k*dim + i) * dXdu1[l][j]));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets1[k][dim], offsets1[l][j], J0(2 * dim + k * dim + i, l*dim + j) * dXdu1[k][i]));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets1[k][dim], offsets1[l][dim], J0(2 * dim + k*dim + i, 2 * dim + l * dim + j) * dXdu1[l][j] * dXdu1[k][i]));
                                        entry_cnt+=4;
                                    }
                            for (int i = 0; i < dim; i++)
                            {
                                for (int j = 0; j < dim; j++)
                                {
                                    entry_K.push_back(Eigen::Triplet<T>(offsets1[k][i], omega_dof[j], J0(k*dim + i, 4 * dim + j)));
                                    entry_K.push_back(Eigen::Triplet<T>(omega_dof[j], offsets1[k][i], J0(4 * dim + j, k*dim + i)));

                                    entry_K.push_back(Eigen::Triplet<T>(offsets1[k][dim], omega_dof[j], J0(2 *dim + k*dim + i, 4 * dim + j) * dXdu1[k][i]));
                                    entry_K.push_back(Eigen::Triplet<T>(omega_dof[j], offsets1[k][dim], J0(4 * dim + j, 2 *dim + k*dim + i) * dXdu1[k][i]));
                                    entry_cnt+=4;
                                }

                                entry_K.push_back(Eigen::Triplet<T>(offsets1[k][i], dof_theta1, J0(k*dim + i, 5 * dim)));
                                entry_K.push_back(Eigen::Triplet<T>(dof_theta1, offsets1[k][i], J0(5 * dim, k*dim + i)));

                                entry_K.push_back(Eigen::Triplet<T>(offsets1[k][dim], dof_theta1, J0(2 * dim + k*dim + i, 5 * dim) * dXdu1[k][i]));
                                entry_K.push_back(Eigen::Triplet<T>(dof_theta1, offsets1[k][dim], J0(5 * dim, 2 * dim + k*dim + i) * dXdu1[k][i]));
                                entry_cnt+=4;
                            }

                        }
                        for (int i = 0; i < dim; i++)
                        {
                            entry_K.push_back(Eigen::Triplet<T>(omega_dof[i], dof_theta1, J0(4 * dim + i, 5 * dim)));
                            entry_K.push_back(Eigen::Triplet<T>(dof_theta1, omega_dof[i], J0(5 * dim, 4 * dim + i)));
                            entry_cnt+=2;
                            for (int j = 0; j < dim; j++)
                            {
                                entry_K.push_back(Eigen::Triplet<T>(omega_dof[i], omega_dof[j], J0(4 * dim + i, 4 * dim + j)));
                                entry_cnt+=1;
                            }
                        }
                        entry_K.push_back(Eigen::Triplet<T>(dof_theta1, dof_theta1, J0(5 * dim, 5 * dim)));
                        entry_cnt+=1;

                        // std::cout << 16 * 16 << " " << entry_cnt << std::endl;
                        // std::getchar();
                    }
                    if (node_k != -1)
                    {
                        computeRodEulerAngleBendingAndTwistEnergyRBSecondHessian(B, kt, 0.0, 
                            referenceTangent1, referenceNormal1,
                            rest_tangent2, rest_normal2, 
                            reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0, J1);
                        
                        J1 *= 0.5;
                        for(int k = 0; k < offsets1.size(); k++)
                        {
                            //dx/dx
                            //dx/dX dX/dx
                            for(int l = 0; l < offsets1.size(); l++)
                                for (int i = 0; i < dim; i++)
                                    for (int j = 0; j < dim; j++)
                                    {
                                        entry_K.push_back(Eigen::Triplet<T>(offsets2[k][i], offsets2[l][j], J1(k*dim + i, l * dim + j)));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets2[k][i], offsets2[l][dim], J1(k*dim + i, 2 * dim + l * dim + j) * dXdu2[l][j]));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets2[l][dim], offsets2[k][i], J1(2 * dim + l * dim + j, k*dim + i) * dXdu2[l][j]));
                                        entry_K.push_back(Eigen::Triplet<T>(offsets2[k][dim], offsets2[l][dim], J1(2 * dim + k*dim + i, 2 * dim + l * dim + j) * dXdu2[l][j] * dXdu2[k][i]));
                                    }
                            for (int i = 0; i < dim; i++)
                            {
                                for (int j = 0; j < dim; j++)
                                {
                                    entry_K.push_back(Eigen::Triplet<T>(offsets2[k][i], omega_dof[j], J1(k*dim + i, 4 * dim + j)));
                                    entry_K.push_back(Eigen::Triplet<T>(omega_dof[j], offsets2[k][i], J1(4 * dim + j, k*dim + i)));
                                    entry_K.push_back(Eigen::Triplet<T>(offsets2[k][dim], omega_dof[j], J1(2 *dim + k*dim + i, 4 * dim + j) * dXdu2[k][i]));
                                    entry_K.push_back(Eigen::Triplet<T>(omega_dof[j], offsets2[k][dim], J1(4 * dim + j, 2 *dim + k*dim + i) * dXdu2[k][i]));
                                }

                                entry_K.push_back(Eigen::Triplet<T>(dof_theta0, offsets2[k][i], J1(5 * dim, k*dim + i)));
                                entry_K.push_back(Eigen::Triplet<T>(offsets2[k][i], dof_theta0, J1(k*dim + i, 5 * dim)));

                                entry_K.push_back(Eigen::Triplet<T>(offsets2[k][dim], dof_theta0, J1(2 * dim + k*dim + i, 5 * dim) * dXdu2[k][i]));
                                entry_K.push_back(Eigen::Triplet<T>(dof_theta0, offsets2[k][dim], J1(5 * dim, 2 * dim + k*dim + i) * dXdu2[k][i]));
                            }

                        }
                        for (int i = 0; i < dim; i++)
                        {
                            entry_K.push_back(Eigen::Triplet<T>(omega_dof[i], dof_theta0, J1(4 * dim + i, 5 * dim)));
                            entry_K.push_back(Eigen::Triplet<T>(dof_theta0, omega_dof[i], J1(5 * dim, 4 * dim + i)));

                            for (int j = 0; j < dim; j++)
                            {
                                entry_K.push_back(Eigen::Triplet<T>(omega_dof[i], omega_dof[j], J1(4 * dim + i, 4 * dim + j)));
                            }
                        }
                        entry_K.push_back(Eigen::Triplet<T>(dof_theta0, dof_theta0, J1(5 * dim, 5 * dim)));
                    }
                    // std::cout << "here here" << std::endl;
                    // std::getchar();
                }
                else
                {
                    if (node_k != -1 && node_j != -1)
                    {
                        // std::cout << node_i << std::endl;
                        std::vector<int> nodes = { node_k, node_i, node_j };
                        std::vector<TV> dXdu = { dXk, dXi, dXj };
                        std::vector<Offset> offsets = { offset_k, offset_i, offset_j };

                        
                        std::vector<int> theta_dofs = {rod->theta_dof_start_offset + edge0, rod->theta_dof_start_offset + edge1};

                        Matrix<T, 20, 20> J;
                    
                        computeRodBendingAndTwistEnergyHessian(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                            referenceNormal2, referenceTangent2,  reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1, J);

                        J *= -1.0;
                        for(int k = 0; k < nodes.size(); k++)
                            for(int l = 0; l < nodes.size(); l++)
                                for(int i = 0; i < dim; i++)
                                    for (int j = 0; j < dim; j++)
                                        {
                                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], offsets[l][j], -J(k*dim + i, l * dim + j)));

                                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], offsets[l][dim], -J(k*dim + i, 3 * dim + l * dim + j) * dXdu[l][j]));
                                            
                                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][dim], offsets[l][j], -J(3 * dim + k * dim + i, l * dim + j) * dXdu[k][i]));

                                            
                                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][dim], 
                                                                                offsets[l][dim], 
                                                                                -J(3 * dim + k*dim + i, 3 * dim + l * dim + j) * dXdu[l][j] * dXdu[k][i]));

                                        }
                        for(int k = 0; k < nodes.size(); k++)
                            for (int j = 0; j < 2; j++)
                                for(int i = 0; i < dim; i++)
                                {
                                    entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], theta_dofs[j], -J(k*dim + i, 18 + j)));
                                    entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][i], -J(18 + j, k * dim + i)));

                                    entry_K.push_back(Eigen::Triplet<T>(offsets[k][dim], theta_dofs[j], -J(3 * dim + k * dim + i, 18 + j) * dXdu[k][i]));
                                    entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][dim], -J(18 + j, 3 * dim + k * dim + i) * dXdu[k][i]));
                                }
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++)
                                entry_K.push_back(Eigen::Triplet<T>(theta_dofs[i], theta_dofs[j], -J(18 + i, 18 + j)));
                    }
                    
                }
                

                cnt ++;
            }
        }
    }
}

template<class T, int dim>
void EoLRodSim<T, dim>::addJointBendingAndTwistingForce(Eigen::Ref<VectorXT> residual)
{
    if constexpr (dim == 3)
    {

        VectorXT residual_cp = residual;
        for (auto& crossing : rod_crossings)
        {
            int node_i = crossing->node_idx;
            
            std::vector<int> rods_involved = crossing->rods_involved;
            std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

            std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
            Vector<T, 3> omega = crossing->omega;
            // Vector<T, 4> omega_acc = crossing->omega_acc;
            Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;
            // Vector<T, 3> omega_acc = crossing->omega_acc;

            int cnt = 0;
            for (int rod_idx : rods_involved)
            {
                Offset offset_i, offset_j, offset_k;
                TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
                
                auto rod = Rods[rod_idx];
                rod->getEntry(node_i, offset_i);
                rod->x(node_i, xi); rod->XdX(node_i, Xi, dXi);

                int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
                int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
                int edge1 = 0, edge0 = 0;
                if (node_j != -1)
                {
                    rod->getEntry(node_j, offset_j);
                    rod->x(node_j, xj);
                    rod->XdX(node_j, Xj, dXj);
                    edge1 = on_rod_idx[rod_idx];
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge1 = 0;
                }
                if (node_k != -1)
                {
                    rod->getEntry(node_k, offset_k);
                    rod->x(node_k, xk);
                    rod->XdX(node_k, Xk, dXk);
                    edge0 = on_rod_idx[rod_idx]-1;
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge0 = rod->numSeg() - 1;

                }

                // std::cout << " node i " << node_i << " node j " << node_j << " node k " << node_k << std::endl;
                // std::cout << " on_rod_idx[rod_idx] " << on_rod_idx[rod_idx] << " " << edge0 << " " << edge1 << std::endl;

                T theta_edge0 = rod->reference_angles[edge0];
                T theta_edge1 = rod->reference_angles[edge1];

                Matrix<T, 2, 2> B = rod->bending_coeffs ;
                T kt = rod->kt;
                T reference_twist_edge0 = rod->reference_twist[edge0];
                T reference_twist_edge1 = rod->reference_twist[edge1];

                TV referenceNormal1 = rod->reference_frame_us[edge0];
                TV referenceNormal2 = rod->reference_frame_us[edge1];

                TV referenceTangent1 = rod->prev_tangents[edge0];
                TV referenceTangent2 = rod->prev_tangents[edge1];       

                TV rest_tangent1 = rod->rest_tangents[edge0];
                TV rest_tangent2 = rod->rest_tangents[edge1];
                
                TV rest_normal1 = rod->rest_normals[edge0];
                TV rest_normal2 = rod->rest_normals[edge1];
                
                if (crossing->is_fixed)
                {
                    Vector<T, 16> F0, F1;
                    F0.setZero(); F1.setZero();
                    // xi is the rigid body
                    
                    if (node_j != -1)
                    {
                        
                        computeRodEulerAngleBendingAndTwistEnergyRBFirstGradient(B, kt, 0.0, 
                            rest_tangent1, rest_normal1, 
                            referenceTangent2, referenceNormal2, 
                            reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1, F0);

                        // std::cout << node_i << " " << node_j << " " << node_k << std::endl;
                        // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
                        // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
                        // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
                        // std::cout << "F0" << std::endl;
                        // std::cout << F0 << std::endl;
                        // std::getchar();

                        F0 *= -0.5;
                        residual.template segment<dim>(offset_i[0]) += F0.template segment<dim>(0);
                        residual.template segment<dim>(offset_j[0]) += F0.template segment<dim>(dim);
                        residual(offset_i[dim]) += F0.template segment<dim>(2*dim).dot(dXi);
                        residual(offset_j[dim]) += F0.template segment<dim>(3*dim).dot(dXj);
                        residual.template segment<dim>(crossing->dof_offset) += F0.template segment<dim>(4*dim);
                        residual[rod->theta_dof_start_offset + edge1] += F0[5*dim];
                    }
                    if (node_k != -1)                    
                    {
                        

                        computeRodEulerAngleBendingAndTwistEnergyRBSecondGradient(B, kt, 0.0, 
                            referenceTangent1, referenceNormal1,
                            rest_tangent2, rest_normal2, 
                            reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0, F1);

                        // std::cout << node_i << " " << node_k << " " << node_k << std::endl;
                        // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
                        // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
                        // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
                        // std::cout << "F1" << std::endl;
                        // std::cout << F1 << std::endl;
                        // std::getchar();
                        
                        F1 *= -0.5;
                        residual.template segment<dim>(offset_i[0]) += F1.template segment<dim>(0);
                        residual.template segment<dim>(offset_k[0]) += F1.template segment<dim>(dim);
                        residual(offset_i[dim]) += F1.template segment<dim>(2*dim).dot(dXi);
                        residual(offset_k[dim]) += F1.template segment<dim>(3*dim).dot(dXk);
                        residual.template segment<dim>(crossing->dof_offset) += F1.template segment<dim>(4*dim);
                        residual[rod->theta_dof_start_offset + edge0] += F1[5*dim];
                    }
                    
                }
                else
                {
                    if (node_k != -1 && node_j != -1)
                    {
                        Vector<T, 20> F;
                        computeRodBendingAndTwistEnergyGradient(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                            referenceNormal2, referenceTangent2,  reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1, F);
                        // std::cout << reference_twist_edge1 << " " << theta_edge0 << " " << theta_edge1 << std::endl;
                        // std::cout << node_i << " " << node_j << " " << node_k << std::endl;
                        // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
                        // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
                        // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
                        // std::cout << "F" << std::endl;
                        // std::cout << F << std::endl;
                        // std::getchar();

                        F *= -1.0;

                        residual.template segment<dim>(offset_k[0]) += F.template segment<dim>(0);
                        residual.template segment<dim>(offset_i[0]) += F.template segment<dim>(dim);
                        residual.template segment<dim>(offset_j[0]) += F.template segment<dim>(dim + dim);

                        residual(offset_k[dim]) += F.template segment<dim>(3*dim).dot(dXk);
                        residual(offset_i[dim]) += F.template segment<dim>(3*dim + dim).dot(dXi);
                        residual(offset_j[dim]) += F.template segment<dim>(3*dim + 2*dim).dot(dXj);

                        
                        residual[rod->theta_dof_start_offset + edge0] += F[18];
                        residual[rod->theta_dof_start_offset + edge1] += F[19];
                        
                    }
                }

                cnt ++;
            }
        }
        if (print_force_mag)
            std::cout << "joint force norm: " << (residual - residual_cp).norm() << std::endl;
    }
}

template<class T, int dim>
T EoLRodSim<T, dim>::addJointBendingAndTwistingEnergy(bool bending, bool twisting)
{
    T energy = 0.0;
    if constexpr (dim == 3)
    {

        for (auto& crossing : rod_crossings)
        {
            int node_i = crossing->node_idx;
            
            std::vector<int> rods_involved = crossing->rods_involved;
            std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

            std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
            Vector<T, 3> omega = crossing->omega;
            Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;
        
            
            int cnt = 0;
            for (int rod_idx : rods_involved)
            {
                Offset offset_i, offset_j, offset_k;
                TV xi = TV::Zero(), xj = TV::Zero(), xk = TV::Zero(), 
                    Xi = TV::Zero(), Xj = TV::Zero(), Xk = TV::Zero();

                auto rod = Rods[rod_idx];
                rod->getEntry(node_i, offset_i);
                rod->x(node_i, xi); rod->X(node_i, Xi);

                int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
                int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
                int edge1 = 0, edge0 = 0;
                
                if (node_j != -1)
                {
                    rod->getEntry(node_j, offset_j);
                    rod->x(node_j, xj);
                    rod->X(node_j, Xj);
                    edge1 = on_rod_idx[rod_idx];
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge1 = 0;
                }
                if (node_k != -1)
                {
                    rod->getEntry(node_k, offset_k);
                    rod->x(node_k, xk);
                    rod->X(node_k, Xk);
                    edge0 = on_rod_idx[rod_idx]-1;
                    if(rod->closed && on_rod_idx[rod_idx] == 0)
                        edge0 = rod->numSeg() - 1;
                }

                T theta_edge0 = rod->reference_angles[edge0];
                T theta_edge1 = rod->reference_angles[edge1];

                Matrix<T, 2, 2> B = bending ? rod->bending_coeffs : Matrix<T, 2, 2>::Zero();
                T kt = twisting ? rod->kt : 0.0;

                T reference_twist_edge0 = rod->reference_twist[edge0];
                T reference_twist_edge1 = rod->reference_twist[edge1];

                TV referenceNormal1 = rod->reference_frame_us[edge0];
                TV referenceNormal2 = rod->reference_frame_us[edge1];

                TV referenceTangent1 = rod->prev_tangents[edge0];
                TV referenceTangent2 = rod->prev_tangents[edge1];       

                TV rest_tangent1 = rod->rest_tangents[edge0];
                TV rest_tangent2 = rod->rest_tangents[edge1];
                TV rest_normal1 = rod->rest_normals[edge0];
                TV rest_normal2 = rod->rest_normals[edge1];

                if (crossing->is_fixed)
                {
                    // xi is the rigid body
                    if (node_j != -1)
                    {
                        T E = 0.5 * computeRodEulerAngleBendingAndTwistEnergyRBFirst(B, kt, 0.0, 
                            rest_tangent1, rest_normal1, 
                            referenceTangent2, referenceNormal2, 
                            reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1);
                        energy += E;
                        
                        // std::cout << E << std::endl;
                        // std::getchar();
                    }
                    if (node_k != -1)
                    {
                        T E = 0.5 * computeRodEulerAngleBendingAndTwistEnergyRBSecond(B, kt, 0.0, 
                            referenceTangent1, referenceNormal1,
                            rest_tangent2, rest_normal2, 
                            reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0);
                       energy += E;
                    }

                    
                }
                else
                {
                    if (node_j != -1 && node_k != -1)
                    {
                        T E = computeRodBendingAndTwistEnergy(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                            referenceNormal2, referenceTangent2, reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1);
                        
                        energy += E;
                    }
                }
                
                cnt++;
            }
        }
        return energy;
    }
    return 0;
}



template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;