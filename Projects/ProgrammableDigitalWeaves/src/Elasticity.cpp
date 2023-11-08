#include "../include/EoLRodSim.h"

template<class T, int dim>
void EoLRodSim<T, dim>::setUniaxialStrain(T theta, T s, TV& strain_dir, TV& ortho_dir)
{
    this->theta = theta;
    pbc_strain_data.clear();
    pbc_strain_data.resize(0);
    strain_dir = TV::Zero();
    if constexpr (dim == 2)
    {
        strain_dir = TV(std::cos(theta), std::sin(theta));
        ortho_dir = TV(-std::sin(theta), std::cos(theta));
    }
    else if constexpr (dim == 3)
    {
        strain_dir = TV(std::cos(theta), std::sin(theta), 0.0);
        ortho_dir = TV(-std::sin(theta), std::cos(theta), 0.0);
    }

    iteratePBCPairs([&](Offset offset_ref_i, Offset offset_ref_j, Offset offset_i, Offset offset_j)
    {
        
        TV Xj = rest_states.template segment<dim>(offset_j[0]);
        TV Xi = rest_states.template segment<dim>(offset_i[0]);
        T Dij = (Xj - Xi).dot(strain_dir);
        T dij = Dij * s;
        // std::cout << Dij << " " << dij << std::endl;
        pbc_strain_data.push_back(std::make_pair(std::make_pair(offset_i, offset_j), std::make_pair(strain_dir, dij)));
    });

}

template<class T, int dim>
void EoLRodSim<T, dim>::setBiaxialStrain(T theta1, T s1, T theta2, T s2, TV& strain_dir, TV& ortho_dir)
{
    // pbc_strain_data.clear();
    // pbc_strain_data.resize(0);
    // if constexpr (dim == 2)
    // {
    //     strain_dir = TV(std::cos(theta1), std::sin(theta1));
    //     ortho_dir = TV(-std::sin(theta1), std::cos(theta1));
    // }
    // iteratePBCReferencePairs([&](int dir_id, int node_i, int node_j){
    //     TV Xj = q0.col(node_j).template segment<dim>(0);
    //     TV Xi = q0.col(node_i).template segment<dim>(0);
    //     if constexpr (dim == 2)
    //     {
    //         T Dij = (Xj - Xi).dot(strain_dir);
    //         T dij = Dij * s1;
    //         // std::cout << Dij << " " << dij << std::endl;
    //         pbc_strain_data.push_back(std::make_pair(IV2(node_i, node_j), std::make_pair(strain_dir, dij)));

    //         Dij = (Xj - Xi).dot(ortho_dir);
    //         dij = Dij * s2;
    //         // std::cout << Dij << " " << dij << std::endl;
    //         pbc_strain_data.push_back(std::make_pair(IV2(node_i, node_j), std::make_pair(ortho_dir, dij)));
    //     }
    // });
}

template<class T, int dim>
void EoLRodSim<T, dim>::setBiaxialStrainWeighted(T theta1, T s1, T theta2, T s2, T w)
{
    // pbc_strain_data.clear();
    // pbc_strain_data.resize(0);
    // TV strain_dir1 = TV::Zero();
    // TV strain_dir2 = TV::Zero();
    // if constexpr (dim == 2)
    // {
    //     strain_dir1 = TV(std::cos(theta1), std::sin(theta1));
    //     strain_dir2 = TV(-std::sin(theta1), std::cos(theta1));
    // }
    // iteratePBCReferencePairs([&](int dir_id, int node_i, int node_j){
    //     TV Xj = q0.col(node_j).template segment<dim>(0);
    //     TV Xi = q0.col(node_i).template segment<dim>(0);
    //     if constexpr (dim == 2)
    //     {
    //         T Dij = (Xj - Xi).dot(strain_dir1);
    //         T dij = Dij * s1;
    //         // std::cout << Dij << " " << dij << std::endl;
    //         pbc_strain_data.push_back(std::make_pair(IV2(node_i, node_j), std::make_pair(strain_dir1, dij)));

    //         Dij = (Xj - Xi).dot(strain_dir2);
    //         dij = Dij * s2;
    //         // std::cout << Dij << " " << dij << std::endl;
    //         pbc_strain_data.push_back(std::make_pair(IV2(node_i, node_j), std::make_pair(strain_dir2, dij)));
    //     }
    // });
}

template<class T, int dim>
void EoLRodSim<T, dim>::computeDeformationGradientUnitCell()
{
    // IV2 ref0 = pbc_ref_unique[0];
    // IV2 ref1 = pbc_ref_unique[1];
    
    // TM x, X;
    // X.col(0) = q0.col(ref0[1]).template segment<dim>(0) - q0.col(ref0[0]).template segment<dim>(0);
    // X.col(1) = q0.col(ref1[1]).template segment<dim>(0) - q0.col(ref1[0]).template segment<dim>(0);
    // x.col(0) = q.col(ref0[1]).template segment<dim>(0) - q.col(ref0[0]).template segment<dim>(0);
    // x.col(1) = q.col(ref1[1]).template segment<dim>(0) - q.col(ref1[0]).template segment<dim>(0);

    // TM F = x * X.inverse();
    
    // std::cout << "F from ref pairs: " << F << std::endl;
    // // std::cout << "F-FT: " << F - F.transpose() << std::endl;
}


// min_F 1/2||F(Xi-Xj) - (xi-xj)||^2
// least square deformation fitting
template<class T, int dim>
void EoLRodSim<T, dim>::fitDeformationGradientUnitCell()
{
    auto computeEnergy = [&](TM _F){
        VectorXT energy(n_rods);
        energy.setZero();
        tbb::parallel_for(0, n_rods, [&](int i){
            TV xi = q.col(rods(0, i)).template segment<dim>(0);
            TV xj = q.col(rods(1, i)).template segment<dim>(0);
            TV Xi = q0.col(rods(0, i)).template segment<dim>(0);
            TV Xj = q0.col(rods(1, i)).template segment<dim>(0);
            energy[i] += 0.5 * (_F * (Xi - Xj) - (xi - xj)).squaredNorm();
        });
        return energy.sum();
    };

    auto computeGradient = [&](TM _F){
        TM dedF = TM::Zero();
        for (int i = 0; i < n_rods; i++)
        {
            TV xi = q.col(rods(0, i)).template segment<dim>(0);
            TV xj = q.col(rods(1, i)).template segment<dim>(0);
            TV Xi = q0.col(rods(0, i)).template segment<dim>(0);
            TV Xj = q0.col(rods(1, i)).template segment<dim>(0);
            dedF += -(_F * (Xi - Xj) - (xi - xj)) * (Xi - Xj).transpose();
        }
        return dedF;
    };


    auto polarDecomposition = [&](TM F, TM& R, TM& RU)
    {
        Eigen::JacobiSVD<TM> svd;
        svd.compute(F, Eigen::ComputeFullU | Eigen::ComputeFullV );
        TM U = svd.matrixU();
        TM V = svd.matrixV();
        TV S = svd.singularValues();
        R = U*V.transpose();
        const auto& SVT = S.asDiagonal() * V.adjoint();
        // from libigl
        if(R.determinant() < 0)
        {
            auto W = V.eval();
            W.col(V.cols()-1) *= -1.;
            R = U*W.transpose();
            RU = W*SVT;
        }
        else
            RU = V*SVT;      
    };


    TM F = TM::Identity();

    while (true)
    {
        TM dedF = computeGradient(F);
        if (dedF.norm() < 1e-5)
            break;
        T E0 = computeEnergy(F);
        T alpha = 1.0;
        while (true)
        {
            TM F_ls = F + alpha * dedF;
            T E1 = computeEnergy(F_ls);
            if (E1 - E0 < 0)
            {
                F = F_ls;
                break;
            }
            alpha *= 0.5;
        }
    }

    TM R, RU;
    polarDecomposition(F, R, RU);

    std::cout << "F fited " <<  F << std::endl;
    // std::cout << "F-F^T: " << F - F.transpose() << std::endl;
    // std::cout << "R: " << R << std::endl;
}

template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;