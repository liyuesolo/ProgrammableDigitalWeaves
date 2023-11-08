#include "../include/EoLRodSim.h"


template<class T, int dim>
void EoLRodSim<T, dim>::derivativeTest()
{
    run_diff_test = true;
    // add_regularizor = false;
    add_stretching = true;
    // add_penalty = false;
    add_bending = true;
    // add_shearing = false;
    add_twisting = true;
    add_rigid_joint = true;
    add_pbc_bending = false;
    add_pbc_twisting = false;
    add_rotation_penalty = false;
    add_pbc = false;
    add_contact_penalty = false;
    add_eularian_reg = false;
    
    deformed_states /= unit;
    VectorXT dq(W.cols());
    // VectorXT dq(W.rows());
    dq.setRandom();
    dq *= 1.0 / dq.norm();
    for (int i = 0; i < dq.rows(); i++)
        dq(i) = dq(i) * 0.5 + 0.5;
    dq *= 0.01;
    // dq(2) += 1.0;
    // dq(3) += 1.0;
    testGradient(dq);
    testHessian(dq);
    testGradient2ndOrderTerm(dq);
    testHessian2ndOrderTerm(dq);
}


template<class T, int dim>
void EoLRodSim<T, dim>::checkMaterialPositionDerivatives()
{
    T epsilon = 1e-6;
    
    for (auto& rod : Rods)
    {
        rod->iterateSegmentsWithOffset([&](int node_i, int node_j, Offset offset_i, Offset offset_j, int rod_idx)
        {
            // std::cout << "node i " << node_i << " node j " << node_j << std::endl;
            TV X1, X0, dX0;
            rod->XdX(node_i, X0, dX0);
            deformed_states[offset_i[dim]] += epsilon;
            rod->X(node_i, X1);
            deformed_states[offset_i[dim]] -= epsilon;
            for (int d = 0; d < dim; d++)
            {
                std::cout << "dXdu: " << (X1[d] - X0[d]) / epsilon << " " << dX0[d] << std::endl;
            }
            // std::getchar();
        });
    }
}


template<class T, int dim>
void EoLRodSim<T, dim>::testGradient2ndOrderTerm(Eigen::Ref<VectorXT> dq)
{
    run_diff_test = true;
    std::cout << "======================== CHECK GRADIENT ========================" << std::endl;
    T epsilon = 1e-6;
    n_dof = W.cols();

    VectorXT gradient(n_dof);
    gradient.setZero();
    computeResidual(gradient, dq);
    
    gradient *= -1;
    T E0 = computeTotalEnergy(dq);
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.001;
    T previous = 0.0;
    for (int i = 0; i < 10; i++)
    {
        T E1 = computeTotalEnergy(dq + dx);
        T dE = E1 - E0;
        
        dE -= gradient.dot(dx);
        // std::cout << "dE " << dE << std::endl;
        if (i > 0)
        {
            std::cout << (previous/dE) << std::endl;
        }
        previous = dE;
        dx *= 0.5;
    }
    run_diff_test = false;
}

template<class T, int dim>
void EoLRodSim<T, dim>::testHessian2ndOrderTerm(Eigen::Ref<VectorXT> dq)
{
    run_diff_test = true;

    std::cout << "===================== check Hessian =====================" << std::endl;
    
    n_dof = W.cols();

    StiffnessMatrix A;
    buildSystemDoFMatrix(dq, A);

    VectorXT f0(n_dof);
    f0.setZero();
    computeResidual(f0, dq);
    f0 *= -1;
    
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    for(int i = 0; i < n_dof; i++) dx[i] += 0.5;
    dx *= 0.001;
    T previous = 0.0;
    for (int i = 0; i < 10; i++)
    {
        
        VectorXT f1(n_dof);
        f1.setZero();
        computeResidual(f1, dq + dx);
        f1 *= -1;
        T df_norm = (f0 + (A * dx) - f1).norm();
        std::cout << "df_norm " << df_norm << std::endl;
        if (i > 0)
        {
            std::cout << (previous/df_norm) << std::endl;
        }
        previous = df_norm;
        dx *= 0.5;
    }
    run_diff_test = false;
}

template<class T, int dim>
void EoLRodSim<T, dim>::testGradient(Eigen::Ref<VectorXT> dq)
{
    // run_diff_test = true;
    // add_regularizor = false;
    // add_stretching = false;
    // add_pbc = false;
    // add_pbc_bending = true;
    // add_bending = false;
    // add_twisting = false;
    // add_rigid_joint = false;
    // add_rotation_penalty = false;
    // add_pbc = false;

    std::cout << "======================== CHECK GRADIENT ========================" << std::endl;
    T epsilon = 1e-6;
    n_dof = W.cols();
    // n_dof = W.rows();
    int n_seg = 0;
    for (auto& rod : Rods)
        n_seg += rod->numSeg();
    int n_crossing = 0;
    for (auto& crossing : rod_crossings)
        n_crossing += 3;
    std::cout << "n crossing " << n_crossing << std::endl;
    std::cout << "n seg " << n_seg << std::endl;
    

    std::cout << W.rows() << " " << W.cols() << std::endl;
    VectorXT gradient(n_dof);
    gradient.setZero();

    computeResidual(gradient, dq);
    // std::cout << gradient << std::endl;
    VectorXT gradient_FD(n_dof);
    gradient_FD.setZero();
    
    // dq.setZero();

    int cnt = 0;
    for(int dof_i = 0; dof_i < n_dof; dof_i++)
    {
        dq(dof_i) += epsilon;
        // std::cout << W * dq << std::endl;
        T E0 = computeTotalEnergy(dq);
        
        dq(dof_i) -= 1.0 * epsilon;
        T E1 = computeTotalEnergy(dq);
        // std::cout << "E1 " << E1 << " E0 " << E0 << std::endl;
        gradient_FD(dof_i) = (E1 - E0) / (1*epsilon);
        if( gradient_FD(dof_i) == 0 && gradient(dof_i) == 0)
            continue;
        // if (std::abs( gradient_FD(d, n_node) - gradient(d, n_node)) < 1e-4)
        //     continue;
        std::cout << " dof " << dof_i << " " << gradient_FD(dof_i) << " " << gradient(dof_i) << std::endl;
        std::getchar();
        cnt++;   
    }
    //  run_diff_test = false;
    // add_regularizor = true;
    // add_stretching = true;
    // add_bending = true;
    // add_twisting = true;
}

template<class T, int dim>
void EoLRodSim<T, dim>::testHessian(Eigen::Ref<VectorXT> dq)
{
    

    // run_diff_test = true;
    // add_regularizor = false;
    // add_stretching = false;
    // add_pbc = false;
    // add_pbc_bending = true;
    // add_bending = false;
    // add_twisting = false;

    std::cout << "======================== CHECK HESSIAN ========================" << std::endl;
    
    int n_seg = 0;
    for (auto& rod : Rods)
        n_seg += rod->numSeg();
    int n_crossing = 0;
    for (auto& crossing : rod_crossings)
        n_crossing += 3;
    std::cout << "n crossing " << n_crossing << std::endl;
    std::cout << "n seg " << n_seg << std::endl;
    std::cout << "rod crosing dof begin " << rod_crossings[0]->reduced_dof_offset << std::endl;
    std::cout << "rod crosing dof end " << rod_crossings.back()->reduced_dof_offset + 3 << std::endl;
    std::cout << "rod theta dof begin " << Rods[0]->theta_reduced_dof_start_offset << std::endl;
    std::cout << "rod theta dof ends " << Rods.back()->theta_reduced_dof_start_offset + Rods.back()->numSeg() << std::endl;
    std::cout << W.rows() << " " << W.cols() << std::endl;

    T epsilon = 1e-6;
    StiffnessMatrix A;
    buildSystemDoFMatrix(dq, A);
    // std::cout << A.coeff(1, 0) << std::endl;
    n_dof = W.cols();
    for(int dof_i = 0; dof_i < n_dof; dof_i++)
    {
        dq(dof_i) += epsilon;
        VectorXT g0(n_dof), g1(n_dof);
        g0.setZero(); g1.setZero();
        computeResidual(g0, dq);

        dq(dof_i) -= 1.0 * epsilon;
        computeResidual(g1, dq);
            
        VectorXT row_FD = (g1 - g0) / (epsilon);

        for(int i = 0; i < n_dof; i++)
        {
            if(A.coeff(dof_i, i) == 0 && row_FD(i) == 0)
                continue;
            // if (std::abs( A.coeff(dof_i, i) - row_FD(i)) < 1e-4)
            //     continue;
            // std::cout << "node i: "  << std::floor(dof_i / T(dof)) << " dof " << dof_i%dof 
            //     << " node j: " << std::floor(i / T(dof)) << " dof " << i%dof 
            //     << " FD: " <<  row_FD(i) << " symbolic: " << A.coeff(i, dof_i) << std::endl;
            std::cout << "H(" << i << ", " << dof_i << ") " << " FD: " <<  row_FD(i) << " symbolic: " << A.coeff(i, dof_i) << std::endl;
            std::getchar();
        }
    }

    // run_diff_test = false;
    // add_regularizor = true;
    // add_stretching = true;
    // add_bending = true;
    // add_twisting = true;
}


template class EoLRodSim<double, 3>;
template class EoLRodSim<double, 2>;