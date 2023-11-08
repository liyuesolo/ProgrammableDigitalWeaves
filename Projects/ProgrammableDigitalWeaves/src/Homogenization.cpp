#include "../include/Homogenization.h"
#include <fstream>
#include <iomanip>

template<class T, int dim>
void Homogenization<T, dim>::testOneSample()
{
    initialize();
    // sim.disable_sliding = false;
    TV strain_dir, ortho_dir;
    
    if (sim.add_pbc)
        sim.setUniaxialStrain(45/180.0 * M_PI, s1, strain_dir, ortho_dir);
    // sim.setUniaxialStrain(0.1885, s1, strain_dir, ortho_dir);
    
    // sim.advanceOneStep();
    // TM2 stress_macro, strain_macro;
    // computeMacroStressStrain(stress_macro, strain_macro);
    // std::cout << stress_macro.norm() << std::endl;
    // TV2 E_nu;
    // materialParametersFromUniaxialStrain(1.885, s1, E_nu);
    // std::cout << "theta: " << 3.3929 << " youngs_modulus " << E_nu(0) << " Poisson Ratio: " << E_nu(1) << std::endl;
    // sim.setUniaxialStrain(0.0/180.0 * M_PI, s1, strain_dir, ortho_dir);
}

template<class T, int dim>
void Homogenization<T, dim>::initialize()
{
    
    sim.print_force_mag = false;
    sim.disable_sliding = true;
    sim.verbose = false;
    // sim.buildPlanePeriodicBCScene3x3Subnodes(8);
    sim.buildSceneFromUnitPatch(18);
    
    // sim.buildPlanePeriodicBCScene3x3();
    
    
    sim.add_penalty = false;

    sim.add_shearing = false;

    
    sim.newton_tol = 1e-6;
        
    sim.k_strain = 1e8;
    
    
    sim.k_yc = 1e8;
    sim.k_pbc = 1e8;
    sim.kr = 1e3;
    
    s1 = 1.1;
    s2 = 1.0;
}

template<class T, int dim>
void Homogenization<T, dim>::computeMacroStressStrain(TM2& stress_marco, TM2& strain_marco)
{

    TV xi = deformed_states.template segment<dim>(sim.pbc_pairs_reference[0].first.first[0]);
    TV xj = deformed_states.template segment<dim>(sim.pbc_pairs_reference[0].first.second[0]);
    TV xk = deformed_states.template segment<dim>(sim.pbc_pairs_reference[1].first.first[0]);
    TV xl = deformed_states.template segment<dim>(sim.pbc_pairs_reference[1].first.second[0]);

    TV Xi = rest_states.template segment<dim>(sim.pbc_pairs_reference[0].first.first[0]);
    TV Xj = rest_states.template segment<dim>(sim.pbc_pairs_reference[0].first.second[0]);
    TV Xk = rest_states.template segment<dim>(sim.pbc_pairs_reference[1].first.first[0]);
    TV Xl = rest_states.template segment<dim>(sim.pbc_pairs_reference[1].first.second[0]);

    TM2 X = TM2::Zero(), x = TM2::Zero();
    X.col(0) = (Xi - Xj).template segment<2>(0);
    X.col(1) = (Xk - Xl).template segment<2>(0);

    x.col(0) = (xi - xj).template segment<2>(0);
    x.col(1) = (xk - xl).template segment<2>(0);


    TM2 F_macro = x * X.inverse();
    
    strain_marco = 0.5 * (F_macro.transpose() + F_macro) - TM2::Identity();

    TM2 R90 = TM2::Zero();
    

    R90.row(0) = TV2(0, -1);
    R90.row(1) = TV2(1, 0);

    TV2 n0 = (R90 * (xj - xi).template segment<2>(0)).normalized(), 
        n1 = (R90 * (xl - xk).template segment<2>(0)).normalized();
    
    VectorXT f(deformed_states.rows());
    f.setZero(); 
    sim.addPBCForce(f);
    f *= -1;
    std::vector<TV> f_bc(2, TV::Zero());

    
    sim.iteratePBCPairsWithDirection([&](int direction, Offset offset_i, Offset offset_j)
    {
        
        TV Xj_bc = rest_states.template segment<dim>(offset_j[0]);
        TV Xi_bc = rest_states.template segment<dim>(offset_i[0]);

        // T length = (Xj_bc - Xi_bc).norm();
        T length = (Xj_bc - Xi_bc).norm();
        Offset bc_node = direction == 0 ? offset_j : offset_i;
        if (std::abs(length) >= 1e-6)
            f_bc[direction].template segment<dim>(0) += f.template segment<dim>(bc_node[0]);
        
    });

    TM2 F_bc = TM2::Zero(), n_bc = TM2::Zero();
    F_bc.col(0) = f_bc[0].template segment<2>(0); 
    F_bc.col(1) = f_bc[1].template segment<2>(0);

    // std::cout << F_bc << std::endl;

    n_bc.col(0) = n1; n_bc.col(1) = n0;

    stress_marco = F_bc * n_bc.inverse();
    
    
}

template<class T, int dim>
void Homogenization<T, dim>::computeYoungsModulusPoissonRatioBatch()
{
    initialize();
    
    int n_angles = 400;
    T cycle = 2. * M_PI;
    std::vector<T> thetas;
    std::vector<T> youngs_moduli;
    std::vector<T> poisson_ratio;
    for (T theta = 0; theta <= cycle; theta += cycle/(T)n_angles)
    {
        T theta6 = std::round( theta * 1e4 ) / 1e4;
        
        thetas.push_back(theta6);
        TV2 E_nu;
        materialParametersFromUniaxialStrain(theta6, s1, E_nu);
        // std::cout << "theta: " << theta / M_PI * 180.0 << " youngs_modulus " << E_nu(0) << " Poisson Ratio: " << E_nu(1) << std::endl;
        std::cout << "theta: " << theta6 << " youngs_modulus " << E_nu(0) << " Poisson Ratio: " << E_nu(1) << std::endl;
        T E = E_nu(0); T nu = E_nu(1);
        youngs_moduli.push_back(E);
        // std::cout << E << " " << youngs_moduli[0] << std::endl;
        poisson_ratio.push_back(nu);
    }
    for(T theta : thetas)
        std::cout << theta << " ";
    std::cout << std::endl;
    for(T E : youngs_moduli)
        std::cout << E << " ";
    std::cout << std::endl;
    for(T nu : poisson_ratio)
        std::cout << nu << " ";
    std::cout << std::endl;
    
}

template<class T, int dim>
void Homogenization<T, dim>::materialParametersFromUniaxialStrain(T theta, T s, TV2& E_nu)
{
    TV strain_dir, ortho_dir;
    sim.setUniaxialStrain(theta, s, strain_dir, ortho_dir);
    
    // sim.setBiaxialStrain(theta, s, theta, 1.0, , strain_dir, ortho_dir);
    sim.advanceOneStep();
    
    E_nu = TV2::Zero();
    TM2 stress, strain;
    computeMacroStressStrain(stress, strain);
    T stretch_in_d = strain_dir.template segment<2>(0).dot(strain * strain_dir.template segment<2>(0));
    
    E_nu(0) = strain_dir.template segment<2>(0).dot(stress * strain_dir.template segment<2>(0)) / stretch_in_d;
    E_nu(1) = -ortho_dir.template segment<2>(0).dot(strain * ortho_dir.template segment<2>(0)) / stretch_in_d;
    sim.resetScene();
    
}


template class Homogenization<double, 3>;
template class Homogenization<double, 2>;   