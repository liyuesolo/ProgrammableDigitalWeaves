#ifndef HOMOGENIZATION_H
#define HOMOGENIZATION_H

#include "EoLRodSim.h"

template<class T, int dim>
class EoLRodSim;

template<class T, int dim>
class Homogenization
{
    using TVDOF = Vector<T, dim+2>;
    using DOFStack = Matrix<T, dim + 2, Eigen::Dynamic>;
    using IV2 = Vector<int, 2>;
    using TV2 = Vector<T, 2>;
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using TM2 = Matrix<T, 2, 2>;
    using CDoF2D = Vector<T, 6>;
    using CHessian2D = Matrix<T, 6, 6>;
    using ComplianceTensor = Matrix<T, 3 * (dim - 1), 3 * (dim - 1)>;
    using TVEntry = Vector<T, 3 * (dim - 1)>;

    using ComplianceTensorFull = Matrix<T, dim * dim, dim * dim>;
    using TVEntryFull = Vector<T, dim * dim>;

    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using Offset = Vector<int, dim + 1>;

public:
    EoLRodSim<T, dim>& sim;
    
    const VectorXT& rest_states = sim.rest_states;
    const VectorXT& deformed_states = sim.deformed_states;

    T s1 = 1.01;
    T s2 = 1.0;

    bool biaxial = false;
    
public:
    Homogenization(EoLRodSim<T, dim>& eol_sim) : sim(eol_sim) {}
    ~Homogenization() {}

    void initialize();
    void testOneSample();

    void marcoMaterialParametersFitting();
    void materialParametersFromUniaxialStrain(T theta, T s, TV2& E_nu);
    void computeYoungsModulusPoissonRatioBatch();

    void computeMacroStressStrain(TM2& stress_marco, TM2& strain_marco);
    
};

#endif