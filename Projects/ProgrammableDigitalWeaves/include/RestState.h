#ifndef CURVATURE_FUNCTION_H
#define CURVATURE_FUNCTION_H

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "HybridC2Curve.h"
// #include "VecMatDef.h"

template<class T, int dim>
class RestState
{
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;

    Vector<T, dim + 1> starting_point, ending_point;

public:
    RestState(Vector<T, dim + 1> q0, Vector<T, dim + 1> q1)
        : starting_point(q0), ending_point(q1) {}
    RestState()
        : starting_point(Vector<T, dim + 1>::Zero()), ending_point(Vector<T, dim + 1>::Ones()) {}
    ~RestState() {}

    virtual T value(T u) {return 0;}
    virtual void gradient(T u, T& dedu) { dedu = 0; }
    virtual void hessian(T u, T& de2du2) { de2du2 = 0; }
    
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h) 
    { 
        X = starting_point.template segment<dim>(0) + 
            (u - starting_point[dim]) / (ending_point[dim] - starting_point[dim]) * (ending_point.template segment<dim>(0) - starting_point.template segment<dim>(0));
        dXdu = (ending_point.template segment<dim>(0) - starting_point.template segment<dim>(0)) / (ending_point[dim] - starting_point[dim]) ;
        d2Xdu2 = TV::Zero();
    }
};

template<class T, int dim>
class LineCurvature : public RestState<T, dim>
{
    using TV = Vector<T, dim>;

public:
    LineCurvature(Vector<T, dim + 1> q0, 
                Vector<T, dim + 1> q1) : RestState<T, dim>(q0, q1) {}
    LineCurvature() : RestState<T, dim>(){}
    
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h);
};

template<class T, int dim>
class DiscreteHybridCurvature : public RestState<T, dim>
{
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;
    HybridC2Curve<T, 2>* curve;
    std::vector<T> data_points_discrete_arc_length;

public:
    DiscreteHybridCurvature() : RestState<T, dim>(){}

    DiscreteHybridCurvature(Vector<T, dim + 1> q0, 
                            Vector<T, dim + 1> q1) : RestState<T, dim>(q0, q1) {}
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h); 

    void setData(HybridC2Curve<T, 2>* c, const std::vector<T>& v)
    {
        curve = c;
        data_points_discrete_arc_length = v;
    }

};

template<class T, int dim>
class PreBendCurvaure : public RestState<T, dim>
{
    T length;
    T theta;
public:
    PreBendCurvaure(T _length, T _theta) : length(_length), theta(_theta) {}
    virtual T value(T u)
    {
        return theta / length;
    }
};

template<class T, int dim>
class CircleCurvature : public RestState<T, dim>
{
private:
    T r;
public:
    CircleCurvature(T _r) : r(_r) {}

    virtual T value(T u);
    virtual void gradient(T u, T& dedu);
    virtual void hessian(T u, T& de2du2);
};

template<class T, int dim>
class SineCurvature : public RestState<T, dim>
{
private:
    T amp;
    T phi;
    T period;

public:
    SineCurvature(T _amp, T _phi, T _period) : amp(_amp), phi(_phi), period(_period) {}
    
    virtual T value(T u);
    //-g -> f
    virtual void gradient(T u, T& dedu);
    // -J -> Hessian
    virtual void hessian(T u, T& de2du2);
};

#endif