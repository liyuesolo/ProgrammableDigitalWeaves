#include "../include/RestState.h"

template<class T, int dim>
void DiscreteHybridCurvature<T, dim>::getMaterialPos(T u, TV& X, 
    TV& dXdu, TV& d2Xdu2, bool g, bool h) 
{
    bool debug = false;

    if (u < 0)
    {
        X = this->starting_point.template segment<dim>(0);
        dXdu = TV::Zero();
        d2Xdu2 = TV::Zero();
        return;
    }

    if (debug)
        std::cout << "getMaterialPos" << std::endl;
    T dx = 1.0;

    int left_node = 0;
    // find which curve the current points lies on
    for (int i = 0; i < curve->data_points.size(); i++)
    {
        T fraction_data_point = i;

        if( fraction_data_point > u || i == curve->data_points.size() - 1)
        {
            left_node = i - 1;
            break;
        }
    }

    int curve_id = 0;
    
    if(left_node > 1)
        curve_id = left_node - 1;

    T t = u - T(left_node * dx);
    
    bool interpolate = false;
    if ( left_node != 0 && left_node != curve->data_points.size() - 2)
        interpolate = true;
    if (curve->data_points.size() == 3 && u >= 1)
    {
        curve_id += 1;
        interpolate= false;
    }
        
    if (debug)
        std::cout << "t " << t << " u " << u <<  " curve idx " << curve_id << " left node " << left_node << " interpolate " << interpolate << std::endl;
    
    // curve->getPosOnCurve(curve_id, t, X, interpolate);

    if constexpr (dim == 3)
    {
        TV2 _X, _dXdu, _d2Xdu2;
        curve->getPosOnCurveWithDerivatives(curve_id, t, _X, _dXdu, _d2Xdu2, g, h, interpolate);
        X = TV(_X[0], _X[1], 0);
        dXdu = TV(_dXdu[0], _dXdu[1], 0);
        d2Xdu2 = TV(_d2Xdu2[0], _d2Xdu2[1], 0);
    }
    else if constexpr (dim == 2)
        curve->getPosOnCurveWithDerivatives(curve_id, t, X, dXdu, d2Xdu2, g, h, interpolate);

    // T shifted = shift ? t * 0.5 + 0.5 : t * 0.5;
    // X *= 0.03;
    if (debug)
    {
        std::cout << u << " " << X.transpose() << std::endl;
        std::cout << "getMaterialPos done" << std::endl;
        // std::getchar();
    }
}
template<class T, int dim>
void LineCurvature<T, dim>::getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h) 
{ 
    T scale = (this->ending_point[dim] - this->starting_point[dim]);
    X = this->starting_point.template segment<dim>(0) + 
        (u - this->starting_point[dim]) / scale
         * (this->ending_point.template segment<dim>(0) - this->starting_point.template segment<dim>(0));
    
    dXdu = (this->ending_point.template segment<dim>(0) - this->starting_point.template segment<dim>(0)) / scale;
    d2Xdu2 = TV::Zero();
}

template<class T, int dim>
T CircleCurvature<T, dim>::value(T u) 
{ 
    // left most / right most boundary
    if (u < 1e-6 || std::abs(u - 2* M_PI * r) < 1e-6)
        return 0;
    T middle_point_arclength = M_PI * r;
    T dis = (u - middle_point_arclength);
    if (dis < -1e-6)
        return 1.0 / r; 
    else if(std::abs(dis) < 1e-6)
        return 0;
    else if(dis > 1e-6)
        return -1.0 / r; 
    else
        std::cout << "undefined curvature for u=" << u << std::endl;
}

template<class T, int dim>
void CircleCurvature<T, dim>::gradient(T u, T& dedu) 
{
     dedu = 0; 
}
template<class T, int dim>
void CircleCurvature<T, dim>::hessian(T u, T& de2du2) 
{ 
    de2du2 = 0; 
}

template<class T, int dim>
T SineCurvature<T, dim>::value(T u) 
{ 
    T t1 = period*period;
    T t3 = period*u;
    T t4 = std::sin(t3);
    T t5 = amp*amp;
    T t7 = std::cos(t3);
    T t8 = t7*t7;
    T t11 = std::pow(t5*t1*t8+1.0,-0.15E1);
    return -amp*t1*t4*t11;
}

//-g -> f
template<class T, int dim>
void SineCurvature<T, dim>::gradient(T u, T& dedu)
{
    T t1 = period*period;
    T t4 = period*u;
    T t5 = std::cos(t4);
    T t6 = amp*amp;
    T t8 = t5*t5;
    T t10 = t6*t1*t8+1.0;
    T t11 = std::pow(t10,-0.15E1);
    T t15 = t1*t1;
    T t18 = std::sin(t4);
    T t19 = t18*t18;
    T t20 = std::pow(t10,-0.25E1);
    dedu = amp*t1*period*t5*t11+0.3E1*t6*amp*t15*period*t19*t20*t5;
}
// -J -> Hessian
template<class T, int dim>
void SineCurvature<T, dim>::hessian(T u, T& de2du2)
{
    T t1 = period*period;
    T t2 = t1*t1;
    T t4 = period*u;
    T t5 = std::sin(t4);
    T t6 = amp*amp;
    T t8 = std::cos(t4);
    T t9 = t8*t8;
    T t11 = t6*t1*t9+1.0;
    T t12 = std::pow(t11,-0.15E1);
    T t17 = t6*amp*t2*t1;
    T t18 = std::pow(t11,-0.25E1);
    T t23 = t6*t6;
    T t25 = t2*t2;
    T t27 = t5*t5;
    T t28 = t27*t5;
    T t29 = std::pow(t11,-0.35E1);
    de2du2 = amp*t2*t5*t12+0.9E1*t17*t9*t18*t5+0.15E2*t23*amp*t25*t28*t29*t9-0.3E1*t17*t28*t18;
}

template class SineCurvature<double, 3>;
template class SineCurvature<double, 2>;   

template class CircleCurvature<double, 3>;
template class CircleCurvature<double, 2>;   

template class LineCurvature<double, 3>;
template class LineCurvature<double, 2>; 

template class DiscreteHybridCurvature<double, 3>;
template class DiscreteHybridCurvature<double, 2>;   
