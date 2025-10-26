#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


#include <utility>


enum class BoundaryType {

    DEFAULT,
    NON_LEAK,
    FREE_EXIT,
    FLUVIAL_H_GIVEN,
    FLUVIAL_FLUX_GIVEN

};


template<BoundaryType BT>
class Boundary;


//Класс границы обычные (задана высота и поток)
template<>
class Boundary<BoundaryType::DEFAULT> {
public:
    Boundary(double h, double q) 
        : h(h)
        , q(q)
    {
    }

    double h;
    double q;
};


//Класс границы непротекания
template<>
class Boundary<BoundaryType::NON_LEAK> {
};



//Класс границы свободного выхода
template<>
class Boundary<BoundaryType::FREE_EXIT> {
public:
    double const_g;
    double const_invariant;
};


//Речное течение дана высота
template<>
class Boundary<BoundaryType::FLUVIAL_H_GIVEN> {
public:
    Boundary(double h) 
        : h(h)
    {
    }

    double h;
};


//Речное течение дан поток
template<>
class Boundary<BoundaryType::FLUVIAL_FLUX_GIVEN> {
public:
    Boundary(double q) 
        : q(q)
    {
    }

    double q;
};


#endif