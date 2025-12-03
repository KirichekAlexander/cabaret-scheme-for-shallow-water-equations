#ifndef AUX_MODULE_H
#define AUX_MODULE_H



#include <functional>
#include "read_matrix.h"


//Синоним для функции двух переменных
using Func2D = double (*)(double, double);

enum class Side {LEFT, RIGHT, TOP, BOTTOM};

/*
Структура ячейки
 ----T----
|         |
|         |
L    c    R
|         |
|         |
 ----B----
*/

//Точка из 3 величин
struct Point3 {
    double u = 0.0;
    double v = 0.0;
    double h = 0.0;
};


struct Vec2
{
    double a;
    double b;

    Vec2 operator+(Vec2 const& other) const;
    Vec2 operator*(double k) const;
    Vec2 operator-(Vec2 const& other) const;
    Vec2 operator/(double k) const;
    Vec2 operator*(Vec2 const& other) const;
    Vec2 operator<(Vec2 const& other) const;
};


Vec2 operator*(double k, Vec2 const& vec);


template<typename T>
struct Vec3{
    T a;
    T b;
    T c;

    Vec3<T> operator+(Vec3<T> const& other) const;
    Vec3<T> operator*(double k) const;
    Vec3<T> operator-(Vec3<T> const& other) const;
    Vec3<T> operator/(double k) const;
    Vec3<T> operator*(Vec3<T> const& other) const;
    Vec3<T> operator<(Vec3<T> const& other) const;
};


template<typename T>
Vec3<T> Vec3<T>::operator+(Vec3<T> const& other) const {
    return Vec3{a + other.a,
                b + other.b, 
                c + other.c};
}


template<typename T>
Vec3<T> Vec3<T>::operator*(double k) const {
    return Vec3<T>{a * k,
                b * k, 
                c * k};
}


template<typename T>
Vec3<T> Vec3<T>::operator-(Vec3<T> const& other) const {
    return Vec3<T>{a - other.a,
                b - other.b,
                c - other.c};
}



template<typename T>
Vec3<T> Vec3<T>::operator/(double k) const {
      return Vec3<T>{a / k,
                  b / k, 
                  c / k};
}


template<typename T>
Vec3<T> Vec3<T>::operator*(Vec3<T> const& other) const {
    return Vec3<T>{a * other.a,
                b * other.b,
                c * other.c};
}



template<typename T>
Vec3<T> operator*(double k, Vec3<T> const& vec) {
    return vec * k;
}


Vec3<Vec2> operator*(Vec3<double> const& vec1, Vec3<Vec2> const& vec2);


Vec3<Vec2> min_invariants(std::initializer_list<Vec3<Vec2>> list);


Vec3<Vec2> max_invariants(std::initializer_list<Vec3<Vec2>> list);


// Синонимы для массивов точек (u, v, h)
using Field1D = std::vector<Point3>;
using Field2D = std::vector<std::vector<Point3>>;


//Структура для просмотра ячейки
struct CellView {
    Point3& L;
    Point3& R;
    Point3& T;
    Point3& B;
    Point3& c;
    double& zL;
    double& zR;
    double& zT;
    double& zB;
    double& zc;
};



#endif