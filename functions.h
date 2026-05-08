#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <utility>
#include "read_matrix.h"


//Ускорение свободного падения
const double g = 1.0;


double z_const(double x);

double u_1(double);


double u_2(double);


double h_1(double);


double h_2(double);


// class Linear_function {
// public:
//     Linear_function(double);
//     double k;
//     double operator()(double);
// };


//Тесты lakes
double z_lakes(double x);
//3.1.1. Lake at rest with an immersed bump

double u_0_311(double x);


double h_0_311(double x);


//3.1.2 Lake at rest with an emerged bump

double u_0_312(double x);

double h_0_312(double x);



//3.1.3 Subcritical flow
double u_0_313(double x);


double h_0_313(double x);


//3.1.5 Transcritical flow with shock
double u_0_315(double x);


double h_0_315(double x);

// 3.2.1 Long channel: 1000 m
double u_0_321(double x);

double h_0_321(double x);

double z_321(double x);

double z_321_derivative(double x);

double h_ex(double x);

double h_ex_derivative(double x);


//4.1.1 Dam break on a wet domain without friction
double u_0_411_left_right(double x);

double h_0_411_left(double x);

double h_0_411_right(double x);

//analytical solution
std::pair<double, double> analytical_solution(double x, double t);

double x_a(double t);

double x_b(double t);

double x_c(double t);


// тесты 2D

//well-balanced
double z2D_1(double x, double y);

double u02D_1(double x, double y);

double v02D_1(double x, double y);

double h02D_1(double x, double y);

double z2Dfear(double x, double y);

double h02D_gauss(double x, double y);

//Тест две волны разряжения
double v02D_2(double x, double y);

double u02D_2(double x, double y);

double h02D_2(double x, double y);


//Численные эксперименты из диплома
//Тест1
double h02D_dip1(double x, double y);


double u02D_dip1(double x, double y);


double v02D_dip1(double x, double y);


//Тест2
double h02D_dip2(double x, double y);


double u02D_dip2(double x, double y);


double v02D_dip2(double x, double y);


double well_balanced_2d_z(double x, double y);

#endif