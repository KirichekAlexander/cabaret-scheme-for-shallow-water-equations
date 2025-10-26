#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <utility>


//Ускорение свободного падения
const double g = 9.81;


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

double h_ex(double x);


//4.1.1 Dam break on a wet domain without friction
double u_0_411_left_right(double x);

double h_0_411_left(double x);

double h_0_411_right(double x);

//analytical solution
std::pair<double, double> analytical_solution(double x, double t);

double x_a(double t);

double x_b(double t);

double x_c(double t);


#endif