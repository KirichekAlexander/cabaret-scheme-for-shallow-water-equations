#include "functions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"


double z_const(double) {
    return 0.0;
}

double u_1(double x) {
    // return 1.0;
    // return -1.0; // волны разряжения //discharge_waves
    // return -0.1; //волны рязряжения //hydrajump_dischargewave
    // return 2.0; // гидродинамические прыжки
    // return 1.0; // гидродинамический прыжок и волна разряжения
    return 0.0;
}


double u_2(double x) {
    return 0.0; // волны разряжения
    // return 0.1; // волны разряжения //hydrajump_dischargewave
    // return -1.0; // гидродинамический прыжки 
    // return -0.5; // гидродинамический прыжок и волна разряжения
    // return 1.0; //discharge_waves
}


double h_1(double x) {
    // return 1.0; // волны разряжения
    // return 0.5; // волны разряжения //hydrajump_dischargewave
    // return 0.6; // гидродинамические прыжки
    // return 1.0; // гидродинамический прыжок и волна разряжения //discharge_waves
    return std::exp(-std::pow(x - 10, 2)) + 3;
    // return 0.5 + 1.0 / 20.0 * x; 
}


double h_2(double x) {
    // return 0.5; //волны разряжения
    // return 1.0; // волны разряжения //hydrajump_dischargewave //discharge_waves
    return 0.5; // гидродинамические прыжки
    // return 0.5; // гидродинамический прыжок и волна разряжения
    // return std::exp(-std::pow(x - 10, 2)) + 3;
    // return 1.0 + 1.0 / 20.0 * x; 
}


// //Реализация линейной функции

// Linear_function::Linear_function(double k) 
//     : k(k)
// {
// }

// double Linear_function::operator()(double x) {
//     return x * k;
// }



//Тесты lakes

double z_lakes(double x) {
    return ((x > 8) and (x < 12) ? 0.2 - 0.05 * (x - 10.0) * (x - 10.0) : 0.0);
}

//3.1.1. Lake at rest with an immersed bump

double u_0_311(double x) {
    return 0.0;
}

double h_0_311(double x) {
    return 0.5 - z_lakes(x);
}

//3.1.2 Lake at rest with an emerged bump

double u_0_312(double x) {
    return 0.0;
}

double h_0_312(double x) {
    return (0.1 > z_lakes(x) ? 0.1 : z_lakes(x)) - z_lakes(x);
}

//3.1.3 Subcritical flow

double u_0_313(double x) {
    return 0.0;
}


double h_0_313(double x) {
    return 2.0 - z_lakes(x);
}

//3.1.5 Transcritical flow with shock
double u_0_315(double x) {
    return 0.0;
}


double h_0_315(double x) {
    return 0.33 - z_lakes(x);
}
//

// 3.2.1 Long channel: 1000 m
double u_0_321(double x) {
    return 0.0;
}


double h_0_321(double x) {
    return h_ex(1000.0);
}


double z_321(double x) {
    return -h_ex(x);
}


double h_ex(double x) {
    return std::pow(4.0 / g, 1.0 / 3.0) * (1.0 + 1.0 / 2.0 * std::exp(-16 * (x / 1000.0 - 1.0 / 2.0) * (x / 1000.0 - 1.0 / 2.0)));  
}
//


//4.1.1 Dam break on a wet domain without friction
double u_0_411_left_right(double x) {
    return 0.0;
}


double h_0_411_left(double x) {
    return 0.005;
}


double h_0_411_right(double x) {
    return 0.001;
}

//analytical solution
std::pair<double, double> analytical_solution(double x, double t) {

    if (x < x_a(t)) {
        return {0.005, 0.0};
    } else if (x_a(t) <= x and x < x_b(t)) {
        return {4 / (9 * g) * (std::sqrt(g * 0.005) - (x - 5.0) / (2.0 * t)) *  (std::sqrt(g * 0.005) - (x - 5.0) / (2.0 * t)) , 2.0 / 3.0 * ((x - 5.0) / t + std::sqrt(g * 0.005))};
    } else if (x_b(t) <= x and x < x_c(t)) {
        return {0.157832 * 0.157832 / g, 2 * (std::sqrt(g * 0.005) - 0.157832)};
    } else {
        return {0.001, 0.0};
    }

}


double x_a(double t) {
    return 5.0 - t * std::sqrt(g * 0.005);
}


double x_b(double t) {
    return 5.0 + t * (2 * std::sqrt(g * 0.005) - 3 * 0.157832);
}


double x_c(double t) {
    return 5.0 + t * (2 * 0.157832 * 0.157832 * (std::sqrt(g * 0.005) - 0.157832) / (0.157832 * 0.157832 - g * 0.001));
}

//



#pragma GCC diagnostic pop