#include <cstdlib>
#include <ctime>
#include "functions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"


double z_const(double) {
    return 0.0;
}

double u_1(double x) {
    return 1.0;
    // return -1.0; // волны разряжения //discharge_waves
    // return -0.1; //волны рязряжения //hydrajump_dischargewave
    // return 2.0; // гидродинамические прыжки
    // return 1.0; // гидродинамический прыжок и волна разряжения
    // return 0.0;
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
    // return std::exp(-std::pow(x - 10, 2)) + 3;
    // return 0.5 + 1.0 / 20.0 * x; 
    return 2.0 - z_lakes(x);
}


double h_2(double x) {
    // return 0.5; //волны разряжения
    // return 1.0; // волны разряжения //hydrajump_dischargewave //discharge_waves
    // return 0.5; // гидродинамические прыжки
    // return 0.5; // гидродинамический прыжок и волна разряжения
    // return std::exp(-std::pow(x - 10, 2)) + 3;
    // return 1.0 + 1.0 / 20.0 * x; 
    return 2.0 - z_lakes(x);
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
    // return h_ex(1000.0);
    return h_ex(x);
}


double z_321(double x) {
    // return -h_ex(x);
    double h = 1.0 / 100.0;
    if (x == 0.0) {
        return 0.0;
    } else {

        int n = static_cast<int>(x / h);
        if (n % 2 == 0) {
            ++n;
        }
        h = x / (n - 1);

        double sum = z_321_derivative(0.0) + z_321_derivative(x);
        for(int i = 1; i < (n - 1); ++i) {

            if (i % 2 == 1) {
                sum += 4 * z_321_derivative(i * h);
            } else {
                sum += 2 * z_321_derivative(i * h);
            }

        }

        return sum * h / 3.0;

    }

}

double z_321_derivative(double x) {
    return (2.25 / (g * std::pow(h_ex(x), 3.0)) - 1.0) * h_ex_derivative(x);
    // tests:
    // return 3 * sqr(x);
}


double h_ex(double x) {
    return std::pow(4.0 / g, 1.0 / 3.0) * (1.0 + 1.0 / 2.0 * std::exp(-16 * sqr(x / 1000.0 - 1.0 / 2.0)));  
}


double h_ex_derivative(double x) {
    return -std::pow(4.0 / g, 1.0 / 3.0) * 2.0 / 125.0 * (x / 1000.0 - 1.0 / 2.0) * std::exp(-16.0 * sqr(x / 1000.0 - 1.0 / 2.0));
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


// тесты 2D

//well-balanced
double z2D_1(double x, double y) {
    return 0.0;
}


double u02D_1(double x, double y) {
    return 0.0;
}


double v02D_1(double x, double y) {
    return 0.0;
}


double h02D_1(double x, double y) {
    return 1.0 - z2Dfear(x, y);
}

double z2Dfear(double x, double y) {
    double s =
        std::sin(35.0*x) +
        0.9*std::sin(33.0*y) +
        0.7*std::sin(29.0*(x+y)) +
        0.6*std::sin(41.0*(x-0.3*y)) +
        0.5*std::sin(37.0*(0.7*x+1.1*y));

    // нормировка "примерно" в [-1,1]
    s /= (1.0 + 0.9 + 0.7 + 0.6 + 0.5);
    s = std::clamp(s, -1.0, 1.0);
    return 0.3 * s; // ~[-0.3, 0.3]
}


double h02D_gauss(double x, double y) {
    return std::exp(-4.5 * (x * x + y * y)) + 1.0;
}


//Тест две волны разряжения
double v02D_2(double x, double y) {
    return 0.0;
}

double u02D_2(double x, double y) {
    return (x > 0.0 ? 0.0 : -3.0);
}


double h02D_2(double x, double y) {
    return (x > 0.0 ? 0.5 : 1.0);
}


//Численные эксперименты из диплома
//Тест1
double h02D_dip1(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r = std::sqrt(x * x + y * y);
    double r0 = 0.03;
    return -(alpha * alpha / (4.0 * beta)) * std::exp(2.0 * beta * (1.0 - (r/r0) * (r/r0))) + 1.0;
}


double u02D_dip1(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r = std::sqrt(x * x + y * y);
    double r0 = 0.03;
    return alpha / r0 * std::exp(beta * (1.0 - (r/r0) * (r/r0))) * y;
}


double v02D_dip1(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r = std::sqrt(x * x + y * y);
    double r0 = 0.03;
    return -alpha / r0 * std::exp(beta * (1.0 - (r/r0) * (r/r0))) * x;
}


//Тест2
double h02D_dip2(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r1 = std::sqrt((x - 0.1) * (x - 0.1) + y * y);
    double r2 = std::sqrt((x + 0.1) * (x + 0.1) + y * y);
    double r0 = 0.03;
    return -(alpha * alpha / (4.0 * beta)) * std::exp(2.0 * beta * (1.0 - (r1/r0) * (r1/r0))) 
           - (alpha * alpha / (4.0 * beta)) * std::exp(2.0 * beta * (1.0 - (r2/r0) * (r2/r0))) + 1.0;
}


double u02D_dip2(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r1 = std::sqrt((x - 0.1) * (x - 0.1) + y * y);
    double r2 = std::sqrt((x + 0.1) * (x + 0.1) + y * y);
    double r0 = 0.03;
    return alpha / r0 * std::exp(beta * (1.0 - (r1/r0) * (r1/r0))) * y - alpha / r0 * std::exp(beta * (1.0 - (r2/r0) * (r2/r0))) * y;
}


double v02D_dip2(double x, double y) {
    double alpha = 0.404;
    double beta = 0.3;
    double r1 = std::sqrt((x - 0.1) * (x - 0.1) + y * y);
    double r2 = std::sqrt((x + 0.1) * (x + 0.1) + y * y);
    double r0 = 0.03;
    return -alpha / r0 * std::exp(beta * (1.0 - (r1/r0) * (r1/r0))) * (x - 0.1) + alpha / r0 * std::exp(beta * (1.0 - (r2/r0) * (r2/r0))) * (x + 0.1);
}


double well_balanced_2d_z(double x, double y) {
    return 0.3 * (std::rand() % 1000) / 1000.0;
}


#pragma GCC diagnostic pop