#ifndef CABARET_SCHEME_2D_H
#define CABARET_SCHEME_2D_H


#include <functional>
#include "read_matrix.h"
#include "functions.h"
#include "boundary_conditions.h"
#include "file_manager.h"


//Синоним для функции двух переменных
using Func2D = double (*)(double, double);


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


struct Vec3 {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    Vec3 operator+(Vec3 const& other) const;
    Vec3 operator*(double k) const;
    Vec3 operator-(Vec3 const& other) const;
    Vec3 operator/(double k) const;
};

//Реализация оператор для Vec3
Vec3 Vec3::operator+(Vec3 const& other) const {
    return Vec3{a + other.a,
                b + other.b, 
                c + other.c};
}


Vec3 Vec3::operator*(double k) const {
    return Vec3{a * k,
                b * k, 
                c * k};
}

Vec3 Vec3::operator-(Vec3 const& other) const {
    return Vec3{a - other.a,
                b - other.b,
                c - other.c};
}

Vec3 Vec3::operator/(double k) const {
      return Vec3{a / k,
                  b / k, 
                  c / k};
}


Vec3 operator*(double k, Vec3 const& vec) {
    return vec * k;
}


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






/*
*Схема КАБАРЕ для уравнений мелкой воды в двумерном случае
*/
class CabaretScheme2D {

public:

    //Конструктор схемы
    explicit CabaretScheme2D(double CFL,
                             int nx, int ny,
                             double l1x, double l2x,
                             double l1y, double l2y,
                             double T,
                             Func2D u0, Func2D v0,
                             Func2D h0, Func2D z0);

    void compute(); // метод вычисления


private:
    // метод строит пространтсвенную сетку
    void build_grid(Row& grid_center, Row& grid_face,
                    int n, double delta,
                    double l1);
    
    void init_data(); // метод инициализации данных

    CellView get_cell(int i, int j); // просмотр ячейки по индексам

    double compute_time_step(); // вычисление временного шага

    Vec3 make_U(Point3 const& p);
    Vec3 make_G(Point3 const& p, double z);
    Vec3 make_H(Point3 const& p, double z);
    Point3 make_Point3_from_U(Vec3 const& U);
    void first_phase(); // первая фаза
    void second_phase(); // вторая фаза
    void third_phase(); // третья фаза

    double CFL;

    int nx; //количество ячеек по x
    int ny; //количество ячеек по y

    //Размеры прямоугольника
    double l1x;
    double l2x;
    double l1y;
    double l2y;
    //

    //шаги сеток по пространству(рассматриваются равномерные сетки)
    double dx;
    double dy;
    //

    double dt; //неравномерный шаг по времени
    double t; // текущее время
    double T; //Время процесса

    //начальные данные
    Func2D u0;
    Func2D v0;
    Func2D h0;
    Func2D z0;
    //

    Row x_center; // x в центрах
    Row x_face_x; // x на вертикальных гранях
    Row y_center; // y в центрах
    Row y_face_y; // y на горизонтальных гранях
 
    Matrix z_center; // z в центрах 
    Matrix z_face_x; // z на вертикальных гранях
    Matrix z_face_y; // z на горизонтальных гранях

    Field2D center; // (u, v, h) в центрах
    Field2D half_step_center; // (u, v, h) в промежуточном временном слою в центрах ячеек
    Field2D face_x; // (u, v, h) на вертикальных гранях
    Field2D face_y; // (u, v, h) на горизонтальных гранях
};


//Реализация CabaretScheme2D
CabaretScheme2D::CabaretScheme2D(double CFL,
                                 int nx, int ny,
                                 double l1x, double l2x,
                                 double l1y, double l2y,
                                 double T,
                                 Func2D u0, Func2D v0,
                                 Func2D h0, Func2D z0) 
    : CFL(CFL)
    , nx(nx)
    , ny(ny)
    , l1x(l1x)
    , l2x(l2x)
    , l1y(l1y)
    , l2y(l2y)
    , dx((l2x - l1x) / nx)
    , dy((l2y- l1y) / ny)
    , t(0.0)
    , T(T)
    , u0(u0)
    , v0(v0)
    , h0(h0)
    , z0(z0)
    , x_center(nx, 0.0)
    , x_face_x(nx + 1, 0.0)
    , y_center(ny, 0.0)
    , y_face_y(ny + 1, 0.0)
    , z_center(nx, ny)
    , z_face_x(nx + 1, ny)
    , z_face_y(nx, ny + 1)
    , center(nx, Field1D(ny))
    , half_step_center(nx, Field1D(ny))
    , face_x(nx + 1, Field1D(ny))
    , face_y(nx, Field1D(ny + 1))
{
    build_grid(x_center, x_face_x, nx, dx, l1x); // вычисление точек по x
    build_grid(y_center, y_face_y, ny, dy, l1y); // вычисление точек по y
    init_data(); // инциализация по начальным данным
}


void CabaretScheme2D::build_grid(Row& grid_center, Row& grid_face,
                                 int n, double delta,
                                 double l1) {

    for(int i = 0; i < n; ++i) {
        grid_center[i] = l1 + (i + 0.5) * delta;
        grid_face[i] = l1 + i * delta;
    }
    grid_face[n] = l1 + n * delta;

}


void CabaretScheme2D::init_data() {

    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {

            z_center[i][j] = z0(x_center[i], y_center[j]);
            z_face_x[i][j] = z0(x_face_x[i], y_center[j]);
            z_face_y[i][j] = z0(x_center[i], y_face_y[j]);

            center[i][j].u = u0(x_center[i], y_center[j]);
            face_x[i][j].u = u0(x_face_x[i], y_center[j]);
            face_y[i][j].u = u0(x_center[i], y_face_y[j]);

            center[i][j].v = v0(x_center[i], y_center[j]);
            face_x[i][j].v = v0(x_face_x[i], y_center[j]);
            face_y[i][j].v = v0(x_center[i], y_face_y[j]);

            center[i][j].h = h0(x_center[i], y_center[j]);
            face_x[i][j].h = h0(x_face_x[i], y_center[j]);
            face_y[i][j].h = h0(x_center[i], y_face_y[j]);

        }
    }

    //заполняю крайнюю верхнюю гортзонтальную грань
    for(int i = 0; i < nx; ++i) {

            z_face_y[i][ny] = z0(x_center[i], y_face_y[ny]);

            face_y[i][ny].u = u0(x_center[i], y_face_y[ny]);

            face_y[i][ny].v = v0(x_center[i], y_face_y[ny]);

            face_y[i][ny].h = h0(x_center[i], y_face_y[ny]);

    }

    //заполняю крайнюю правую вертикальную грань
    for(int i = 0; i < ny; ++i) {

            z_face_x[nx][i] = z0(x_face_x[nx], y_center[i]);

            face_x[nx][i].u = u0(x_face_x[nx], y_center[i]);

            face_x[nx][i].v = v0(x_face_x[nx], y_center[i]);

            face_x[nx][i].h = h0(x_face_x[nx], y_center[i]);
            
    }

}


CellView CabaretScheme2D::get_cell(int i, int j) {
    return CellView{face_x[i][j],
                    face_x[i + 1][j],
                    face_y[i][j + 1],
                    face_y[i][j],
                    center[i][j],
                    z_face_x[i][j],
                    z_face_x[i + 1][j],
                    z_face_y[i][j + 1],
                    z_face_y[i][j],
                    z_center[i][j]};
}


void CabaretScheme2D::compute() {

    while(t != T) {

        dt = compute_time_step();

        first_phase();
        second_phase();
        third_phase();


    }

}


double CabaretScheme2D::compute_time_step() {

    //Изначально шаг - максимальное число double
    double res = std::numeric_limits<double>::max();

    //Проход по всем центрам ячеек
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {

            Point3 const& c = center[i][j];
            res = std::min({res, CFL * dx / (std::abs(c.u) + std::sqrt(g * c.h)),
                            CFL * dy / (std::abs(c.v) + std::sqrt(g * c.h))});

        }
    }

    return res;

}


Vec3 CabaretScheme2D::make_U(Point3 const& p) {
    return Vec3{p.h,        // h
                p.h * p.u,  // hu
                p.h * p.v}; // hv
}


Vec3 CabaretScheme2D::make_G(Point3 const& p, double z) {
    return Vec3{p.h * p.u,                                                     // hu
                p.h * sqr(p.u) + g * (sqr(p.h + z) - 2 * z * (p.h + z)) * 0.5, // hu^2 + g(H^2 - 2zH) / 2
                p.h * p.v * p.u};                                              // hvu
}


Vec3 CabaretScheme2D::make_H(Point3 const& p, double z) {
    return Vec3{p.h * p.v,                                                      // hv
                p.h * p.v * p.u,                                                // hvu
                p.h * sqr(p.v) + g * (sqr(p.h + z) - 2 * z * (p.h + z)) * 0.5}; // hv^2 + g(H^2 - 2zH) / 2
}


Point3 CabaretScheme2D::make_Point3_from_U(Vec3 const& U) {
    return Point3{U.b / U.a, // u
                  U.c / U.a, // v
                  U.a};      // h
}


void CabaretScheme2D::first_phase() {

    //Проход по всем ячейкам
    for(int i = 0 ; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            //Решаем первую фазу в векторном виде
            CellView cell = get_cell(i, j);
            Point3& half_step_c = half_step_center[i][j];
            Vec3 U_c = make_U(cell.c);
            Vec3 G_R = make_G(cell.R, cell.zR);
            Vec3 G_L = make_G(cell.L, cell.zL);
            Vec3 H_T = make_H(cell.T, cell.zT);
            Vec3 H_B = make_H(cell.B, cell.zB);
            Vec3 U_half_step_c = U_c + 0.5 * dt * ((G_L - G_R) / dx + (H_B - H_T) / dy);
            half_step_c = make_Point3_from_U(U_half_step_c);
            //
            
        }
    }

}


void CabaretScheme2D::third_phase() {

    //Проход по всем ячейкам
    for(int i = 0 ; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            //Решаем третью фазу в векторном виде
            CellView cell = get_cell(i, j);
            Point3& half_step_c = half_step_center[i][j];
            Vec3 U_half_step_c = make_U(half_step_c);
            Vec3 G_R = make_G(cell.R, cell.zR);
            Vec3 G_L = make_G(cell.L, cell.zL);
            Vec3 H_T = make_H(cell.T, cell.zT);
            Vec3 H_B = make_H(cell.B, cell.zB);
            Vec3 U_c = U_half_step_c + 0.5 * dt * ((G_L - G_R) / dx + (H_B - H_T) / dy);
            cell.c = make_Point3_from_U(U_c);
            //
            
        }
    }
}


#endif