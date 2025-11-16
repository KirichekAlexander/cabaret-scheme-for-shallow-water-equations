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


// Синонимы для массивов точек (u, v, h)
using Field1D = std::vector<Point3>;
using Field2D = std::vector<std::vector<Point3>>;


//Структура для просмотра ячейки
struct CellView {
    Point3 &L;
    Point3 &R;
    Point3 &T;
    Point3 &B;
    Point3 &c;
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

    //метод вычисления
    void compute();


private:
    //метод строит пространтсвенную сетку
    void build_grid(Row& grid_center, Row& grid_face,
                    int n, double delta,
                    double l1);
    
    //метод инициализации данных
    void init_data();

    double compute_time_step();

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


void CabaretScheme2D::compute() {



}


#endif