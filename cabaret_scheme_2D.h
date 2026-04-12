#ifndef CABARET_SCHEME_2D_H
#define CABARET_SCHEME_2D_H


#include <functional>
#include "read_matrix.h"
#include "functions.h"
#include "boundary_conditions.h"
#include "file_manager.h"



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
                             Func2D h0, Func2D z0,
                             std::string const& file_name);

    void compute(); // метод вычисления


private:
    // метод строит пространтсвенную сетку
    void build_grid(Row& grid_center, Row& grid_face,
                    int n, double delta,
                    double l1);
    
    void init_data(); // метод инициализации данных
    // заполнение выбранных полей
    void init_fields(int i1, int i2,
                     int j1, int j2,
                     Matrix& z, Field2D& points,
                     Row const& x, Row const& y);
    void init_boundaries(Side side); // инициализация границ 

    // создать ячейку по заданным массивам
    CellView make_cell(Field2D& center_field,
                       Field2D& face_x_field,
                       Field2D& face_y_field,
                       Matrix&  z_center_field,
                       Matrix&  z_face_x_field,
                       Matrix&  z_face_y_field,
                       int i, int j);

    CellView get_cell(int i, int j); // просмотр ячейки текущего временного слоя по индексам
    CellView get_next_cell(int i, int j); // просмотр ячейки следующего временного слоя, при этом центр на половинном слое, по индексам

    double compute_time_step(); // вычисление временного шага

    Vec3r make_U(Point3 const& p);
    Vec3r make_G(Point3 const& p, double z, double hL, double hR);
    Vec3r make_H(Point3 const& p, double z, double hT, double hB);
    Point3 make_Point3_from_U(Vec3r const& U);
    void first_phase(); // первая фаза

    // функции второй фазы
    void second_phase(); // вторая фаза
    void compute_boundaries(Side side); // вычисление границ
    Point3 make_Point3_from_invariants_x(Inv3 const& I, double z);
    Point3 make_Point3_from_invariants_y(Inv3 const& I, double z);
    Inv3 choose_invariants(Vec3r const& char_speeds_sum, Inv3 const& I_minus, Inv3 const& I_plus);
    Inv3 compute_invariants(int i, int j, Side side); // вычисление инвариантов
    Inv3 make_invariants_x(Point3 const& p, double z, double h_half_step); // вычисление инвариантов вдоль оси x
    Inv3 make_invariants_y(Point3 const& p, double z, double h_half_step); // вычисление инвариантов вдоль оси y
    Vec3r make_char_speeds_x(Point3 const& p); // вычисление хар-их скоростей вдоль оси x
    Vec3r make_char_speeds_y(Point3 const& p); // вычисление хар-их скоростей вдоль оси y
    //

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
    Field2D next_face_x; // (u, v, h) на вертикальных гранях на новом временном слою
    Field2D next_face_y; // (u, v, h) на горизонтальных гранях на новом временном слою


    std::string file_name;
    FileManager2D file_manager;
    int num_layer;

};


//Реализация CabaretScheme2D
CabaretScheme2D::CabaretScheme2D(double CFL,
                                 int nx, int ny,
                                 double l1x, double l2x,
                                 double l1y, double l2y,
                                 double T,
                                 Func2D u0, Func2D v0,
                                 Func2D h0, Func2D z0,
                                 std::string const& file_name) 
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
    , next_face_x(nx + 1, Field1D(ny))
    , next_face_y(nx, Field1D(ny + 1))
    , file_name(file_name)
    , file_manager(nx, ny)
    , num_layer(0)
{
    build_grid(x_center, x_face_x, nx, dx, l1x); // вычисление точек по x
    build_grid(y_center, y_face_y, ny, dy, l1y); // вычисление точек по y
    init_data(); // инциализация по начальным данным
    file_manager.save_layer(file_name, 0, num_layer, t, x_center, y_center, z_center, center);
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

    //Инициализация внутренних точек
    init_fields(0, nx, 0, ny, z_center, center, x_center, y_center); // инициализация центральных точек
    init_fields(1, nx, 0, ny, z_face_x, face_x, x_face_x, y_center); // инициализация вертикальных граней
    init_fields(0, nx, 1, ny, z_face_y, face_y, x_center, y_face_y); // инициализация горизонтальных граней
    //


    // инициализация границ (приоритет к граничным данным)
    init_boundaries(Side::LEFT);
    init_boundaries(Side::RIGHT);
    init_boundaries(Side::TOP);
    init_boundaries(Side::BOTTOM);

}


void CabaretScheme2D::init_fields(int i1, int i2,
                                  int j1, int j2,
                                  Matrix& z, Field2D& points,
                                  Row const& x, Row const& y) {

    for(int i = i1; i < i2; ++i) {
        for(int j = j1; j < j2; ++j) {
            z[i][j] = z0(x[i], y[j]);
            points[i][j].u = u0(x[i], y[j]);
            points[i][j].v = v0(x[i], y[j]);
            // points[i][j].h = h0(x[i], y[j]);
            points[i][j].h = 1.0 - z[i][j]; // well-balanced init
        }
    }

}


void CabaretScheme2D::init_boundaries(Side side) {

    //пока что реализация только граничных условий непротекания
    switch (side)
    {

    case Side::LEFT: // крайняя левая вертикальная грань

        for(int j = 0; j < ny; ++j) {

            z_face_x[0][j] = z0(x_face_x[0], y_center[j]);
            
            face_x[0][j].u = 0.0;
            face_x[0][j].v = v0(x_face_x[0], y_center[j]);
            // face_x[0][j].h = h0(x_face_x[0], y_center[j]);
            face_x[0][j].h = 1.0 - z_face_x[0][j]; // well-balanced init

        }    

        break;
    
    case Side::RIGHT: // карйняя правая вертикальная грань

        for(int j = 0; j < ny; ++j) {

            z_face_x[nx][j] = z0(x_face_x[nx], y_center[j]);
            
            face_x[nx][j].u = 0.0;
            face_x[nx][j].v = v0(x_face_x[nx], y_center[j]);
            // face_x[nx][j].h = h0(x_face_x[nx], y_center[j]);
            face_x[nx][j].h = 1.0 - z_face_x[nx][j]; // well-balanced init

        }    

        break;

    case Side::TOP: // крайняя верхняя горизонтальная грань

        for(int i = 0; i < nx; ++i) {

            z_face_y[i][ny] = z0(x_center[i], y_face_y[ny]);
            
            face_y[i][ny].u = u0(x_center[i], y_face_y[ny]);
            face_y[i][ny].v = 0.0;
            // face_y[i][ny].h = h0(x_center[i], y_face_y[ny]);
            face_y[i][ny].h = 1.0 - z_face_y[i][ny]; // well-balanced init

        }    

        break;

    case Side::BOTTOM: // крайняя нижняя грань

    for(int i = 0; i < nx; ++i) {

            z_face_y[i][0] = z0(x_center[i], y_face_y[0]);
            
            face_y[i][0].u = u0(x_center[i], y_face_y[0]);
            face_y[i][0].v = 0.0;
            // face_y[i][0].h = h0(x_center[i], y_face_y[0]);
            face_y[i][0].h = 1.0 - z_face_y[i][0]; // well-balanced init

        }    

        break;

    }

}


CellView CabaretScheme2D::make_cell(Field2D& center_field,
                                    Field2D& face_x_field,
                                    Field2D& face_y_field,
                                    Matrix&  z_center_field,
                                    Matrix&  z_face_x_field,
                                    Matrix&  z_face_y_field,
                                    int i, int j) {

    return CellView{face_x_field[i][j],
                    face_x_field[i + 1][j],
                    face_y_field[i][j + 1],
                    face_y_field[i][j],
                    center_field[i][j],
                    z_face_x_field[i][j],
                    z_face_x_field[i + 1][j],
                    z_face_y_field[i][j + 1],
                    z_face_y_field[i][j],
                    z_center_field[i][j]};                              
                    
}



CellView CabaretScheme2D::get_cell(int i, int j) {
    return make_cell(center, face_x, face_y, z_center, z_face_x, z_face_y, i, j);
}


CellView CabaretScheme2D::get_next_cell(int i, int j) {
    return make_cell(half_step_center, next_face_x, next_face_y, z_center, z_face_x, z_face_y, i, j);
}


void CabaretScheme2D::compute() {

    while(t < T) {

        dt = compute_time_step();

        first_phase();
        second_phase();
        third_phase();

        ++num_layer;
        std::cout << "t: " << t << " layer: " << num_layer << std::endl;
        t += dt;

        file_manager.save_layer(file_name, 0, num_layer, t, x_center, y_center, z_center, center);


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


Vec3r CabaretScheme2D::make_U(Point3 const& p) {
    return Vec3r{p.h,        // h
                p.h * p.u,  // hu
                p.h * p.v}; // hv
}


Vec3r CabaretScheme2D::make_G(Point3 const& p, double z, double hL, double hR) {
    return Vec3r{p.h * p.u,                                                     // hu
                p.h * sqr(p.u) + 0.5 * g * (hL + hR) * (p.h + z) ,             // hu^2 + g(H^2 - 2zH) / 2
                p.h * p.v * p.u};                                              // hvu
}


Vec3r CabaretScheme2D::make_H(Point3 const& p, double z, double hT, double hB) {
    return Vec3r{p.h * p.v,                                                      // hv
                p.h * p.v * p.u,                                                // hvu
                p.h * sqr(p.v) + 0.5 * g * (hT + hB) * (p.h + z)};              // hv^2 + g(H^2 - 2zH) / 2
}



Point3 CabaretScheme2D::make_Point3_from_U(Vec3r const& U) {
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
            Vec3r U_c = make_U(cell.c);
            Vec3r G_R = make_G(cell.R, cell.zR, cell.L.h, cell.R.h);
            Vec3r G_L = make_G(cell.L, cell.zL, cell.L.h, cell.R.h);
            Vec3r H_T = make_H(cell.T, cell.zT, cell.T.h, cell.B.h);
            Vec3r H_B = make_H(cell.B, cell.zB, cell.T.h, cell.B.h);
            Vec3r U_half_step_c = U_c + 0.5 * dt * ((G_L - G_R) / dx + (H_B - H_T) / dy);
            half_step_c = make_Point3_from_U(U_half_step_c);
            //
            
        }
    }

}


void CabaretScheme2D::second_phase() {
    //Вычисление границ
    compute_boundaries(Side::LEFT);
    compute_boundaries(Side::RIGHT);
    compute_boundaries(Side::TOP);
    compute_boundaries(Side::BOTTOM);
    //

    //Вычисление потоковых величин на вертикальных гранях
    for(int i = 1; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
            Inv3 invariants_left_cell = compute_invariants(i - 1, j, Side::LEFT);
            Inv3 invariants_right_cell = compute_invariants(i, j, Side::RIGHT);
            Vec3r char_speeds_center_half_step_left = make_char_speeds_x(half_step_center[i - 1][j]);
            Vec3r char_speeds_center_half_step_right = make_char_speeds_x(half_step_center[i][j]);
            Vec3r sum_char_speeds = char_speeds_center_half_step_left + char_speeds_center_half_step_right;
            Inv3 invariants = choose_invariants(sum_char_speeds, invariants_left_cell, invariants_right_cell);
            CellView cell_next = get_next_cell(i, j);
            cell_next.L = make_Point3_from_invariants_x(invariants, cell_next.zL);
        }
    }

    //Вычисление потоковых величин на горизонтальных гранях
    for(int j = 1; j < ny; ++j) {
        for(int i = 0; i < nx; ++i) {
            Inv3 invariants_bottom_cell = compute_invariants(i, j - 1, Side::BOTTOM);
            Inv3 invariants_top_cell = compute_invariants(i, j, Side::TOP);
            Vec3r char_speeds_center_half_step_left = make_char_speeds_y(half_step_center[i][j - 1]);
            Vec3r char_speeds_center_half_step_right = make_char_speeds_y(half_step_center[i][j]);
            Vec3r sum_char_speeds = char_speeds_center_half_step_left + char_speeds_center_half_step_right;
            Inv3 invariants = choose_invariants(sum_char_speeds, invariants_bottom_cell, invariants_top_cell);
            CellView cell_next = get_next_cell(i, j);
            cell_next.B = make_Point3_from_invariants_y(invariants, cell_next.zB);
        }
    }

    std::swap(face_x, next_face_x);
    std::swap(face_y, next_face_y);
    
}


void CabaretScheme2D::compute_boundaries(Side side) {

    //пока что реализация только граничных условий непротекания
    switch (side)
    {

    case Side::LEFT: // крайняя левая вертикальная грань

        for(int j = 0; j < ny; ++j) {
            
            next_face_x[0][j].u = 0.0; // нормальная скорость ноль

            Inv3 invariants_right_cell = compute_invariants(0, j, Side::RIGHT); // вычисление приходящих инвариантов
            Vec3r char_speeds_center_half_step_right = make_char_speeds_x(half_step_center[0][j]); // вычисление хар-их скоростей
            next_face_x[0][j].h = (-invariants_right_cell.b.I) * invariants_right_cell.a.c / g - z_face_x[0][j]; // h = I_2x^2/(4g)
            if (char_speeds_center_half_step_right.c < 0) {
                next_face_x[0][j].v = invariants_right_cell.c.I;
            //Иначе берём инвариант с предыдущего временного слоя
            } else {
                next_face_x[0][j].v = half_step_center[0][j].v;
            }

        }    

        break;
    
    case Side::RIGHT: // карйняя правая вертикальная грань

        for(int j = 0; j < ny; ++j) {

            next_face_x[nx][j].u = 0.0; // нормальная скорость ноль

            Inv3 invariants_left_cell = compute_invariants(nx - 1, j, Side::LEFT); // вычисление приходящих инвариантов
            Vec3r char_speeds_center_half_step_left = make_char_speeds_x(half_step_center[nx - 1][j]); // вычисление хар-их скоростей
            next_face_x[nx][j].h = invariants_left_cell.a.I * invariants_left_cell.a.c / g - z_face_x[nx][j]; // h = I_1x^2/(4g)
            if (char_speeds_center_half_step_left.c > 0) {
                next_face_x[nx][j].v = invariants_left_cell.c.I;
            //Иначе берём инвариант с предыдущего временного слоя
            } else {
                next_face_x[nx][j].v = half_step_center[nx - 1][j].v;
            }

        }    

        break;

    case Side::TOP: // крайняя верхняя горизонтальная грань

        for(int i = 0; i < nx; ++i) {

            next_face_y[i][ny].v = 0.0; // нормальная скорость ноль

            Inv3 invariants_bottom_cell = compute_invariants(i, ny - 1, Side::BOTTOM); // вычисление приходящих инвариантов
            Vec3r char_speeds_center_half_step_bottom = make_char_speeds_y(half_step_center[i][ny - 1]); // вычисление хар-их скоростей
            next_face_y[i][ny].h = invariants_bottom_cell.a.I * invariants_bottom_cell.a.c / g - z_face_y[i][ny]; // h = I_1y^2/(4g)
            if (char_speeds_center_half_step_bottom.c > 0) {
                next_face_y[i][ny].u = invariants_bottom_cell.c.I;
            //Иначе берём инвариант с предыдущего временного слоя
            } else {
                next_face_y[i][ny].u = half_step_center[i][ny - 1].u;
            }

        }    

        break;

    case Side::BOTTOM: // крайняя нижняя грань

        for(int i = 0; i < nx; ++i) {

            next_face_y[i][0].v = 0.0; // нормальная скорость ноль

            Inv3 invariants_top_cell = compute_invariants(i, 0, Side::TOP); // вычисление приходящих инвариантов
            Vec3r char_speeds_center_half_step_top = make_char_speeds_y(half_step_center[i][0]); // вычисление хар-их скоростей
            next_face_y[i][0].h = (-invariants_top_cell.b.I) * invariants_top_cell.b.c / g - z_face_y[i][0]; // h = I_1y^2/(4g)
            if (char_speeds_center_half_step_top.c < 0) {
                next_face_y[i][0].u = invariants_top_cell.c.I;
            //Иначе берём инвариант с предыдущего временного слоя
            } else {
                next_face_y[i][0].u = half_step_center[i][0].u;
            }

        }      

        break;

    }

}



Point3 CabaretScheme2D::make_Point3_from_invariants_x(Inv3 const& I, double z){
    double H = (I.a.I - I.b.I) / (g * (1 / I.a.c + 1 / I.b.c));
    return Point3{I.a.I - g / I.a.c * H,
                  I.c.I,
                  H - z};
}


Point3 CabaretScheme2D::make_Point3_from_invariants_y(Inv3 const& I, double z) {
    double H = (I.a.I - I.b.I) / (g * (1 / I.a.c + 1 / I.b.c));
    return Point3{I.c.I,
                  I.a.I - g / I.a.c * H,
                  H - z};
}


Inv3 CabaretScheme2D::choose_invariants(Vec3r const& char_speeds_sum, Inv3 const& I_minus, Inv3 const& I_plus) {
    return Inv3{(char_speeds_sum.a > 0 ? I_minus.a : I_plus.a),
                (char_speeds_sum.b > 0 ? I_minus.b : I_plus.b),
                (char_speeds_sum.c > 0 ? I_minus.c : I_plus.c)};
}


Inv3 CabaretScheme2D::compute_invariants(int i, int j, Side side) {

    CellView cell = get_cell(i, j);
    Inv3 I_stream_n;      // потоковое состояние на слое n (из стороны side)
    Inv3 I_center_half;   // центровое состояние на слое n+1/2
    Inv3 I_center_n;      // центровое на слое n
    Inv3 I_opp_stream_n;  // потоковое с противоположной стороны на слое n


    // Вычисление инвариантов вдоль оси x
    if (side == Side::LEFT or side == Side::RIGHT) {
        Point3 const& p = (side == Side::LEFT ? cell.L : cell.R);
        double z = (side == Side::LEFT ? cell.zL : cell.zR);
        I_stream_n = make_invariants_x(p, z, half_step_center[i][j].h);
        I_center_half = make_invariants_x(half_step_center[i][j], cell.zc, half_step_center[i][j].h);
        I_center_n = make_invariants_x(cell.c, cell.zc, half_step_center[i][j].h);
        Point3 const& p_other = (side == Side::LEFT ? cell.R : cell.L);
        z = (side == Side::LEFT ? cell.zR : cell.zL);
        I_opp_stream_n = make_invariants_x(p_other, z, half_step_center[i][j].h);

    // Иначе вычисляем инварианты вдоль оси y
    } else {
        Point3 const& p = (side == Side::BOTTOM ? cell.B : cell.T);
        double z = (side == Side::BOTTOM ? cell.zB : cell.zT);
        I_stream_n = make_invariants_y(p, z, half_step_center[i][j].h);
        I_center_half = make_invariants_y(half_step_center[i][j], cell.zc, half_step_center[i][j].h);
        I_center_n = make_invariants_y(cell.c, cell.zc, half_step_center[i][j].h);
        Point3 const& p_other = (side == Side::BOTTOM ? cell.T : cell.B);
        z = (side == Side::BOTTOM ? cell.zT : cell.zB);
        I_opp_stream_n = make_invariants_y(p_other, z, half_step_center[i][j].h);
    }

    // Инвариант приходящий на новый временной слой
    Inv3 I_stream_next = 2 * I_center_half - I_stream_n;
    // Нужно использовать принцип максимумов для инвариантов
    Inv3 G = (I_center_half - I_center_n) / (0.5 * dt) +
             (side == Side::LEFT or side == Side::RIGHT ? make_char_speeds_x(half_step_center[i][j]) / dx :
              make_char_speeds_y(half_step_center[i][j]) / dy) * 
              (side == Side::LEFT or side == Side::BOTTOM ? I_opp_stream_n - I_stream_n : I_stream_n - I_opp_stream_n);
    Inv3 m = min_invariants({I_stream_n, I_center_n, I_opp_stream_n}) + dt * G;
    Inv3 M = max_invariants({I_stream_n, I_center_n, I_opp_stream_n}) + dt * G;
    //
    I_stream_next = min_invariants({max_invariants({m, I_stream_next}), M});
    return I_stream_next;

}


Inv3 CabaretScheme2D::make_invariants_x(Point3 const& p, double z, double h_half_step) {
    double c = std::sqrt(g * h_half_step);
    return Inv3{InvC{p.u + g / c * (z + p.h), c},
                InvC{p.u - g / c * (z + p.h), c},
                InvC{p.v, c}};
}


Inv3 CabaretScheme2D::make_invariants_y(Point3 const& p, double z, double h_half_step) {
    double c = std::sqrt(g * h_half_step);
    return Inv3{InvC{p.v + g / c * (z + p.h), c},
                InvC{p.v - g / c * (z + p.h), c},
                InvC{p.u, c}};
}


Vec3r CabaretScheme2D::make_char_speeds_x(Point3 const& p) {
    return Vec3r {p.u + std::sqrt(g * p.h),
                  p.u - std::sqrt(g * p.h),
                  p.u};
}


Vec3r CabaretScheme2D::make_char_speeds_y(Point3 const& p) {
    return Vec3r {p.v + std::sqrt(g * p.h),
                  p.v - std::sqrt(g * p.h),
                  p.v};
}



void CabaretScheme2D::third_phase() {

    //Проход по всем ячейкам
    for(int i = 0 ; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            //Решаем третью фазу в векторном виде
            CellView cell = get_cell(i, j);
            Point3& half_step_c = half_step_center[i][j];
            Vec3r U_half_step_c = make_U(half_step_c);
            Vec3r G_R = make_G(cell.R, cell.zR, cell.L.h, cell.R.h);
            Vec3r G_L = make_G(cell.L, cell.zL, cell.L.h, cell.R.h);
            Vec3r H_T = make_H(cell.T, cell.zT, cell.T.h, cell.B.h);
            Vec3r H_B = make_H(cell.B, cell.zB, cell.T.h, cell.B.h);
            Vec3r U_c = U_half_step_c + 0.5 * dt * ((G_L - G_R) / dx + (H_B - H_T) / dy);
            cell.c = make_Point3_from_U(U_c);
            //
            
        }
    }
}


#endif