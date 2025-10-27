#ifndef CABARET_SCHEME_H
#define CABARET_SCHEME_H


#include <functional>
#include "read_matrix.h"
#include "functions.h"
#include "boundary_conditions.h"
#include "file_manager.h"

#include "TECIO.h"
#include "TECXXX.h"


//Тип сетки равномерная или неравномерная
enum Type_grid {EVEN, UNEVEN};


//Тип инварианта из левой или правой ячейки
enum Type_invariant {LEFT_CELL, RIGHT_CELL};


//Тип поведения автомодельного решения
enum Type_behavior_automodel_solution {DISCHARGE_WAVE, HYDRAULIC_JUMP};


//Тип потоковой величины
enum Type_value {H, U};


/*
*Схема Кабаре для уравнений мелкой воды над неровным дном
*для несжимаемой идеальной жидкости с одной точкой разрыва для граничных условий непротекания:
 * - - - * - - - * - - - * - - - *
 |             / | \             |
 |           /   |   \           |
 |         /     |     \         |
 |       *       |       *       |
 |     /         |         \     |
 |   /           |           \   |
 | /             |             \ |
 * - - - * - - - * - - - * - - - *
*/
template<BoundaryType LeftBT, BoundaryType RightBT>
class Cabaret_scheme {
public:


    /*
    *Инициализация данных, на вход вычислениям: тип сетки по x (по t сетка неравномерная и строится по специальному правилу), CFL, точка разрыва,
    *u_1, u_2, h_1, h_2 (начальные функции сокрости и глубины), z(x) (поверхность стационарного дна),начальная точка по x, конечная точка по x,
    *начальная точка по t, конечная точка по t, количество точек между начальной точкой
    *и точкой разрыва, количество
    *точек между разрывом и конечной точкой.
    */
   Cabaret_scheme(Type_grid type_grid, Boundary<LeftBT> left_boundary, Boundary<RightBT> right_boundary
    , double CFL, double break_pt, double (* u_1)(double), double (* u_2)(double), double (* h_1)(double), double (* h_2)(double)
    , double (* z)(double), double start_x_pt, double end_x_pt, double start_t_pt, double end_t_pt, int cnt_x_pts_before_break_pt
    , int cnt_x_pts_after_break_pt, std::string file_name, bool continuity, std::pair<double, double> (* analytical_solution)(double, double));

    
    //Запуск вычислений
    void compute();


    //Обработка случаев гидродинамических прыжков
    void process_hydraulic_jump(double&, double&, double, double, double, double);


    //Обработка случаев волн рарзряжения
    void process_discharge_wave(double&, double&, double, double, double, double, double, double, double, bool);


    //Автомодельное решение задачи о распаде разрыва
    double automodeling_solution(double, double);


    //Возвращение конкретной потоковой величины
    Row& stream_value(Type_value);


    //Возвращение конкретной консервативной величины
    Row& conservative_value(Type_value);


private:
    
    //Вычисление сетки по x
    void compute_x_grid();

    //Вычисление равномерной сетки, на входе: вектор, начальный индекс, количество точек, начальная точка, конечная точка
    void compute_even_grid(Row&, int, int, double, double);
    

    //Фазы

    //Здесь происходят основные вычисления(запуски фаз)
    void phases();

    //инициализация границ
    void init_boundaries(Row& stream_u, Row& stream_h);

    
    //Вычисление начальных потоковых и консервативных значений для u и h
    void compute_start_streaming_and_conservative_values(Row&, Row&, Row&, Row&);

    
    //Вычисление текущего шага по t
    void compute_tau(int);


    //Первая фаза    
    void first_phase(int);


    void compute_boundaries(int);
    
    //Вторая фаза    
    void second_phase(int);

    
    //Вычисление инварианта, возвращает пару инвариант и коэфф g/c
    std::pair<double, double> compute_invariant(Type_invariant, int, int);


    //Третья фаза
    void third_phase(int);

    //


    void save_error(std::string const&, int, Row&, Row&);

    //Начальные данные

    //Тип сетки
    Type_grid type_grid;


    //Граничные условия
    Boundary<LeftBT> left_boundary;
    Boundary<RightBT> right_boundary;

    
    //CFL
    double CFL;


    //Точка разрыва
    double break_pt;

    
    //Условие Коши до точки разрыва для u
    double (* u_1)(double);

    
    //Условие Коши после точки разрыва для u
    double (* u_2)(double);

    
    //Условие Коши до точки разрыва для h
    double (* h_1)(double);

    
    //Условие Коши после точки разрыва для h
    double (* h_2)(double);

    
    //Функция дна
    double (* z)(double);

    
    //Начальная точка по x
    double start_x_pt;


    //Конечная точка по x
    double end_x_pt;


    //Начальная точка по t
    double start_t_pt;


    //Конечная точка по t
    double end_t_pt;


    //Количество точек по x до точки разрыва, не учитывая левую граничную точку
    int cnt_x_pts_before_break_pt;

    
    //Количество точек по x после точки разрыва, не учитывая правую граничную точку
    int cnt_x_pts_after_break_pt;

    
    //Общее количество точек по x
    int cnt_x_pts;


    //Индекс точки разрыва
    int idx_break_pt;


    //Количество точек по времени
    int cnt_t_pts;

    //


    //Сетки

    //Сетка по x
    Row x_grid;


    Row conservative_x_grid;


    //Сетка по t
    Row t_grid;

    
    //Сетка по потокам u
    Row stream_u_grid;

    
    //Сетка по потокам h
    Row stream_h_grid;


    //Новый временной слой u
    Row new_stream_u_grid;


    //Новый временной слой h
    Row new_stream_h_grid;

    //Сетка по консервативным u
    Row conservative_u_grid;

    
    //Сетка по консервативным h
    Row conservative_h_grid;

    
    //Центральные консервативные u
    Row center_conservative_u_grid;


    //Центральные консервативные h
    Row center_conservative_h_grid;

    
    //z в узлах
    Row stream_z_grid;

    
    //z в полуцелых точках
    Row conservative_z_grid;


    //сеткам по шагам по x
    Row h_grid;


    //сетка по шагам по t
    Row tau_grid;


    //Константный G на левой границе
    double const_g_left_border;


    //Константный G на правой границе
    double const_g_right_border;


    //Константный инвариант на левой границе
    double const_invariant_left_border;


    //Константный инвариант на правой границе
    double const_invariant_right_border;

    //

    
    //Имя файла куда записываются данные
    std::string file_name;
    FileManager file_manager;

    //Решаем непрерывную задачу?
    bool continuity;


    //Величины нужные для автомодельного решения
    
    //u в центральной зоне
    double automodel_u;


    //h в центральной зоне
    double automodel_h;

    
    //u в левом постоянном течении
    double automodel_u_1;


    //u в правом потстоянном течении
    double automodel_u_2;


    //h в левом постоянном течении
    double automodel_h_1;


    //h в правом постоянном течении
    double automodel_h_2;


    //c в центральной зоне
    double automodel_c;


    //константный инвариант слева
    double automodel_const_invariant_left;


    //констатный инвариант слева
    double automodel_const_invariant_right;

    
    //c в левом постоянном течении
    double automodel_c_1;


    //c в правом постоянном течении
    double automodel_c_2;


    double automodel_S_1;


    double automodel_S_2;


    Row automodel_conservative_h_grid;


    Row automodel_conservative_u_grid;


    Row automodel_stream_h_grid;


    Row automodel_stream_u_grid;


    //Функция аналитического решения
            /*  h       u
                |       | */
    std::pair<double, double> (* analytical_solution)(double, double);

    //Сетки для аналитического решения
    Row analytical_solution_conservative_h_grid;


    Row analytical_solution_conservative_u_grid;


    Row analytical_solution_stream_h_grid;


    Row analytical_solution_stream_u_grid;


    //Функции нужные для вычисления c в центральной зоне
    double phi_k(double, int);


    double F(double);


    double derivative_phi_k(int);


    double derivative_F();

};




//Реализация класса схемы Кабаре


template<BoundaryType LeftBT, BoundaryType RightBT>
Cabaret_scheme<LeftBT, RightBT>::Cabaret_scheme(Type_grid type_grid, Boundary<LeftBT> left_boundary
    , Boundary<RightBT> right_boundary, double CFL, double break_pt, double (* u_1)(double), double (* u_2)(double)
    , double (* h_1)(double), double (* h_2)(double), double (* z)(double), double start_x_pt, double end_x_pt, double start_t_pt
    , double end_t_pt, int cnt_x_pts_before_break_pt, int cnt_x_pts_after_break_pt, std::string file_name, bool continuity
    , std::pair<double, double> (* analytical_solution)(double, double))
    //Инициализация всего
        : type_grid(type_grid)
        , left_boundary(left_boundary)
        , right_boundary(right_boundary)
        , CFL(CFL)
        , break_pt(break_pt)
        , u_1(u_1)
        , u_2(u_2)
        , h_1(h_1)
        , h_2(h_2)
        , z(z)
        , start_x_pt(start_x_pt)
        , end_x_pt(end_x_pt)
        , start_t_pt(start_t_pt)
        , end_t_pt(end_t_pt)
        , cnt_x_pts_before_break_pt(cnt_x_pts_before_break_pt)
        , cnt_x_pts_after_break_pt(cnt_x_pts_after_break_pt)
        , cnt_x_pts(3 + cnt_x_pts_after_break_pt + cnt_x_pts_before_break_pt)
        , idx_break_pt(cnt_x_pts_before_break_pt + 1)
        , cnt_t_pts(0)
        , x_grid(cnt_x_pts)
        , conservative_x_grid(cnt_x_pts - 1)
        , t_grid(1, start_t_pt)
        , stream_u_grid(cnt_x_pts)
        , stream_h_grid(cnt_x_pts)
        , new_stream_u_grid(cnt_x_pts)
        , new_stream_h_grid(cnt_x_pts)
        , conservative_u_grid(cnt_x_pts - 1)
        , conservative_h_grid(cnt_x_pts - 1)
        , center_conservative_u_grid(cnt_x_pts - 1)
        , center_conservative_h_grid(cnt_x_pts - 1)
        , stream_z_grid(cnt_x_pts)
        , conservative_z_grid(cnt_x_pts - 1)
        , h_grid(cnt_x_pts - 1) 
        , const_g_left_border(g / std::sqrt(g * h_1(start_x_pt)))
        , const_g_right_border(g / std::sqrt(g * h_2(end_x_pt)))
        , const_invariant_left_border(u_1(start_x_pt) + const_g_left_border * (h_1(start_x_pt) + z(start_x_pt)))
        , const_invariant_right_border(u_2(end_x_pt) - const_g_right_border * (h_2(end_x_pt) + z(end_x_pt)))
        , file_name(file_name)
        , continuity(continuity)
        , automodel_u(0.0)
        , automodel_h(0.0)
        , automodel_u_1(0.0)
        , automodel_u_2(0.0)
        , automodel_h_1(0.0)
        , automodel_h_2(0.0)
        , automodel_c(0.0)
        , automodel_const_invariant_left(0.0)
        , automodel_const_invariant_right(0.0)
        , automodel_c_1(0.0)
        , automodel_c_2(0.0)
        , automodel_S_1(0.0)
        , automodel_S_2(0.0)
        , automodel_conservative_h_grid(cnt_x_pts - 1)
        , automodel_conservative_u_grid(cnt_x_pts - 1)
        , automodel_stream_h_grid(cnt_x_pts)
        , automodel_stream_u_grid(cnt_x_pts)
        , analytical_solution(analytical_solution)
    //
{ 
}

template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute() {

    //Вычисление сетки по x
    compute_x_grid();

    try {

        //Запуск фаз
        phases();

    } catch(std::exception const& err) {

        //Можем поймать сверхзвуковое течение
        std::cout << err.what() << std::endl;

    }

}

template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute_x_grid() {

    if (type_grid == EVEN) {
        
        //Вычисление равномерной сетки по x слева от разрыва
        compute_even_grid(x_grid, 0, cnt_x_pts_before_break_pt + 2, start_x_pt, break_pt);
        
        //Вычисление равномерной сетки по x справа от разрыва
        compute_even_grid(x_grid, cnt_x_pts_before_break_pt + 1, cnt_x_pts_after_break_pt + 2, break_pt, end_x_pt);

    } else {

        // Пока без неравномерных сеток

    }

    //Вычисление шагов сетки по x
    for (int i = 0; i < cnt_x_pts - 1; ++i) {
        h_grid[i] = x_grid[i + 1] - x_grid[i];
        conservative_x_grid[i] = x_grid[i] + h_grid[i] / 2.0;
    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute_even_grid(Row& grid, int start_idx, int cnt_pts, double start_pt, double end_pt) {

    //Вычисление равномерной сетки
    for (int i = 0; i < cnt_pts; ++i) {

        double t = static_cast<double>(i) / (cnt_pts - 1);
        grid[start_idx + i] = start_pt + t * (end_pt - start_pt);

    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::phases() {
    std::cout << "START COMPUTING" << std::endl;
    //Вычисляем начальные данные
    compute_start_streaming_and_conservative_values(stream_u_grid, stream_h_grid, conservative_u_grid, conservative_h_grid);
    if (analytical_solution) {

        analytical_solution_conservative_h_grid = conservative_h_grid;
        analytical_solution_conservative_u_grid = conservative_u_grid;
        analytical_solution_stream_h_grid = stream_h_grid;
        analytical_solution_stream_u_grid = stream_u_grid;

    }

    //Начальный индекс по времени нулевой
    int t_idx = 0;

    if (!file_name.empty()) {

        file_manager.init_file(file_name + ".plt", 0, cnt_x_pts - 1);
        file_manager.save_layer(t_grid[0], conservative_x_grid, conservative_h_grid, conservative_u_grid, conservative_z_grid);

    }


    //Цикл по всему времени
    while(t_grid[t_idx] != end_t_pt) {

        std::cout << "TIME LAYER: " << t_grid[t_idx] << std::endl;

        //На каждом шагу считаем шаг по времени
        compute_tau(t_idx);

        //Первая фаза
        first_phase(t_idx);

        //Вторая фаза
        second_phase(t_idx);

        //Третья фаза
        third_phase(t_idx);

        ++t_idx;

        if (tau_grid[t_idx - 1] != 0.0 and !file_name.empty()) {
            file_manager.save_layer(t_grid[t_idx], conservative_x_grid, conservative_h_grid, conservative_u_grid, conservative_z_grid);
        }

    }

    if (!file_name.empty()) {
        file_manager.end_file();
    }

    if (analytical_solution) {
        
        if (!file_name.empty()) {

            file_manager.init_file(file_name + "_analytical.plt", 0, cnt_x_pts - 1);
            file_manager.save_layer(t_grid[0], conservative_x_grid, analytical_solution_conservative_h_grid, analytical_solution_conservative_u_grid
                                  , conservative_z_grid);

        }
        std::pair<double, double> analytical_solution_values;
        cnt_t_pts = t_grid.size();
        for(int t_idx = 1; t_idx < cnt_t_pts; ++t_idx) {

            for(int i = 0; i < (cnt_x_pts - 1); ++i) {

                analytical_solution_values = analytical_solution(x_grid[i], t_grid[t_idx]);
                analytical_solution_stream_h_grid[i] = analytical_solution_values.first;
                analytical_solution_stream_u_grid[i] = analytical_solution_values.second;

                analytical_solution_values = analytical_solution(conservative_x_grid[i], t_grid[t_idx]);
                analytical_solution_conservative_h_grid[i] = analytical_solution_values.first;
                analytical_solution_conservative_u_grid[i] = analytical_solution_values.second;

            }
            analytical_solution_values = analytical_solution(x_grid[cnt_x_pts - 1], t_grid[t_idx]);
            analytical_solution_stream_h_grid[cnt_x_pts - 1] = analytical_solution_values.first;
            analytical_solution_stream_u_grid[cnt_x_pts - 1] = analytical_solution_values.second;
        
            if (tau_grid[t_idx - 1] != 0.0 and !file_name.empty()) {
                file_manager.save_layer(t_grid[t_idx], conservative_x_grid, analytical_solution_conservative_h_grid
                                      , analytical_solution_conservative_u_grid, conservative_z_grid);

            }

        }

        if (!file_name.empty()) {
            file_manager.end_file();
        }

    }

    std::cout << "END COMPUTING" << std::endl;

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::init_boundaries(Row& stream_u, Row& stream_h) {

    if constexpr (LeftBT == BoundaryType::DEFAULT) {

        stream_h[0] = left_boundary.h;
        stream_u[0] = left_boundary.q / left_boundary.h;

    } else if constexpr (LeftBT == BoundaryType::NON_LEAK) {

        stream_u[0] = 0.0;
        stream_h[0] = h_1(x_grid[0]);

    } else if constexpr (LeftBT == BoundaryType::FREE_EXIT) {

        stream_u[0] = u_1(x_grid[0]);
        stream_h[0] = h_1(x_grid[0]);
        left_boundary.const_g = g / std::sqrt(g * stream_h[0]);
        left_boundary.const_invariant = stream_u[0] + left_boundary.const_g * (stream_h[0] + z(x_grid[0]));

    } else if constexpr (LeftBT == BoundaryType::FLUVIAL_H_GIVEN) {

        stream_h[0] = left_boundary.h;
        stream_u[0] = u_1(x_grid[0]);

    } else if constexpr (LeftBT == BoundaryType::FLUVIAL_FLUX_GIVEN) {

        stream_h[0] = h_1(x_grid[0]);
        stream_u[0] = left_boundary.q / stream_h[0];

    }


    if constexpr(RightBT == BoundaryType::DEFAULT) {

        stream_h[cnt_x_pts - 1] = right_boundary.h;
        stream_u[cnt_x_pts - 1] = -right_boundary.q / right_boundary.h;

    } else if constexpr (RightBT == BoundaryType::NON_LEAK) {

        stream_u[cnt_x_pts - 1] = 0.0;
        stream_h[cnt_x_pts - 1] = h_2(x_grid[cnt_x_pts - 1]);

    } else if constexpr (RightBT == BoundaryType::FREE_EXIT) {

        stream_u[cnt_x_pts - 1] = u_2(x_grid[cnt_x_pts - 1]);
        stream_h[cnt_x_pts - 1] = h_2(x_grid[cnt_x_pts - 1]);
        right_boundary.const_g = g / std::sqrt(g * stream_h[cnt_x_pts - 1]);
        right_boundary.const_invariant = stream_u[cnt_x_pts - 1] - right_boundary.const_g * (stream_h[cnt_x_pts - 1] + z(x_grid[cnt_x_pts - 1]));

    } else if constexpr (RightBT == BoundaryType::FLUVIAL_H_GIVEN) {

        stream_h[cnt_x_pts - 1] = right_boundary.h;
        stream_u[cnt_x_pts - 1] = u_2(x_grid[cnt_x_pts - 1]);

    } else if constexpr (RightBT == BoundaryType::FLUVIAL_FLUX_GIVEN) {

        stream_h[cnt_x_pts - 1] = h_2(x_grid[cnt_x_pts - 1]);
        stream_u[cnt_x_pts - 1] = -right_boundary.q / stream_h[cnt_x_pts - 1];

    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute_start_streaming_and_conservative_values(Row& stream_u, Row& stream_h, Row& conservative_u
                                                                                                                 , Row& conservative_h) {

    //Инициализация на краях
    init_boundaries(stream_u, stream_h);
    stream_z_grid[0] = z(x_grid[0]);
    stream_z_grid[cnt_x_pts - 1]= z(x_grid[cnt_x_pts - 1]);


    for (int i = 0; i < (cnt_x_pts - 1); ++i) {
        
        //Инициализация консервативных величин
        double half_h_i = h_grid[i] / 2.0;

        conservative_u[i] = (i < (cnt_x_pts_before_break_pt + 1) ? u_1(x_grid[i] + half_h_i) : u_2(x_grid[i] + half_h_i));
        conservative_h[i] = (i < (cnt_x_pts_before_break_pt + 1) ? h_1(x_grid[i] + half_h_i) : h_2(x_grid[i] + half_h_i));
        conservative_z_grid[i] = z(x_grid[i] + half_h_i);


        //Инициализация потоковых величин
        if (i > 0) {

            stream_u[i] = (i < idx_break_pt ? u_1(x_grid[i]) : 
                (i > idx_break_pt ? u_2(x_grid[i]) : (continuity ? u_1(x_grid[i]) : (conservative_u[i - 1] + conservative_u[i]) / 2.0)));

            stream_h[i] = (i < idx_break_pt ? h_1(x_grid[i]) : 
                (i > idx_break_pt ? h_2(x_grid[i]) : (continuity ? h_1(x_grid[i]) : (conservative_h[i - 1] + conservative_h[i]) / 2.0)));

            stream_z_grid[i] = z(x_grid[i]);

        }

    }

    std::cout << "START VALUES INITIALIZED" << std::endl;

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute_tau(int t_idx) {

    //Инициализация шага по времени максимальным значением double
    double tau = std::numeric_limits<double>::max();

    //Выбор шага по времени по специальному правилу
    for(int i = 0; i < (cnt_x_pts - 1); ++i) {

        // double cur_c = std::sqrt(g * stream_h_grid[t_idx][i]);
        double cur_c = std::sqrt(g * conservative_h_grid[i]);

        if (std::abs(conservative_u_grid[i]) > cur_c) {

            file_manager.end_file();    
            throw std::runtime_error("Сверхзвуковое течение в compute tau");   

        }

        // double cur_num = CFL * h_grid[i] / (std::abs(stream_u_grid[t_idx][i]) + cur_c); //!!! ВОЗМОЖНО ШАГ НУЖНО ОПРЕДЕЛЯТЬ ПО КОНСЕРВАТИВНЫМ ВЕЛИЧИНАМ
        double cur_num = CFL * h_grid[i] / (std::abs(conservative_u_grid[i]) + cur_c); 

        tau = std::min({cur_num, tau});

    }

    if (t_idx < 200) {
        tau = 0.0;
    }

    if ((end_t_pt - t_grid[t_idx]) < tau) {

        tau = end_t_pt - t_grid[t_idx];
        t_grid.push_back(end_t_pt);

    } else {

        t_grid.push_back(t_grid[t_idx] + tau);

    }

    tau_grid.push_back(tau);
}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::first_phase(int t_idx) {

    for(int i = 0; i < (cnt_x_pts - 1); ++i) {

        //Вычисление h в центре ячейки
        center_conservative_h_grid[i] = conservative_h_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * stream_u_grid[i] -
            stream_h_grid[i + 1] * stream_u_grid[i + 1]) / h_grid[i];

        //Вычисление u в центре ячейки
        center_conservative_u_grid[i] = (conservative_h_grid[i] * conservative_u_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * 
            sqr(stream_u_grid[i]) - stream_h_grid[i + 1] * sqr(stream_u_grid[i + 1]) + g * (stream_h_grid[i] + 
            stream_h_grid[i + 1]) / 2.0 * (stream_z_grid[i] + stream_h_grid[i] - stream_z_grid[i + 1] - stream_h_grid[i + 1])) /
            h_grid[i]) / center_conservative_h_grid[i];


        // center_conservative_u_grid[i] = (conservative_h_grid[i] * conservative_u_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * 
        //     sqr(stream_u_grid[i]) - stream_h_grid[i + 1] * sqr(stream_u_grid[i + 1]) + g * conservative_h_grid[i] * (stream_z_grid[i] + stream_h_grid[i] - stream_z_grid[i + 1] - stream_h_grid[i + 1])) /
        //     h_grid[i]) / center_conservative_h_grid[i];

    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::compute_boundaries(int t_idx) {

    std::pair<double, double> right_invariant_and_g = compute_invariant(RIGHT_CELL, t_idx, 0);
    if constexpr(LeftBT == BoundaryType::DEFAULT) {

        new_stream_h_grid[0] = left_boundary.h;
        new_stream_u_grid[0] = left_boundary.q / left_boundary.h;

    } else if constexpr (LeftBT == BoundaryType::NON_LEAK) {

        new_stream_u_grid[0] = 0.0;
        new_stream_h_grid[0] = -right_invariant_and_g.first / right_invariant_and_g.second - stream_z_grid[0];


    } else if constexpr (LeftBT == BoundaryType::FREE_EXIT) {

        new_stream_u_grid[0] = (left_boundary.const_invariant * right_invariant_and_g.second + right_invariant_and_g.first * 
            left_boundary.const_g) / (left_boundary.const_g + right_invariant_and_g.second);

        new_stream_h_grid[0] = (left_boundary.const_invariant - right_invariant_and_g.first) / (left_boundary.const_g + 
            right_invariant_and_g.second) - stream_z_grid[0];

    } else if constexpr (LeftBT == BoundaryType::FLUVIAL_H_GIVEN) {

        new_stream_h_grid[0] = left_boundary.h;
        new_stream_u_grid[0] = right_invariant_and_g.first + right_invariant_and_g.second * (left_boundary.h + stream_z_grid[0]);

    } else if constexpr (LeftBT == BoundaryType::FLUVIAL_FLUX_GIVEN) {

        new_stream_h_grid[0] = (-(right_invariant_and_g.first + stream_z_grid[0] * right_invariant_and_g.second) 
                                + std::sqrt(sqr(right_invariant_and_g.first + stream_z_grid[0] * right_invariant_and_g.second) 
                                + 4 * right_invariant_and_g.second * left_boundary.q)) 
                                / (2 * right_invariant_and_g.second);
        new_stream_u_grid[0] = left_boundary.q / new_stream_h_grid[0];

    }


    std::pair<double, double> left_invariant_and_g = compute_invariant(LEFT_CELL, t_idx, cnt_x_pts - 2);
    if constexpr(RightBT == BoundaryType::DEFAULT) {

        new_stream_h_grid[cnt_x_pts - 1] = right_boundary.h;
        new_stream_u_grid[cnt_x_pts - 1] = -right_boundary.q / right_boundary.h;

    } else if constexpr (RightBT == BoundaryType::NON_LEAK) {

        new_stream_u_grid[cnt_x_pts - 1] = 0.0;
        new_stream_h_grid[cnt_x_pts - 1] = left_invariant_and_g.first / left_invariant_and_g.second - stream_z_grid[cnt_x_pts - 1];

    } else if constexpr (RightBT == BoundaryType::FREE_EXIT) {

        new_stream_u_grid[cnt_x_pts - 1] = (left_invariant_and_g.first * right_boundary.const_g + right_boundary.const_invariant * 
            left_invariant_and_g.second) / (left_invariant_and_g.second + right_boundary.const_g);

        new_stream_h_grid[cnt_x_pts - 1] = (left_invariant_and_g.first - right_boundary.const_invariant) / (left_invariant_and_g.second + 
            right_boundary.const_g) - stream_z_grid[cnt_x_pts - 1];

    } else if constexpr(RightBT == BoundaryType::FLUVIAL_H_GIVEN) {

        new_stream_h_grid[cnt_x_pts - 1] = right_boundary.h;
        new_stream_u_grid[cnt_x_pts - 1] = left_invariant_and_g.first - left_invariant_and_g.second * (right_boundary.h + stream_z_grid[cnt_x_pts - 1]);

    } else if constexpr (RightBT == BoundaryType::FLUVIAL_FLUX_GIVEN) {

        new_stream_h_grid[cnt_x_pts - 1] = (-(left_invariant_and_g.second * stream_z_grid[cnt_x_pts - 1] - left_invariant_and_g.first) 
                                + std::sqrt(sqr(left_invariant_and_g.second * stream_z_grid[cnt_x_pts - 1] - left_invariant_and_g.first) 
                                - 4 * left_invariant_and_g.second * right_boundary.q)) 
                                / (2 * left_invariant_and_g.second);
        new_stream_u_grid[cnt_x_pts - 1] = -right_boundary.q / new_stream_h_grid[cnt_x_pts - 1];

    }
    
}




template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::second_phase(int t_idx) {

    //Вычисление границ
    compute_boundaries(t_idx);

    //Вычисления между краями
    for(int i = 1; i < (cnt_x_pts - 1); ++i) {

            //Инвариант левой ячейки
            std::pair<double, double> left_invariant_and_g = compute_invariant(LEFT_CELL, t_idx, i - 1);
            
            //Инвариант правой ячейки
            std::pair<double, double> right_invariant_and_g = compute_invariant(RIGHT_CELL, t_idx, i);

            //Вычисление потоковых величин
            new_stream_u_grid[i] = (left_invariant_and_g.first * right_invariant_and_g.second + right_invariant_and_g.first * 
                left_invariant_and_g.second) / (left_invariant_and_g.second + right_invariant_and_g.second);

            new_stream_h_grid[i] = (left_invariant_and_g.first - right_invariant_and_g.first) / (left_invariant_and_g.second + 
                right_invariant_and_g.second) - stream_z_grid[i];
            //

    }

    //Сохранение в текущие временные слои
    std::swap(new_stream_u_grid, stream_u_grid);
    std::swap(new_stream_h_grid, stream_h_grid);

}


template<BoundaryType LeftBT, BoundaryType RightBT>
std::pair<double, double> Cabaret_scheme<LeftBT, RightBT>::compute_invariant(Type_invariant type_invariant, int t_idx, int cell_idx) {


    //Вычисление инварианта
    double cur_c = std::sqrt(g * center_conservative_h_grid[cell_idx]);

    //Проверка на дозвуковое течение
    if (std::abs(center_conservative_u_grid[cell_idx]) > cur_c) {

        file_manager.end_file();
        throw std::runtime_error("Сверхзвуковое течение в compute invariant");

    }

    double G = g / cur_c;
    int num_koeff = (type_invariant == LEFT_CELL ? 1 : -1);

    double cell_left_invariant = stream_u_grid[cell_idx] + num_koeff * G * (stream_z_grid[cell_idx] + 
        stream_h_grid[cell_idx]);

    double cell_conservative_invariant = conservative_u_grid[cell_idx] + num_koeff * G * (conservative_z_grid[cell_idx] +
        conservative_h_grid[cell_idx]);

    double cell_center_conservative_invariant = center_conservative_u_grid[cell_idx] + num_koeff * G * (conservative_z_grid[cell_idx] + 
        center_conservative_h_grid[cell_idx]);

    double cell_right_invariant = stream_u_grid[cell_idx + 1] + num_koeff * G * (stream_z_grid[cell_idx + 1] + 
        stream_h_grid[cell_idx + 1]);

    double p = 0.0;
    double invariant = (2.0 * cell_center_conservative_invariant - (1 - p) * (type_invariant == LEFT_CELL ? cell_left_invariant : cell_right_invariant)) / (1 + p);
    
    double additional_term = num_koeff * tau_grid[t_idx] * center_conservative_u_grid[cell_idx] * G * (stream_z_grid[cell_idx + 1] -
        stream_z_grid[cell_idx]) / h_grid[cell_idx];

    double min_invariant = std::min({cell_left_invariant, cell_conservative_invariant, cell_right_invariant}) + additional_term;
    double max_invariant = std::max({cell_left_invariant, cell_conservative_invariant, cell_right_invariant}) + additional_term;

    if (invariant > max_invariant) {

        invariant = max_invariant;

    } else if (invariant < min_invariant) {

        invariant = min_invariant;

    }
    //
    return {invariant, G};

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::third_phase(int t_idx) {

    for(int i = 0; i < (cnt_x_pts - 1); ++i) {

        //Вычисление h на новом временном слое в центре отрезка
        conservative_h_grid[i] = center_conservative_h_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * stream_u_grid[i] -
            stream_h_grid[i + 1] * stream_u_grid[i + 1]) / h_grid[i];

        // Вычисление u на новом временном слое в центре отрезка
        conservative_u_grid[i] = (center_conservative_h_grid[i] * center_conservative_u_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * 
            sqr(stream_u_grid[i]) - stream_h_grid[i + 1] * sqr(stream_u_grid[i + 1]) + g * (stream_h_grid[i] + 
            stream_h_grid[i + 1]) / 2.0 * (stream_z_grid[i] + stream_h_grid[i] - stream_z_grid[i + 1] - stream_h_grid[i + 1])) / 
            h_grid[i]) / conservative_h_grid[i];



        // conservative_u_grid[i] = (center_conservative_h_grid[i] * center_conservative_u_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * 
        //     sqr(stream_u_grid[i]) - stream_h_grid[i + 1] * sqr(stream_u_grid[i + 1]) + g * conservative_h_grid[i] * (stream_z_grid[i] + stream_h_grid[i] - stream_z_grid[i + 1] - stream_h_grid[i + 1])) / 
        //     h_grid[i]) / conservative_h_grid[i];


        // conservative_u_grid[i] = (center_conservative_h_grid[i] * center_conservative_u_grid[i] + tau_grid[t_idx] / 2.0 * (stream_h_grid[i] * 
        //     sqr(stream_u_grid[i]) - stream_h_grid[i + 1] * sqr(stream_u_grid[i + 1]) + g * 0.5 * (sqr(stream_h_grid[i]) - sqr(stream_h_grid[i + 1]))) / 
        //     h_grid[i]) / conservative_h_grid[i];

    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::process_hydraulic_jump(double& cur_h, double& cur_u, double cur_automodel_h, double cur_automodel_u,
    double x_pt_or_hydraulic_jump_border, double x_pt_or_hydraulic_jump_border_1) {

    if (x_pt_or_hydraulic_jump_border < x_pt_or_hydraulic_jump_border_1) {

        cur_h = cur_automodel_h;
        cur_u = cur_automodel_u;

    } else {

        cur_h = automodel_h;
        cur_u = automodel_u;

    }
}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::process_discharge_wave(double& cur_h, double& cur_u, double cur_automodel_h_left, double cur_automodel_u_left,
    double cur_automodel_h_right, double cur_automodel_u_right, double x_pt, double discharge_wave_border_left,
    double discharge_wave_border_right, bool is_left) {
    
    if (x_pt < discharge_wave_border_left) {

        cur_h = cur_automodel_h_left;
        cur_u = cur_automodel_u_left;

    } else if (x_pt > discharge_wave_border_right) {

        cur_h = cur_automodel_h_right;
        cur_u = cur_automodel_u_right;

    } else {

        cur_u = cur_automodel_u_left + (cur_automodel_u_right - cur_automodel_u_left) *
            (x_pt - discharge_wave_border_left) / (discharge_wave_border_right - discharge_wave_border_left);

        cur_h = sqr((is_left ? automodel_const_invariant_left - cur_u : automodel_const_invariant_right - cur_u)) / (4.0 * g);
    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
double Cabaret_scheme<LeftBT, RightBT>::automodeling_solution(double start_cmp_seg, double end_cmp_seg) {
    
    //Начальная инициализация с_0
    automodel_u_1 = u_1(start_x_pt);
    automodel_u_2 = u_2(end_x_pt);
    automodel_h_1 = h_1(start_x_pt);
    automodel_h_2 = h_2(end_x_pt);
    automodel_c_1 = std::sqrt(g * automodel_h_1);
    automodel_c_2 = std::sqrt(g * automodel_h_2);
    automodel_c = 0.25 * (automodel_u_1 - automodel_u_2) + 0.5 * (automodel_c_1 + automodel_c_2);
    double k = (z(end_x_pt) - (z(start_x_pt))) / (end_x_pt - start_x_pt);

    //Итерации метода Ньютона
    double eps = 1e-10;
    while (true) {

        double automodel_c_prev = automodel_c;
        automodel_S_1 = automodel_c / automodel_c_1;
        automodel_S_2 = automodel_c / automodel_c_2;
        automodel_c = automodel_c_prev - (F(automodel_c_prev) - (automodel_u_1 - automodel_u_2)) / (derivative_F());

        if (std::abs(automodel_c - automodel_c_prev) < eps) {

            break;

        }

    }

    //Поиск h и u 
    automodel_h = sqr(automodel_c) / g;
    automodel_S_1 = automodel_c / automodel_c_1;
    automodel_S_2 = automodel_c / automodel_c_2;
    automodel_u = 0.5 * (automodel_u_1 + automodel_u_2 + phi_k(automodel_c, 2) - phi_k(automodel_c, 1));
    automodel_const_invariant_left = automodel_u + 2.0 * automodel_c;
    automodel_const_invariant_right = automodel_u - 2.0 * automodel_c;
    std::cout << "Автомодельные переменные: " << "\n" << "c = " << automodel_c << "\n" << "u = " << automodel_u << "\n" << "h = " << automodel_h << "\n";

    Type_behavior_automodel_solution type_left_behavior = (automodel_c > automodel_c_1 ? HYDRAULIC_JUMP : DISCHARGE_WAVE);
    Type_behavior_automodel_solution type_right_behavior = (automodel_c > automodel_c_2 ? HYDRAULIC_JUMP : DISCHARGE_WAVE);

    //Переменные для гидродинамических прыжков
    double D_1 = automodel_u_1 - automodel_c * std::sqrt(1 + sqr(automodel_S_1)) / std::sqrt(2);
    double D_2 = automodel_u_2 + automodel_c * std::sqrt(1 + sqr(automodel_S_2)) / std::sqrt(2);

    //Переменные для волн разряжения
    double D_1_left = automodel_u_1 - automodel_c_1;
    double special_c_1 = automodel_c_1 + (automodel_u_1 - automodel_u) / 2.0;
    double D_1_right = automodel_u - special_c_1;
    double D_2_right = automodel_u_2 + automodel_c_2;
    double special_c_2 = automodel_c_2 - (automodel_u_2 - automodel_u) / 2.0;
    double D_2_left = automodel_u + special_c_2;

    compute_start_streaming_and_conservative_values(automodel_stream_u_grid, automodel_stream_h_grid,
        automodel_conservative_u_grid, automodel_conservative_h_grid);


    if (!file_name.empty()) {
        file_manager.init_file(file_name + "_automodel.plt", 0, cnt_x_pts - 1);
        file_manager.save_layer(t_grid[0], conservative_x_grid, automodel_conservative_h_grid, automodel_conservative_u_grid
                              , conservative_z_grid);
    }

    int t_idx = 1;
    std::cout << type_left_behavior << " " << type_right_behavior << std::endl;

    //Разбор всевозможных случаев
    while(t_idx != cnt_t_pts) {
        //Слагаемое за счёт перехода в другую систему координат
        double additional_term_x = g * k * sqr(t_grid[t_idx] - start_t_pt) / 2.0;
        double additional_term_u = g * k * (t_grid[t_idx] - start_t_pt);

        //Границы для гидродинамических прыжков
        // double hydraulic_jump_left_border = (D_1) * t_grid[t_idx] + break_pt + additional_term_x;
        // double hydraulic_jump_right_border = (D_2) * t_grid[t_idx] + break_pt + additional_term_x;
        double hydraulic_jump_left_border = (D_1) * t_grid[t_idx] + break_pt;
        double hydraulic_jump_right_border = (D_2) * t_grid[t_idx] + break_pt;

        //Граница для волн разряжения
        // double discharge_wave_left_border_first = (D_1_left) * (t_grid[t_idx] - start_t_pt) + break_pt + additional_term_x;
        // double discharge_wave_left_border_second = (D_1_right) * (t_grid[t_idx] - start_t_pt) + break_pt + additional_term_x;
        // double discharge_wave_right_border_first = (D_2_left) * (t_grid[t_idx] - start_t_pt) + break_pt + additional_term_x;
        // double discharge_wave_right_border_second = (D_2_right) * (t_grid[t_idx] - start_t_pt) + break_pt + additional_term_x;
        double discharge_wave_left_border_first = (D_1_left) * (t_grid[t_idx] - start_t_pt) + break_pt;
        double discharge_wave_left_border_second = (D_1_right) * (t_grid[t_idx] - start_t_pt) + break_pt;
        double discharge_wave_right_border_first = (D_2_left) * (t_grid[t_idx] - start_t_pt) + break_pt;
        double discharge_wave_right_border_second = (D_2_right) * (t_grid[t_idx] - start_t_pt) + break_pt;


        for(int i = 0; i < cnt_x_pts; ++i) {
            if (i < (cnt_x_pts_before_break_pt + 1)) { //Левое расположение
                  
                if (type_left_behavior == HYDRAULIC_JUMP) { // Гидродинамический прыжок слева
                    //Потоковые величины
                    process_hydraulic_jump(
                        automodel_stream_h_grid[i], automodel_stream_u_grid[i], automodel_h_1, automodel_u_1,
                        x_grid[i] + additional_term_x, hydraulic_jump_left_border);

                    //Консервативные величины
                    process_hydraulic_jump(
                        automodel_conservative_h_grid[i], automodel_conservative_u_grid[i], automodel_h_1, automodel_u_1,
                        conservative_x_grid[i] + additional_term_x, hydraulic_jump_left_border);    
                    

                } else { // Волна разряжения слева

                    // //Потоковые величины
                    process_discharge_wave(
                        automodel_stream_h_grid[i], automodel_stream_u_grid[i], automodel_h_1, automodel_u_1, automodel_h, automodel_u,
                        x_grid[i] + additional_term_x, discharge_wave_left_border_first, discharge_wave_left_border_second, true);

                    //Консервативные величины    
                    process_discharge_wave(
                        automodel_conservative_h_grid[i], automodel_conservative_u_grid[i], automodel_h_1, automodel_u_1, automodel_h, automodel_u,
                        conservative_x_grid[i] + additional_term_x, discharge_wave_left_border_first, discharge_wave_left_border_second, true);

                }
                automodel_conservative_u_grid[i] -= additional_term_u;

            } else if (i == (cnt_x_pts_before_break_pt + 1)) { // Центральное расположение

                automodel_stream_h_grid[i] = automodel_h;
                automodel_stream_u_grid[i] = automodel_u;

            } else { // Правое расположения

                if (type_right_behavior == HYDRAULIC_JUMP) { // Гидродинамический прыжок справа

                    //Потоковые величины
                    process_hydraulic_jump(
                        automodel_stream_h_grid[i], automodel_stream_u_grid[i], automodel_h_2, automodel_u_2,
                        hydraulic_jump_right_border, x_grid[i] + additional_term_x);

                    //Консервативные величины
                    process_hydraulic_jump(
                        automodel_conservative_h_grid[i - 1], automodel_conservative_u_grid[i - 1], automodel_h_2, automodel_u_2,
                        hydraulic_jump_right_border ,conservative_x_grid[i - 1] + additional_term_x);    

                    

                } else { // Волна разряжения справа
                    //Потоковые величины
                    process_discharge_wave(
                        automodel_stream_h_grid[i], automodel_stream_u_grid[i], automodel_h, automodel_u, automodel_h_2, automodel_u_2,
                        x_grid[i] + additional_term_x, discharge_wave_right_border_first, discharge_wave_right_border_second, false);

                    //Консервативные величины    
                    process_discharge_wave(
                        automodel_conservative_h_grid[i - 1], automodel_conservative_u_grid[i - 1], automodel_h, automodel_u, automodel_h_2, automodel_u_2,
                        conservative_x_grid[i - 1] + additional_term_x, discharge_wave_right_border_first, discharge_wave_right_border_second, false);

                }
                automodel_conservative_u_grid[i - 1] -= additional_term_u;

            }

            automodel_stream_u_grid[i] -= additional_term_u;

        }

        if (tau_grid[t_idx - 1] != 0.0 and !file_name.empty()) {
            file_manager.save_layer(t_grid[t_idx], conservative_x_grid, automodel_conservative_h_grid, automodel_conservative_u_grid
                                  , conservative_z_grid);
        }

        ++t_idx;

    }

    if (!file_name.empty()) {
        file_manager.end_file();
    }


    double err = 0.0;
    Row err_vec = conservative_h_grid - automodel_conservative_h_grid;
    save_error(file_name + "_err.dat", cnt_x_pts - 1, err_vec, conservative_x_grid);
    int sz = err_vec.size();
    // std::cout << "sz=" << sz << std::endl;
    for(int i = 0; i < sz; ++i) {

        if ((x_grid[i] + h_grid[i] / 2.0) > start_cmp_seg and (x_grid[i] + h_grid[i] / 2.0) < end_cmp_seg) {

            if (std::abs(err_vec[i]) > err) {

                err = std::abs(err_vec[i]);

            }

        }
        // std::cout << "("<< err_vec[i] << ", " << x_grid[i] + h_grid[i] / 2.0 << ")," << std::endl; 

    }
    std::cout << cnt_x_pts << " " << err << std::endl;
    return err;

}


template<BoundaryType LeftBT, BoundaryType RightBT>
double Cabaret_scheme<LeftBT, RightBT>::phi_k(double c, int k) {

    double c_k = (k == 1 ? automodel_c_1 : automodel_c_2);
    double S_k = (k == 1 ? automodel_S_1 : automodel_S_2);

    if (S_k > 1) {

        return (c - c_k) * (S_k + 1.0) * std::sqrt(1.0 + 1.0 / sqr(S_k)) / std::sqrt(2.0);

    }

    return 2.0 * (c - c_k);

}


template<BoundaryType LeftBT, BoundaryType RightBT>
double Cabaret_scheme<LeftBT, RightBT>::F(double c) {

    return phi_k(c, 1) + phi_k(c, 2);

}


template<BoundaryType LeftBT, BoundaryType RightBT>
double Cabaret_scheme<LeftBT, RightBT>::derivative_phi_k(int k) {

    double S_k = (k == 1 ? automodel_S_1 : automodel_S_2);

    if (S_k > 1) {

        return (2.0 * sqr(S_k) + 1.0 + 1.0 / sqr(S_k)) / (std::sqrt(2.0) * S_k * std::sqrt(1.0 + 1.0 / sqr(S_k))); 

    }

    return 2.0;

}


template<BoundaryType LeftBT, BoundaryType RightBT>
double Cabaret_scheme<LeftBT, RightBT>::derivative_F() {

    return derivative_phi_k(1) + derivative_phi_k(2);

}


template<BoundaryType LeftBT, BoundaryType RightBT>
void Cabaret_scheme<LeftBT, RightBT>::save_error(std::string const& file_name, int cnt_pts_j, Row& err_vec, Row& x_grid) {
    std::ofstream file_out;
    file_out.open(file_name);
    file_out << "VARIABLES = \"x\", \"y\"" << std::endl;

    for(int i = 0; i < cnt_pts_j; ++i) {

        file_out << x_grid[i] << ' ' << err_vec[i] << std::endl;

    }

    file_out.close();

}


template<BoundaryType LeftBT, BoundaryType RightBT>
Row& Cabaret_scheme<LeftBT, RightBT>::stream_value(Type_value _type_stream_value) {

    switch (_type_stream_value) {

    case H:
        return stream_h_grid;

    case U:
        return stream_u_grid;

    }

}


template<BoundaryType LeftBT, BoundaryType RightBT>
Row& Cabaret_scheme<LeftBT, RightBT>::conservative_value(Type_value _type_conservative_value) {

    switch (_type_conservative_value) {

    case H:
        return conservative_h_grid;

    case U:
        return conservative_u_grid;

    }

}



//


#endif
