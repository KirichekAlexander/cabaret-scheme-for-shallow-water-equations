#include "cabaret_scheme.h"


int
main() {
    // // Lakes tests
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 12.5;
    // double (* u_0)(double) = u_0_313;
    // double (* h_0)(double) = h_0_313;
    // double (* z)(double) = z_lakes;
    // // //3.1.1
    // // // Boundary<BoundaryType::DEFAULT> left_boundary(0.5, 0.0);
    // // // Boundary<BoundaryType::DEFAULT> right_boundary(left_boundary);
    // // //3.1.3
    // Boundary<BoundaryType::FLUVIAL_FLUX_GIVEN> left_boundary(3.5); // с 4.0 и 4.42 был сверхзвук
    // Boundary<BoundaryType::FLUVIAL_H_GIVEN> right_boundary(2.0);
    // // //3.1.5
    // // // Boundary<BoundaryType::FLUVIAL_FLUX_GIVEN> left_boundary(0.18); // с 4.0 и 4.42 был сверхзвук
    // // // Boundary<BoundaryType::FLUVIAL_H_GIVEN> right_boundary(0.33);
    // double start_x_pt = 0.0;
    // double end_x_pt = 25.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 1000.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 200;
    // int cnt_x_pts_after_break_pt = 200;
    // bool continuity = true;

    // std::string file_name = "./data/3.1.3.lake";
    // Cabaret_scheme<BoundaryType::FLUVIAL_FLUX_GIVEN, BoundaryType::FLUVIAL_H_GIVEN> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0, u_0, h_0, h_0, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt, file_name, continuity, nullptr);
    // cabaret_scheme_1.compute();
    // std::vector<double> h_cons = cabaret_scheme_1.conservative_value(H);
    // double c = 3.5 * 3.5 / (8 * g) + 2.0;
    // for(int i = 0; i < h_cons.size(); ++i) {
    //     std::cout << "eq=" << sqr(h_cons[i]) * h_cons[i] + (z(i * (25.0 / 402.0) + (25.0 / 804.0)) - c) * sqr(h_cons[i]) + 3.5 * 3.5 / (2 * g) << std::endl;
    // }
    


    // Dam break tests
    Type_grid type_grid = EVEN;
    double CFL = 0.3;
    double break_pt = 5.0;
    double (* u_0_left)(double) = u_0_411_left_right;
    double (* u_0_right)(double) = u_0_411_left_right;
    double (* h_0_left)(double) = h_0_411_left;
    double (* h_0_right)(double) = h_0_411_right;
    double (* z)(double) = z_uneven;
    //4.1.1.
    Boundary<BoundaryType::FREE_EXIT> left_boundary;
    Boundary<BoundaryType::FREE_EXIT> right_boundary;

    double start_x_pt = 0.0;
    double end_x_pt = 10.0;
    double start_t_pt = 0.0;
    double end_t_pt = 0.25;
    // double start_cmp_seg = 6.0;
    // double end_cmp_seg = 14.0;
    int cnt_x_pts_before_break_pt = 50;
    int cnt_x_pts_after_break_pt = 50;
    bool continuity = false;

    std::string file_name = "./data/rim_unev";
    Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_1(type_grid
        , left_boundary
        , right_boundary
        , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
        , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt, file_name, continuity, analytical_solution);
    double err1 = cabaret_scheme_1.compute();
    cabaret_scheme_1.automodeling_solution(0.0, 10.0);

    // Cabaret_scheme<BoundaryType::DEFAULT, BoundaryType::DEFAULT> cabaret_scheme_2(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 2 + 1, cnt_x_pts_after_break_pt * 2 + 1, file_name, continuity, analytical_solution);
    // double err2 = cabaret_scheme_2.compute();
    // std::cout << "Ошибки на аналитическом решении: обычная сетка: " << err1 << " , на дробленной: " << err2 << std::endl 
    //           << "Порядок аппркосимации: " << std::log2(err1 / err2) << std::endl;
    


    // //Тесты вкусно и точка
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 500.0;
    // double (* u_0_left)(double) = u_0_321;
    // double (* u_0_right)(double) = u_0_321;
    // double (* h_0_left)(double) = h_0_321;
    // double (* h_0_right)(double) = h_0_321;
    // double (* z)(double) = z_321;
    // //3.2.1
    // Boundary<BoundaryType::FLUVIAL_FLUX_GIVEN> left_boundary(1.5);
    // Boundary<BoundaryType::FLUVIAL_H_GIVEN> right_boundary(h_ex(1000));

    // double start_x_pt = 0.0;
    // double end_x_pt = 1000.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 30000.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 1000;
    // int cnt_x_pts_after_break_pt = 1000;
    // bool continuity = true;

    // std::string file_name = "./data/3.2.1_hex";
    // Cabaret_scheme<BoundaryType::FLUVIAL_FLUX_GIVEN, BoundaryType::FLUVIAL_H_GIVEN> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt, file_name, continuity, nullptr);
    // cabaret_scheme_1.compute();
    //


    // //Проверка порядка аппроксимации(распад гаусса)
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 10.0;
    // double (* u_0_left)(double) = u_1;
    // double (* u_0_right)(double) = u_2;
    // double (* h_0_left)(double) = h_1;
    // double (* h_0_right)(double) = h_2;
    // double (* z)(double) = z_const;
    // Boundary<BoundaryType::NON_LEAK> left_boundary;
    // Boundary<BoundaryType::NON_LEAK> right_boundary;

    // double start_x_pt = 0.0;
    // double end_x_pt = 20.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 1.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 512;
    // int cnt_x_pts_after_break_pt = 512;
    // bool continuity = true;

    // std::string file_name = "./data/approx_order";
    // Cabaret_scheme<BoundaryType::NON_LEAK, BoundaryType::NON_LEAK> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_1.compute();
    // Row value_cabaret_scheme = cabaret_scheme_1.conservative_value(H);



    // Cabaret_scheme<BoundaryType::NON_LEAK, BoundaryType::NON_LEAK> cabaret_scheme_2(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 2 + 1, cnt_x_pts_after_break_pt * 2 + 1
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_2.compute();
    // Row value_cabaret_scheme_2 = cabaret_scheme_2.conservative_value(H);
    // Row interp_value_cabaret_scheme_2 = Row(value_cabaret_scheme.size());
    // interpolate_half_sum_row(value_cabaret_scheme_2, interp_value_cabaret_scheme_2);



    
    // Cabaret_scheme<BoundaryType::NON_LEAK, BoundaryType::NON_LEAK> cabaret_scheme_4(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 4 + 3, cnt_x_pts_after_break_pt * 4 + 3
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_4.compute();
    // Row value_cabaret_scheme_4 = cabaret_scheme_4.conservative_value(H);
    // Row interp_value_cabaret_scheme_4 = Row(value_cabaret_scheme_2.size());
    // interpolate_half_sum_row(value_cabaret_scheme_4, interp_value_cabaret_scheme_4);
    // double h1 = (end_x_pt - start_x_pt) / (cnt_x_pts_before_break_pt + cnt_x_pts_after_break_pt + 2);
    // double h2 = h1 * 0.5;
    // Row err1 = value_cabaret_scheme - interp_value_cabaret_scheme_2;
    // Row err2 = value_cabaret_scheme_2 - interp_value_cabaret_scheme_4;
    // save_vec(file_name + "_conservative/err_h.dat", "err", err1, cabaret_scheme_1.get_conservative_x_grid());
    // save_vec(file_name + "_conservative/err_h0.5.dat", "err", err2, cabaret_scheme_2.get_conservative_x_grid());
    // double norm_1 = l2_norm(err1, h1);
    // double norm_2 = l2_norm(err2, h2);
    // double approx_order = std::log2(norm_1 / norm_2);
    // std::cout << "Порядок аппроксимации схемы: " << approx_order << std::endl << norm_1 << ' ' << norm_2 << std::endl;
    


    //Проверка порядка аппроксимации(две волны разряжения)
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 0.0;
    // double (* u_0_left)(double) = u_1;
    // double (* u_0_right)(double) = u_2;
    // double (* h_0_left)(double) = h_1;
    // double (* h_0_right)(double) = h_2;
    // double (* z)(double) = z_const;
    // Boundary<BoundaryType::FREE_EXIT> left_boundary;
    // Boundary<BoundaryType::FREE_EXIT> right_boundary;

    // double start_x_pt = -10.0;
    // double end_x_pt = 10.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 1.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 100;
    // int cnt_x_pts_after_break_pt = 100;
    // bool continuity = false;

    // std::string file_name = "./data/approx_order";
    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_1.compute();
    // double automodel_err_1 = cabaret_scheme_1.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme = cabaret_scheme_1.conservative_value(H);



    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_2(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 2 + 1, cnt_x_pts_after_break_pt * 2 + 1
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_2.compute();
    // double automodel_err_2 = cabaret_scheme_2.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme_2 = cabaret_scheme_2.conservative_value(H);
    // Row interp_value_cabaret_scheme_2 = Row(value_cabaret_scheme.size());
    // interpolate_half_sum_row(value_cabaret_scheme_2, interp_value_cabaret_scheme_2);



    
    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_4(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 4 + 3, cnt_x_pts_after_break_pt * 4 + 3
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_4.compute();
    // double automodel_err_4 = cabaret_scheme_4.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme_4 = cabaret_scheme_4.conservative_value(H);
    // Row interp_value_cabaret_scheme_4 = Row(value_cabaret_scheme_2.size());
    // interpolate_half_sum_row(value_cabaret_scheme_4, interp_value_cabaret_scheme_4);
    // double h1 = (end_x_pt - start_x_pt) / (cnt_x_pts_before_break_pt + cnt_x_pts_after_break_pt + 2);
    // double h2 = h1 * 0.5;
    // Row err1 = value_cabaret_scheme - interp_value_cabaret_scheme_2;
    // Row err2 = value_cabaret_scheme_2 - interp_value_cabaret_scheme_4;
    // save_vec(file_name + "_conservative/err_h.dat", "err", err1, cabaret_scheme_1.get_conservative_x_grid());
    // save_vec(file_name + "_conservative/err_h0.5.dat", "err", err2, cabaret_scheme_2.get_conservative_x_grid());
    // double norm_1 = l2_norm(err1, h1);
    // double norm_2 = l2_norm(err2, h2);
    // double approx_order = std::log2(norm_1 / norm_2);
    // std::cout << "Порядок аппроксимации схемы: " << approx_order << std::endl << norm_1 << ' ' << norm_2 << std::endl;
    // std::cout << "Порядок аппроксимации на обычной и дробленной сетке: " << std::log2(automodel_err_1 / automodel_err_2) 
    //           << std::endl << automodel_err_1 << ' ' << automodel_err_2 << std::endl; 
    // std::cout << "Порядок аппроксимации на дробленной и сильно дробленной сетке: " << std::log2(automodel_err_2 / automodel_err_4) 
    //           << std::endl << automodel_err_2 << ' ' << automodel_err_4 << std::endl; 
    //

    // //Тест одна волна разряжения(one discharge wave)
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 0.0;
    // double (* z)(double) = z_const;
    // Boundary<BoundaryType::FREE_EXIT> left_boundary;
    // Boundary<BoundaryType::FREE_EXIT> right_boundary;

    // double start_x_pt = -10.0;
    // double end_x_pt = 10.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 1.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 100;
    // int cnt_x_pts_after_break_pt = 100;
    // bool continuity = false; // флаг разрыва начальных данных в узле разрыва храним полусумму величин в начальный момент времени

    // std::string file_name = "./data/one_discharge_wave_1";
    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, one_wave_u, one_wave_u, one_wave_h, one_wave_h, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_1.compute();
    // double automodel_err_1 = cabaret_scheme_1.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme = cabaret_scheme_1.conservative_value(H);

    // file_name =  "./data/one_discharge_wave_2";
    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_2(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, one_wave_u, one_wave_u, one_wave_h, one_wave_h, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 2 + 1, cnt_x_pts_after_break_pt * 2 + 1
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_2.compute();
    // double automodel_err_2 = cabaret_scheme_2.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme_2 = cabaret_scheme_2.conservative_value(H);
    // Row interp_value_cabaret_scheme_2 = Row(value_cabaret_scheme.size());
    // interpolate_half_sum_row(value_cabaret_scheme_2, interp_value_cabaret_scheme_2);

    // file_name =  "./data/one_discharge_wave_4";
    // Cabaret_scheme<BoundaryType::FREE_EXIT, BoundaryType::FREE_EXIT> cabaret_scheme_4(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, one_wave_u, one_wave_u, one_wave_h, one_wave_h, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt * 4 + 3, cnt_x_pts_after_break_pt * 4 + 3
    //     , file_name, continuity, nullptr);
    // cabaret_scheme_4.compute();
    // double automodel_err_4 = cabaret_scheme_4.automodeling_solution(-10.0, 10.0);
    // Row value_cabaret_scheme_4 = cabaret_scheme_4.conservative_value(H);
    // Row interp_value_cabaret_scheme_4 = Row(value_cabaret_scheme_2.size());
    // interpolate_half_sum_row(value_cabaret_scheme_4, interp_value_cabaret_scheme_4);
    // double h1 = (end_x_pt - start_x_pt) / (cnt_x_pts_before_break_pt + cnt_x_pts_after_break_pt + 2);
    // double h2 = h1 * 0.5;
    // Row err1 = value_cabaret_scheme - interp_value_cabaret_scheme_2;
    // Row err2 = value_cabaret_scheme_2 - interp_value_cabaret_scheme_4;
    // save_vec("./data/one_wave_conservative_err_h.dat", "err", err1, cabaret_scheme_1.get_conservative_x_grid());
    // save_vec("./data/one_wave_conservative_err_h0.5.dat", "err", err2, cabaret_scheme_2.get_conservative_x_grid());
    // double norm_1 = l2_norm(err1, h1);
    // double norm_2 = l2_norm(err2, h2);
    // double approx_order = std::log2(norm_1 / norm_2);
    // std::cout << "Порядок аппроксимации схемы: " << approx_order << std::endl;
    // std::cout << "Порядок аппроксимации на обычной и дробленной сетке сравнение с точным решением: " << std::log2(automodel_err_1 / automodel_err_2) << std::endl; 
    // std::cout << "Порядок аппроксимации на дробленной и сильно дробленной сетке сравнение с точным решением: " << std::log2(automodel_err_2 / automodel_err_4) << std::endl; 


    // // Тест разрыв над наклонным дном
    // Type_grid type_grid = EVEN;
    // double CFL = 0.3;
    // double break_pt = 5.0;
    // double (* u_0_left)(double) = u_0_411_left_right;
    // double (* u_0_right)(double) = u_0_411_left_right;
    // double (* h_0_left)(double) = h_0_411_left;
    // double (* h_0_right)(double) = h_0_411_right;
    // double (* z)(double) = z_const;
    // //4.1.1.
    // Boundary<BoundaryType::DEFAULT> left_boundary(5.0, 0.0);
    // Boundary<BoundaryType::DEFAULT> right_boundary(1.0, 0.0);

    // double start_x_pt = 0.0;
    // double end_x_pt = 10.0;
    // double start_t_pt = 0.0;
    // double end_t_pt = 6.0;
    // // double start_cmp_seg = 6.0;
    // // double end_cmp_seg = 14.0;
    // int cnt_x_pts_before_break_pt = 50;
    // int cnt_x_pts_after_break_pt = 50;
    // bool continuity = false;

    // std::string file_name = "./data/4.1.1__2z_uneven";
    // Cabaret_scheme<BoundaryType::DEFAULT, BoundaryType::DEFAULT> cabaret_scheme_1(type_grid
    //     , left_boundary
    //     , right_boundary
    //     , CFL, break_pt, u_0_left, u_0_right, h_0_left, h_0_right, z
    //     , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt, file_name, continuity, nullptr);
    // double err1 = cabaret_scheme_1.compute();
    // cabaret_scheme_1.automodeling_solution(0.0, 10.0);

    // std::string file_name = "./data/hydrajumps";
    // Cabaret_scheme cabaret_scheme_1(type_grid, NON_LEAK, CFL, break_pt
    //         , u_1, u_2, h_1, h_2, linear_function
    //         , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
    //         , file_name);
    // cabaret_scheme_1.compute();
    // cabaret_scheme_1.automodeling_solution(start_cmp_seg, end_cmp_seg);

//     file_name = "";
//     cnt_x_pts_before_break_pt = 200;
//     cnt_x_pts_after_break_pt = 200;
//     Cabaret_scheme cabaret_scheme_2(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_2.compute();

//     cnt_x_pts_before_break_pt = 400;
//     cnt_x_pts_after_break_pt = 400;
//     Cabaret_scheme cabaret_scheme_3(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_3.compute();
    
//     cnt_x_pts_before_break_pt = 600;
//     cnt_x_pts_after_break_pt = 600;
//     Cabaret_scheme cabaret_scheme_4(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_4.compute();
    
//     cnt_x_pts_before_break_pt = 800;
//     cnt_x_pts_after_break_pt = 800;
//     Cabaret_scheme cabaret_scheme_5(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_5.compute();
    
//     cnt_x_pts_before_break_pt = 1000;
//     cnt_x_pts_after_break_pt = 1000;
//     Cabaret_scheme cabaret_scheme_6(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_6.compute();

//     cnt_x_pts_before_break_pt = 1500;
//     cnt_x_pts_after_break_pt = 1500;
//     Cabaret_scheme cabaret_scheme_7(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_7.compute();

//     cnt_x_pts_before_break_pt = 2000;
//     cnt_x_pts_after_break_pt = 2000;
//     Cabaret_scheme cabaret_scheme_8(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_8.compute();

//     cnt_x_pts_before_break_pt = 2500;
//     cnt_x_pts_after_break_pt = 2500;
//     Cabaret_scheme cabaret_scheme_9(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_9.compute();

//     cnt_x_pts_before_break_pt = 3000;
//     cnt_x_pts_after_break_pt = 3000;
//     Cabaret_scheme cabaret_scheme_10(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_2, h_1, h_2, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//             , file_name);
//     cabaret_scheme_10.compute();

//     cabaret_scheme_2.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_3.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_4.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_5.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_6.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_7.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_8.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_9.automodeling_solution(start_cmp_seg, end_cmp_seg);
//     cabaret_scheme_10.automodeling_solution(start_cmp_seg, end_cmp_seg);

//         std::string file_name = "./data/output1_free_exit";
//     Cabaret_scheme cabaret_scheme_1(type_grid, NON_LEAK, CFL, break_pt
//         , u_1, u_1, h_1, h_1, linear_function
//         , start_x_pt, end_x_pt, start_t_pt, end_t_pt, cnt_x_pts_before_break_pt, cnt_x_pts_after_break_pt
//         , file_name);
// //     cabaret_scheme.compute();
//     cabaret_scheme_1.compute();
// //     Запуск вычислений
// //     cabaret_scheme.compute();
//     Row value_cabaret_scheme = cabaret_scheme_1.conservative_value(H);



//     file_name = "./data/output2.dat";
//     Cabaret_scheme cabaret_scheme_2(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_1, h_1, h_1, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, 2 * cnt_x_pts_before_break_pt + 1, 2 * cnt_x_pts_after_break_pt + 1
//             , file_name);

//     cabaret_scheme_2.compute();
//     Row value_cabaret_scheme_2 = cabaret_scheme_2.conservative_value(H);
//     Row interp_value_cabaret_scheme_2 = Row(value_cabaret_scheme.size());
//     interpolate_half_sum_row(value_cabaret_scheme_2, interp_value_cabaret_scheme_2);



    
//     file_name = "./data/output3.dat";
//     Cabaret_scheme cabaret_scheme_4(type_grid, NON_LEAK, CFL, break_pt
//             , u_1, u_1, h_1, h_1, linear_function
//             , start_x_pt, end_x_pt, start_t_pt, end_t_pt, 4 * cnt_x_pts_before_break_pt + 3, 4 * cnt_x_pts_after_break_pt + 3
//             , file_name);

//     cabaret_scheme_4.compute();
//     Row value_cabaret_scheme_4 = cabaret_scheme_4.conservative_value(H);
//     Row interp_value_cabaret_scheme_4 = Row(value_cabaret_scheme_2.size());
//     interpolate_half_sum_row(value_cabaret_scheme_4, interp_value_cabaret_scheme_4);

//     double norm_1 = max_elem(value_cabaret_scheme - interp_value_cabaret_scheme_2);
//     double norm_2 = max_elem(value_cabaret_scheme_2 - interp_value_cabaret_scheme_4);
//     double approx_order = std::log2(norm_1 / norm_2);
//     std::cout << "Порядок аппроксимации схемы: " << approx_order << std::endl << norm_1 << ' ' << norm_2 << std::endl;
    // cabaret_scheme.compute();
    // cabaret_scheme_2.compute();

    // double start_cmp_seg = 6.0;
    // double end_cmp_seg = 14.0;
    // double err_1 = cabaret_scheme.automodeling_solution(start_cmp_seg, end_cmp_seg);
    // double err_2 = cabaret_scheme_2.automodeling_solution(start_cmp_seg, end_cmp_seg);
    // double approx_order_with_automodel_solution = std::log2(err_1 / err_2);
    // std::cout << "Порядок аппроксимации схемы с использованием аналитического решения: " << approx_order_with_automodel_solution << std::endl;

}