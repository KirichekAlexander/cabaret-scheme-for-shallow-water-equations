#include "file_manager.h"


//Реализация файлового менеджера
FileManager::FileManager(int cnt_pts) 
    //Переменные для TECINI112
    : title("Shallow Water Solution")
    , variables("x y h u")
    , scratch_dir(".")
    , file_type(0)
    , v_is_double(1)

    //Переменные для TECZNE112
    , zone_title("t=")
    , zone_type(0)
    , jmax(2)
    , kmax(1)
    , strand_id(1)
    , parent_zone(0)
    , is_block(1)
    , num_face_connections(0)
    , face_neighbor_mode(0)
    , total_num_face_nodes(0)
    , num_connected_boundary_faces(0)
    , total_num_boundary_connections(0)
    , passive_var_list(nullptr)
    , value_location({1, 1, 1, 1})
    , share_var_from_zone({0, 0, 0, 0})
    , share_connectivity_from_zone(0)

    //Переменные для TECDAT112
    , x_values(cnt_pts, 0.0)
    , y_values(cnt_pts, 0.0)
    , h_values(cnt_pts, 0.0)
    , u_values(cnt_pts, 0.0)
{
}


//Сохранение временного слоя
void FileManager::save_layer(std::string fn, INTEGER4 d, int im, int num_layer, double t, std::vector<double>& x
                           , std::vector<double>& h, std::vector<double>& u, std::vector<double>& z) {
    
    //очищение файлов
    if (num_layer == 0) {

        std::filesystem::remove_all(fn);
        std::filesystem::create_directory(fn);

    }

    imax = im;
    num_points = 2 * imax;
    debug = d;
    TECINI112(const_cast<char*>(title.c_str())
            , const_cast<char*>(variables.c_str())
            , const_cast<char*>((fn + "/" + std::to_string(num_layer) + ".plt").c_str())
            , const_cast<char*>(scratch_dir.c_str())
            , &file_type
            , &debug
            , &v_is_double
    );

    solution_time = t;
    icellmax = jcellmax = kcellmax = 0;
    TECZNE112(const_cast<char*>((zone_title + std::to_string(solution_time)).c_str())
            , &zone_type
            , &imax
            , &jmax
            , &kmax 
            , &icellmax
            , &jcellmax
            , &kcellmax
            , &solution_time
            , &strand_id
            , &parent_zone
            , &is_block
            , &num_face_connections
            , &face_neighbor_mode
            , &total_num_face_nodes
            , &num_connected_boundary_faces
            , &total_num_boundary_connections
            , passive_var_list
            , value_location.data()
            , share_var_from_zone.data()
            , &share_connectivity_from_zone
    );

    for(int i = 0; i < imax; ++i) {

        x_values[i] = x_values[i + imax] = x[i];
        y_values[i] = z[i];
        y_values[i + imax] = z[i] + h[i];
        h_values[i] = h_values[i + imax] = h[i];
        u_values[i] = u_values[i + imax] = u[i];
        
    }

    TECDAT112(&num_points, x_values.data(), &v_is_double);
    TECDAT112(&num_points, y_values.data(), &v_is_double);
    TECDAT112(&num_points, h_values.data(), &v_is_double);
    TECDAT112(&num_points, u_values.data(), &v_is_double);

    TECEND112();
}



FileManager2D::FileManager2D(int nx_, int ny_)
    : nx(nx_)
    , ny(ny_)

    // TECINI112
    , title("Shallow Water Solution 2D")
    , variables("x y z H h u v")
    , scratch_dir(".")
    , file_type(0)
    , v_is_double(1)

    // TECZNE112
    , zone_title_prefix("t=")
    , zone_type(0)              // Ordered
    , imax(nx)
    , jmax(ny)
    , kmax(2)                   // <-- ВАЖНО: два слоя (дно + поверхность)
    , icellmax(0)
    , jcellmax(0)
    , kcellmax(0)
    , solution_time(0.0)
    , strand_id(1)
    , parent_zone(0)
    , is_block(1)               // Block data
    , num_face_connections(0)
    , face_neighbor_mode(0)
    , total_num_face_nodes(0)
    , num_connected_boundary_faces(0)
    , total_num_boundary_connections(0)
    , passive_var_list(nullptr)
    , value_location({1,1,1,1,1,1,1})
    , share_var_from_zone({0,0,0,0,0,0,0})
    , share_connectivity_from_zone(0)

    // TECDAT112
    , num_points(static_cast<INTEGER4>(nx * ny * 2))  // <-- ВАЖНО
    , x_values(num_points, 0.0)
    , y_values(num_points, 0.0)
    , z_values(num_points, 0.0)
    , H_values(num_points, 0.0)
    , h_values(num_points, 0.0)
    , u_values(num_points, 0.0)
    , v_values(num_points, 0.0)
{
}

void FileManager2D::save_layer(const std::string& folder,
                               INTEGER4 debug,
                               int layer,
                               double t,
                               Row const& x_center,
                               Row const& y_center,
                               Matrix const& z_center,
                               Field2D const& center)
{
    if (layer == 0) {
        std::filesystem::remove_all(folder);
        std::filesystem::create_directory(folder);
    }

    std::string file_name = folder + "/" + std::to_string(layer) + ".plt";

    TECINI112(const_cast<char*>(title.c_str()),
              const_cast<char*>(variables.c_str()),
              const_cast<char*>(file_name.c_str()),
              const_cast<char*>(scratch_dir.c_str()),
              &file_type,
              &debug,
              &v_is_double);

    solution_time = t;
    imax = nx;
    jmax = ny;
    kmax = 2;
    icellmax = jcellmax = kcellmax = 0;

    std::string zone_title = zone_title_prefix + std::to_string(solution_time);

    TECZNE112(const_cast<char*>(zone_title.c_str()),
              &zone_type,
              &imax,
              &jmax,
              &kmax,
              &icellmax,
              &jcellmax,
              &kcellmax,
              &solution_time,
              &strand_id,
              &parent_zone,
              &is_block,
              &num_face_connections,
              &face_neighbor_mode,
              &total_num_face_nodes,
              &num_connected_boundary_faces,
              &total_num_boundary_connections,
              passive_var_list,
              value_location.data(),
              share_var_from_zone.data(),
              &share_connectivity_from_zone);

    // k=0 -> дно, k=1 -> свободная поверхность
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int idx = i + nx * (j + ny * k);

                const double x = x_center[i];
                const double y = y_center[j];

                const double zb = z_center[i][j];
                const double h  = center[i][j].h;
                const double u  = center[i][j].u;
                const double v  = center[i][j].v;

                const double zw = zb + h;

                x_values[idx] = x;
                y_values[idx] = y;

                if (k == 0) {
                    // слой дна
                    z_values[idx] = zb;
                    H_values[idx] = zb;     // можно так
                    h_values[idx] = 0.0;    // глубина на дне не рисуется
                    u_values[idx] = u;      // можно оставить u,v как у поверхности (удобно для раскраски)
                    v_values[idx] = v;
                } else {
                    // слой свободной поверхности
                    z_values[idx] = zw;
                    H_values[idx] = zw;     // поверхность
                    h_values[idx] = h;
                    u_values[idx] = u;
                    v_values[idx] = v;
                }
            }
        }
    }

    TECDAT112(&num_points, x_values.data(), &v_is_double);
    TECDAT112(&num_points, y_values.data(), &v_is_double);
    TECDAT112(&num_points, z_values.data(), &v_is_double);
    TECDAT112(&num_points, H_values.data(), &v_is_double);
    TECDAT112(&num_points, h_values.data(), &v_is_double);
    TECDAT112(&num_points, u_values.data(), &v_is_double);
    TECDAT112(&num_points, v_values.data(), &v_is_double);

    TECEND112();
}
