#include "file_manager.h"


//Реализация файлового менеджера
FileManager::FileManager() 
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
{
}


//Инициализация файла
void FileManager::init_file(std::string fn, INTEGER4 d, int im) {

    file_name = fn;
    debug = d;
    
    imax = im;
    num_points = 2 * imax;
    x_values = y_values = h_values = u_values = z_values = std::vector<double>(num_points, 0.0);

    TECINI112(const_cast<char*>(title.c_str())
            , const_cast<char*>(variables.c_str())
            , const_cast<char*>(file_name.c_str())
            , const_cast<char*>(scratch_dir.c_str())
            , &file_type
            , &debug
            , &v_is_double
    );

}


//Сохранение временного слоя
void FileManager::save_layer(double t, std::vector<double>& x, std::vector<double>& h, std::vector<double>& u, std::vector<double>& z) {

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

}


//Закрытие файла
void FileManager::end_file() {
    TECEND112();
}