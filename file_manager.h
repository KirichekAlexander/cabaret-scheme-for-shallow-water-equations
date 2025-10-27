#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H

#include <iostream>
#include <string>
#include "TECIO.h"
#include "TECXXX.h"


class FileManager {

public:

    FileManager();
    void init_file(std::string fn, INTEGER4 d, int imax);
    void save_layer(double t, std::vector<double>& x, std::vector<double>& h, std::vector<double>& u, std::vector<double>& z);
    void end_file();



private:

    //Переменные для TECINI112
    std::string title;
    std::string variables;
    std::string file_name;
    std::string scratch_dir;
    INTEGER4 file_type;
    INTEGER4 debug;
    INTEGER4 v_is_double;

    //Переменные для TECZNE112
    std::string zone_title;
    INTEGER4 zone_type;
    INTEGER4 imax;             
    INTEGER4 jmax;     
    INTEGER4 kmax;

    INTEGER4 icellmax;
    INTEGER4 jcellmax; 
    INTEGER4 kcellmax;

    double solution_time;
    INTEGER4 strand_id;
    INTEGER4 parent_zone;
    INTEGER4 is_block;        
    INTEGER4 num_face_connections;
    INTEGER4 face_neighbor_mode;
    INTEGER4 total_num_face_nodes;
    INTEGER4 num_connected_boundary_faces;
    INTEGER4 total_num_boundary_connections;
    INTEGER4* passive_var_list;
    std::vector<INTEGER4> value_location; 
    std::vector<INTEGER4> share_var_from_zone;
    INTEGER4 share_connectivity_from_zone;

    //Переменные для TECDAT112
    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> h_values;
    std::vector<double> u_values;
    std::vector<double> z_values;
    INTEGER4 num_points;

};



#endif