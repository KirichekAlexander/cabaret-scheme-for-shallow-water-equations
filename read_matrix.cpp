#include "read_matrix.h"


//Matrix реализация
Matrix::Matrix() 
    : _rows(0)
    , _cols(0)
{
}


Matrix::Matrix(Row const& row) 
    : matrix(1, Row(row.size(), 0.0))
    , _rows(1)
    , _cols(row.size())
{
    for(size_t i = 0; i < _cols; ++i) {
        matrix[0][i] = row[i];
    }
}


Matrix::Matrix(size_t n, size_t m, int type_matrix) 
    : matrix(n, Row(m, 0.0))
    , _rows(n)
    , _cols(m)
{
    if (type_matrix != 0) {
        for(size_t i = 0; i < _rows; ++i) {
            if (i < _cols) {
                matrix[i][i] = 1.0;
            }
        }
    }
}


Matrix::Matrix(Matrix const& _matrix, Type_copy_matrix transposition) 
    : _rows((transposition == TRANSPOSITION ? _matrix._cols : _matrix._rows))
    , _cols((transposition == TRANSPOSITION ? _matrix._rows : _matrix._cols))
{
    if (transposition == TRANSPOSITION) {
        matrix = std::vector<Row>(_rows, Row(_cols, 0.0));
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j) {
                matrix[i][j] = _matrix[j][i];
            }
        }
    } else {
        matrix = _matrix.matrix;
    }
}


size_t Matrix::rows() const {
    return _rows;
}


size_t Matrix::cols() const {
    return _cols;
}


void Matrix::read_column(Column& col, size_t idx) const {
    size_t col_size = col.size();
    for(size_t i = 0; i < col_size; ++i) {
        col[i] = matrix[i][idx];
    }
}


void Matrix::read_matrix(char const* file_path) {
    matrix.clear();
    if (file_path != nullptr) {
        std::ifstream file_input;
        file_input.open(file_path);
        file_input >> *this;
        file_input.close();
    } else {
        std::cin >> *this;
    }
}


void Matrix::print_matrix(char const* file_path) const {
    if (file_path != nullptr) {
        std::ofstream file_output;
        file_output.open(file_path);
        file_output << *this;
        file_output.close();
    } else {
        std::cout << *this;
    }
}


double Matrix::matrix_norm() const {
    double max_norm = 0.0;
    for (size_t i = 0; i < _cols; ++i) {
        double col_sum = 0.0;
        for (size_t j = 0; j < _rows; ++j) {
            col_sum += std::abs((*this)[j][i]);
        }
        max_norm = std::max(max_norm, col_sum);
    }
    return max_norm;
}


double Matrix::max_elem() const{
    double max_elem = 0.0;
    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < _cols; ++j) {
            if (std::abs(matrix[i][j]) > max_elem) {
                max_elem = std::abs(matrix[i][j]);
            }
        }
    }
    return max_elem;
}


void Matrix::push_back_row(Row const& row) {
    matrix.push_back(row);
    ++_rows;
}


void Matrix::clear() {
    matrix.clear();
    _rows = 0;
    _cols = 0;
}


Matrix Matrix::operator*(double num) const {
    Matrix res = Matrix(*this);
    for(size_t i = 0; i < _rows; ++i) {
        for(size_t j = 0; j < _cols; ++j) {
            res[i][j] *= num;
        }
    }
    return res;
}


Matrix Matrix::operator-(Matrix const& other) const {
    if (_rows != other._rows or _cols != other._cols) {
        throw std::invalid_argument("not equal shapes");
    }
    Matrix res = Matrix(*this);
    for(size_t i = 0; i < _rows; ++i) {
        res[i] -= other[i]; 
    }
    return res;
}


Matrix Matrix::operator+(Matrix const& other) const {
    if (_rows != other._rows or _cols != other._cols) {
        throw std::invalid_argument("not equal shapes");
    }
    Matrix res = Matrix(*this);
    for(size_t i = 0; i < _rows; ++i) {
        res[i] += other[i]; 
    }
    return res;
}


Matrix& Matrix::operator=(Matrix const& other) {
    if (this != &other) {
        matrix = other.matrix;
        _rows = other._rows;
        _cols = other._cols;
    }
    return *this;
}


Row& Matrix::operator[](size_t idx) {
    if (idx >= _rows) {
        throw std::invalid_argument("idx is not exist");
    }
    return matrix[idx];
}


Row const& Matrix::operator[](size_t idx) const {
    if (idx >= _rows) {
        throw std::invalid_argument("idx is not exist");
    }
    return matrix[idx];
}


Matrix Matrix::operator*(Matrix const& m2) const {
    if (_cols != m2._rows) {
        throw std::invalid_argument("cant mul matrix");
    }
    Matrix res = Matrix(_rows, m2._cols);
    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < m2._cols; ++j) {
            for (size_t k = 0; k < _cols; ++k) {
                res[i][j] += (*this)[i][k] * m2[k][j];
            }
        }
    }
    return res;
}


//


//Реализация других операторов и функций
std::ostream& operator<<(std::ostream& out, Matrix const& matrix) {
    for(size_t i = 0; i < matrix._rows; ++i) {
        for (size_t j = 0; j < matrix._cols; ++j) {
            out << std::setprecision(15) << matrix[i][j]; 
            out << (j != matrix._cols - 1 ? "," : "");
        }
        out << (i != matrix._rows - 1 ? "\n" : "");
    }
    return out;
}


std::ostream& operator<<(std::ostream& out, Row const& row) {
    size_t cnt = row.size();
    for(size_t i = 0; i < cnt; ++i) {
        out << std::setprecision(15) << row[i] << (i != cnt - 1 ? "," : "");
    }
    return out << "\n";
}


std::ostream& operator<<(std::ostream& out, std::vector<int> const& row) {
    size_t cnt = row.size();
    for(size_t i = 0; i < cnt; ++i) {
        out << row[i] << (i != cnt - 1 ? "," : "");
    }
    return out << "\n";
}


double operator*(Row const& vec1, Row const& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("not equal shapes vec1 vec2");
    }
    double res = 0;
    size_t vec_size = vec1.size();
    for(size_t i = 0; i < vec_size; ++i) {
        res += vec1[i]*vec2[i];
    }
    return res;
}


Row operator-(Row const& vec1, Row const& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("not equal shapes vec1 vec2");
    }
    Row res = vec1;
    size_t vec_size = vec1.size();
    for(size_t i = 0; i < vec_size; ++i) {
        res[i] -= vec2[i];
    }
    return res;
}


Row operator+(Row const& vec1, Row const& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("not equal shapes vec1 vec2");
    }
    Row res = vec1;
    size_t vec_size = vec1.size();
    for(size_t i = 0; i < vec_size; ++i) {
        res[i] += vec2[i];
    }
    return res;
}


Row& operator+=(Row& vec1, Row const& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("not equal shapes vec1 vec2");
    }
    size_t vec_size = vec1.size();
    for(size_t i = 0; i < vec_size; ++i) {
        vec1[i] += vec2[i];
    }
    return vec1;
}


Row& operator-=(Row& vec1, Row const& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("not equal shapes vec1 vec2");
    }
    size_t vec_size = vec1.size();
    for(size_t i = 0; i < vec_size; ++i) {
        vec1[i] -= vec2[i];
    }
    return vec1;
}


Row& operator*=(Row& vec, double num) {
    size_t vec_size = vec.size();
    for(size_t i = 0; i < vec_size; ++i) {
        vec[i] *= num;
    }
    return vec;
}


Row& operator/=(Row& vec, double num) {
    size_t vec_size = vec.size();
    for(size_t i = 0; i < vec_size; ++i) {
        vec[i] /= num;
    }
    return vec;
}


Row operator*(Row const& vec, double num) {
    Row res = vec;
    size_t vec_size = vec.size();
    for(size_t i = 0; i < vec_size; ++i) {
        res[i] *= num;
    }
    return res;
}


std::istream& operator >>(std::istream& in, Matrix& matrix) {
    std::string cur_str;
    while(std::getline(in, cur_str)) {
        if (!cur_str.empty()) {
            std::istringstream sin(cur_str);
            Row cur_row;
            while(std::getline(sin, cur_str, ',')) {
                cur_row.push_back(std::stod(cur_str));
            }
            matrix.matrix.push_back(cur_row);
        }
    }
    matrix._rows = matrix.matrix.size();
    matrix._cols = (matrix._rows == 0 ? 0 : matrix[0].size());
    return in;
}

//


//Реализация функции sqr

double sqr(double value) {

    return value * value;

}


//Реализация функции интерполяции 

void interpolate_row(Row const& row_sender, Row& row_reciever) {

    int sz = row_reciever.size();

    for (int i = 0; i < sz; ++i) {

        row_reciever[i] = row_sender[i * 2];

    }
}


//Реализация другой функции интерполяции

void interpolate_half_sum_row(Row const& row_sender, Row& row_reciever) {

    int sz = row_reciever.size();

    for(int i = 0; i < sz; ++i) {

        row_reciever[i] = (row_sender[2 * i] + row_sender[2 * i + 1]) / 2.0;

    }

} 


//Реализация функции максимального элемента вектора

double max_elem(Row const& row) { 

    double max_num = 0.0;
    int sz = row.size();

    for(int i = 0; i < sz; ++i) {

        if (std::abs(row[i]) > max_num) {
            
            max_num = std::abs(row[i]);

        }

    }

    return max_num;

}