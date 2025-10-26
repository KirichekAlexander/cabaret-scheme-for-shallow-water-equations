#ifndef READ_MATRIX_H
#define READ_MATRIX_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <stdexcept>
#include <utility>
#include <iomanip>
#include <chrono>
#include <limits>


//Ключи для копирования матрицы
enum Type_copy_matrix {TRANSPOSITION, NON_TRANSPOSITION};


//Синонимы для строки и столбца
using Row = std::vector<double>;
using Column = std::vector<double>;


//Мой класс матрицы
class Matrix {
public:
    //Дефолт конструктор, создаёт пустую матрицу
    Matrix();
    //Конструктор от вектора
    Matrix(Row const&);
    //Копирование матрицы с или без транспонирования(по дефолту без транспонирования)
    Matrix(Matrix const&, Type_copy_matrix transposition = NON_TRANSPOSITION);
    //создаёт нулевую матрицу с заданными значениями
    Matrix(size_t, size_t, int type_matrix = 0);

    //возвращает количество столбцов
    size_t cols() const;
    //возвращает количество строк
    size_t rows() const;
    //прочитать столбец по индексу
    void read_column(Column&, size_t) const;
    //прочитать матрицу из файла или просто из входного потока по умолчанию
    void read_matrix(char const* file_path = nullptr);
    //записать матрицу в файл или просто в выходной  поток по умолчанию
    void print_matrix(char const* file_path = nullptr) const;
    //матричная норма(максимальная сумма элементов столбца, для вектор стобцов используется такая же норма)
    double matrix_norm() const;
    //другая матричная норма(максимальный элемент матрицы по модулю)
    double max_elem() const;
    //Добавление строки в матрицу
    void push_back_row(Row const&);

    //очистить матрицу
    void clear();
    //оператор умножения матрицы на число
    Matrix operator*(double) const;
    //оператор "-" для матриц
    Matrix operator-(Matrix const&) const;
    //оператор "+" для матриц
    Matrix operator+(Matrix const&) const;
    //оператор "=" для матриц
    Matrix& operator=(Matrix const&);
    //операторы взятия строки по индексу
    Row& operator[](size_t idx);
    Row const& operator[](size_t idx) const;
    //оператор "*" для матриц
    Matrix operator*(Matrix const&) const;
    friend std::istream& operator>>(std::istream&, Matrix&);
    friend std::ostream& operator<<(std::ostream&, Matrix const&);
private:
    std::vector<Row> matrix;
    size_t _rows;
    size_t _cols;
};


//Скалярное умножение вектор строк
double operator*(Row const&, Row const&);


//Оператор "-" для вектор строк
Row operator-(Row const&, Row const&);


//Оператор "+" для вектор строк
Row operator+(Row const&, Row const&);


//Оператор "+=" для вектор строк
Row& operator+=(Row&, Row const&);


//Оператор "-=" для вектор строк
Row& operator-=(Row&, Row const&);


//Оператор "*=" для вектор строк на число
Row& operator*=(Row&, double); 


//Оператор "/=" для вектор строк на число
Row& operator/=(Row&, double); 


//Оператор "*" для вектор строк на число
Row operator*(Row const&, double);


//Вывод матрицы
std::ostream& operator<<(std::ostream& ,Matrix const&);


//Чтение матрицы
std::istream& operator>>(std::istream&, Matrix&);


//Вывод векторов
std::ostream& operator<<(std::ostream&, Row const&);


//Вывод вектора целых чисел
std::ostream& operator<<(std::ostream& out, std::vector<int> const&);

//Функция возвращающая квадрат вещественного числа
double sqr(double);


//Функция интерполирущая матрицу к меньшей матрице, на входе матрица, которую интерполирует и матрица, куда идёт запись
void interpolate_row(Row const&, Row&);


//Функция для интерполяции вектора к вектору через полусумму
void interpolate_half_sum_row(Row const&, Row&);


//Функиця возвращаю максимальный элемент вектора
double max_elem(Row const&);




#endif