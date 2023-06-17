#include <iostream>
#include <vector>
#include <string>
#include <fstream>

template<class T>
class Matrix {
private:
    size_t N, M;
    std::vector<std::vector<T>> matrix;
public:
    explicit Matrix(size_t n, size_t m, std::vector<std::vector<T>> new_matrix) : N(n), M(m), matrix(new_matrix) {};

    Matrix(size_t n, size_t m) : N(n), M(m) {
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    Matrix() {
        std::cin >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    Matrix(const std::string &filename) {
        std::ifstream file_in(filename);
        file_in >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                file_in >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    size_t get_lines_number() {
        return N;
    }

    size_t get_colomns_number() {
        return M;
    }

    std::vector<T> &operator[](size_t i) {
        return matrix[i];
    }

    static Matrix<T> &zero_matrix(size_t n, size_t m) {
        std::vector<T> B(m, 0);
        std::vector<std::vector<T>> A(n, B);
        Matrix<T>* Zereos = new Matrix(n, m, A);
        return *Zereos;
    }

    static Matrix<T> &E_matrix(size_t n){
        std::vector<T> B(n, 0);
        std::vector<std::vector<T>> A(n, B);
        for(int i = 0; i < n; ++i){
            A[i][i] = 1;
        }
        Matrix<T>* E = new Matrix(n, n, A);
        return *E;
    }

};


template<class T>
std::ostream &operator<<(std::ostream &out, Matrix<T> A) {
    for (size_t i = 0; i < A.get_lines_number(); ++i) {
        for (size_t j = 0; j < A.get_colomns_number(); ++j) {
            out << A[i][j] << " ";
        }
        out << "\n";
    }
    return out;
}

template<class T>
Matrix<T> &operator*(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_lines_number()) {
        std::cerr << "Error: Unable to perfom matrix multiplication";
        exit(EXIT_FAILURE);
    }
    std::vector<std::vector<T>> C;
    for (size_t i = 0; i < A.get_lines_number(); ++i) {
        std::vector<T> new_line;
        for (size_t j = 0; j < B.get_colomns_number(); ++j) {
            T new_element = 0;
            for (size_t k = 0; k < A.get_colomns_number(); ++k) {
                new_element += A[i][k] * B[k][j];
            }
            new_line.push_back(new_element);
        }
        C.push_back(new_line);
    }
    Matrix<T>* new_matrix = new Matrix<T>(A.get_lines_number(), B.get_colomns_number(), C);
    return *new_matrix;
}


int main() {
    std::cout << Matrix<int>::zero_matrix(2, 2);
    std::cout << Matrix<int>::E_matrix(3);
    Matrix<int> a;
    Matrix<int> b;

    std::cout << a << "\n" << b << "\n" << a*b;
    return 0;
}
