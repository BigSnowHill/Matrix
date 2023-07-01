#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

template<class T>
class Matrix {
private:
    size_t N, M;
    std::vector<std::vector<T>> matrix;
public:
    // Cunstractor if all elemenets are known
    explicit Matrix<T>(size_t n, size_t m, std::vector<std::vector<T>> new_matrix) : N(n), M(m), matrix(new_matrix) {};

    // Cunstractor if only size is known, take elemenets from console
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

    // Cunstractor take size and elements from console
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

    // Copy Constractor
//    Matrix(Matrix<T> A) : N(A.get_lines_number()), M(A.get_colomns_number()) {
//        std::vector<T> line(N, 0);
//        std::vector<std::vector<T>> lines(M, line);
//        for (int i = 0; i < N; ++i) {
//            for (int j = 0; j < M; ++j) {
//                lines[i][j] = A[i][j];
//            }
//        }
//        this->matrix = lines;
//    }

    // File Constractor
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

    void operator*=(T number) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                matrix[i][j] *= number;
            }
        }
    }

    static Matrix<T> zero_matrix(size_t n, size_t m) {
        std::vector<T> B(m, 0);
        std::vector<std::vector<T>> A(n, B);
        Matrix<T>* Zereos = new Matrix(n, m, A);
        return *Zereos;
    }

    static Matrix<T> E_matrix(size_t n) {
        std::vector<T> B(n, 0);
        std::vector<std::vector<T>> A(n, B);
        for (int i = 0; i < n; ++i) {
            A[i][i] = 1;
        }
        Matrix<T> *E = new Matrix(n, n, A);
        return *E;
    }

};

template<class T>
T determinant(Matrix<T>);

template<class T>
Matrix<T> Transporation(Matrix<T> A) {
    if (A.get_lines_number() != A.get_colomns_number()) {
        std::cerr << "Impossible to make operation for non square matrix.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<T> line(A.get_colomns_number(), 0);
    std::vector<std::vector<T>> lines(A.get_lines_number(), line);
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            lines[j][i] = A[i][j];
        }
    }
    Matrix<T>* B = new Matrix<T>(A.get_lines_number(), A.get_colomns_number(), lines);
    return *B;
}

template<class T>
T Alg_Ad(Matrix<T> A, size_t x, size_t y) {
    if (A.get_lines_number() != A.get_colomns_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    if (A.get_lines_number() == 1) {
        return A[0][0];
    }
    std::vector<std::vector<T>> B_v;
    for (int i = 0; i < A.get_lines_number(); ++i) {
        if (i == x) continue;
        std::vector<T> line;
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            if (j == y) continue;
            line.push_back(A[i][j]);
        }
        B_v.push_back(line);
    }
    Matrix<T> B(A.get_lines_number() - 1, A.get_colomns_number() - 1, B_v);
    return pow(-1, x + y) * determinant<T>(B);
}

template<class T>
T determinant(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    T det = 0;
    for (int i = 0; i < A.get_colomns_number(); ++i) {
        det += A[0][i] * Alg_Ad<T>(A, 0, i);
    }
    return det;
}

template<class T>
Matrix<double> & operator!(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<T> line(A.get_colomns_number(), 0);
    std::vector<std::vector<T>> lines(A.get_lines_number(), line);
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            lines[i][j] = Alg_Ad<T>(A, i, j);
        }
    }
    Matrix<T>* C = new Matrix<T>(A.get_lines_number(), A.get_colomns_number(), lines);
    *(C) *= (1 / determinant(A));
    Matrix<T>* B =  new Matrix<T>(Transporation(*(C)));
    return *(B);
}


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
        std::cerr << "Error: impossible to multiply matrices.\n";
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
    Matrix<T> *new_matrix = new Matrix<T>(A.get_lines_number(), B.get_colomns_number(), C);
    return *new_matrix;
}

template<class T>
Matrix<T> &operator*(Matrix<T> A, T c)
{
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    std::vector<T> line(m, 0);
    std::vector<std::vector<T>> lines(n, line);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            lines[i][j] = A[i][j] * c;
        }
    }
    Matrix<T> *Ac = new Matrix<T>(n, m, lines);
    return *Ac;
}

template<class T>
Matrix<T> &operator*(T c, Matrix<T> A) {
    return A * c;
}

template<class T>
Matrix<T> &operator+(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        std::cerr << "Error: impossible to add matrices.";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    std::vector<T> line(m, 0);
    std::vector<std::vector<T>> lines(n, line);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            lines[i][j] = A[i][j] + B[i][j];
        }
    }
    Matrix<T> *C = new Matrix<T>(n, m, lines);
    return *C;
}

template<class T>
Matrix<T> &operator-(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        std::cerr << "Error: impossible to subtract matrices.";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    std::vector<T> line(m, 0);
    std::vector<std::vector<T>> lines(n, line);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            lines[i][j] = A[i][j] - B[i][j];
        }
    }
    Matrix<T> *C = new Matrix<T>(n, m, lines);
    return *C;
}


template<class T>
bool operator==(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        return false;
    }
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < B.get_colomns_number(); ++j) {
            if (A[i][j] != B[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
bool operator!=(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        return true;
    }
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < B.get_colomns_number(); ++j) {
            if (A[i][j] == B[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
bool operator==(Matrix<T> &A, T c) {
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    if (c == 0) {
        Matrix<T> B = Matrix<T>::zero_matrix(n, m);
        return A == B;
    }
    if (n != m) {
        std::cerr << "Impossible to compare not square matrix and scalar.";
        exit(EXIT_FAILURE);
    }
    Matrix<T> B = Matrix<T>::E_matrix(n);
    B *= c;
    return A == B;
}


template<class T>
bool operator==(T c, Matrix<T> &A) {
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    if (c == 0) {
        Matrix<T> B = Matrix<T>::zero_matrix(n, m);
        return A == B;
    }
    if (n != m) {
        std::cerr << "Impossible to compare not square matrix and scalar.";
        exit(EXIT_FAILURE);
    }
    Matrix<T> B = Matrix<T>::E_matrix(n);
    B *= c;
    return A == B;
}