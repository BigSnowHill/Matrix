#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>
#include <functional>
#include <future>


size_t block_size = 2;

template<typename T>
class Matrix {
private:
    size_t N, M;
    std::vector <std::vector<T>> matrix;
public:
    // Cunstractor if all elemenets are known
    explicit Matrix(size_t n, size_t m, std::vector <std::vector<T>> new_matrix) : N(n), M(m), matrix(new_matrix) {};

    // Cunstractor if only size is known, take elemenets from console
    Matrix(size_t n, size_t m) : N(n), M(m) {
        for (size_t i = 0; i < N; ++i) {
            std::vector <T> elem_vector;
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
            std::vector <T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    // Copy Constractor
    Matrix(Matrix<T> &A) : N(A.get_lines_number()), M(A.get_colomns_number()) {
        std::vector <T> line(N, 0);
        std::vector <std::vector<T>> lines(M, line);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                lines[i][j] = A[i][j];
            }
        }
        matrix = lines;
    }

    // File Constractor
    Matrix(const std::string &filename) {
        std::ifstream file_in(filename);
        file_in >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector <T> elem_vector;
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

    std::vector <T> &operator[](size_t i) {
        return matrix[i];
    }


    static Matrix<T> E_matrix(size_t n) {
        Matrix<T> M(n, n);
        std::function<void(size_t, size_t)> func;

        func = [&](size_t start, size_t end) {
            if (end - start <= block_size) {
                for (size_t i = start; i < end; ++i) {
                    if (i / n == i % n)
                        M[i / n][i % n] = 1;
                }
                return;
            }

            size_t mid = (start + end) / 2;
            std::future<void> async_func1 = std::async(func, start, mid);
            std::future<void> async_func2 = std::async(func, mid, end);

            async_func1.wait();
            async_func2.wait();

        };

        func(0, n * n);

        return M;
    }

    static Matrix<T> zero_matrix(size_t n, size_t m) {
        std::vector <T> B(m, 0);
        std::vector <std::vector<T>> A(n, B);
        Matrix<T> Zereos = new Matrix(n, m, A);
        return Zereos;
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
    std::vector <T> line(A.get_colomns_number(), 0);
    std::vector <std::vector<T>> lines(A.get_lines_number(), line);
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            lines[j][i] = A[i][j];
        }
    }
    Matrix<T> *B = new Matrix<T>(A.get_lines_number(), A.get_colomns_number(), lines);
    return *B;
}

template<class T>
T Minor(Matrix<T> A, size_t I, size_t J) {
    if (A.get_lines_number() != A.get_colomns_number()) {
        std::cerr << "Impossible to make operation for non square matrix.\n";
        exit(EXIT_FAILURE);
    }
    size_t k = 0;
    Matrix<T> A(Matrix<T>::zero_matrix(A.get_lines_number(), A.get_colomns_number()));
    for (size_t i = 0; i < n; ++i) {
        size_t h = 0;
        for (size_t j = 0; j < n; ++j) {
            if (i != I && j != J) {
                B[k][h] = A[i][j];
                h++;
            }
        }
        if (i != I)
            k++;
    }
    return determinant(B);

}

template<class T>
T Alg_Ad(Matrix<T> A, size_t i, size_t j) {
    return pow(-1, i + j) * Minor_async(A, i, j);
}

template<class T>
T determinant(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    if (n == 1) return data[0][0];

    T det = 0;
    std::queue <std::future<T>> sub_det;
    auto func = [&](size_t i) { return Alg_Ad(A, 0, i); };
    for (size_t i = 0; i < m; ++i) { sub_det.push(std::async(func, i)); }
    for (size_t i = 0; i < m; ++i) {
        det += data[0][i] * sub_det.front().get();
        sub_det.pop();
    }
    return det;
}


template<class T>
Matrix<double> operator!(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number();
    size_t start = 0;
    size_t end = n * n;
    std::future <T> detMatrix_f = std::async(std::launch::async, [&](Matrix<T> A) { return determinant(A); });
    Matrix<double> Res(Matrix<double>::zero_matrix(n, n));
    std::function<void(size_t, size_t)> func;
    func = [&](size_t start, size_t end) {
        if (end - start <= block_size) {
            for (size_t i = start; i < end; ++i) {
                Res[i / n][i % n] = Alg_Ad(A, i / n, i % n);
            }
            return ;
        }
        size_t mid = (start + end) / 2;
        std::future<void> async_func1 = std::async(func, start, mid);
        std::future<void> async_func2 = std::async(func, mid, end);

        async_func1.wait();
        async_func2.wait();
    };

    func(start, end);
    std::future <Matrix<double>> transRes_f = std::async(std::launch::async, [&Res]() { return Transporation(Res); });
    T detMatrix = detMatrix_f.get();
    if (detMatrix == 0) {
        std::cerr << "Impossible to find inverse for matrix with determinant equal 0.\n";
        exit(EXIT_FAILURE);
    }
    Matrix<double> transRes = transRes_f.get();
    return transRes * (1.0 / detMatrix);
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
Matrix<T> operator*(Matrix<T> &N, Matrix<T> &M){
    if (N.get_colomns_number() != M.get_lines_number()) {
        std::cerr << "Error: impossible to multiply matrices.\n";
        exit(EXIT_FAILURE);
    }
    size_t n_N = N.get_lines_number(), m_N = N.get_colomns_number();
    size_t n_M = M.get_lines_number(), m_M = M.get_colomns_number();
    size_t start = 0;
    size_t end = n_N * m_M;
    Matrix<T> res(n_N, m_M);
    std::function<T(size_t, size_t)> el_func;
    std::function<void(size_t, size_t)> func;
    el_func = [&](size_t line, size_t column){
        T sum = 0;
        for (size_t i = 0; i < m_N; ++i){
            sum += N[line][i] * M[i][column];
        }
        return sum;
    };
    func = [&](size_t start, size_t end){
        if (end - start <= block_size){
            for (size_t i = start; i < end; ++i){
                res[i / m_M][i % m_M] = el_func(i / m_M, i % m_M);
            }
            return ;
        }
        size_t mid = (end + start) / 2;
        std::future<void> async_func1 = std::async(func, start, mid);
        std::future<void> async_func2 = std::async(func, mid, end);
        async_func1.wait();
        async_func2.wait();
    };
    func(start, end);
    return res;
}

template<class T>
Matrix<T> operator*(Matrix<T> A, T c) {
    size_t start = 0;
    size_t end = A.get_lines_number() * A.get_colomns_number();
    std::function<void(size_t, size_t)> func;
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    func = [&](size_t start, size_t end) {
        if (end - start <= block_size) {
            for (size_t i = start; i < end; ++i) {
                A[i / n][i % m] *= c;
            }
            return;
        }
        size_t mid = (start + end) / 2;
        std::future<void> async_func1 = std::async(func, start, mid);
        std::future<void> async_func2 = std::async(func, mid, end);
        async_func1.wait();
        async_func2.wait();
    };
    func(start, end);
    Matrix<T>* B = new Matrix<T>(A);
    return *B;
}

template<class T>
Matrix<T> operator*( T c, Matrix<T> A){
    return A * c;
}

template<class T>
Matrix<T> operator+(Matrix<T> &N, Matrix<T> &M) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        std::cerr << "Error: impossible to add matrices.";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    size_t start = 0;
    size_t end = n * m;
    std::function<void(size_t, size_t)> func;
    func = [&](size_t start, size_t end) {
        if (end - start <= block_size) {
            for (size_t i = start; i < end; ++i) {
                N[i / n][i % m] += M[i / n][i % m];
            }
            return;
        }
        size_t mid = (start + end) / 2;
        std::future<void> async_func1 = std::async(func, start, mid);
        std::future<void> async_func2 = std::async(func, mid, end);
        async_func1.wait();
        async_func2.wait();
    };
    func(start, end);
    Matrix<T>* A = new Matrix<T>(N);
    return *A;
}

template<class T>
Matrix<T> operator-(Matrix<T> &N, Matrix<T> &M) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        std::cerr << "Error: impossible to add matrices.";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    size_t start = 0;
    size_t end = n * m;
    std::function<void(size_t, size_t)> func;
    func = [&](size_t start, size_t end) {
        if (end - start <= block_size) {
            for (size_t i = start; i < end; ++i) {
                N[i / n][i % m] -= M[i / n][i % m];
            }
            return;
        }
        size_t mid = (start + end) / 2;
        std::future<void> async_func1 = std::async(func, start, mid);
        std::future<void> async_func2 = std::async(func, mid, end);
        async_func1.wait();
        async_func2.wait();
    };
    func(start, end);
    Matrix<T>* A = new Matrix<T>(N);
    return *A;
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
