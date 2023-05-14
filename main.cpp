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
    explicit Matrix(size_t n, size_t m, std::vector<std::vector<T>> new_matrix): N(n), M(m), matrix(new_matrix){};

    Matrix(size_t n, size_t m) : N(n), M(m){
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for(size_t j = 0; j < M; ++j){
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    Matrix(){
        std::cin >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for(size_t j = 0; j < M; ++j){
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    Matrix(const std::string& filename){
        std::ifstream file_in(filename);
        file_in >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector<T> elem_vector;
            for(size_t j = 0; j < M; ++j){
                T element;
                file_in >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

};


int main() {
    Matrix<int> a;
    return 0;
}
