// LinAlgLib.h

#include <vector>
#include <string>

#ifndef LINALGLIB_H
#define LINALGLIB_H

void getFileName(int argc, char **argv, std::string& inpFile,
                 std::string& outFile);

std::vector<std::vector<double>> readMatrix(const std::string fileName); 

std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>>
                                                      matrix);

std::vector<std::vector<double>> MatSubs(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b);

std::vector<std::vector<double>> MatAdd(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b);

std::vector<std::vector<double>> MatMul(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b);

std::vector<double> VecAdd(const std::vector<double> a,
                          const std::vector<double> b);

std::vector<double> VecSust(const std::vector<double> a,
                          const std::vector<double> b);

std::vector<double> MatVecMul(const std::vector<std::vector<double>> a,
                          const std::vector<double> b);

std::vector<std::vector<double>> Identity(const size_t matSize);

std::vector<std::vector<double>> randomMatrix(const size_t matSize, 
                                              const double inf, const double sup);

std::vector<std::vector<double>> randomSymMatrix(const size_t matSize, 
                                              const double inf, const double sup);

void JacobiRotation(const std::vector<std::vector<double>>,
                    std::vector<std::vector<double>> &,
                    std::vector<std::vector<double>> &);

void printMatrix(const std::vector<std::vector<double>> A);

#endif
