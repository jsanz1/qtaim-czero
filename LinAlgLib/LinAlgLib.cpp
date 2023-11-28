// Lineal Algebra Libraries
#include <iostream>
#include <iomanip>
#include <fstream>  // ifstream, open(), close()
#include <sstream>  // istringstream
#include <random>
#include "LinAlgLib.h"

void getFileName(int argc, char **argv, std::string& inpFile,
                 std::string& outFile) {
    if(argc == 3){
         inpFile = argv[1];
         outFile = argv[2];
    }          
    else if(argc < 2){
         std::cout << "Input filename: " << std::endl;
         std::cin >> inpFile;
    }          
    else       
         inpFile = argv[1];
}    

std::vector<std::vector<double>> readMatrix(const std::string nombreArchivo) {
    std::ifstream Archivo;
    Archivo.open(nombreArchivo, std::fstream::in);
    std::vector<double> Vector {};
    double valor {};
    std::string linea {};
    std::vector<std::vector<double>> matriz;

    if(Archivo.is_open()){
         while(!Archivo.eof()){
              getline (Archivo, linea);
              std::istringstream is(linea);
              if (linea != "") {
                   while(is >> valor){
                        Vector.push_back(valor);
                   }
                   matriz.push_back(Vector);
                   Vector.clear();
              }
         }
    }
    else
         std::cout << "Cannot open file." << std::endl;

    Archivo.close(); 

    return matriz;
}

std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>> a){
               
//  Declaring and defining the vector with transpose elements
    std::vector<std::vector<double>> transposed(a.size(),std::vector<double>(a.size())); 
               
    for(size_t i=0; i < transposed.size(); ++i){
         for(size_t j=0; j < transposed.size(); ++j){
              transposed[i][j] = a[j][i];
         }     
    }          
               
    return transposed;
}   

std::vector<std::vector<double>> MatSubs(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b){
               
    std::vector<std::vector<double>> c(a.size(),std::vector<double>(a.size()));
               
    for(size_t i = 0; i < a.size(); ++i){
         for(size_t j = 0; j < a.size(); ++j){
              c[i][j] = a[i][j]-b[i][j];
         }     
    }          
               
    return c;  
}

std::vector<std::vector<double>> MatAdd(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b){
               
    std::vector<std::vector<double>> c(a.size(),std::vector<double>(a.size()));
               
    for(size_t i = 0; i < a.size(); ++i){
         for(size_t j = 0; j < a.size(); ++j){
              c[i][j] = a[i][j]+b[i][j];
         }     
    }          
               
    return c;  
}

std::vector<std::vector<double>> MatMul(const std::vector<std::vector<double>> a,
                          const std::vector<std::vector<double>> b){
               
    // Only for symmetrical matrices
    std::vector<std::vector<double>> c(a.size(),std::vector<double>(a.size()));
               
    for(size_t i = 0; i < a.size(); ++i){
      for(size_t j = 0; j < a.size(); ++j){
         for(size_t k = 0; k < a.size(); ++k){
                   c[i][j] += a[i][k]*b[k][j];
              }
         }     
    }                         
    return c;  
}  

std::vector<double> VecAdd(const std::vector<double> a,
                          const std::vector<double> b){

    std::vector<double> c(a.size());

    for(size_t i = 0; i < a.size(); ++i){
              c[i] = a[i]+b[i];
         }
    
    return c;
}

std::vector<double> VecSust(const std::vector<double> a,
                          const std::vector<double> b){

    std::vector<double> c(a.size());

    for(size_t i = 0; i < a.size(); ++i){
              c[i] = a[i] - b[i];
         }
    return c;
}

std::vector<double> MatVecMul(const std::vector<std::vector<double>> a,
                          const std::vector<double> b){

    // Only for symmetrical matrices
    std::vector<double> c(b.size());
    for(size_t i = 0; i < a.size(); ++i){
       for(size_t j = 0; j < b.size(); ++j){
                   c[i] += a[i][j]*b[j];
              }
         }    
    return c;
}

// Create a new ramdon elements Matrix to test
std::vector<std::vector<double>> randomMatrix(const size_t matSize, 
                                              const double inf, const double sup) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> random_real(inf, sup);
               
    std::vector<std::vector<double>> matriz(matSize,std::vector<double>(matSize));
 
    for (unsigned i = 0; i < matriz.size(); ++i) {
         for (unsigned j = 0; j < matriz[0].size(); ++j) {
              matriz[i][j] = random_real(mt);
         }     
    }          
    
    return matriz;
}  

                                                                        
// Create a new ramdon elements for a Symmetric Matrix to test
std::vector<std::vector<double>> randomSymMatrix(const size_t matSize, 
                                              const double inf, const double sup) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> random_real(inf, sup);
               
    std::vector<std::vector<double>> matriz(matSize,std::vector<double>(matSize));
 
    for (unsigned i = 0; i < matriz.size(); ++i) {
         for (unsigned j = i; j < matriz[0].size(); ++j) {
              matriz[i][j] = random_real(mt);
              matriz[j][i] = matriz[i][j];
         }     
    }          
    
    return matriz;
}  

// Creates a Identity Matriz of size "matSize"
std::vector<std::vector<double>> Identity(const size_t matSize){ 

    std::vector<std::vector<double>> a(matSize,std::vector<double>(matSize));
    // a = std::vector<std::vector<double> >(matSize, std::vector<double>(matSize, 0.));
    for (size_t i = 0; i < matSize; ++i){
         a[i][i] = 1.;
    }

    return a;

}

// Print Matrix elements function
void printMatrix(const std::vector<std::vector<double>> A){
    for(size_t i = 0; i < A.size(); ++i){
         for(size_t j = 0; j < A[i].size(); ++j){
              // 6 decimal + 8 = 14
              std::cout << std::scientific << std::setw(14) 
                        << std::setprecision(6) << A[i][j];
         }     
         std::cout << std::endl;
    }          
} 
  
//--------------------------------------------------------------------//
// Makes Jacobi rotation to diagonalize a matriz
// Diagonalize "matrix" to store in eigenValue and eigenVector 
// matrices pass by reference (matrices that already exist)
//--------------------------------------------------------------------//
void JacobiRotation(const std::vector<std::vector<double>> matriz,
                    std::vector<std::vector<double>> & eigenValue,
                    std::vector<std::vector<double>> & eigenVector){
    int const MAXITER = 1000;
    double const epsilon {1e-12};
    eigenValue = matriz;
    eigenVector = Identity(eigenValue.size());     
    for (int iter = 0; iter < MAXITER; ++iter) {
    // safe guard
         if (iter == MAXITER - 1) {
              std::cout << "MAXITER exceeded" << std::endl;
              return;
         }
    // look for largest non diagonal element
    // first, zero out diagonal elements
    std::vector<std::vector<double> > eigenZero = eigenValue;
    for (size_t i = 0; i < matriz.size(); ++i)
         eigenZero[i][i] = 0.;
    //printMatrix(eigenZero);     
    double nonDiag {};
    int iMax, jMax;
    for (size_t i = 0; i < eigenZero[i].size()-1; ++i) {
         for (size_t j = i+1; j < eigenZero[i].size(); ++j) {
              if (std::abs(eigenZero[i][j]) > nonDiag) {         
                   nonDiag = std::abs(eigenZero[i][j]);
                   iMax = i;
                   jMax = j;
              }
         }     
    }
    // begin diagonalizing
    double theta;
    if (std::abs(eigenValue[iMax][jMax]) > epsilon) {
     //    if (std::abs(eigenValue[iMax][iMax]-eigenValue[jMax][jMax]) < epsilon) 
     //         theta = 0.7853981633974483;
     //    else {
              theta = atan((2.*eigenValue[iMax][jMax])/
                      (eigenValue[iMax][iMax]-eigenValue[jMax][jMax]))/2.;
       //  }
         // create U
         std::vector<std::vector<double>> U = Identity(eigenValue.size());     
         U[iMax][iMax] = cos(theta);
         U[jMax][jMax] = U[iMax][iMax];
         U[iMax][jMax] = -sin(theta);
         U[jMax][iMax] = -U[iMax][jMax];
         // diagonalizing
         eigenValue = MatMul(Transpose(U),MatMul(eigenValue,U));
         eigenVector = MatMul(eigenVector,U);
    } // if (std::abs(eigenValue[iMax][jMax]) > epsilon) 
    else {
         // std::cout << "In " << iter << " steps" << std::endl;
         // Sort eigenvalues and eigenvectors 0,1,2
          for (size_t i = 0; i < eigenValue.size(); ++i) {
               auto iter = i;
               double valor = eigenValue[i][i];
               for (size_t j = i+1; j < eigenValue.size(); ++j) {
                    if (eigenValue[j][j] > valor) {
                         iter = j;
                         valor = eigenValue[j][j];
                    }
               }
               if (iter != i) {
                    eigenValue[iter][iter] = eigenValue[i][i];
                    eigenValue[i][i] = valor;
                    for (size_t j = 0; j < eigenValue.size(); ++j) {
                         valor = eigenVector[j][i];
                         eigenVector[j][i] = eigenVector[j][iter];
                         eigenVector[j][iter] = valor;
                    }
               }
          }
         // Sort eigenvalues and eigenvetors 3,4,5 
     /*       for (size_t i = 3; i < 6; ++i) {
               auto iter = i;
               double valor = eigenValue[i][i];
               for (size_t j = i+1; j < 6; ++j) {
                    if (eigenValue[j][j] > valor) {
                         iter = j;
                         valor = eigenValue[j][j];
                    }
               }
               if (iter != i) {
                    eigenValue[iter][iter] = eigenValue[i][i];
                    eigenValue[i][i] = valor;
                    for (size_t j = 3; j < 6; ++j) {
                         valor = eigenVector[j][i];
                         eigenVector[j][i] = eigenVector[j][iter];
                         eigenVector[j][iter] = valor;
                    }
               }
          }

         // std::cout << "eigenValue" << std::endl;
         // printMatrix(eigenValue);
         // std::cout << "eigenVector" << std::endl;
         // printMatrix(eigenVector);
         // std::cout << "Vec^T A Vec = Val" << std::endl;
         // printMatrix(MatMul(Transpose(eigenVector),MatMul(matriz,eigenVector)));
       */ return;
    }
    } // for (int iter = 0; iter < MAXITER; ++iter) 
}

