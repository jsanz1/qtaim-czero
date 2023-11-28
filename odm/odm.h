#ifndef ODM_H
#define ODM_H


/*
Clase con las definiciones de Chi, Psi, Rho y sus operadores diferenciales.
Los vectores son de tres dimensiones únicamente.
*/


#include <iostream>
#include <vector>
#include "../WaveFunction/WaveFunction.h"

class Chi {
public:
    // Methods
    double powC(double, int);
    double div(int,double);
    std::vector<double> chi(const int &, 
                            const std::vector<std::vector<double>>&,
                            const std::array<std::array<int,3>,56>&,
                            const std::vector<double>&,                 
                            const std::vector<double>&,                   
                            const std::vector<double>&,                    
                            const std::vector<double>&);                    
    // constructor 
    Chi(){};
    // destructor
    ~Chi(){};
};

class Psi : public Chi {
public:
    // Method
    std::vector<double> psi (const int &, const std::vector<double>&, 
                             const std::vector<std::vector<double>>&);
    // constructor
    Psi (){};
    // destructor
    ~Psi (){};
};

// RHO (density)
class Rho : public Psi {                                              
public:                                                                 
    // Method                                                           
    void rho (double &, const std::vector<double>&,
                  WaveFunction&);                                       
    // constructor                                                      
    Rho (){};                                                         
    // destructor                                                       
    ~Rho (){};                                                        
};       


class Gradient : public Psi {
public:      
    // Methods                                          
    std::vector<double> chi_d (const int &,                                
                               const std::vector<std::vector<double>>&,
                               const std::array<std::array<int,3>,56>&,
                               const std::vector<double>&,                 
                               const std::vector<double>&,                 
                               const std::vector<double>&, 
                               const std::vector<double>&); 
                
    void gradient (std::vector<double> &, const std::vector<double>&, 
                   WaveFunction&); 
                                      
    // constructor                                                      
    Gradient (){};                                                         
    // destructor                                                       
    ~Gradient (){};     
};

class Hessian : public Gradient {                                         
public:
    // Methods                                                                 
    std::vector<double> chi_d2 (const int &, const int &,                            
                                const std::vector<std::vector<double>>&, 
                                const std::array<std::array<int,3>,56>&,
                                const std::vector<double>&,              
                                const std::vector<double>&,              
                                const std::vector<double>&,              
                                const std::vector<double>&);             
    void hessian (std::vector<std::vector<double>>&, 
                  const std::vector<double>&,           
                  WaveFunction&);  
                     
    void laplacian (double &, const std::vector<std::vector<double>>&);
     
    // constructor                                                      
    Hessian (){};                                                      
    // destructor                                                       
    ~Hessian (){};                                                     
}; 

class Methods : public Hessian {
public:
    void gradTrayectory(const int &npoints,                             
                        const double &direction,                        
                        std::ofstream &print,                           
                        std::vector<double> &r,                         
                        WaveFunction &molecule);                        
                                                                        
    void gradField(WaveFunction&);    

};

#endif
