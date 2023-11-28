// WaveFunction.h
// WaveFunction Base Class
#include <iostream>
#include <iomanip>
#include <array> 
#include <vector>

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

class WaveFunction {
public:
    // Constructor
    WaveFunction(){};
    // Destructor
    ~WaveFunction(){};

    // Variables
    int nMO_;
    int nPrim_;
    int nCentre_;
    std::vector< std::vector<double> > Coord_;
    std::vector< std::vector<double> > coeffMO_;
    std::vector<double> Centre_;
    std::vector<double> Alpha_;
    std::vector<double> Type_;
    std::vector<double> Occ_;

    // Angular Cartersian 
    // Gamess, mcpgrd.src file, line 1933  
    std::array<std::array<int,3>,56> primType_ = {{
          {0,0,0},
         {1,0,0}, {0,1,0}, {0,0,1},
         {2,0,0}, {0,2,0}, {0,0,2},
         {1,1,0}, {1,0,1}, {0,1,1},
         {3,0,0}, {0,3,0}, {0,0,3},
         {1,2,0}, {2,1,0}, {2,0,1},
         {1,0,2}, {0,1,2}, {0,2,1},
         {1,1,1},
         {4,0,0}, {0,4,0}, {0,0,4},
         {3,1,0}, {3,0,1}, {1,3,0},
         {0,3,1}, {1,0,3}, {0,1,3},
         {2,2,0}, {2,0,2}, {0,2,2},
         {2,1,1}, {1,2,1}, {1,1,2},
         {5,0,0}, {0,5,0}, {0,0,5},
         {4,1,0}, {0,4,1}, {1,0,4},
         {0,1,4}, {1,4,0}, {1,2,2},
         {1,0,4}, {4,1,0}, {2,3,0},
         {2,1,2}, {0,5,0}, {0,3,2},
         {0,1,4}, {4,0,1}, {2,2,1},
         {2,0,3}, {0,4,1}, {0,2,3}
        }};
    
    // Method to assign values                                          
    void Read(int &argc, char *argv[]);          

    // Method to show atom and coordinates
    void showAtoms();
       // Method to show wfn file info
    void printInfo();
};
#endif

