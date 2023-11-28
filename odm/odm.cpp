//---------------------------------//
// Creación de ODM para la tesis EJAS 
// 10/02/23  
// Revisar documentación necesaria de la QTAIM
// (Bader) y archivo anexo
//---------------------------------//

#include <vector>
#include <cmath>
#include <fstream>
#include "../LinAlgLib/LinAlgLib.h"
#include "../WaveFunction/WaveFunction.h"
#include "odm.h"



//  Replace pow for a ~ 0 
double Chi::powC(double a, int b){                          
//  For (0)^(-) // change to -10
         if ( fabs(a) < 1e-20 && b < 0 ) {                                           
             return 0.;                                                     
         }  
//  For (0)^(0)  // < 1e-10  (not making sense)                                                       
         if ( fabs(a) < 1e-20 && b == 0 ) {                              
              return 1.;                                                
         } 
//  For (0.0)^(0)  // < 1e-10  (not making sense)                                                       
         if ( a == 0.0 && b == 0 ) {
              return 1.;
         }
         if ( a == 0.0 && b < 0 ) {
              return 0.;
         }
         if ( a == 0.0 && b > 0 ) {
              return 0.;
         }
// For (n)^(0) 
         if ( b == 0 ) {
              return 1.;
         }
         // Added
         if ( a > 0 && b < 0 ) {
              return div(1,pow(a,(-1)*b));
         }    
         else  return pow(a,b);                                               
}      

//  Replace a/b for exeption 0                                               
double Chi::div(int a, double b){                                      
         if ( fabs(b) < 1e-20 ) {                              
              return 0.0;                                                
         }           
         if ( a == 0 && fabs(b) < 1e-20 ) {
              return 0.0;
         }
         else  return a/b;                                         
} 

// Chi (primitive)   
std::vector<double> Chi::chi(const int &nPrim,                                                 
                             const std::vector<std::vector<double>>& Coord,   
                             const std::array<std::array<int,3>,56>& primType,
                             const std::vector<double>& Centre,                 
                             const std::vector<double>& Alpha,                 
                             const std::vector<double>& Type,                      
                             const std::vector<double>& r){

    // Evaluate chi value
    std::vector<double> chi(nPrim);
         
    for(int i=0; i < nPrim; ++i){
         chi[i]= powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                 powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                 powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])*
                 exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                + powC(r[2] - Coord[Centre[i]-1][2], 2)));           
     }
   
   // Asignar valor de chi al dato miembro
    return chi;
}                     

// Orbital Molecular
std::vector<double> Psi::psi (const int& nMO, const std::vector<double>& chi, 
                              const std::vector<std::vector<double>>& coeffMO){ 

    std::vector<double> MO(nMO);

    // For each Molecular orbital                                  
    for (int j = 0; j < nMO; ++j){                                 
         for (size_t i = 0; i < chi.size(); ++i) {                         
              MO[j] = MO[j] + coeffMO[j][i]*chi[i]; 
            // std::cout << MO[j] << std::endl; 
              }                                                         
         }    
         return MO;         
        
}

// Rho                                                               
void Rho::rho(double &rho_, const std::vector<double>& r,        
               WaveFunction& molecule){                            
                                                                        
    // Data from wfn                                                    
    int& nPrim = molecule.nPrim_;                                       
    int& nMO = molecule.nMO_;                                           
    std::vector< std::vector<double> >& Coord = molecule.Coord_;        
    std::array<std::array<int,3>,56>& primType = molecule.primType_;    
    std::vector<double>& Centre = molecule.Centre_;                     
    std::vector<double>& Alpha = molecule.Alpha_;                       
    std::vector<double>& Type = molecule.Type_;                         
    std::vector<double>& Occ = molecule.Occ_;                           
    std::vector<std::vector<double>>& coeffMO = molecule.coeffMO_;      
                                                                        
    // Chi evaluation                                                  
    std::vector<double> chi_r = chi(nPrim,Coord,primType,Centre,Alpha,Type,r);                                     
    // Psi evaluation
    std::vector<double> psi_r = psi(nMO,chi_r,coeffMO);                 

    // Rho evaluation
    for (int j = 0; j < nMO; ++j){                                      
         rho_ += Occ[j]*psi_r[j]*psi_r[j];                           
    }                                                                   
}

// Chi_1st_derivative 
std::vector<double> Gradient::chi_d(const int& d,
                                    const std::vector<std::vector<double>>& Coord,
                                    const std::array<std::array<int,3>,56>& primType,
                                    const std::vector<double>& Centre,          
                                    const std::vector<double>& Alpha,           
                                    const std::vector<double>& Type,            
                                    const std::vector<double>& r){ 
         
    std::vector<double> chi_d(Alpha.size());                 

    // Redefinición del gradiente de acuerdo a la ecuación [1.a],
    // archivo: resumen_errores_2021.pdf

    switch (d){  // d -> 0 para x, 1 para y, 2 para z
      case 0:  for (size_t i=0; i < Alpha.size(); ++i){   
                   chi_d[i]=
                       powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                       powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])* 
                       exp((-Alpha[i])*( powC(r[0] - Coord[Centre[i]-1][0], 2)
                                       + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                       + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                       (primType[Type[i]-1][0]*
                        powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]-1)
                        - 2*Alpha[i]*
                        powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]+1));
                }  
       break;   
              
         case 1:  for (size_t i=0; i < Alpha.size(); ++i){ 
                       chi_d[i]= 
                       powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                       powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])*
                       exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                      + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                      + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                       (primType[Type[i]-1][1]*
                        powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]-1)
                        - 2*Alpha[i]*
                          powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]+1));
                  }
         break;
              
         case 2:  for (size_t i=0; i < Alpha.size(); ++i){ 
                       chi_d[i]= 
                       powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                       powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                       exp((-Alpha[i])*( powC(r[0] - Coord[Centre[i]-1][0], 2)
                                       + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                       + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                      (primType[Type[i]-1][2]*
                       powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]-1) 
                       - 2*Alpha[i]*
                         powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]+1));

                  }
         break;
         }

         return chi_d;                                                         
}         

// Gradient 
void Gradient::gradient(std::vector<double> & gradient_, 
                	    const std::vector<double>& r,
                        WaveFunction& molecule){                          
                                                                      
    // Data from wfn                                                    
    int& nPrim = molecule.nPrim_;                                        
    int& nMO = molecule.nMO_;                                            
    std::vector<std::vector<double>>& Coord = molecule.Coord_;         
    std::array<std::array<int,3>,56>& primType = molecule.primType_;   
    std::vector<double>& Centre = molecule.Centre_;                      
    std::vector<double>& Alpha = molecule.Alpha_;                        
    std::vector<double>& Type = molecule.Type_;                          
    std::vector<double>& Occ = molecule.Occ_;                            
    std::vector<std::vector<double>>& coeffMO = molecule.coeffMO_;       
    
    // Resize grad vector                                                 
    gradient_.resize(3);                      
     
    //De r, evaluar Chi and Psi                                                 
    std::vector<double> chi_r  = chi(nPrim,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> psi_r =  psi(nMO,chi_r,coeffMO);                 
    
    // Derivatives for chi_dr (0,1,2 -> x,y,z)
    std::vector<double> chi_dx  = chi_d(0,Coord,primType,Centre,Alpha,Type,r);  
    std::vector<double> chi_dy  = chi_d(1,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dz  = chi_d(2,Coord,primType,Centre,Alpha,Type,r);
                                                                        
    std::vector<double> psi_dx =  psi(nMO,chi_dx,coeffMO);                 
    std::vector<double> psi_dy =  psi(nMO,chi_dy,coeffMO);  
    std::vector<double> psi_dz =  psi(nMO,chi_dz,coeffMO);  
                                                                                                          
    for (int j = 0; j < nMO; ++j){
         // Revisar definiciones del gradiente de la densidad 
         gradient_[0]  += Occ[j]*psi_dx[j]*psi_r[j]; 
         gradient_[1]  += Occ[j]*psi_dy[j]*psi_r[j]; 
         gradient_[2]  += Occ[j]*psi_dz[j]*psi_r[j];                           
    }                             
}

// Modificaciones para el gradiente 
 
// Método del campo gradiente
void Methods::gradTrayectory(const int &npoints, 
                              const double &direction, 
                              std::ofstream &print,
                              std::vector<double> &r, 
                              WaveFunction &molecule){

    // Constant Variables
    const float epsilon = 10e-12;                                       
    const float MaxStepSize = 0.2;                                      
    const float MinStepSize = 0.05;                                     
                                                                        
    //  Gradient vector and dg                                       
    float dg;                                           

    // Begin Trayectory                                                     
    for (auto i = 0; i < npoints; ++i) {                                
         std::vector<double> grad;
         std::vector<std::vector<double>> hess;
         double lap;
         gradient(grad,r,molecule);

         //  Magnitudes of gradients                                    
         dg = sqrt(grad[0]*grad[0]+                             
                   grad[1]*grad[1]+                             
                   grad[2]*grad[2]);                            
                                                                        
//  Check to see if dg is smaller than epsilon                          
         if((dg <= epsilon)){ 
            break;                         
            std::cout << "Ya llegamos" << std::endl;                  
            std::cout << dg << std::endl;               
         }                                                     

//  Check to see if dg is larger than maxStepSize              
         //  Verify with dg, dgp before and after                       
         if(dg > MaxStepSize){                                          
            grad[0] *= MaxStepSize/dg;                            
            grad[1] *= MaxStepSize/dg;                            
            grad[2] *= MaxStepSize/dg;                            
         }                                                              
                                                                        
//  Check to see if dg is larger than minStepSize                       
         //  Verify with dg, dgp before and after                       
         if(dg < MinStepSize){                                          
            grad[0] *= MinStepSize/dg;                            
            grad[1] *= MinStepSize/dg;                            
            grad[2] *= MinStepSize/dg;                            
         }                                                              
                                                                        
// Check to see if we're out of area of interest                        
// Update point                                                
// Next step is to add this to a vector (push_back)            
         r[0] += direction*grad[0];                                 
         r[1] += direction*grad[1];                                 
         r[2] += direction*grad[2];                                 
                
         hessian(hess, r, molecule);
         laplacian(lap, hess);

    print << std::fixed  << std::setprecision(10) 
          << std::setw(14) 
          << r[1]    << " "                                               
          << r[2]    << " "
          << -grad[1] << " "                                               
          << -grad[2] << " "
          << lap     << std::endl;

// Trayectory saved in print ofstream object         
    }
}       

// gradField modificado para obtener un campo gradiente
// en el plano yz

void Methods::gradField(WaveFunction &molecule){

    // Parameters
    int nPoints{100};
    std::vector<double> direction{0.25,-0.25};
    const double cirleAngle_{360.0}, radius{1.0};   // 360 Degrees (a circle)        
    const double nSubSeeds{100.};           // Number of seeeds around    
    const int nSubSeeds_int{100};
    const double PI{3.141592653};      // Pi value       
    // Object File to print
    std::ofstream print{"gradientFieldh2o.dat", std::ios::out | std::ios::app};

    // Values needed
    double angleDeg = cirleAngle_/nSubSeeds; // Angle in Degrees          
    double angleRad = angleDeg * PI/180.;   // Angle conversion to       
                                           // radians          
    // Starting at A nuclei
    // Declaring variables                                              
    double y_A{0.0}, z_A{0.0};        // y & z will vary and radius    
    double x_0 {1.0};                 // Diferentes valores de x 

    // For nucleus A as the centre point (y_0,z_0)
    double y_A0{molecule.Coord_[0][1]}, z_A0{molecule.Coord_[0][2]}; 
       
    // print << "Estamos en A" << std::endl;  // Comentario de control

    // For cicle: to evaluate for each angle and obtain y, z subseeds
    for (int i=0; i < nSubSeeds_int; ++i){                                    

        y_A  = y_A0 + radius * cos(i*angleRad);     // Evaluation of y and z     
        z_A  = z_A0 + radius * sin(i*angleRad);     // centered in (y0, z0)      
                                                                        
        std::vector<double> r{x_0, y_A, z_A};
  
         for (size_t i=0; i < direction.size(); ++i){
              gradTrayectory(nPoints,direction[i],print,r,molecule);
         }                                                           
    }                      
  
    // Starting at B nuclei                                             
    // Declaring variables                                              
    double y_B{0.0}, z_B{0.0};        // y & z will vary and radius    
    // Diferentes valores de x                                  
    
    // For nucleus A as the centre point (z_0,zp_0)                     
    double y_B0{molecule.Coord_[1][1]}, z_B0{molecule.Coord_[1][2]};   
                                        
    // print << "Estamos en B " << std::endl;

    // For cicle: to evaluate for each angle and obtain y, z subseeds
    for (int i=0; i < nSubSeeds_int; ++i){                                  
                                                                        
        y_B  = y_B0 + radius * cos(i*angleRad);     // Evaluation of x and y     
        z_B  = z_B0 + radius * sin(i*angleRad);     // centered in (x0, y0)      
                                                                       
        std::vector<double> r{x_0,y_B,z_B};                                   
                                                                        
        gradTrayectory(nPoints,direction[0],print,r,molecule);              
        //gradTrayectory(nPoints,directionDown,print,r,rp,molecule);              
                                                                        
    }                      

    print.close();
}

// Chi second derivative                                                       
std::vector<double> Hessian::chi_d2(const int& d1,const int & d2,
                                    const std::vector<std::vector<double>>& Coord,
                                    const std::array<std::array<int,3>,56>& primType,
                                    const std::vector<double>& Centre,  
                                    const std::vector<double>& Alpha,   
                                    const std::vector<double>& Type,    
                                    const std::vector<double>& r){      
                                                                        
    int d{d1+d2};
    std::vector<double> chi_d2(Alpha.size());                              
                                                                        
    switch (d){                                                    
         case 0:  for (size_t i=0; i < Alpha.size(); ++i){           
                      chi_d2[i]= exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                                + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                                + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                                 powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                                 powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])*
                                (powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]-2)*
                                 ((primType[Type[i]-1][0]*primType[Type[i]-1][0])-primType[Type[i]-1][0])
                                 + 4*Alpha[i]*Alpha[i]*
                                   powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]+2) 
                                 - Alpha[i]*(4*primType[Type[i]-1][0]+2)*
                                   powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]));
                      
         }                                                
         break;                                                    

         case 1:  for (size_t i=0; i < Alpha.size(); ++i){                
                      chi_d2[i]= exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                                + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                                + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                                 powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])*
                                 (primType[Type[i]-1][0]*
                                  powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]-1)
                                  - 2*Alpha[i]*
                                    powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]+1))*
                                 (primType[Type[i]-1][1]*
                                  powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]-1)
                                  - 2*Alpha[i]*
                                    powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]+1));
                   } 
         break;
                                                                        
         case 2:  if(d1==d2){
                    for (size_t i=0; i < Alpha.size(); ++i){           
                       chi_d2[i]= exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                                 + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                                 + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                                  powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                                  powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2])*
                                 (powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]-2)*
                                  ((primType[Type[i]-1][1]*primType[Type[i]-1][1])-primType[Type[i]-1][1]) 
                                  + 4*Alpha[i]*Alpha[i]*
                                    powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]+2)
                                  - Alpha[i]*(4*primType[Type[i]-1][1]+2)*
                                    powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]));
                    }    
                  }
                  else {
                    for (size_t i=0; i < Alpha.size(); ++i){
                        chi_d2[i]= exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                                  + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                                  + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                                   powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                                 (primType[Type[i]-1][0]*
                                  powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]-1)
                                  - 2*Alpha[i]*
                                    powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0]+1))*
                                 (primType[Type[i]-1][2]*
                                  powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]-1)
                                  - 2*Alpha[i]*
                                    powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]+1));
                       } 
                  } 
         break;                                                    
         
         case 3:
              for (size_t i=0; i < Alpha.size(); ++i){                
                   chi_d2[i] =  exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                               + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                               + powC(r[2] - Coord[Centre[i]-1][2], 2)))*
                                 powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                                 (primType[Type[i]-1][1]*
                                  powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]-1)
                                  - 2*Alpha[i]*
                                    powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1]+1))*
                                 (primType[Type[i]-1][2]*
                                  powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]-1)
                                  - 2*Alpha[i]*
                                    powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]+1));
                  }  
         break;

         case 4:  for (size_t i=0; i < Alpha.size(); ++i){
                        chi_d2[i]= exp((-Alpha[i])*(powC(r[0] - Coord[Centre[i]-1][0], 2)
                                                  + powC(r[1] - Coord[Centre[i]-1][1], 2)
                                                  + powC(r[2] - Coord[Centre[i]-1][2], 2)))*             
                                   powC(r[0]-Coord[Centre[i]-1][0], primType[Type[i]-1][0])*
                                   powC(r[1]-Coord[Centre[i]-1][1], primType[Type[i]-1][1])*
                                  (powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]-2)*
                                   ((primType[Type[i]-1][2]*primType[Type[i]-1][2])-primType[Type[i]-1][2])
                                   + 4*Alpha[i]*Alpha[i]*
                                     powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]+2)
                                   - Alpha[i]*(4*primType[Type[i]-1][2]+2)*
                                     powC(r[2]-Coord[Centre[i]-1][2], primType[Type[i]-1][2]));
                     
                  }
         break;                                                    
         }                                                              
                                                                        
         return chi_d2;                                                  
}                         

void Hessian::hessian(std::vector<std::vector<double>>& Hess, 
                       const std::vector<double>& r,
                       WaveFunction& molecule){       
                                                                        
    // Data from wfn                                                    
    int& nPrim = molecule.nPrim_;                                        
    int& nMO = molecule.nMO_;                                            
    std::vector< std::vector<double> >& Coord = molecule.Coord_;  
    std::array<std::array<int,3>,56>& primType = molecule.primType_;          
    std::vector<double>& Centre = molecule.Centre_;                      
    std::vector<double>& Alpha = molecule.Alpha_;                        
    std::vector<double>& Type = molecule.Type_;                          
    std::vector<double>& Occ = molecule.Occ_;                            
    std::vector<std::vector<double>>& coeffMO = molecule.coeffMO_;       
                                                                        
    // From r & rp, Chi and Psi                                                 
    std::vector<double> chi_r  = chi(nPrim,Coord,primType,Centre,Alpha,Type,r);  
    std::vector<double> psi_r =  psi(nMO,chi_r,coeffMO);                
        
    // Chi_d2 for r (dx,dy,dz)                                                                
    std::vector<double> chi_dx2   = chi_d2(0,0,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dxdy  = chi_d2(0,1,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dxdz  = chi_d2(0,2,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dy2   = chi_d2(1,1,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dydz  = chi_d2(1,2,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dz2   = chi_d2(2,2,Coord,primType,Centre,Alpha,Type,r);
 
    std::vector<double> chi_dx  = chi_d(0,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dy  = chi_d(1,Coord,primType,Centre,Alpha,Type,r);
    std::vector<double> chi_dz  = chi_d(2,Coord,primType,Centre,Alpha,Type,r);
    
    // Psi_d2 for r (dx,dy,dz)                         
    std::vector<double> psi_dx2 =  psi(nMO,chi_dx2,coeffMO);              
    std::vector<double> psi_dxdy = psi(nMO,chi_dxdy,coeffMO);              
    std::vector<double> psi_dxdz = psi(nMO,chi_dxdz,coeffMO);              
    std::vector<double> psi_dy2  = psi(nMO,chi_dy2,coeffMO);                                                        
    std::vector<double> psi_dydz = psi(nMO,chi_dydz,coeffMO);
    std::vector<double> psi_dz2  = psi(nMO,chi_dz2,coeffMO); 

    // For crossed terms, r & rp
    std::vector<double> psi_dx =  psi(nMO,chi_dx,coeffMO);              
    std::vector<double> psi_dy =  psi(nMO,chi_dy,coeffMO);              
    std::vector<double> psi_dz =  psi(nMO,chi_dz,coeffMO);              
 

    // Sum & evaluation
    // Hessian Matrix declaration (ab initio with 0)
    // std::vector< std::vector<double> > Hess(6,std::vector<double>(6));
       Hess.resize(3,std::vector<double>(3));

    // Hessian Matrix Hess_ from Apuntes A [pág. 2] 
    // Hess_[0][0] = dx2   Hess_[0][1] = dxdy   Hess_[0][2] = dxdz
    // Hess_[1][0] = dxdy  Hess_[1][1] = dy2    Hess_[1][2] = dydz 
    // Hess_[2][0] = dxdz  Hess_[2][1] = dydz   Hess_[2][2] = dz2     
    
    for (int j=0; j < nMO; ++j){                                       
         // r 
         Hess[0][0] += Occ[j]*psi_r[j]*psi_dx2[j];                        
         Hess[0][1] += Occ[j]*psi_r[j]*psi_dxdy[j];
         Hess[0][2] += Occ[j]*psi_r[j]*psi_dxdz[j];
         Hess[1][1] += Occ[j]*psi_r[j]*psi_dy2[j]; 
         Hess[1][2] += Occ[j]*psi_r[j]*psi_dydz[j];                   
         Hess[2][2] += Occ[j]*psi_r[j]*psi_dz2[j];          
    
        // Simetric terms
        Hess[1][0] = Hess[0][1];
        Hess[2][0] = Hess[0][2];
        Hess[2][1] = Hess[1][2];
    }
}

void Hessian::laplacian(double& laplacian_, 
                         const std::vector<std::vector<double>>& hess){
    
    double lap_diagElements{0.0};

    for (size_t i=0; i < hess[0].size() ; ++i){                                       
       lap_diagElements += hess[i][i]; 
    }
 
    laplacian_ = lap_diagElements; 
    
}






