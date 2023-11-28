// WaveFunction_class.cpp
// WaveFunction Base class 
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include "WaveFunction.h"
#include "../LinAlgLib/LinAlgLib.h"

// enum necesary for Read
enum Section{  
    CENTRE = 1,
    TYPE = 2,  
    EXPONENTS = 3,
    MO = 4,    
    ELSE= 5    
};                  

// enumSection necesary for Read
Section enumSection(const std::string& readSection){
    if (readSection == "CENTRE") return CENTRE;
    if (readSection == "TYPE") return TYPE;
    if (readSection == "EXPONENTS") return EXPONENTS;
    if (readSection == "MO") return MO;     
    else return ELSE;
}

//------------------------------------------------------------------//
//  Open wfn 
//------------------------------------------------------------------//
void WaveFunction::Read(int &argc, char *argv[]){
  //  Get file name from command line 
    std::string fileName;
    if(argc > 1){
        fileName = argv[1];
    }
    else {
         std::cout << "Please enter file name: " << std::endl;
         std::cin  >> fileName;
    }

    //  Read wfn                                          
    std::ifstream wfnFile{fileName, std::ios::in};
    if(wfnFile.is_open()) {
         // Read MO, nPrim_ and nCenter
         wfnFile.ignore(256,'\n');
         wfnFile.seekg(8, wfnFile.cur) >> nMO_;
         wfnFile.seekg(14, wfnFile.cur) >> nPrim_;
         wfnFile.seekg(11, wfnFile.cur) >> nCentre_;
         wfnFile.ignore(6,'\n');

         // CHECK::::::::::::::::::::::::::::::::::::::
         //std::cout << "nMo" << nMO_ << std::endl;
         //std::cout << "nPrim" << nPrim_ << std::endl;
         //std::cout << "nCentre" << nCentre_ << std::endl;

         //::::::::::::::::::::::::::::::::::::::::::::
         
         // Read Coord         
         std::string wfnLineC;
         std::getline(wfnFile, wfnLineC); 

 // HERE I NEED TO CHANGE:
 // Save the atomic symbols in a vector
 // To show on screen and choose them from the command line
 // ....(work on this latter)...

         for (int i=0; i < nCentre_; ++i){
              std::vector<double> xyz(3);
              std::getline(wfnFile, wfnLineC);
              std::istringstream thisLine(wfnLineC);
              thisLine.seekg(24, thisLine.cur)
                  >> xyz[0] >> xyz[1] >> xyz[2];
              Coord_.push_back(xyz);
         }   

         //std::cout << "x " << Coord_[1][0] << std::endl;
         //std::cout << "y " << Coord_[1][1] << std::endl;
         //std::cout << "z " << Coord_[1][2] << std::endl;

         // Read the Rest
         std::string wfnLine;        

         // Begin reading 
         while (std::getline(wfnFile, wfnLine)){
              int centreInt, typeInt;  
              double nOcc;
              std::string categoria, alphaVal, coefMoVal;
              std::istringstream thisLine(wfnLine);
              thisLine.seekg(0, thisLine.cur) >> categoria;
               
              switch (enumSection(categoria)){

              case ELSE:
                   thisLine.ignore(256,'\n');
                   break;
               
              case CENTRE:
                   thisLine.ignore(16, thisLine.cur);
                   while (thisLine >> centreInt){
                        Centre_.push_back(centreInt);
                   }
                   break;
              
               case TYPE:
                    thisLine.ignore(14, thisLine.cur);
                    while (thisLine >> typeInt){
                        Type_.push_back(typeInt);
                    }
                    break;
               
              case EXPONENTS:
                   thisLine.ignore(2, thisLine.cur);
                   while (thisLine >> alphaVal){
                        std::string rD ("D");
                        alphaVal.replace(alphaVal.find(rD), rD.length(), "e");
                        std::istringstream alphaNumbers(alphaVal); 
                        double alphaNum;    
                        alphaNumbers >> alphaNum;
                        Alpha_.push_back(alphaNum);
                   }
                   break;
               
              case MO:
                   // MO Occupation number
                   thisLine.seekg(34, thisLine.cur) >> nOcc;
                   Occ_.push_back(nOcc);

                   // MO coefficients
                   std::vector<double> coeffMO;
                   double nLinesMO = ceil(static_cast<double>(nPrim_)/5);
                   for (auto i = 0; i < nLinesMO; ++i){
                        std::string wfnLineMO;
                        std::string coeffMo;
                        std::getline(wfnFile, wfnLineMO);
                        std::istringstream actualLine(wfnLineMO);
               
                        while (actualLine >> coeffMo){
                             // String to find†
                             std::string rD("D");   
                             // Find and replace
                             coeffMo.replace(coeffMo.find(rD),rD.length(), "e");
                             std::istringstream coeffLine(coeffMo);
                             double coeffNum;
                             coeffLine >> coeffNum;
                             coeffMO.push_back(coeffNum);
                        }
                   }  
                   coeffMO_.push_back(coeffMO);
                   break;
              } // switch (enumSection(c))
         } // while (std::getline(wfnFile, wfnLine))
    } // if(wfnFile.is_open())

    else {     
         std::cerr << "Can't open file \""
                   << fileName << "\"."<< std::endl;
         exit(EXIT_FAILURE);
    }  

// Close the input file   
    wfnFile.close(); 
}

//------------------------------------------------------------------//
//  Show atom information to select two atoms
//------------------------------------------------------------------//
void WaveFunction::showAtoms(){

   std::cout << "Total number of Atoms: " << nCentre_ << "\n"
             << "N Atom - Coordinates: " << std::endl;
   for (int i=0; i < nCentre_; ++i){
            std::cout << i+1 << "  - " << std::fixed << std::setprecision(10)
                                       << Coord_[i][0] << ", "  
                                       << Coord_[i][1] << ", "
                                       << Coord_[i][2] << std::endl;
   }
}

void WaveFunction::printInfo(){
   
   
   //  File to write data
   std::ofstream infoFile{"info_wfn.dat", std::ios::out | std::ios::app};
   // Write on file
   infoFile << "-----------------------------------------------------------" << "\n"
            << "nAtoms = " << nCentre_ << "  nMO = " << nMO_ 
                                     << "  nPrim = " << nPrim_ << "\n"
            << "-----------------------------------------------------------" << std::endl; 
   infoFile << "CENTRE ASSIGNMENTS VECTOR " << std::endl;
            for (size_t i=0; i < Centre_.size() ; ++i){
                 infoFile << Centre_[i]  << ",";
            } infoFile << std::endl;
   infoFile << "-----------------------------------------------------------" << "\n"
            << "TYPE ASSIGNMENTS VECTOR " << std::endl;
            for (size_t i=0; i < Type_.size() ; ++i){
                 infoFile  << Type_[i]  << ",";
            } infoFile << std::endl;
   infoFile << "-----------------------------------------------------------" << "\n"
            << "ATOM COORDINATES " << std::endl;
            for(size_t i = 0; i < Coord_.size(); ++i){
               for(size_t j = 0; j < Coord_[i].size(); ++j){
                 // 12 decimals + 8 = 20
                 infoFile  << std::scientific << std::setw(20)
                           << std::setprecision(12) << Coord_[i][j] << ", ";
               }
              infoFile << std::endl;
            }
   infoFile << "-----------------------------------------------------------" << "\n"
            << "ALPHA EXPONENTS VECTOR " << std::endl;
            for (size_t i=0; i < Alpha_.size() ; ++i){
                 infoFile << Alpha_[i]  << ",";
            } infoFile << std::endl;
   infoFile << "-----------------------------------------------------------" << "\n"
            << "OCCUPATION NUMBERS VECTOR " << std::endl;
            for (size_t i=0; i < Occ_.size() ; ++i){
                 infoFile << Occ_[i]  << ",";
            } infoFile << std::endl;
   infoFile << "-----------------------------------------------------------" << "\n"
            << "COEFFICIENTS OM " << std::endl;
            for(size_t i = 0; i < coeffMO_.size(); ++i){
               infoFile << "{";
               for(size_t j = 0; j < coeffMO_[i].size(); ++j){
                 // 12 decimals + 8 = 20
                 infoFile << std::scientific << std::setw(20)
                          << std::setprecision(12) << coeffMO_[i][j] << ",";
               }
               infoFile << "}, " << std::endl;
            }
   infoFile << "-----------------------------------------------------------" << std::endl;
  
}

