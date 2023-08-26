//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial Program for spinoidal decomposition using Finite Difference method.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 17.02.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=SpinoidalDecomposition.cpp
// make run
//________________________________________________________________________________________________

#include "General.h"

using namespace std::chrono;

#define TECPLOT
//#define VISIT

const int Nx = 512, Ny = 512, SaveInterval = 2000;
const double dt = 1.0E-3, Tf = 250.0;
const double Lx = Nx, Ly = Ny;

const double M = 1.0, kappa = 0.5, alpha = 1.0;

double **C, **D, Laplacian, mu, dCdx, dCdy, dx, dy, t = 0.0;
int ip, im, jp, jm, TimeStep = 0;

char *s;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Allocate ()
{
    Allocate(C,Nx,Ny);
    Allocate(D,Nx,Ny);
    Allocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(C,Nx,Ny);
    Deallocate(D,Nx,Ny);
    Deallocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            C[i][j] = 0.4 + 0.02*(0.5 - RandomNumber());
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double ComputeEnergy ( double **C )
{
    double energy = 0.0;
    
    for (int i{}; i < Nx; ++i)
    {
        ip = (i == (Nx-1) ? 0 : i+1);
        im = (i == 0 ? (Nx-1) : i-1);
        
        for (int j{}; j < Ny; ++j)
        {
            jp = (j == (Ny-1) ? 0 : j+1);
            jm = (j == 0 ? (Ny-1) : j-1);
            
            dCdx = 0.5*(C[ip][j] - C[im][j])/dx;
            dCdy = 0.5*(C[i][jp] - C[i][jm])/dy;
            
            energy += alpha*C[i][j]*C[i][j]*(1.0 - C[i][j])*(1.0 - C[i][j]) + 0.5*kappa*(dCdx*dCdx + dCdy*dCdy);
        }
    }
    
    return energy*dx*dy;
}

//________________________________________________________________________________________________
// Main code
//________________________________________________________________________________________________
int main()
{
    std::cout.flags( std::ios::dec | std::ios::fixed );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Allocate();
    
    dx = Lx/Nx;
    dy = Ly/Ny;
    
    Initialize();
    
    // Update concentration
    while (t < Tf)
    {
        std::cout << "Time step = " << TimeStep << ", Current time = " << t << ", Energy = " << ComputeEnergy(C) << std::endl;
        
        // Step 1: Compute D
        for (int i{}; i < Nx; ++i)
        {
            ip = (i == (Nx-1) ? 0 : i+1);
            im = (i == 0 ? (Nx-1) : i-1);
            
            for (int j{}; j < Ny; ++j)
            {
                jp = (j == (Ny-1) ? 0 : j+1);
                jm = (j == 0 ? (Ny-1) : j-1);
                
                Laplacian = (C[ip][j] - 2.0*C[i][j] + C[im][j])/(dx*dx) + (C[i][jp] - 2.0*C[i][j] + C[i][jm])/(dy*dy);
                
                mu = 2.0*alpha*C[i][j]*(1.0 - C[i][j])*(1.0 - 2.0*C[i][j]);
                
                D[i][j] = mu - kappa*Laplacian;
            }
        }
        
        // Step 2: Update C
        for (int i = 0; i < Nx; ++i)
        {
            ip = (i == (Nx-1) ? 0 : i+1);
            im = (i == 0 ? (Nx-1) : i-1);
            
            for (int j = 0; j < Ny; ++j)
            {
                jp = (j == (Ny-1) ? 0 : j+1);
                jm = (j == 0 ? (Ny-1) : j-1);
                
                Laplacian = (D[ip][j] - 2.0*D[i][j] + D[im][j])/(dx*dx) + (D[i][jp] - 2.0*D[i][j] + D[i][jm])/(dy*dy);
                
                C[i][j] += dt*M*Laplacian;
                
                if (C[i][j] >= 0.99999) C[i][j] = 0.99999;
                if (C[i][j] <= 0.00001) C[i][j] = 0.00001;
            }
        }
        
        // Step 3: Write file
        if (TimeStep % SaveInterval == 0)
        {
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d.tec",TimeStep/SaveInterval);
            
            std::ofstream Output(s,std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"C\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    Output << i*dx/Lx << "\t" << j*dy/Ly << "\t" << C[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%04d.vtk",TimeStep/SaveInterval);
            
            std::ofstream Output(s,std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Spinoidal decomposition" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << i*dx/Lx << "\t" << j*dy/Ly << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS C float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << C[i][j] << std::endl;
            
            Output.close();
            #endif
        }
        
        TimeStep++;
        
        t += dt;
        
        CheckNaN(C,Nx,Ny);
    }
    
    Deallocate();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
