//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial program for simulating dendritic solidification.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 27.01.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=DendriticGrowth.cpp
// make run
//________________________________________________________________________________________________

#include "General.h"

using namespace std::chrono;

#define TECPLOT
//#define VISIT

const int Nx = 512, Ny = 512, SaveInterval = 50;
const double Lx = 9.0, Ly = 9.0, dx = Lx/Nx, dy = Ly/Ny;
const double xmin = 0.0, xmax = xmin + Lx, ymin = 0.0, ymax = ymin + Ly;
const double dt = 5.0E-5, Tf = 0.4;

const double tau = 3.0E-4, epsilonb = 1.0E-2, kappa = 1.8, delta = 2.0E-2;
const double beta = 6.0, alpha = 0.9, Gamma = 10.0, Teq = 1.0, Theta0 = 0.2, r0 = 5.0;

double **phi, **phiold, **T, **Told, **Epsilon, **Epsilon_theta;
double Laplacian, Theta, phix, phiy, term0, term1, term2, term3, phixp, phixm, phiyp, phiym, m, t = 0.0;
int ip, im, jp, jm, TimeStep = 0;

char *s;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Allocate ()
{
    Allocate(phi,Nx,Ny);
    Allocate(T,Nx,Ny);
    Allocate(phiold,Nx,Ny);
    Allocate(Told,Nx,Ny);
    Allocate(Epsilon,Nx,Ny);
    Allocate(Epsilon_theta,Nx,Ny);
    
    Allocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(phi,Nx,Ny);
    Deallocate(T,Nx,Ny);
    Deallocate(phiold,Nx,Ny);
    Deallocate(Told,Nx,Ny);
    Deallocate(Epsilon,Nx,Ny);
    Deallocate(Epsilon_theta,Nx,Ny);
    
    Deallocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            double r = (i-Nx*0.5)*(i-Nx*0.5) + (j-Ny*0.5)*(j-Ny*0.5);
            
            phi[i][j] = (r < r0 ? 1.0 : 0.0);
            
            T[i][j] = 0.0;
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double PhiX ( const int i, 
              const int j )
{
    ip = (i == (Nx-1) ? 0 : i+1);
    im = (i == 0 ? (Nx-1) : i-1);
    
    return 0.5*(phiold[ip][j]-phiold[im][j])/dx;
}

//________________________________________________________________________________________________
// 
//_______________________________________________________________________________________________
double PhiY ( const int i, 
              const int j )
{
    jp = (j == (Ny-1) ? 0 : j+1);
    jm = (j == 0 ? (Ny-1) : j-1);
    
    return 0.5*(phiold[i][jp]-phiold[i][jm])/dy;
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
    
    Initialize();
    
    // Update concentration
    while (t < Tf)
    {
        std::cout << "Time step = " << TimeStep << ", Current time = " << t << std::endl;
        
        // Step 1
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                phiold[i][j] = phi[i][j];
                Told[i][j] = T[i][j];
            }
        }
        
        // Step 2: Compute Epsilon and (d Epsilon/d theta)
        for (int i{}; i < Nx; ++i)
        {
            ip = (i == (Nx-1) ? 0 : i+1);
            im = (i == 0 ? (Nx-1) : i-1);
            
            for (int j{}; j < Ny; ++j)
            {
                jp = (j == (Ny-1) ? 0 : j+1);
                jm = (j == 0 ? (Ny-1) : j-1);
                
                phix = 0.5*(phiold[ip][j] - phiold[im][j])/dx;
                phiy = 0.5*(phiold[i][jp] - phiold[i][jm])/dy;
                
                Theta = atan2(phiy,phix);
                
                Epsilon[i][j] = epsilonb + epsilonb*delta*cos(beta*(Theta-Theta0));
                Epsilon_theta[i][j] = -epsilonb*beta*delta*sin(beta*(Theta-Theta0));
            }
        }
        
        // Step 3: Update phase-field parameter and temperature
        for (int i{}; i < Nx; ++i)
        {
            ip = (i == (Nx-1) ? 0 : i+1);
            im = (i == 0 ? (Nx-1) : i-1);
            
            for (int j{}; j < Ny; ++j)
            {
                jp = (j == (Ny-1) ? 0 : j+1);
                jm = (j == 0 ? (Ny-1) : j-1);
                
                // Update phase-field parameter
                Laplacian = (phiold[ip][j] - 2.0*phiold[i][j] + phiold[im][j])/(dx*dx) + (phiold[i][jp] - 2.0*phiold[i][j] + phiold[i][jm])/(dy*dy);
                m = alpha*atan(Gamma*(Teq-Told[i][j]))/PI;
                
                term0 = 0.5*(Epsilon[i][jp]*Epsilon_theta[i][jp]*PhiX(i,jp) - Epsilon[i][jm]*Epsilon_theta[i][jm]*PhiX(i,jm))/dy;
                term1 = 0.5*(Epsilon[ip][j]*Epsilon_theta[ip][j]*PhiY(ip,j) - Epsilon[im][j]*Epsilon_theta[im][j]*PhiY(im,j))/dx;
                term2 = Epsilon[i][j]*Epsilon[i][j]*Laplacian;
                term3 = phiold[i][j]*(1.0-phiold[i][j])*(phiold[i][j]-0.5+m);
                
                phi[i][j] += dt*(term0 - term1 + term2 + term3)/tau;
                
                // Update temperature
                Laplacian = (Told[ip][j] - 2.0*Told[i][j] + Told[im][j])/(dx*dx) + (Told[i][jp] - 2.0*Told[i][j] + Told[i][jm])/(dy*dy);
                
                T[i][j] += dt*Laplacian + kappa*(phi[i][j]-phiold[i][j]);
            }
        }
        
        // Step 4: Write file
        if (TimeStep % SaveInterval == 0)
        {
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d.tec",TimeStep/SaveInterval);
            
            std::ofstream Output(s,std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"phi\", \"T\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    Output << i*dx << "\t" << j*dy << "\t" << phi[i][j] << "\t" << T[i][j] << std::endl;
            
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
            Output << "Dendritic growth" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << i*dx << "\t" << j*dy << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << phi[i][j] << std::endl;
            
            Output.close();
            #endif
        }
        
        TimeStep++;
        
        t += dt;
        
        CheckNaN(phi,Nx,Ny);
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
