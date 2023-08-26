//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial Program to simulate viscous fingering using PseudoSpectral method.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 23.01.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=ViscousFingering.cpp
// make run
//________________________________________________________________________________________________

#include "General.h"
#include "Fourier_Quadrature.h"

using namespace std::chrono;

#define TECPLOT
//#define VISIT

const int Nx = 512, Ny = 512, SaveInterval = 200;
const double Pe = 3000.0, R = 2.5, ratio = 1.0;
const double Lx = Pe, Ly = Pe*ratio, dx = Lx/Nx, dy = Ly/Ny;
const double Cr = 0.0, Cl = 1.0, epsilon = 0.01, xr = Lx*0.5, xl = Lx*0.125;
const double xmin = 0.0, xmax = xmin + Lx, ymin = 0.0, ymax = ymin + Ly;
const double dt = 0.05, Tf = 1000.0;

double **C, **Cbar, **Cx, **Cy, t = 0.0;
double **Si, **J, **Jbar, **N, **SiY, **SiX, **Jm;
char *s;

double x, y, cx, cy, K;

int kx, ky, TimeStep = 0;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Allocate ()
{
    Allocate(C,2*Nx,2*Ny);
    Allocate(Cbar,2*Nx,2*Ny);
    Allocate(Si,2*Nx,2*Ny);
    Allocate(J,2*Nx,2*Ny);
    Allocate(Jbar,2*Nx,2*Ny);
    Allocate(N,2*Nx,2*Ny);
    
    Allocate(Cx,Nx,Ny);
    Allocate(Cy,Nx,Ny);
    Allocate(SiX,Nx,Ny);
    Allocate(SiY,Nx,Ny);
    
    Allocate(Jm,Nx,2*Ny);
    Allocate(s,50);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(C,2*Nx,2*Ny);
    Deallocate(Cbar,2*Nx,2*Ny);
    Deallocate(Si,2*Nx,2*Ny);
    Deallocate(J,2*Nx,2*Ny);
    Deallocate(Jbar,2*Nx,2*Ny);
    Deallocate(N,2*Nx,2*Ny);
    
    Deallocate(Cx,Nx,Ny);
    Deallocate(Cy,Nx,Ny);
    Deallocate(SiX,Nx,Ny);
    Deallocate(SiY,Nx,Ny);
    
    Deallocate(Jm,Nx,2*Ny);
    Deallocate(s,50);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    // Initialize concentration field.
    double random, wn, delta = (Nx < 300 ? 20000.0 : 30000.0);
    
    for (int n{}; n < 5; ++n)
    {
        random = RandomNumber();
        wn = (n+1.0)*10.0;
        
        for (int j{}; j < Ny; ++j)
        {
            y = j*dy;
            
            for (int i{}; i < Nx; ++i)
                C[i][j] += 0.004*(random-0.5)*cos(wn*PI*(y-ymin)/Ly);
        }
    }
    
    for (int i{}; i < Nx; ++i)
    {
        x = i*dx + Lx/2.5;
        
        for (int j{}; j < Ny; ++j)
        {
            C[i][j] = (-0.71 + (x-xmin)/Lx - 4.0/sqrt(delta) + C[i][j])*sqrt(delta);
            C[i][j] = 0.5*(1.0-Erf(C[i][j]));
        }
    }
    
    // Modification over initial concentration field to maintain periodicity
    for (int i{}; i < Nx; ++i)
    {
        x = i*dx;
        
        for (int j{}; j < Ny; ++j)
            C[i][j] *= 0.5*(1.0+tanh((x-xl)/(Lx*0.005)));
    }
    
    // Initialize streamfunction
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            Si[i][j] = 0.0;
    
    // Calculate SiX and SiY
    FirstDerivativeFourierX(SiX,Si,Nx,Ny,Lx);
    FirstDerivativeFourierY(SiY,Si,Nx,Ny,Ly);
    
    // Calculate Cx and Cy
    FirstDerivativeFourierX(Cx,C,Nx,Ny,Lx);
    FirstDerivativeFourierY(Cy,C,Nx,Ny,Ly);
    
    // Calculate J
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            J[i][j] = SiY[i][j]*Cx[i][j] - SiX[i][j]*Cy[i][j];
    
    FourierTransformXY(J,1,Nx,Ny);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < 2*Ny; ++j)
            Jm[i][j] = J[i][j];
}

//________________________________________________________________________________________________
// Wave numbers
//________________________________________________________________________________________________
void WaveNumbers ( int &kx, 
                   int &ky,  
                   const int i, 
                   const int j )
{
    // Shifted location of the waves on spectral space
    // 0 1 2 3 ...... (N/2-2) (N/2-1) (-N/2) -(N/2-1) -(N/2-2) ....... -3 -2 -1
    
    kx = (i < Nx/2 ? i : i-Nx);
    ky = (j < Ny/2 ? j : j-Ny);
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
    
    CreateFFTPlan(Maximum(Nx,Ny));
    Allocate();
    
    Initialize();
    
    // Update concentration
    while (t < Tf)
    {
        std::cout << "Time step = " << TimeStep << ", Current time = " << t << std::endl;
        
        // Step 1: Provisional values of concentration (Cbar)
        FourierTransformXY(C,1,Nx,Ny);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                // Shifted location of the waves on spectral space
                // 0 1 2 3 ...... (N/2-2) (N/2-1) (-N/2) -(N/2-1) -(N/2-2) ....... -3 -2 -1
                
                kx = (i < Nx/2 ? i : i-Nx);
                ky = (j < Ny/2 ? j : j-Ny);
                
                cx = 2.0*PI*kx/Lx;
                cy = 2.0*PI*ky/Ly;
                
                K = cx*cx + cy*cy;
                
                Cbar[i][2*j] = C[i][2*j] + dt*(1.5*J[i][2*j] - Jm[i][2*j]);
                Cbar[i][2*j+1] = C[i][2*j+1] + dt*(1.5*J[i][2*j+1] - Jm[i][2*j+1]);
                
                // Partial update of spectral coefficients of C
                C[i][2*j] -= 0.5*dt*K*(Cbar[i][2*j] + C[i][2*j]);
                C[i][2*j+1] -= 0.5*dt*K*(Cbar[i][2*j+1] + C[i][2*j+1]);
                
                // Partial spectral coefficients of rhs of Poisson equation for Si
                K = -R*(((kx == 0) && (ky == 0)) ? 0.0 : 1.0/K)*cy;
                
                Si[i][2*j]   = -K*Cbar[i][2*j+1];
                Si[i][2*j+1] =  K*Cbar[i][2*j];
            }
        }
        
        // Cbar spectral to physical
        FourierTransformXY(Cbar,-1,Nx,Ny);
        
        // Cx and Cy
        FirstDerivativeFourierX(Cx,Cbar,Nx,Ny,Lx);
        FirstDerivativeFourierY(Cy,Cbar,Nx,Ny,Ly);
        
        // N
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                N[i][j] = SiX[i][j]*Cx[i][j] + SiY[i][j]*Cy[i][j];
        
        FourierTransformXY(N,1,Nx,Ny);
        
        // Step 2: Spectral coefficients of rhs of Poisson equation for Si
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                // Shifted location of the waves on spectral space
                // 0 1 2 3 ...... (N/2-2) (N/2-1) (-N/2) -(N/2-1) -(N/2-2) ....... -3 -2 -1
                
                kx = (i < Nx/2 ? i : i-Nx);
                ky = (j < Ny/2 ? j : j-Ny);
                
                cx = 2.0*PI*kx/Lx;
                cy = 2.0*PI*ky/Ly;
                
                K = (((kx == 0) && (ky == 0)) ? 0.0 : 1.0/(cx*cx + cy*cy));
                
                Si[i][2*j]   -= R*K*N[i][2*j];
                Si[i][2*j+1] -= R*K*N[i][2*j+1];
            }
        }
        
        // Physical values of Si
        FourierTransformXY(Si,-1,Nx,Ny);
        
        // SiY and SiX
        FirstDerivativeFourierX(SiX,Si,Nx,Ny,Lx);
        FirstDerivativeFourierY(SiY,Si,Nx,Ny,Ly);
        
        // Jbar
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                Jbar[i][j] = SiY[i][j]*Cx[i][j] - SiX[i][j]*Cy[i][j];
        
        // Spectral coefficients of Jbar
        FourierTransformXY(Jbar,1,Nx,Ny);
        
        // Step 3: Full update of spectral coefficients of C
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                C[i][2*j] -= 0.5*dt*(Jbar[i][2*j] + J[i][2*j]);
                C[i][2*j+1] -= 0.5*dt*(Jbar[i][2*j+1] + J[i][2*j+1]);
            }
        }
        
        // C spectral to physical
        FourierTransformXY(C,-1,Nx,Ny);
        
        // Calculate Cx and Cy
        FirstDerivativeFourierX(Cx,C,Nx,Ny,Lx);
        FirstDerivativeFourierY(Cy,C,Nx,Ny,Ly);
        
        // Jm = J
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < 2*Ny; ++j)
                Jm[i][j] = J[i][j];
        
        // Calculate J
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                J[i][j] = SiY[i][j]*Cx[i][j] - SiX[i][j]*Cy[i][j];
        
        // Spectral coefficients of J
        FourierTransformXY(J,1,Nx,Ny);
        
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
            
            int Nx0 = ceil(0.135*Nx);
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"C\", \"U\", \"V\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << (Nx-Nx0) << std::endl;
            
            for (int i{Nx0}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    Output << i*dx/Pe << "\t" << j*dy/Pe << "\t" << C[i][j] << "\t" << SiY[i][j] << "\t" << -SiX[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%04d.vtk",TimeStep/SaveInterval);
            
            std::ofstream Output(s,std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            int Nx0 = ceil(0.135*Nx);
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Saffman-Taylor instability" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << (Nx-Nx0) << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << (Nx-Nx0)*Ny << " FLOAT" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{Nx0}; i < Nx; ++i)
                    Output << i*dx/Pe << "\t" << j*dy/Pe << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << (Nx-Nx0)*Ny << std::endl;
            Output << "SCALARS Concentration float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{Nx0}; i < Nx; ++i)
                    Output << C[i][j] << std::endl;
            
            Output.close();
            #endif
        }
        
        TimeStep++;
        
        t += dt;
        
        CheckNaN(C,Nx,Ny);
    }
    
    Deallocate();
    DestroyFFTPlan(Maximum(Nx,Ny));
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
