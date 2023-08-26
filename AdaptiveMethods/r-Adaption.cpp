//_______________________________________________________________________________
//_______________________________________________________________________________
// Serial implementation of adaptive mesh redistribution.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 27.01.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=r-Adaption.cpp
// make run
//_______________________________________________________________________________

#include "General.h"

using namespace std::chrono;

double t = 1.0, frequency = 0.5;

//_______________________________________________________________________________
// Lambda : N
// Q      : N x N
//_______________________________________________________________________________
void EigenDD ( double *Lambda,
               double **Q,
               const int N )
{
    for (int j{}; j < N; ++j)
        Lambda[j] = -4.0*sin((j+1.0)*PI/(2.0*(N+1.0)))*sin((j+1.0)*PI/(2.0*(N+1.0)));
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Q[i][j] = sqrt(2.0/(N+1.0))*sin(PI*(i+1.0)*(j+1.0)/(N+1.0));
}

//_______________________________________________________________________________
// x      : N
// RHS    : N-2
// Lambda : N-2
// Q      : (N-2) x (N-2)
//_______________________________________________________________________________
void HelmholtzSolver ( double *x, 
                       double *rhs, 
                       double *Lambda, 
                       double **Q,
                       double *temporary,  
                       const int N,
                       const double lambda )
{
    // Modify right-hand side using the given boundary conditions
    rhs[0] -= x[0];
    rhs[N-3] -= x[N-1];
    
    // Modify right-hand side
    for (int i{}; i < (N-2); ++i)
    {
        temporary[i] = 0.0;
        
        for (int k{}; k < (N-2); ++k)
            temporary[i] += rhs[k]*Q[k][i];
    }
    
    for (int i{}; i < (N-2); ++i)
        rhs[i] = temporary[i];
    
    // Calculate modified solution
    for (int i{1}; i < (N-1); ++i)
        rhs[i-1] = (Absolute(Lambda[i-1]-lambda) < 1.0E-16 ? 0.0 : rhs[i-1]/(Lambda[i-1]-lambda));
    
    // Calculate the final solution
    for (int i{}; i < (N-2); ++i)
    {
        temporary[i] = 0.0;
        
        for (int k{}; k < (N-2); ++k)
            temporary[i] += rhs[k]*Q[i][k];
    }
    
    for (int i{1}; i < (N-1); ++i)
        x[i] = temporary[i-1];
}

//_______________________________________________________________________________
// X       : Nx x Ny
// RHS     : (Nx-2) x (Ny-2)
// Lambdax : Nx-2
// Lambday : Ny-2
// Qx      : (Nx-2) x (Nx-2)
// Qy      : (Ny-2) x (Ny-2)
//_______________________________________________________________________________
void HelmholtzSolver ( double **X, 
                       double **RHS, 
                       double *Lambdax, 
                       double **Qx, 
                       double *Lambday, 
                       double **Qy, 
                       double **Temporary, 
                       const int Nx, 
                       const int Ny, 
                       const double lambda )
{
    // Modify right-hand side using the given boundary conditions
    for (int i{1}; i < (Nx-1); ++i)
    {
        RHS[i-1][0] -= X[i][0];
        RHS[i-1][Ny-3] -= X[i][Ny-1];
    }
    
    for (int j{1}; j < (Ny-1); ++j)
    {
        RHS[0][j-1] -= X[0][j];
        RHS[Nx-3][j-1] -= X[Nx-1][j];
    }
    
    // Modify right-hand side
    for (int i{}; i < (Nx-2); ++i)
    {
        for (int j{}; j < (Ny-2); ++j)
        {
            Temporary[i][j] = 0.0;
            
            for (int k{}; k < (Ny-2); ++k)
                Temporary[i][j] += RHS[i][k]*Qy[k][j];
        }
    }
    
    for (int i{}; i < (Nx-2); ++i)
    {
        for (int j{}; j < (Ny-2); ++j)
        {
            RHS[i][j] = 0.0;
            
            for (int k{}; k < (Nx-2); ++k)
                RHS[i][j] += Qx[k][i]*Temporary[k][j];
        }
    }
    
    // Calculate modified solution
    for (int i{1}; i < (Nx-1); ++i)
        for (int j{1}; j < (Ny-1); ++j)
            RHS[i-1][j-1] = (Absolute(Lambdax[i-1]+Lambday[j-1]-lambda) < 1.0E-16 ? 0.0 : RHS[i-1][j-1]/(Lambdax[i-1]+Lambday[j-1]-lambda));
    
    // Calculate the final solution
    for (int i{}; i < (Nx-2); ++i)
    {
        for (int j{}; j < (Ny-2); ++j)
        {
            Temporary[i][j] = 0.0;
            
            for (int k{}; k < (Ny-2); ++k)
                Temporary[i][j] += RHS[i][k]*Qy[j][k];
        }
    }
    
    for (int i{1}; i < (Nx-1); ++i)
    {
        for (int j{1}; j < (Ny-1); ++j)
        {
            X[i][j] = 0.0;
            
            for (int k{}; k < (Nx-2); ++k)
                X[i][j] += Qx[i-1][k]*Temporary[k][j-1];
        }
    }
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void RightTerm ( double *rhs, 
                 double *omega, 
                 double *x, 
                 const double OmegaMax, 
                 const double lambda, 
                 const int N )
{
    for (int i{1}; i < (N-1); ++i)
    {
        rhs[i-1] = (x[i+1]-2.0*x[i]+x[i-1]) - lambda*x[i];
        
        rhs[i-1] += (i == 1 ? (x[i]-x[i-1]) : 0.5*(x[i]-x[i-2]))*omega[i-1]/OmegaMax;
        rhs[i-1] -= (i == (N-2) ? (x[i+1]-x[i]) : 0.5*(x[i+2]-x[i]))*omega[i+1]/OmegaMax;
    }
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void RightTerm ( double **RHS, 
                 double **Omega, 
                 double **X, 
                 const double OmegaMax, 
                 const double lambda, 
                 const int Nx, 
                 const int Ny )
{
    for (int i{1}; i < (Nx-1); ++i)
    {
        for (int j{1}; j < (Ny-1); ++j)
        {
            RHS[i-1][j-1] = (X[i+1][j]-2.0*X[i][j]+X[i-1][j]) + (X[i][j+1]-2.0*X[i][j]+X[i][j-1]) - lambda*X[i][j];
            
            RHS[i-1][j-1] += (i == 1 ? (X[i][j]-X[i-1][j]) : 0.5*(X[i][j]-X[i-2][j]))*Omega[i-1][j]/OmegaMax;
            RHS[i-1][j-1] -= (i == (Nx-2) ? (X[i+1][j]-X[i][j]) : 0.5*(X[i+2][j]-X[i][j]))*Omega[i+1][j]/OmegaMax;
            RHS[i-1][j-1] += (j == 1 ? (X[i][j]-X[i][j-1]) : 0.5*(X[i][j]-X[i][j-2]))*Omega[i][j-1]/OmegaMax;
            RHS[i-1][j-1] -= (j == (Ny-2) ? (X[i][j+1]-X[i][j]) : 0.5*(X[i][j+2]-X[i][j]))*Omega[i][j+1]/OmegaMax;
        }
    }
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void MonitorFunction ( double **Phi,
                       double **Omega, 
                       double **Omeganew, 
                       double &OmegaMax, 
                       const int Nx, 
                       const int Ny )
{
    const double beta = 20.0;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            Omega[i][j] = 1.0;
            Omega[i][j] += beta*beta*(i == 0 ? (Phi[1][j]-Phi[0][j])*(Phi[1][j]-Phi[0][j]) : (i == (Nx-1) ? (Phi[Nx-1][j]-Phi[Nx-2][j])*(Phi[Nx-1][j]-Phi[Nx-2][j]) : 0.25*(Phi[i+1][j]-Phi[i-1][j])*(Phi[i+1][j]-Phi[i-1][j])));
            Omega[i][j] += beta*beta*(j == 0 ? (Phi[i][1]-Phi[i][0])*(Phi[i][1]-Phi[i][0]) : (j == (Ny-1) ? (Phi[i][Ny-1]-Phi[i][Ny-2])*(Phi[i][Ny-1]-Phi[i][Ny-2]) : 0.25*(Phi[i][j+1]-Phi[i][j-1])*(Phi[i][j+1]-Phi[i][j-1])));
            
            Omega[i][j] = sqrt(Omega[i][j]);
        }
    }
    
    // Smoothing
    for (int smoothing{}; smoothing < 4; ++smoothing)
    {
        for (int i{1}; i < (Nx-1); ++i)
        {
            Omeganew[i][0] = 0.25*(2.0*Omega[i][0] + Omega[i+1][0] + Omega[i-1][0]);
            Omeganew[i][Ny-1] = 0.25*(2.0*Omega[i][Ny-1] + Omega[i+1][Ny-1] + Omega[i-1][Ny-1]);
        }
        
        for (int j{1}; j < (Ny-1); ++j)
        {
            Omeganew[0][j] = 0.25*(2.0*Omega[0][j] + Omega[0][j+1] + Omega[0][j-1]);
            Omeganew[Nx-1][j] = 0.25*(2.0*Omega[Nx-1][j] + Omega[Nx-1][j+1] + Omega[Nx-1][j-1]);
        }
        
        for (int i{1}; i < (Nx-1); ++i)
            for (int j{1}; j < (Ny-1); ++j)
                Omeganew[i][j] = 0.0625*(4.0*Omega[i][j] + 2.0*(Omega[i+1][j]+Omega[i-1][j]+Omega[i][j+1]+Omega[i][j-1]) + (Omega[i+1][j+1]+Omega[i-1][j+1]+Omega[i+1][j-1]+Omega[i-1][j-1]));
        
        for (int i{1}; i < (Nx-1); ++i)
        {
            Omega[i][0] = Omeganew[i][0];
            Omega[i][Ny-1] = Omeganew[i][Ny-1];
        }
        
        for (int j{1}; j < (Ny-1); ++j)
        {
            Omega[0][j] = Omeganew[0][j];
            Omega[Nx-1][j] = Omeganew[Nx-1][j];
        }
        
        for (int i{1}; i < (Nx-1); ++i)
            for (int j{1}; j < (Ny-1); ++j)
                Omega[i][j] = Omeganew[i][j];
    }
    
    // Determine maximum of the monitor function
    OmegaMax = 0.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            OmegaMax = Maximum(OmegaMax,Omega[i][j]);
}

//_______________________________________________________________________________
// Non-Premixed flame
//_______________________________________________________________________________
double ExactFunction ( double x, 
                       double y )
{
    // Input parameters.
    double Pe = 40.0, VelocityRatio = 1.0, alpha = 0.1;
    double xi = (3.2/7), yi = 3.2*VelocityRatio, xib = 0, yib = 0.5*yi;
    
    double z1 = 0.0, z2 = 0.0, a, b, term;
    
    for (int n{}; n < 20; ++n)
    {
        a = 0.25*Pe*Pe + n*n*PI*PI;
        b = 2.0*PI*frequency*Pe;
        
        SquareRoot(a,b);
        
        a = 0.5*Pe - a;
        b *= -1.0;
        
        if (n == 0)
            term = ((1.0-alpha)*xib-alpha*yib)*exp(a*x)*cos(b*x+2.0*PI*frequency*t)*cos(n*PI*y);
        else
            term = (-2.0*sin(n*PI*alpha)*(xib+yib))*Pe*exp(a*x)*cos(n*PI*y)*((Pe-a)*cos(b*x+2.0*PI*frequency*t) - b*sin(b*x+2.0*PI*frequency*t))/(n*PI*((Pe-a)*(Pe-a)+b*b));
        
        z2 += term;
    }
    
    for (int n{}; n < 20; ++n)
    {
        a = 0.5*Pe - sqrt(0.25*Pe*Pe + n*n*PI*PI);
        
        if (n == 0)
            term = ((1.0-alpha)*xi-alpha*yi)*exp(a*x)*cos(n*PI*y);
        else
            term = ((-2.0*sin(n*PI*alpha)*(xi+yi))/(n*PI*(1.0-(a/Pe))))*exp(a*x)*cos(n*PI*y);
        
        z1 += term;
    }
    
    return (z1+z2);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main()
{
    std::cout.flags( std::ios::dec | std::ios::fixed );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    const int Nx = 150, Ny = 30;
    const double xi = 0.0, xo = 10.0, yi = -1.0, yo = 1.0;
    double dt, OmegaMax, lambda;
    
    double **Phi, ***Xn, **Omega, **Omeganew, *omega;
    double ***X, ***Xnew, *x, **RHS, *rhs, **Temporary, *temporary;
    double *Lambdax, **Qx, *Lambday, **Qy;
    
    char *s;
    
    // Allocate memory
    Allocate(Phi,Nx,Ny);
    
    Allocate(Xn,2,Nx,Ny);
    
    Allocate(Omega,Nx,Ny);
    Allocate(Omeganew,Nx,Ny);
    Allocate(omega,Maximum(Nx,Ny));
    
    Allocate(X,2,Nx,Ny);
    Allocate(Xnew,2,Nx,Ny);
    Allocate(x,Maximum(Nx,Ny));
    Allocate(RHS,Nx-2,Ny-2);
    Allocate(rhs,Maximum(Nx,Ny)-2);
    Allocate(Temporary,Nx-2,Ny-2);
    Allocate(temporary,Maximum(Nx,Ny)-2);
    
    Allocate(Lambdax,Nx-2);
    Allocate(Qx,Nx-2,Nx-2);
    Allocate(Lambday,Ny-2);
    Allocate(Qy,Ny-2,Ny-2);
    
    Allocate(s,200);
    
    // Compute necessary matrix decompositions (eigenvalues and eigenvectors)
    EigenDD(Lambdax,Qx,Nx-2);
    EigenDD(Lambday,Qy,Ny-2);
    
    // Main algorithm starts here
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            Xn[0][i][j] = xi + i*(xo-xi)/(Nx-1.0);
            Xn[1][i][j] = yi + j*(yo-yi)/(Ny-1.0);
            
            X[0][i][j] = Xn[0][i][j];
            X[1][i][j] = Xn[1][i][j];
        }
    }
    
    double resolution = Minimum((xo-xi)/(Nx-1.0),(yo-yi)/(Ny-1.0))*0.05;
    
    int TimeSteps = 100;
    
    for (int timestep{}; timestep < TimeSteps; ++timestep)
    {
        t = timestep/(frequency*TimeSteps);
        
        int updatesteps = 0;
        
        bool flag = true;
        
        while ( flag )
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < Ny; ++j)
                {
                    Xnew[0][i][j] = X[0][i][j];
                    Xnew[1][i][j] = X[1][i][j];
                }
            }
            
            // Update values on the present coordinates
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    Phi[i][j] = tanh(ExactFunction(X[0][i][j],X[1][i][j])/0.01);
            
            // Monitor function
            MonitorFunction(Phi,Omega,Omeganew,OmegaMax,Nx,Ny);
            
            dt = 0.5/OmegaMax;
            lambda = 1.0/(dt*OmegaMax);
            
            // Adapt mesh on the top edge
            for (int i{}; i < Nx; ++i)
            {
                x[i] = X[0][i][Ny-1];
                omega[i] = Omega[i][Ny-1];
            }
            
            RightTerm(rhs,omega,x,OmegaMax,lambda,Nx);
            
            HelmholtzSolver(x,rhs,Lambdax,Qx,temporary,Nx,lambda);
            
            for (int i{1}; i < (Nx-1); ++i)
                X[0][i][Ny-1] = x[i];
            
            // Adapt mesh on the bottom edge
            for (int i{}; i < Nx; ++i)
            {
                x[i] = X[0][i][0];
                omega[i] = Omega[i][0];
            }
            
            RightTerm(rhs,omega,x,OmegaMax,lambda,Nx);
            
            HelmholtzSolver(x,rhs,Lambdax,Qx,temporary,Nx,lambda);
            
            for (int i{1}; i < (Nx-1); ++i)
                X[0][i][0] = x[i];
            
            // Adapt mesh on the left edge
            for (int j{}; j < Ny; ++j)
            {
                x[j] = X[1][0][j];
                omega[j] = Omega[0][j];
            }
            
            RightTerm(rhs,omega,x,OmegaMax,lambda,Ny);
            
            HelmholtzSolver(x,rhs,Lambday,Qy,temporary,Ny,lambda);
            
            for (int j{1}; j < (Ny-1); ++j)
                X[1][0][j] = x[j];
            
            // Adapt mesh on the right edge
            for (int j{}; j < Ny; ++j)
            {
                x[j] = X[1][Nx-1][j];
                omega[j] = Omega[Nx-1][j];
            }
            
            RightTerm(rhs,omega,x,OmegaMax,lambda,Ny);
            
            HelmholtzSolver(x,rhs,Lambday,Qy,temporary,Ny,lambda);
            
            for (int j{1}; j < (Ny-1); ++j)
                X[1][Nx-1][j] = x[j];
            
            // Adapt interior mesh
            RightTerm(RHS,Omega,X[0],OmegaMax,lambda,Nx,Ny);
            
            HelmholtzSolver(X[0],RHS,Lambdax,Qx,Lambday,Qy,Temporary,Nx,Ny,lambda);
            
            RightTerm(RHS,Omega,X[1],OmegaMax,lambda,Nx,Ny);
            
            HelmholtzSolver(X[1],RHS,Lambdax,Qx,Lambday,Qy,Temporary,Nx,Ny,lambda);
            
            double error = 0.0;
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    error = Maximum(error,Absolute(X[0][i][j]-Xnew[0][i][j]),Absolute(X[1][i][j]-Xnew[1][i][j]));
            
            std::cout << "Iteration number : " << updatesteps << ", Error = " << (error/resolution-1.0) << std::endl;
            
            if (error < resolution)
                flag = false;
            else
                updatesteps++;
        }
        
        std::cout << "Number of iterations to converge = " << updatesteps << std::endl;
        std::cout << std::endl << "Initial condition has been successfully adapted!" << std::endl;
        
        // Write file
        sprintf(s,"Output/Redistribution/Field-%d.tec",timestep);
        
        std::ofstream TecplotWrite(s,std::ios::out);
        TecplotWrite.flags(std::ios::dec);
        TecplotWrite.precision(8);
        
        if ( !TecplotWrite )
            ErrorMessage("Output file couldnot be opened!");
        
        TecplotWrite << "TITLE = \"Physical plane\"\nVariables = \"X\",\"Y\",\"Phi\",\"Omega\"" << std::endl;
        TecplotWrite << "Zone I = " << Ny << ", J = " << Nx << ", DATAPACKING=POINT" << std::endl;
        
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                TecplotWrite << X[0][i][j] << "\t" << X[1][i][j] << "\t" << Phi[i][j] << "\t" << Omega[i][j] << std::endl;
        
        TecplotWrite.close();
        
        // Write file for function plot
        sprintf(s,"Output/Flame/Field-%d.tec",timestep);
        
        const int Nxh = 500, Nyh = 100;
        
        std::ofstream TecplotWriteFlame(s,std::ios::out);
        TecplotWriteFlame.flags(std::ios::dec);
        TecplotWriteFlame.precision(8);
        
        if ( !TecplotWriteFlame )
            ErrorMessage("Output file couldnot be opened!");
        
        TecplotWriteFlame << "TITLE = \"Physical plane\"\nVariables = \"X\",\"Y\",\"Phi\"" << std::endl;
        TecplotWriteFlame << "Zone I = " << Nyh << ", J = " << Nxh << ", DATAPACKING=POINT" << std::endl;
        
        for (int i{}; i < Nxh; ++i)
        {
            for (int j{}; j < Nyh; ++j)
            {
                double x = xi + i*(xo-xi)/(Nxh-1.0);
                double y = yi + j*(yo-yi)/(Nyh-1.0);
                
                TecplotWriteFlame << x << "\t" << y << "\t" << ExactFunction(x,y) << std::endl;
            }
        }
        
        TecplotWriteFlame.close();
        
        std::cout << "Time step = " << timestep << std::endl;
    }
    
    // Deallocate memory
    Deallocate(Phi,Nx,Ny);
    
    Deallocate(Xn,2,Nx,Ny);
    
    Deallocate(Omega,Nx,Ny);
    Deallocate(Omeganew,Nx,Ny);
    Deallocate(omega,Maximum(Nx,Ny));
    
    Deallocate(X,2,Nx,Ny);
    Deallocate(Xnew,2,Nx,Ny);
    Deallocate(x,Maximum(Nx,Ny));
    Deallocate(RHS,Nx-2,Ny-2);
    Deallocate(rhs,Maximum(Nx,Ny)-2);
    Deallocate(Temporary,Nx-2,Ny-2);
    Deallocate(temporary,Maximum(Nx,Ny)-2);
    
    Deallocate(Lambdax,Nx-2);
    Deallocate(Qx,Nx-2,Nx-2);
    Deallocate(Lambday,Ny-2);
    Deallocate(Qy,Ny-2,Ny-2);
    
    Deallocate(s,200);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
