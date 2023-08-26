//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial implementation of p-adaption.
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
// make TARGET=p-Adaption.cpp
// make run
//________________________________________________________________________________________________

#include "General.h"

using namespace std::chrono;

double t = 1.0, frequency = 0.5;

//________________________________________________________________________________________________
// Non-Premixed flame
//________________________________________________________________________________________________
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

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
bool IsAdaptionNeeded ( const int i, 
                        const int j, 
                        const double xi, 
                        const double yi, 
                        const double dx, 
                        const double dy )
{
    double x0 = xi+i*dx, y0 = yi+j*dy, x1 = x0+dx, y1 = y0+dy;
    
    double Weight0 = Absolute((tanh(ExactFunction(x1,0.5*(y0+y1))/0.01)-tanh(ExactFunction(x0,0.5*(y0+y1))/0.01))/(x1-x0));
    double Weight1 = Absolute((tanh(ExactFunction(0.5*(x0+x1),y1)/0.01)-tanh(ExactFunction(0.5*(x0+x1),y0)/0.01))/(y1-y0));
    
    double CellWeight = Maximum(Weight0,Weight1);
    
    return (CellWeight > 0.2 ? true : false);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main()
{
    std::cout.flags(std::ios::dec);
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    const int Nx = 150, Ny = 30, MaxLevel = 1;
    const double xi = 0.0, xo = 10.0, yi = -1.0, yo = 1.0;
    const double dx = (xo-xi)/(Nx-1.0), dy = (yo-yi)/(Ny-1.0);
    
    char *s;
    
    // Allocate memory
    Allocate(s,200);
    
    int TimeSteps = 100;
    
    for (int timestep{}; timestep < TimeSteps; ++timestep)
    {
        t = timestep/(frequency*TimeSteps);
        
        // Write file for p-adaption
        sprintf(s,"Output/p-Adaption/Field-%d.tec",timestep);
        
        std::ofstream TecplotWrite(s,std::ios::out);
        TecplotWrite.flags(std::ios::dec);
        TecplotWrite.precision(8);
        
        if ( !TecplotWrite )
            ErrorMessage("Output file couldnot be opened.");
        
        TecplotWrite << "TITLE = \"Mesh refinement\"\nVariables = \"X\",\"Y\",\"Degree\"" << std::endl;
        TecplotWrite << "Zone N = " << 4*(Nx-1)*(Ny-1) << ", E = " << (Nx-1)*(Ny-1) << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        for (int i{}; i < (Nx-1); ++i)
        {
            for (int j{}; j < (Ny-1); ++j)
            {
                double x0 = xi+i*dx;
                double y0 = yi+j*dy;
                double value = (IsAdaptionNeeded(i,j,xi,yi,dx,dy) ? 2.0 : 1.0);
                
                TecplotWrite << x0 << "\t" << y0 << "\t" << value << std::endl;
                TecplotWrite << x0+dx << "\t" << y0 << "\t" << value << std::endl;
                TecplotWrite << x0+dx << "\t" << y0+dy << "\t" << value << std::endl;
                TecplotWrite << x0 << "\t" << y0+dy << "\t" << value << std::endl;
            }
        }
        
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                TecplotWrite << ((Ny-1)*i+j)*4 + 1 << "\t" << ((Ny-1)*i+j)*4 + 2 << "\t" << ((Ny-1)*i+j)*4 + 3 << "\t" << ((Ny-1)*i+j)*4 + 4 << std::endl;
        
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
    
    Deallocate(s,200);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
