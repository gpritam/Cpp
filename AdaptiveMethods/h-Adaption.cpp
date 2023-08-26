//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial implementation of adaptive mesh refinement.
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
// make TARGET=h-Adaption.cpp
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
// vertices : 4x2
//________________________________________________________________________________________________
struct Node
{
    int level;
    
    double **vertices;
    
    struct Node *NextNode0, *NextNode1, *NextNode2, *NextNode3;
    struct Node *PreviousNode;
};

//________________________________________________________________________________________________
// Allocate the memory
//________________________________________________________________________________________________
void Allocate ( struct Node **B,
                const int i,
                const int j,
                const double dx,
                const double dy,
                const double xi, 
                const double yi )
{
    struct Node *Temporary;
    
    if ( (Temporary = (struct Node *)malloc(sizeof(struct Node)) ) == nullptr )
        ErrorMessage("Error: Memory cannot be allocated!");
    
    Temporary->level = 0;
    
    Allocate(Temporary->vertices,4,2);
    
    Temporary->vertices[0][0] = xi + i*dx;
    Temporary->vertices[0][1] = yi + j*dy;
    
    Temporary->vertices[1][0] = Temporary->vertices[0][0] + dx;
    Temporary->vertices[1][1] = Temporary->vertices[0][1];
    
    Temporary->vertices[2][0] = Temporary->vertices[0][0] + dx;
    Temporary->vertices[2][1] = Temporary->vertices[0][1] + dy;
    
    Temporary->vertices[3][0] = Temporary->vertices[0][0];
    Temporary->vertices[3][1] = Temporary->vertices[0][1] + dy;
    
    Temporary->NextNode0 = nullptr;
    Temporary->NextNode1 = nullptr;
    Temporary->NextNode2 = nullptr;
    Temporary->NextNode3 = nullptr;
    
    Temporary->PreviousNode = nullptr;
    
    *B = Temporary;
}

//________________________________________________________________________________________________
// Deallocate the memory
//________________________________________________________________________________________________
void Deallocate ( struct Node **B )
{
    struct Node *Temporary = *B;
    
    if (Temporary != nullptr)
    {
        Deallocate(&(Temporary->NextNode0));
        Deallocate(&(Temporary->NextNode1));
        Deallocate(&(Temporary->NextNode2));
        Deallocate(&(Temporary->NextNode3));
        
        if ((*B)->PreviousNode != nullptr)
        {
            Deallocate(Temporary->vertices,4,2);
            
            free(Temporary);
        }
        else
        {
            Temporary->NextNode0 = nullptr;
            Temporary->NextNode1 = nullptr;
            Temporary->NextNode2 = nullptr;
            Temporary->NextNode3 = nullptr;
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
bool IsRefinementNeeded ( struct Node *B,
                          const int MaxLevel )
{
    if (B->NextNode0 == nullptr)
    {
        double x0 = B->vertices[0][0], y0 = B->vertices[0][1], x1 = B->vertices[2][0], y1 = B->vertices[2][1];
        
        double Weight0 = Absolute((tanh(ExactFunction(x1,0.5*(y0+y1))/0.01)-tanh(ExactFunction(x0,0.5*(y0+y1))/0.01))/(x1-x0));
        double Weight1 = Absolute((tanh(ExactFunction(0.5*(x0+x1),y1)/0.01)-tanh(ExactFunction(0.5*(x0+x1),y0)/0.01))/(y1-y0));
        
        double CellWeight = Maximum(Weight0,Weight1);
        
        if ( (CellWeight > 0.2) && (B->level < MaxLevel) )
            return true;
        else
            return false;
    }
    else
        return false;
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void AdaptiveMeshRefinement ( struct Node **B,
                              const int MaxLevel )
{
    if (IsRefinementNeeded(*B,MaxLevel))
    {
        struct Node *Temporary;
        
        // Create first child node
        if ( (Temporary = (struct Node *)malloc(sizeof(struct Node)) ) == nullptr )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        Allocate(Temporary->vertices,4,2);
        
        Temporary->level = (*B)->level+1;
        
        Temporary->vertices[0][0] = (*B)->vertices[0][0];
        Temporary->vertices[0][1] = (*B)->vertices[0][1];
        
        Temporary->vertices[1][0] = 0.5*((*B)->vertices[0][0]+(*B)->vertices[1][0]);
        Temporary->vertices[1][1] = 0.5*((*B)->vertices[0][1]+(*B)->vertices[1][1]);
        
        Temporary->vertices[2][0] = 0.25*((*B)->vertices[0][0]+(*B)->vertices[1][0]+(*B)->vertices[2][0]+(*B)->vertices[3][0]);
        Temporary->vertices[2][1] = 0.25*((*B)->vertices[0][1]+(*B)->vertices[1][1]+(*B)->vertices[2][1]+(*B)->vertices[3][1]);
        
        Temporary->vertices[3][0] = 0.5*((*B)->vertices[0][0]+(*B)->vertices[3][0]);
        Temporary->vertices[3][1] = 0.5*((*B)->vertices[0][1]+(*B)->vertices[3][1]);
        
        Temporary->NextNode0 = nullptr;
        Temporary->NextNode1 = nullptr;
        Temporary->NextNode2 = nullptr;
        Temporary->NextNode3 = nullptr;
        
        Temporary->PreviousNode = *B;
        
        (*B)->NextNode0 = Temporary;
        
        // Create second child node
        if ( (Temporary = (struct Node *)malloc(sizeof(struct Node)) ) == nullptr )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        Allocate(Temporary->vertices,4,2);
        
        Temporary->level = (*B)->level+1;
        
        Temporary->vertices[0][0] = 0.5*((*B)->vertices[0][0]+(*B)->vertices[1][0]);
        Temporary->vertices[0][1] = 0.5*((*B)->vertices[0][1]+(*B)->vertices[1][1]);
        
        Temporary->vertices[1][0] = (*B)->vertices[1][0];
        Temporary->vertices[1][1] = (*B)->vertices[1][1];
        
        Temporary->vertices[2][0] = 0.5*((*B)->vertices[1][0]+(*B)->vertices[2][0]);
        Temporary->vertices[2][1] = 0.5*((*B)->vertices[1][1]+(*B)->vertices[2][1]);
        
        Temporary->vertices[3][0] = 0.25*((*B)->vertices[0][0]+(*B)->vertices[1][0]+(*B)->vertices[2][0]+(*B)->vertices[3][0]);
        Temporary->vertices[3][1] = 0.25*((*B)->vertices[0][1]+(*B)->vertices[1][1]+(*B)->vertices[2][1]+(*B)->vertices[3][1]);
        
        Temporary->NextNode0 = nullptr;
        Temporary->NextNode1 = nullptr;
        Temporary->NextNode2 = nullptr;
        Temporary->NextNode3 = nullptr;
        
        Temporary->PreviousNode = *B;
        
        (*B)->NextNode1 = Temporary;
        
        // Create third child node
        if ( (Temporary = (struct Node *)malloc(sizeof(struct Node)) ) == nullptr )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        Allocate(Temporary->vertices,4,2);
        
        Temporary->level = (*B)->level+1;
        
        Temporary->vertices[0][0] = 0.25*((*B)->vertices[0][0]+(*B)->vertices[1][0]+(*B)->vertices[2][0]+(*B)->vertices[3][0]);
        Temporary->vertices[0][1] = 0.25*((*B)->vertices[0][1]+(*B)->vertices[1][1]+(*B)->vertices[2][1]+(*B)->vertices[3][1]);
        
        Temporary->vertices[1][0] = 0.5*((*B)->vertices[1][0]+(*B)->vertices[2][0]);
        Temporary->vertices[1][1] = 0.5*((*B)->vertices[1][1]+(*B)->vertices[2][1]);
        
        Temporary->vertices[2][0] = (*B)->vertices[2][0];
        Temporary->vertices[2][1] = (*B)->vertices[2][1];
        
        Temporary->vertices[3][0] = 0.5*((*B)->vertices[3][0]+(*B)->vertices[2][0]);
        Temporary->vertices[3][1] = 0.5*((*B)->vertices[3][1]+(*B)->vertices[2][1]);
        
        Temporary->NextNode0 = nullptr;
        Temporary->NextNode1 = nullptr;
        Temporary->NextNode2 = nullptr;
        Temporary->NextNode3 = nullptr;
        
        Temporary->PreviousNode = *B;
        
        (*B)->NextNode2 = Temporary;
        
        // Create fourth child node
        if ( (Temporary = (struct Node *)malloc(sizeof(struct Node)) ) == nullptr )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        Allocate(Temporary->vertices,4,2);
        
        Temporary->level = (*B)->level+1;
        
        Temporary->vertices[0][0] = 0.5*((*B)->vertices[0][0]+(*B)->vertices[3][0]);
        Temporary->vertices[0][1] = 0.5*((*B)->vertices[0][1]+(*B)->vertices[3][1]);
        
        Temporary->vertices[1][0] = 0.25*((*B)->vertices[0][0]+(*B)->vertices[1][0]+(*B)->vertices[2][0]+(*B)->vertices[3][0]);
        Temporary->vertices[1][1] = 0.25*((*B)->vertices[0][1]+(*B)->vertices[1][1]+(*B)->vertices[2][1]+(*B)->vertices[3][1]);
        
        Temporary->vertices[2][0] = 0.5*((*B)->vertices[2][0]+(*B)->vertices[3][0]);
        Temporary->vertices[2][1] = 0.5*((*B)->vertices[2][1]+(*B)->vertices[3][1]);
        
        Temporary->vertices[3][0] = (*B)->vertices[3][0];
        Temporary->vertices[3][1] = (*B)->vertices[3][1];
        
        Temporary->NextNode0 = nullptr;
        Temporary->NextNode1 = nullptr;
        Temporary->NextNode2 = nullptr;
        Temporary->NextNode3 = nullptr;
        
        Temporary->PreviousNode = *B;
        
        (*B)->NextNode3 = Temporary;
        
        // Check again all the newly created elements
        AdaptiveMeshRefinement(&((*B)->NextNode0),MaxLevel);
        AdaptiveMeshRefinement(&((*B)->NextNode1),MaxLevel);
        AdaptiveMeshRefinement(&((*B)->NextNode2),MaxLevel);
        AdaptiveMeshRefinement(&((*B)->NextNode3),MaxLevel);
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CountElement ( struct Node *B, 
                    int &Ne )
{
    if (B != nullptr)
    {
        if (B->NextNode0 == nullptr)
            Ne++;
        else
        {
            CountElement(B->NextNode0,Ne);
            CountElement(B->NextNode1,Ne);
            CountElement(B->NextNode2,Ne);
            CountElement(B->NextNode3,Ne);
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void WriteVertices ( struct Node *B,
                     char *s )
{
    if (B != nullptr)
    {
        if (B->NextNode0 == nullptr)
        {
            std::ofstream TecplotWrite(s,std::ios::app);
            TecplotWrite.flags(std::ios::dec);
            TecplotWrite.precision(8);
            
            if ( !TecplotWrite )
                ErrorMessage("Output file couldnot be opened.");
            
            TecplotWrite << B->vertices[0][0] << "\t" << B->vertices[0][1] << std::endl;
            TecplotWrite << B->vertices[1][0] << "\t" << B->vertices[1][1] << std::endl;
            TecplotWrite << B->vertices[2][0] << "\t" << B->vertices[2][1] << std::endl;
            TecplotWrite << B->vertices[3][0] << "\t" << B->vertices[3][1] << std::endl;
            
            TecplotWrite.close();
        }
        else
        {
            WriteVertices(B->NextNode0,s);
            WriteVertices(B->NextNode1,s);
            WriteVertices(B->NextNode2,s);
            WriteVertices(B->NextNode3,s);
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void WriteCells ( struct Node *B,
                  int &vertexnumber,
                  char *s )
{
    if (B != nullptr)
    {
        if (B->NextNode0 == nullptr)
        {
            std::ofstream TecplotWrite(s,std::ios::app);
            TecplotWrite.flags(std::ios::dec);
            TecplotWrite.precision(8);
            
            if ( !TecplotWrite )
                ErrorMessage("Output file couldnot be opened.");
            
            TecplotWrite << vertexnumber+1 << "\t" << vertexnumber+2 << "\t" << vertexnumber+3 << "\t" << vertexnumber+4 << std::endl;
            
            vertexnumber += 4;
            
            TecplotWrite.close();
        }
        else
        {
            WriteCells(B->NextNode0,vertexnumber,s);
            WriteCells(B->NextNode1,vertexnumber,s);
            WriteCells(B->NextNode2,vertexnumber,s);
            WriteCells(B->NextNode3,vertexnumber,s);
        }
    }
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
    
    struct Node *A[Nx-1][Ny-1];
    
    // Allocate memory
    Allocate(s,200);
    
    for (int i{}; i < (Nx-1); ++i)
        for (int j{}; j < (Ny-1); ++j)
            Allocate(&A[i][j],i,j,dx,dy,xi,yi);
    
    int TimeSteps = 100;
    
    for (int timestep{}; timestep < TimeSteps; ++timestep)
    {
        t = timestep/(frequency*TimeSteps);
        
        // Adaptive mesh refinement
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                AdaptiveMeshRefinement(&A[i][j],MaxLevel);
        
        // Write file
        sprintf(s,"Output/Refinement/Field-%d.tec",timestep);
        
        int Ne = 0, vertexnumber = 0;
        
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                CountElement(A[i][j],Ne);
        
        std::ofstream TecplotWrite(s,std::ios::out);
        TecplotWrite.flags(std::ios::dec);
        TecplotWrite.precision(8);
        
        if ( !TecplotWrite )
            ErrorMessage("Output file couldnot be opened.");
        
        TecplotWrite << "TITLE = \"Mesh refinement\"\nVariables = \"X\",\"Y\"" << std::endl;
        TecplotWrite << "Zone N = " << 4*Ne << ", E = " << Ne << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        TecplotWrite.close();
        
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                WriteVertices(A[i][j],s);
        
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                WriteCells(A[i][j],vertexnumber,s);
        
        for (int i{}; i < (Nx-1); ++i)
            for (int j{}; j < (Ny-1); ++j)
                Deallocate(&A[i][j]);
        
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
    Deallocate(s,200);
    
    for (int i{}; i < (Nx-1); ++i)
    {
        for (int j{}; j < (Ny-1); ++j)
        {
            Deallocate((A[i][j])->vertices,4,2);
            
            free(A[i][j]);
        }
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
