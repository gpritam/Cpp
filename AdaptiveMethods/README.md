## Adaptive meshing techniques

Adaptive meshing is a powerful technique used in computational simulations to improve the accuracy and efficiency of numerical solutions. It involves dynamically modifying the mesh, which is a discretized representation of the computational domain, to better capture the solution features. There are several methods available for adaptive meshing, including adaptive mesh redistribution, adaptive mesh refinement, and p-adaption. Those three adaptive methodologies are demonstrated by adapting the underlying grid to simulate laminar, non-premixed flame as below.
A laminar, non-premixed flame is a type of flame where the fuel and oxidizer are not mixed prior to combustion. These flames exhibit distinct layers of fuel and oxidizer, with combustion occurring at the flame front. They find applications in various combustion systems and are extensively studied in combustion research. Accurate modeling techniques are required to simulate and analyze the behavior of laminar, non-premixed flames.


| ![Simulation of laminar, non-premixed flame](https://github.com/gpritam/Cpp/blob/main/AdaptiveMethods/Output/Flame/Flame.gif?raw=true) | ![Adaptive mesh redistribution](https://github.com/gpritam/Cpp/blob/main/AdaptiveMethods/Output/Redistribution/Redistribution.gif?raw=true) | ![Adaptive mesh refinement](https://github.com/gpritam/Cpp/blob/main/AdaptiveMethods/Output/Refinement/Refinement.gif?raw=true) | ![p-Refinement](https://github.com/gpritam/Cpp/blob/main/AdaptiveMethods/Output/p-Adaption/p-Adaption.gif?raw=true) |


## Adaptive Mesh Redistribution or, r-Adaption

Adaptive mesh redistribution is a technique that aims to redistribute the existing mesh nodes in a way that improves the mesh quality and resolves solution features more effectively. It involves moving the mesh nodes based on certain criteria, such as solution gradients or error indicators. The goal is to achieve a more uniform distribution of mesh nodes, which can lead to better accuracy and convergence of the numerical solution. This technique is particularly useful when the initial mesh is not well-suited for capturing the solution features or when the solution undergoes significant changes during the simulation. To build and run this code, invoke the following commands

```
make TARGET=r-Adaption.cpp
make run
```
Upon successful completion of code execution, all the output data files in .tec format can be found in the Output/Redistribution directory. The output of the program can be visualized using Tecplot. At first, images (in .jpeg format) are generated in the same directory by running the script Movie-Redistribution.mcr in the Tecplot environment. Note that the |MFBD| variable in the .mcr script needs to be changed for the successful execution of the script. This particular variable denotes the path of the .tec data files. To generate the final movie file (in .gif format), invoke the following command

```
convert -resize 100% -delay 5 -loop 0 *.jpeg Redistribution.gif
```
The last command requires imagemagick to be installed.



## Adaptive Mesh Refinement or, h-Adaption

Adaptive mesh refinement, on the other hand, focuses on selectively refining or coarsening the mesh in regions of interest. Instead of redistributing the existing mesh nodes, this technique adds or removes mesh nodes in specific areas to achieve a higher resolution where it is needed the most. The refinement criteria are typically based on solution gradients, error estimates, or other indicators of solution features. By adaptively refining the mesh, the computational resources can be concentrated on regions that require more accuracy, while coarser meshes are used in areas with less variation. This approach can significantly reduce the computational cost while maintaining the desired level of accuracy.
To build and run this code, invoke the following commands

```
make TARGET=h-Adaption.cpp
make run
```
Upon successful completion of code execution, all the output data files in .tec format can be found in the Output/Refinement directory. The output of the program can be visualized using Tecplot. At first, images (in .jpeg format) are generated in the same directory by running the script Movie-Refinement.mcr in the Tecplot environment. Note that the |MFBD| variable in the .mcr script needs to be changed for the successful execution of the script. This particular variable denotes the path of the .tec data files. To generate the final movie file (in .gif format), invoke the following command

```
convert -resize 100% -delay 5 -loop 0 *.jpeg Refinement.gif
```
The last command requires imagemagick to be installed.

## p-Adaption

p-Adaption, also known as polynomial degree adaption, is a technique that involves changing the order of the polynomial basis functions used to represent the solution within each element of the mesh. Instead of uniformly using the same polynomial degree throughout the domain, p-adaption allows for varying the degree locally based on the solution requirements. Higher polynomial degrees provide more accurate representations of the solution, but they also increase the computational cost. By adaptively adjusting the polynomial degree, p-adaption can achieve a balance between accuracy and computational efficiency. This technique is particularly useful when the solution exhibits sharp gradients or discontinuities. To build and run this code, invoke the following commands

```
make TARGET=p-Adaption.cpp
make run
```
Upon successful completion of code execution, all the output data files in .tec format can be found in the Output/p-Adaption directory. The output of the program can be visualized using Tecplot. At first, images (in .jpeg format) are generated in the same directory by running the script Movie-p-Adaption.mcr in the Tecplot environment. Note that the |MFBD| variable in the .mcr script needs to be changed for the successful execution of the script. This particular variable denotes the path of the .tec data files. To generate the final movie file (in .gif format), invoke the following command

```
convert -resize 100% -delay 5 -loop 0 *.jpeg p-Adaption.gif
```
The last command requires imagemagick to be installed.

In summary, adaptive mesh redistribution focuses on redistributing the existing mesh nodes to improve mesh quality, adaptive mesh refinement selectively adds or removes mesh nodes to achieve higher resolution in regions of interest, and p-adaption adjusts the polynomial degree locally to balance accuracy and computational cost. These techniques can be used individually or in combination to enhance the accuracy and efficiency of numerical simulations. The choice of technique depends on the specific problem and the desired trade-off between accuracy and computational cost.

## References

1. H. D. Ceniceros and T. Y. Hou An Efficient Dynamically Adaptive Mesh for Potentially Singular Solutions.
 Journal of Computational Physics 172, pp. 609-639 (2001)
