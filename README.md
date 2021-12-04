# GenSDFRaytracer
An implementation of a Raytracer, using generalized SDFs (signed distance functions). 


Source code, executables and a makefile can be found in the src directory. 

SDF.h contains structs and function declarations.

renderFuncs.c contains functions needed to determine calculations of distance 
functions and gradients, taking parameters of the desired shape. Other functions
necessary to compute the shading and output pixel values can be found as well. 

shapes.c contains the SDFs themselves in their canonical format, taking only a
position and light ray orientation, returning a single distance (these functions
often call functions in renderFuncs.c). 

render.c contains the main function, defining the "camera" and shape, and rendering
the image. The SDFs and gradient functions chosen here are those defined in shapes.c.


The bash executable ./execute.sh (pixel width) (pixel height) (output filename)
can be used to compile and run the code. The result is a rendered image.
