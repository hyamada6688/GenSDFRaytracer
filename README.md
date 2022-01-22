# GenSDFRaytracer
An implementation of a Raytracer, using generalized SDFs (signed distance functions),
written in C, with images being rendered with Python PIL.


Source code, executables and a makefile can be found in the src directory. 

SDF.h contains struct and function declarations.

renderFuncs.c contains functions needed to determine calculations of distance 
functions and gradients, taking parameters of the desired shape. Other functions
necessary to compute the shading and output pixel values can be found as well. 

shapes.c contains the SDFs themselves in their canonical format: taking only a
position and light ray orientation, returning a single distance (these functions
often call functions in renderFuncs.c). 

render.c contains the main function, defining positions for the "camera" and shape,
and rendering the image. The SDFs and gradient functions chosen here are those defined 
in shapes.c. To render multiple shapes in the same image, a linked list is used. To define
multiple shapes, set the "next" atttribute to another "node":

struct node *shape2 = malloc(sizeof(struct node));
.
.
.  //initialize the shape with an RGB color and SDF, gradient chosen from shapes.c
shape1->next = shape2;    //set shape2 to be the "next" of shape1


The bash executable ./execute.sh (pixel width) (pixel height) can be used to compile and
run the code. The result is a temporary file of a rendered image, whose bitmap is saved 
in data/small.csv. 
  
Some sample renders can be found in the Renders folder.
