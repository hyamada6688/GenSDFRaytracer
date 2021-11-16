Tracer: Render.o RenderFuncs.o shapes.o
	gcc -o Tracer Render.o shapes.o RenderFuncs.o

RenderFuncs.o: renderFuncs.c SDF.h
	gcc -c -g renderFuncs.c

shapes.o: renderFuncs.c shapes.c SDF.h
	gcc -c -g shapes.c

Render.o: render.c SDF.h
	gcc -c -g render.c

