Tracer: Render.o RenderFuncs.o SDFs.o
	gcc -o Tracer Render.o SDFs.o RenderFuncs.o

RenderFuncs.o: ../src/RenderFuncs.c ../src/SDF.h
	gcc -c -g src/RenderFuncs.c

SDFs.o: ../src/RenderFuncs.c ../src/SDFs.c ../src/SDF.h
	gcc -c -g src/SDFs.c

Render.o: ../src/Render.c ../src/SDF.h
	gcc -c -g ../src/Render.c

