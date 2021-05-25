#include <stdio.h>
#include <stdlib.h>
#include "SDF.h"
//#include <GLUS.h>

/*struct Camera *newCam(struct vec3 *focus, struct vec3 *corner[4]);

struct vec3 *newVec(double x, double y, double z);

double dot3(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecsubtract(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecmult(double scalar, struct vec3 *vector);

struct vec3 *normalize(struct vec3 *OG);

double genSDFsphere(struct vec3 *center, double r, struct vec3 *v, struct vec3 *p);*/

int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp);

double genSDFbox(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *v, struct vec3 *p);

struct vec3 *boxGrad(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p);

/*int main(int argc, char *argv[]){
	
	struct vec3 *corners[4]; // = malloc(sizeof(struct vec3)*4);
	corners[0] = newVec(5.0, -1.0, 1.5);
	corners[1] = newVec(5.0, 1.0, 1.5);
	corners[2] = newVec(5.375, 1.0, 0.5);
	corners[3] = newVec(5.375, -1.0, 0.5);
	
	struct vec3 *focus = newVec(8.0, 0.0, 1.0);

	
	struct Camera *camera = newCam(focus, corners);
	struct vec3 *center = newVec(0.0, 0.0, 2.0);	

	Render(2.0, center, camera, 100, 50);
}*/
/*double firstCirc(struct vec3 *v, struct vec3 *p){
        struct vec3 *center = 
	struct vec3 *a = vecsubtract(p, center);
        double result;
        if (pow(dot3(a, v), 2.0) + pow(r, 2.0) - dot3(p, p) < 0.0) return 1111.1;
        else result = -dot3(v, a) + pow((double) pow((dot3(a, v)), 2.0) + pow(r, 2.0) - dot3(p, p), 0.5);
        //printf("%f\n", result);
        if (result > 0.0) return result;
        else return 1111.1;
}*/

double circleSDF(struct vec3 *v, struct vec3 *p){
        return genSDFsphere(newVec(0.0, 0.0, 1.7), 1.7, v, p);
}


struct vec3 *circleSDFGrad(struct vec3 *p){
	return normalize(vecsubtract(p, newVec(0.0, 0.0, 1.0)));
}

double boxSDF(struct vec3 *v, struct vec3 *p){
	struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0)); 
	struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
	struct vec3 *n3 = normalize(newVec(0.0,0.0,1.0));
	struct vec3 *corn = newVec(0.0, 0.0, 1.7);
	double dims[3] = {3.0, 3.0, 1.0};
	struct vec3 *norms[3];
	norms[0] = n1;
	norms[1] = n2;
	norms[2] = n3;
	double result = genSDFbox(norms, corn, dims, v, p);
	free(n1);
	free(n2);
	free(n3);
	//free(dims);
	free(corn);
	return result;
}

struct vec3 *boxSDFgrad(struct vec3 *p){
	struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = normalize(newVec(0.0,0.0,1.0));
        struct vec3 *corn = newVec(0.0, 0.0, 1.7);
        double dims[3] = {3.0, 3.0, 1.0};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

/*double boxGrad(struct vec3 *p){
	if (dot3(p->x, newVec(1.1) <)
}*/

int main(int argc, char *argv[]){
	struct vec3 *focus = newVec(8.0, 0.0, 2.5);
        struct vec3 *light = newVec(10.0, 10.0, 10.0);


	struct vec3 *corners[4];// = malloc(sizeof(struct vec3)*4);
        corners[0] = newVec(4.7, focus->y - 1, focus->z + 0.5);
        corners[1] = newVec(4.7, focus->y + 1, focus->z + 0.5);
        corners[2] = newVec(5.75, focus->y + 1, focus->z - 0.5);
        corners[3] = newVec(5.75, focus->y -1, focus->z - 0.5);
        
	struct Node *head = malloc(sizeof(struct Node));
	//struct Node *head;
	head->next = 0;
	head->shape = newShape(newVec(0, 255, 0));
	//head->shape->SDF = boxGrad;
	head->shape->SDF = boxSDF;
	head->shape->gradSDF = boxSDFgrad;
	//head->shape->gradSDF = circleSDFGrad;
	//head->gradSDF = firstCirc;
        
	struct Shape *sphere = newShape(newVec(0, 255, 0));
	sphere->SDF = circleSDF;
	sphere->gradSDF = circleSDFGrad;
	newNode(head, sphere);
		

	struct Camera *camera = newCam(focus, corners, light);
        //struct vec3 *center = newVec(0.0, 0.0, 2.0);
	FILE *fp = fopen(argv[3], "wt");
	//printf("HiHi\n");
	//printf("%f\n", boxSDF(normalize(newVec(-8.0, 0.0, -0.95)), focus));
	Render(head, camera, atoi(argv[1]), atoi(argv[2]), fp);
	fclose(fp);
	free(head);
	free(focus);
	free(light);	
	return 0;
}

