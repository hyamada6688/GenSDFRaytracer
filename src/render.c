#include <stdio.h>
#include <stdlib.h>
#include "SDF.h"
#include "shapes.h"
//#include "SDFs.c"
//#include <GLUS.h>

/*struct vec3 *newVec(double x, double y, double z);

struct Camera *newCam(struct vec3 *focus, struct vec3 *corner[4], struct vec3 *light);

int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp);
*/
/*double circleSDF(struct vec3 *v, struct vec3 *p);


struct vec3 *circleSDFGrad(struct vec3 *p);

double stemSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *stemgrad(struct vec3 *p);

double boxSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *boxSDFgrad(struct vec3 *p);


double planeSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *planeSDFGrad(struct vec3 *p);

double leg1SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg1SDFgrad(struct vec3 *p);

double leg2SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg2SDFgrad(struct vec3 *p);

double leg3SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg3SDFgrad(struct vec3 *p);

double leg4SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg4SDFgrad(struct vec3 *p);

double roundBoxSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *roundboxSDFgrad(struct vec3 *p);

double cylinder(struct vec3 *v, struct vec3 *p);

struct vec3 *cylinderGrad(struct vec3 *p);

struct vec3 *prismSDFgrad(struct vec3 *p);

double prismgenSDF(struct vec3 *v, struct vec3 *p);
*/
int main(int argc, char *argv[]){
	printf("%s, %s, %s\n", argv[1],argv[2],argv[3]);
	struct vec3 *focus = newVec(8.0, 0.0, 2.5);
        //struct vec3 *light = newVec(5.0, 3.5 , 3.5);
	struct vec3 *light = newVec(8,-1,7);

	struct vec3 *corners[4];// = malloc(sizeof(struct vec3)*4);
        corners[0] = newVec(focus->x - 3, focus->y - 1, focus->z -0.5);
        corners[1] = newVec(focus->x - 3, focus->y + 1, focus->z - 0.5);
        corners[2] = newVec(focus->x - 2, focus->y + 1, focus->z - 1.5);
        corners[3] = newVec(focus->x - 2, focus->y -1, focus->z - 1.5);
        
	struct Node *head = malloc(sizeof(struct Node));
	//struct Node *head;
/*	head->next = 0;
	head->shape = newShape(newVec(120, 180, 250));
	//head->shape->SDF = boxGrad;
	head->shape->SDF = roundBoxSDF;
	head->shape->gradSDF = roundboxSDFgrad;
	//head->shape->gradSDF = circleSDFGrad;
        //head->shape->SDF = circleSDF;
*/	
	//struct Shape *box = newShape(newVec(200, 200, 200));
	head->shape = newShape(newVec(120, 180, 250));
	head->shape->SDF = prismgenSDF;
	head->shape->gradSDF = prismSDFgrad;
	//newNode(head, box);
	//head->shape = box;
	head->next = 0; 
	
	/*struct Shape *prism = newShape(newVec(255, 240, 220));
	prism->SDF = prismgenSDF;
        prism->gradSDF = prismSDFgrad;
	head->shape = prism;
	head->next = 0;
	*/
/*	struct Shape *sphere = newShape(newVec(255, 0, 100));
	sphere->SDF = circleSDF;
	sphere->gradSDF = circleSDFGrad;
	newNode(head, sphere);
	
	struct Shape *stem = newShape(newVec(0, 255, 70));
	stem->SDF = stemSDF;
	stem->gradSDF = stemgrad;
	newNode(head, stem);
		
	struct Shape *floor = newShape(newVec(255, 155, 51));
	floor->SDF = planeSDF;
	floor->gradSDF = planeSDFGrad;
	newNode(head, floor);	

	struct Shape *leg1 = newShape(newVec(255, 240, 220));
	leg1->SDF = leg1SDF;
	leg1->gradSDF = leg1SDFgrad;
	newNode(head, leg1);

	struct Shape *leg2 = newShape(newVec(255, 240, 220));
        leg2->SDF = leg2SDF;
        leg2->gradSDF = leg2SDFgrad;
        newNode(head, leg2);

	struct Shape *leg3 = newShape(newVec(255, 240, 220));
        leg3->SDF = leg3SDF;
        leg3->gradSDF = leg3SDFgrad;
        newNode(head, leg3);

	struct Shape *leg4 = newShape(newVec(255, 240, 220));
        leg4->SDF = leg4SDF;
        leg4->gradSDF = leg4SDFgrad;
        newNode(head, leg4);*/
	

	struct Camera *camera = newCam(focus, corners, light);
        //struct vec3 *center = newVec(0.0, 0.0, 2.0);
	FILE *fp = fopen(argv[3], "wt");
	printf("HiHi\n");
	//printf("%f\n", boxSDF(normalize(newVec(-8.0, 0.0, -0.95)), focus));
	Render(head, camera, atoi(argv[1]), atoi(argv[2]), fp);
	fclose(fp);
	free(head);
	free(focus);
	free(light);	
	return 0;
}

