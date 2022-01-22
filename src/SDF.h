#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Camera{
	struct vec3 *focus;
	struct vec3 *corner[4];
	struct vec3 *light;
};

struct Node{
	struct Node *next;
	struct Shape *shape;
};


//typedef struct vec3 *(*r3r3)(struct Sphere *shape, struct vec3 *v, struct vec3 *pos);

//typedef double (*r3r1)(struct Sphere *shape, struct vec3 *pos);

typedef struct vec3 *(*r3r3)(struct vec3 *);

typedef double (*r3r1)(struct vec3 *, struct vec3 *);

typedef struct contact *(*evalfunc)(struct vec3 *, struct vec3 *);

struct Shape {
        struct vec3 *color;
	r3r3 gradSDF;
	r3r1 SDF;
	evalfunc eval;
};

/*struct Shape{
	struct vec3 *color;
        r3r3 gradSDF;
        r3r1 SDF;
	union {
		struct Sphere sphereFacts;
		struct Plane planeFacts;
		struct Cone coneFacts;
		struct Cylinder cylFacts;
		struct Prism prismFacts;
		struct Box boxFacts;
		struct 
		
	};
	double roundRad = 0;
}*/

struct Sphere {
	struct Shape base;
	double radius;
	struct vec3 *center;
	r3r1 SDF;
	struct vec3 *color;
};

struct vec3{
	double x;
	double y;
	double z;
};

struct Plane {
	struct Shape base;
	struct vec3 *normal;
	struct vec3 *pointOn;
};

struct Box {
	struct vec3 *corn;
	struct vec3 *norms[3];
	double dims[3];	
};

struct Cone {
	struct vec3 *a;
	struct vec3 *b;
	double h;
	double theta; 
};

struct Cylinder {
	struct vec3 *a;
        struct vec3 *b;
        double h;
};

struct Prism{
	struct vec3 *b;
	int faces;
	struct vdc3 *n0;
	double width;
};


struct SD {
	struct vec3 *pos;
	double dist;
	r3r3 gradient;
	struct vec3 *color;
	struct vec3 *grad;
};

struct contact{
	struct vec3 *p;
	double dist;
	struct vec3 *grad;
};

double max(double arr[]);

int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp);
struct Camera *newCam(struct vec3 *focus, struct vec3 *corner[4], struct vec3 *light);

struct vec3 *newVec(double x, double y, double z);

double dot3(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecsubtract(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecmult(double scalar, struct vec3 *vector);

struct vec3 *vecadd(struct vec3 *v1, struct vec3 *v2);

int vecEq(struct vec3 *v1, struct vec3 *v2);

struct vec3 *normalize(struct vec3 *OG);

double vecnorm(struct vec3 *OG);

struct vec3 *antiComponent(struct vec3 *pb, struct vec3 *a);

struct vec3 *component(struct vec3 *pb, struct vec3 *a);

struct vec3 *shade(struct Node *head, struct Camera *camera, struct SD *sd);

struct SD *calcSDF(struct Node *head, struct vec3 *v, struct vec3 *p);

void newNode(struct Node *head, struct Shape *shape);

struct Plane *newPlane(struct vec3 *normal, struct vec3 *pointon);

struct Sphere *newSphere(double radius, struct vec3 *center);

struct Shape *newShape(struct vec3 *color);




double genSDFsphere(struct vec3 *center, double r, struct vec3 *v, struct vec3 *p);

double genSDFplane(struct vec3 *normal, struct vec3 *pointOn, struct vec3 *v, struct vec3 *p);


//struct vec3 *shade(struct Node *head, struct Camera *camera, struct vec3 *v, struct vec3 *pos, struct SD *sd);

double genSDFbox(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *v, struct vec3 *p);

struct vec3 *boxGrad(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p);

struct vec3 *RoundboxGrad(double r, struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p);

double cylSDF(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a);

struct vec3 *cylgrad(struct vec3 *p, struct vec3 *b, struct vec3 *a);

struct contact *evalCylinder(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a, double h, struct vec3 *b);

double coneSDF(struct vec3 *v, struct vec3 *p, struct vec3 *b, struct vec3 *a, double theta);

struct vec3 *coneGrad(struct vec3 *p, struct vec3 *b, struct vec3 *a, double t);

double prismSDF(struct vec3 *base, struct vec3 *a, struct vec3 *n0, int faces, double height, double facewidth, struct vec3 *v, struct vec3 *p);

struct vec3 *prismGrad(struct vec3 *base, struct vec3 *a, struct vec3 *n0, int faces, double height, double facewidth, struct vec3 *p);

double maxim(double arr[], int size);

int valCompare(double a, double b);
