#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double circleSDF(struct vec3 *v, struct vec3 *p);


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

double cylinder1(struct vec3 *v, struct vec3 *p);

struct vec3 *cylinderGrad1(struct vec3 *p);

double cylinder2(struct vec3 *v, struct vec3 *p);

struct vec3 *cylinderGrad2(struct vec3 *p);

struct vec3 *prismSDFgrad(struct vec3 *p);

double prismgenSDF(struct vec3 *v, struct vec3 *p);

double coneGenSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *coneSDFGrad(struct vec3 *p);

struct contact *cylEval(struct vec3 *v, struct vec3 *p);

double cylinder(struct vec3 *v, struct vec3 *p);

struct vec3 *cylinderGrad(struct vec3 *p);

double infCyl(struct vec3 *v, struct vec3 *p);

struct vec3 *infCylGrad(struct vec3 *p);

double infcone(struct vec3 *v, struct vec3 *p);

struct vec3 *infconeGrad(struct vec3 *p);

double iceCreamCone(struct vec3 *v, struct vec3 *p);

struct vec3 *iceCreamConeGrad(struct vec3 *p);

double iceCream(struct vec3 *v, struct vec3 *p);

struct vec3 *iceCreamGrad(struct vec3 *p);

