#include <stdio.h>
#include <math.h>

float deg2rad(float deg) {
    return deg * M_PI / 180.0;
}

void dotsToVector(const float a[3], const float b[3], float* result) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

void nVector(const float a[3], const float b[3], float* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

float nVectorMagnitude(const float a[3], const float b[3]){
  float nvec[3];
  nVector(a, b, nvec);

  return sqrt(pow(nvec[0], 2) + pow(nvec[1], 2) + pow(nvec[2], 2));
}

int checkDot(const float dot[3], const float vec1[3], const float vec2[3], const float vec3[3], const float vec4[3], float a){
  float vec12[3];
  float vec23[3];
  dotsToVector(vec1, vec2, vec12);
  dotsToVector(vec2, vec3, vec23);

  float nvec[3];
  nVector(vec12, vec23, nvec);

  float constantT = (nvec[0] * (vec1[0] - 2 * a) + nvec[1] * vec1[1] + nvec[2] * vec1[2]) / (nvec[0] * (2 * a - dot[0]) - nvec[1] * dot[1] - nvec[2] * dot[2]);

  float dotOnPlane[3] = {2 * a + (2 * a - dot[0]) * constantT, -dot[1] * constantT, -dot[2] * constantT};

  if(constantT > 0.0f || constantT < -1.0f) return 1;

  float vec1andDot[3];
  float vec2andDot[3];
  float vec3andDot[3];
  float vec4andDot[3];

  dotsToVector(vec1, dotOnPlane, vec1andDot);
  dotsToVector(vec2, dotOnPlane, vec2andDot);
  dotsToVector(vec3, dotOnPlane, vec3andDot);
  dotsToVector(vec4, dotOnPlane, vec4andDot);

  float p1 = nVectorMagnitude(vec1andDot, vec2andDot) / 2.0f;
  float p2 = nVectorMagnitude(vec2andDot, vec3andDot) / 2.0f;
  float p3 = nVectorMagnitude(vec3andDot, vec4andDot) / 2.0f;
  float p4 = nVectorMagnitude(vec4andDot, vec1andDot) / 2.0f;
  
  if(fabsf(p1 + p2 + p3 + p4 - powf(a, 2)) >= 1e-3f) return 1;

  return 0;
}

void rotateVector(const float vector[3], const float angles[3], float* result){
  result[0] = cos(deg2rad(angles[2])) * (vector[0] * cos(deg2rad(angles[1])) + (vector[1] * sin(deg2rad(angles[0])) + vector[2] * cos(deg2rad(angles[0]))) * sin(deg2rad(angles[1]))) - sin(deg2rad(angles[2])) * (vector[1] * cos(deg2rad(angles[0])) - vector[2] * sin(deg2rad(angles[0])));
  result[1] = cos(deg2rad(angles[2])) * (vector[1] * cos(deg2rad(angles[0])) - vector[2] * sin(deg2rad(angles[0]))) + sin(deg2rad(angles[2])) * (vector[0] * cos(deg2rad(angles[1])) + (vector[1] * sin(deg2rad(angles[0])) + vector[2] * cos(deg2rad(angles[0]))) * sin(deg2rad(angles[1])));
  result[2] = cos(deg2rad(angles[1])) * (vector[1] * sin(deg2rad(angles[0])) + vector[2] * cos(deg2rad(angles[0]))) - vector[0] * sin(deg2rad(angles[1]));
}

int checkField(int fieldX, int fieldY, float xpA, float ypA, float xpB, float ypB, float xpC, float ypC, float xpD, float ypD){
  int tacnostAB = 0;
  int tacnostBC = 0;
  int tacnostCD = 0;
  int tacnostDA = 0;

  float kAB;
  float kBC;
  float kCD;
  float kDA;

  int isVerticalAB = 0;
  int isVerticalBC = 0;
  int isVerticalCD = 0;
  int isVerticalDA = 0;

  if(xpA != xpB)
    kAB = (ypA - ypB)/(xpA - xpB);
  else
    isVerticalAB = 1; 

  if(xpB != xpC)
    kBC = (ypB - ypC)/(xpB - xpC);
  else
    isVerticalBC = 1;

  if(xpC != xpD)
    kCD = (ypC - ypD)/(xpC - xpD);
  else
    isVerticalCD = 1;

  if(xpD != xpA)
    kDA = (ypD - ypA)/(xpD - xpA);
  else
    isVerticalDA = 1;

    //AB duz
  if(fieldY <= fmaxf(ypA, ypB) && fieldY >= fminf(ypA, ypB) && (isVerticalAB ? xpA : (fieldY - ypA + kAB * xpA)/kAB) >= fieldX)
    tacnostAB = 1;

    //BC duz
  if(fieldY <= fmaxf(ypB, ypC) && fieldY >= fminf(ypB, ypC) && (isVerticalBC ? xpB : (fieldY - ypB + kBC * xpB)/kBC) >= fieldX)
    tacnostBC = 1;

    //CD duz
  if(fieldY <= fmaxf(ypC, ypD) && fieldY >= fminf(ypC, ypD) && (isVerticalCD ? xpC : (fieldY - ypC + kCD * xpC)/kCD) >= fieldX)
    tacnostCD = 1;

    //DA duz
  if(fieldY <= fmaxf(ypD, ypA) && fieldY >= fminf(ypD, ypA) && (isVerticalDA ? xpD : (fieldY - ypD + kDA * xpD)/kDA) >= fieldX)
    tacnostDA = 1;
  
  if((tacnostAB + tacnostBC + tacnostCD + tacnostDA) % 2 == 1)
    return 1;

  return 0;
}

void render(const float angles[3], float a){
  float root = a / 2.0;
  float resolution = a * 3 / 2;

  float vecA[3];
  float vecB[3];
  float vecC[3];
  float vecD[3];
  float vecE[3];
  float vecF[3];
  float vecG[3];
  float vecH[3];

  rotateVector((float[]){root, root, root}, angles, vecA);
  rotateVector((float[]){root, -root, root}, angles, vecB);
  rotateVector((float[]){root, -root, -root}, angles, vecC);
  rotateVector((float[]){root, root, -root}, angles, vecD);
  rotateVector((float[]){-root, root, root}, angles, vecE);
  rotateVector((float[]){-root, -root, root}, angles, vecF);
  rotateVector((float[]){-root, -root, -root}, angles, vecG);
  rotateVector((float[]){-root, root, -root}, angles, vecH);

  int prikazivanjeABCD = 0;
  if(checkDot(vecC, vecA, vecB, vecF, vecE, a) && checkDot(vecD, vecA, vecB, vecF, vecE, a) && checkDot(vecA, vecB, vecC, vecG, vecF, a) && checkDot(vecD, vecB, vecC, vecG, vecF, a) && checkDot(vecB, vecA, vecD, vecH, vecE, a) && checkDot(vecC, vecA, vecD, vecH, vecE, a) && checkDot(vecA, vecC, vecD, vecH, vecG, a) && checkDot(vecB, vecC, vecD, vecH, vecG, a) && checkDot(vecA, vecE, vecF, vecG, vecH, a) && checkDot(vecB, vecE, vecF, vecG, vecH, a) && checkDot(vecC, vecE, vecF, vecG, vecH, a) && checkDot(vecD, vecE, vecF, vecG, vecH, a)) prikazivanjeABCD = 1;
  
  int prikazivanjeABFE = 0;
  if(checkDot(vecF, vecA, vecB, vecC, vecD, a) && checkDot(vecE, vecA, vecB, vecC, vecD, a) && checkDot(vecA, vecB, vecC, vecG, vecF, a) && checkDot(vecE, vecB, vecC, vecG, vecF, a) && checkDot(vecB, vecA, vecD, vecH, vecE, a) && checkDot(vecF, vecA, vecD, vecH, vecE, a) && checkDot(vecA, vecC, vecD, vecH, vecG, a) && checkDot(vecB, vecC, vecD, vecH, vecG, a) && checkDot(vecF, vecC, vecD, vecH, vecG, a) && checkDot(vecE, vecC, vecD, vecH, vecG, a) && checkDot(vecA, vecE, vecF, vecG, vecH, a) && checkDot(vecB, vecE, vecF, vecG, vecH, a)) prikazivanjeABFE = 1;

  int prikazivanjeBCGF = 0;
  if(checkDot(vecF, vecA, vecB, vecC, vecD, a) && checkDot(vecG, vecA, vecB, vecC, vecD, a) && checkDot(vecC, vecA, vecB, vecF, vecE, a) && checkDot(vecG, vecA, vecB, vecF, vecE, a) && checkDot(vecB, vecA, vecD, vecH, vecE, a) && checkDot(vecC, vecA, vecD, vecH, vecE, a) && checkDot(vecG, vecA, vecD, vecH, vecE, a) && checkDot(vecF, vecA, vecD, vecH, vecE, a) && checkDot(vecB, vecC, vecD, vecH, vecG, a) && checkDot(vecF, vecC, vecD, vecH, vecG, a) && checkDot(vecB, vecE, vecF, vecG, vecH, a) && checkDot(vecC, vecE, vecF, vecG, vecH, a)) prikazivanjeBCGF = 1;

  int prikazivanjeADHE = 0;
  if(checkDot(vecA, vecB, vecC, vecG, vecF, a) && checkDot(vecD, vecB, vecC, vecG, vecF, a) && checkDot(vecD, vecA, vecB, vecF, vecE, a) && checkDot(vecH, vecA, vecB, vecC, vecD, a) && checkDot(vecH, vecA, vecB, vecF, vecE, a) && checkDot(vecH, vecB, vecC, vecG, vecF, a) && checkDot(vecE, vecA, vecB, vecC, vecD, a) && checkDot(vecE, vecB, vecC, vecG, vecF, a) && checkDot(vecA, vecC, vecD, vecH, vecG, a) && checkDot(vecE, vecC, vecD, vecH, vecG, a) && checkDot(vecA, vecE, vecF, vecG, vecH, a) && checkDot(vecD, vecE, vecF, vecG, vecH, a)) prikazivanjeADHE = 1;

  int prikazivanjeCDHG = 0;
  if(checkDot(vecC, vecA, vecB, vecF, vecE, a) && checkDot(vecC, vecA, vecD, vecH, vecE, a) && checkDot(vecD, vecA, vecB, vecF, vecE, a) && checkDot(vecD, vecB, vecC, vecG, vecF, a) && checkDot(vecH, vecA, vecB, vecC, vecD, a) && checkDot(vecH, vecA, vecB, vecF, vecE, a) && checkDot(vecH, vecB, vecC, vecG, vecF, a) && checkDot(vecG, vecA, vecB, vecC, vecD, a) && checkDot(vecG, vecA, vecB, vecF, vecE, a) && checkDot(vecG, vecA, vecD, vecH, vecE, a) && checkDot(vecC, vecE, vecF, vecG, vecH, a) && checkDot(vecD, vecE, vecF, vecG, vecH, a)) prikazivanjeCDHG = 1;

  int prikazivanjeEFGH = 0;
  if(checkDot(vecE, vecA, vecB, vecC, vecD, a) && checkDot(vecE, vecB, vecC, vecG, vecF, a) && checkDot(vecE, vecC, vecD, vecH, vecG, a) && checkDot(vecF, vecA, vecB, vecC, vecD, a) && checkDot(vecF, vecA, vecD, vecH, vecE, a) && checkDot(vecF, vecC, vecD, vecH, vecG, a) && checkDot(vecG, vecA, vecB, vecC, vecD, a) && checkDot(vecG, vecA, vecB, vecF, vecE, a) && checkDot(vecG, vecA, vecD, vecH, vecE, a) && checkDot(vecH, vecA, vecB, vecC, vecD, a) && checkDot(vecH, vecA, vecB, vecF, vecE, a) && checkDot(vecH, vecB, vecC, vecG, vecF, a)) prikazivanjeEFGH = 1;

  float xpA = (vecA[1] * a) / (2 * a - vecA[0]);
  float ypA = (vecA[2] * a) / (2 * a - vecA[0]);
 
  float xpB = (vecB[1] * a) / (2 * a - vecB[0]);
  float ypB = (vecB[2] * a) / (2 * a - vecB[0]);

  float xpC = (vecC[1] * a) / (2 * a - vecC[0]);
  float ypC = (vecC[2] * a) / (2 * a - vecC[0]);

  float xpD = (vecD[1] * a) / (2 * a - vecD[0]);
  float ypD = (vecD[2] * a) / (2 * a - vecD[0]);

  float xpE = (vecE[1] * a) / (2 * a - vecE[0]);
  float ypE = (vecE[2] * a) / (2 * a - vecE[0]);

  float xpF = (vecF[1] * a) / (2 * a - vecF[0]);
  float ypF = (vecF[2] * a) / (2 * a - vecF[0]);

  float xpG = (vecG[1] * a) / (2 * a - vecG[0]);
  float ypG = (vecG[2] * a) / (2 * a - vecG[0]);

  float xpH = (vecH[1] * a) / (2 * a - vecH[0]);
  float ypH = (vecH[2] * a) / (2 * a - vecH[0]);

 for(int i = 0; i < resolution / 2; i++){
    for(int j = 0; j < resolution; j++){
      int fieldX = j - resolution / 2;
      int fieldY = resolution / 2 - i * 2;

      if(j == 0)
        printf("\n     ");

      if(prikazivanjeABCD && checkField(fieldX, fieldY, xpA, ypA, xpB, ypB, xpC, ypC, xpD, ypD))
        printf("O");
      else if(prikazivanjeABFE && checkField(fieldX, fieldY, xpA, ypA, xpB, ypB, xpF, ypF, xpE, ypE))
        printf("*");
      else if(prikazivanjeBCGF && checkField(fieldX, fieldY, xpB, ypB, xpC, ypC, xpG, ypG, xpF, ypF))
        printf("!");
      else if(prikazivanjeADHE && checkField(fieldX, fieldY, xpA, ypA, xpD, ypD, xpH, ypH, xpE, ypE))
        printf("~");
      else if(prikazivanjeCDHG && checkField(fieldX, fieldY, xpC, ypC, xpD, ypD, xpH, ypH, xpG, ypG))
        printf("^");
      else if(prikazivanjeEFGH && checkField(fieldX, fieldY, xpE, ypE, xpF, ypF, xpG, ypG, xpH, ypH))
        printf(".");
      else
        printf(" ");
    }
  }
  printf("\n\n");
}

int main(){
  printf("Unesi stranicu kocke (preporuceno vise od 30): ");
  int a;
  scanf("%d", &a);
  if(a > 130)
    a = 130;

  printf("Unesi ugao X: ");
  float fi;
  scanf("%f", &fi);
  render((float[]){fi, 0.0, 0.0}, a);   
  
  printf("Unesi ugao Y: ");
  float teta;
  scanf("%f", &teta);
  render((float[]){fi, teta, 0.0}, a);   

  printf("Unesi ugao Z: ");
  float psi;
  scanf("%f", &psi);
  render((float[]){fi, teta, psi}, a);   

  return 0;
}
