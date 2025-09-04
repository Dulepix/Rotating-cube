#include <stdio.h>
#include <math.h>

#define DEG2RAD (M_PI / 180.0)

int main() {
    printf("This program will rotate a square as you like (Exit code is '2006')\n");
    printf("Enter the square's side length: ");

    int length = 40;
    scanf("%d", &length);
    if(length > 100)
      length = 100;

    float root = (length / 2.0f) * sqrt(2.0f);
    int resolution = ceil(root) * 2;

  while(1){

    int angle;
    printf("\nUnesi ugao: ");
    scanf("%d", &angle);
    
    if(angle == 2006)
      break;

    float xA = root * cos((45 + angle) * DEG2RAD);
    float yA = root * sin((45 + angle) * DEG2RAD);
    float xB = root * cos((135 + angle) * DEG2RAD);
    float yB = root * sin((135 + angle) * DEG2RAD);
    float xC = root * cos((225 + angle) * DEG2RAD);
    float yC = root * sin((225 + angle) * DEG2RAD);
    float xD = root * cos((315 + angle) * DEG2RAD);
    float yD = root * sin((315 + angle) * DEG2RAD);

    for(int i = 0; i < resolution / 2; i++){
      for(int j = 0; j < resolution; j++){
        int fieldX = j - resolution / 2;
        int fieldY = resolution / 2 - i * 2;
        
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

        if(xA != xB)
          kAB = (yA - yB)/(xA - xB);
        else
          isVerticalAB = 1; 

        if(xB != xC)
          kBC = (yB - yC)/(xB - xC);
        else
          isVerticalBC = 1;

        if(xC != xD)
          kCD = (yC - yD)/(xC - xD);
        else
          isVerticalCD = 1;

        if(xD != xA)
          kDA = (yD - yA)/(xD - xA);
        else
          isVerticalDA = 1;

        //AB duz
        if(fieldY <= fmaxf(yA, yB) && fieldY >= fminf(yA, yB) && (isVerticalAB ? xA : (fieldY - yA + kAB * xA)/kAB) >= fieldX)
          tacnostAB = 1;

        //BC duz
        if(fieldY <= fmaxf(yB, yC) && fieldY >= fminf(yB, yC) && (isVerticalBC ? xB : (fieldY - yB + kBC * xB)/kBC) >= fieldX)
          tacnostBC = 1;

        //CD duz
        if(fieldY <= fmaxf(yC, yD) && fieldY >= fminf(yC, yD) && (isVerticalCD ? xC : (fieldY - yC + kCD * xC)/kCD) >= fieldX)
          tacnostCD = 1;

        //DA duz
        if(fieldY <= fmaxf(yD, yA) && fieldY >= fminf(yD, yA) && (isVerticalDA ? xD : (fieldY - yD + kDA * xD)/kDA) >= fieldX)
          tacnostDA = 1;

        if(j == 0)
          printf("\n     ");

        if((tacnostAB + tacnostBC + tacnostCD + tacnostDA) % 2 == 1)
          printf("O");
        else
          printf(" ");

      }
    }
    printf("\n\n");
  }
  return 0;
}
