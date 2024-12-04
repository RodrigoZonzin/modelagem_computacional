#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

float dt  = 0.01;
float days = 100.0;
float dx = 0.05;
float dy = 0.05;
int x_size = 50;
int y_size = 1;

float x_diffusion(float*** u, int x, int y, int t){
    float dif = 0;
    if (x == 0){
        dif = (u[t][x+1][y] - u[t][x][y]) /(dx*dx);
    }
    else if (x == (x_size -1)){
        dif = (u[t][x-1][y] - u[t][x][y])/(dx*dx);
    }
    else{
        dif = (u[t][x+1][y] - 2*u[t][x][y] + u[t][x-1][y])/(dx*dx); 	
    }
    return dif;
}

void imprimeDados(FILE* arq, float*** u, int tn){
    for (int x = 0; x < x_size; x++){
        fprintf(arq, "%f ", u[tn][x][0]);
    }
    fprintf(arq,"\n");      
    fprintf(arq,"--");
    fprintf(arq,"\n");
}

int main(int argc, char* argv[]) {
    float r = 0.25,
        a = 0.05, 
        m = 0.15,
        Dh = 0.1,
        Dp = 0.05;

    int numIteracoesGravadas = 10;    
    int interval = (int)days/numIteracoesGravadas;

    float ***h = (float***)malloc(2*sizeof(float**));
    float ***p = (float***)malloc(2*sizeof(float**));

    for(int t=0; t < 2; t++) {
        h[t] = (float**)malloc(x_size*sizeof(float*));
        p[t] = (float**)malloc(x_size*sizeof(float*));

        for(int x=0; x < x_size; x++){
            h[t][x] = (float*)malloc(y_size*sizeof(float));
            p[t][x] = (float*)malloc(y_size*sizeof(float));

            for(int y=0; y < y_size; y++) {
                h[t][x][y] =  0;
                p[t][x][y] =  0;
            }
        }
    }

    h[0][12][0] = 30.0;
    p[0][30][0] = 5.0;
    //write space file
    FILE* espaco = fopen("x.dat", "w");
    float inc = 0.0;
    for (int x=0; x < x_size; x++){
        fprintf(espaco, "%.1f\n", inc);
        inc += (float)dx;
    }
    fclose(espaco);

    espaco = fopen("y.dat", "w");
    inc = 0.0;
    for (int y = 0; y < y_size; y++){
        fprintf(espaco, "%.1f\n", inc);
        inc += (float)dy;
    }
    fclose(espaco);
    //write time file
    float time = 0.0;
    FILE* tempo = fopen("t.dat", "w");
    for (int t=0; t <= (int)days; t++){
        if ((t % interval) == 0){
            fprintf(tempo, "%.0f \n", time);
        }
        time = time + (float)dt;//(float)dt;
    }
    fclose(tempo);

    FILE* hfile = fopen("presa.dat", "w");
    FILE* pfile = fopen("predador.dat", "w");

    FILE* htempo = fopen("presa_t.dat", "w");
    FILE* ptempo = fopen("predador_t.dat", "w");

    int t = 0;
    int tn = 1;
    float totalPresas = 0, totalPredadores = 0;

    #pragma omp parallel for num_threads(6) 
    for(int x=0; x < x_size; x++){
        for(int y=0; y < y_size; y++)
        {
            totalPresas += h[t][x][y];
            totalPredadores += p[t][x][y];
        }
    }

    imprimeDados(hfile,(float***)h,t);
    imprimeDados(pfile,(float***)p,t);

    for(int step = 1; step < days+1; step++){
        if(step % 2 == 0){
            tn = 0;
            t = 1;
        }
        else{
            tn = 1;
            t = 0;
        }
        totalPresas = 0;
        totalPredadores = 0;
        #pragma omp parallel for num_threads(6)
        for (int x=0; x < x_size; x++){
            for (int y=0; y < y_size; y++){
                //EDP presa
                h[tn][x][y] = ( r*h[t][x][y] - a*h[t][x][y]*p[t][x][y] + Dh*x_diffusion(h,x,y,t) )*dt + h[t][x][y];

                //EDP Predador
                p[tn][x][y] = ( a*h[t][x][y]*p[t][x][y] - m*p[t][x][y] + Dp*x_diffusion(p,x,y,t) )*dt + p[t][x][y];

                #pragma omp critical
                {
                    totalPresas += h[tn][x][y];
                    totalPredadores += p[tn][x][y];
                }
            }
        }

        if((step % interval) == 0){
            imprimeDados(hfile,(float***)h,tn);
            imprimeDados(pfile,(float***)p,tn);
            fprintf(htempo, "%f\n", totalPresas);
            fprintf(ptempo, "%f\n", totalPredadores);
        }
    }

    for(int t=0; t < 2; t++){
        for(int x=0; x < x_size; x++){
            free(h[t][x]);
            free(p[t][x]);
        }
        free(h[t]);
        free(p[t]);
    }
    free(h);
    free(p);

    fclose(hfile);
    fclose(pfile);
    fclose(htempo);
    fclose(ptempo);
    return 0;
}
