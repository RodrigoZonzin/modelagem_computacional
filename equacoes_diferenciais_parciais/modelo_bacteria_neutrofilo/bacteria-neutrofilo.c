#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#define N_FIELDS 3
#define T_FINAL 100.0
#define dt 0.0001

int y_size = 1;

float r = 0.2,
    lambda_nb = 0.05,
    D_b = 0.1,
    beta_ch = 0.8,
    mu_ch = 0.4,
    D_ch = 0.8,
    pmin = 0.01,
    pmax = 1.0,
    keq = 5.0,
    m_n = 0.5,
    D_n = 0.25,
    NBlood = 500.0, 
    v_n = 0.; 

//(-1, 0) faz plotar em todos os pontos do dominio
float source_points[4] = {10.0, 15.0, 30.0, 32.0};

enum CellType {
    B,
    N,
    Ch
};

typedef struct TAxis {
    float min;
    float max; 
    float dim;
    int size;
    float delta;
    float *values;
    int *indexes;
} Axis; 

typedef struct TField1D {
    enum CellType cell_type;
    float *cur_values;
    float *next_values;
    Axis x_axis;
    Axis y_axis;
    char name[100];
    float *source_points;
    float *sink_points;
} Field1D;

float new_precision(float v, int i){
    return floorf(powf(10,i)*v)/powf(10,i);
}

int get_index(Field1D *field, int i, int j){
    return i;// + j*(x_size);
}

float get_current_value(Field1D *field, int i, int j){
    return field->cur_values[get_index(field, i, j)];
}

float get_next_value(Field1D *field, int i, int j){
    return field->next_values[get_index(field, i, j)];
}

void set_current_value(Field1D *field, int i, int j, float v){
    field->cur_values[get_index(field, i, j)] = v;
}

void set_next_value(Field1D *field, int i, int j, float v){
    field->next_values[get_index(field, i, j)] = v;
}

float chemotaxis(Field1D *p_field, Field1D *ch_field, int x, int y){
    float resx = 0, resy = 0, flux_left = 0, flux_right = 0;
    float p = get_current_value(p_field, x, y);
    float ch = get_current_value(ch_field, x, y);
    float dx = p_field->x_axis.delta, dy = 0;

    if (x > 0){
        float pm1 = get_current_value(p_field, x - 1, y);
        float chm1 = get_current_value(ch_field, x - 1, y);
        if(ch -  chm1 > 0)        
            flux_left = -((ch - chm1) * pm1)/dx;
        else        
            flux_left = -((ch - chm1) * p)/dx;
    }
    if(x < p_field->x_axis.size - 1){
        float pp1 = get_current_value(p_field, x + 1, y);
        float chp1 = get_current_value(ch_field, x + 1, y);
        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dx;         
        else        
            flux_right = (chp1 - ch)*pp1 /dx;        
    }
    resx = (flux_left + flux_right)/dx;

    flux_left = 0, flux_right = 0;
    if (y > 0) {
        
        float pm1 = get_current_value(p_field, x, y - 1);
        float chm1 = get_current_value(ch_field, x, y - 1);

        if(ch - chm1 > 0)        
            flux_left = -(ch - chm1)*pm1 /dy;         
        else
            flux_left = -(ch - chm1)*p /dy;
    }
    if(y < y_size-1) {
        
        float pp1 = get_current_value(p_field, x, y + 1);
        float chp1 = get_current_value(ch_field, x, y + 1);

        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dy; 
        else 
            flux_right = (chp1 - ch)*pp1 /dy; 
    }
    resy = (flux_left + flux_right)/dy;

    //return resx + resy; 
    return resx;
}

float diffusion(Field1D *field, int x, int y){
    float ux = get_current_value(field, x, y), uxp1, uxm1, diff_x, diff_y;
    float dx = field->x_axis.delta;
    float dy = 0;

    if (x == 0) {
        uxp1 = get_current_value(field, x + 1, y);
        diff_x = (uxp1 - ux) /(dx*dx);
    }
    else if (x == (field->x_axis.size - 1)) {        
        uxm1 = get_current_value(field, x - 1, y);
        diff_x = (uxm1 - ux)/(dx*dx);
    }
    else {
        uxp1 = get_current_value(field, x + 1, y);   
        uxm1 = get_current_value(field, x - 1, y);
        diff_x = (uxp1 - 2*ux + uxm1)/(dx*dx);
    }
    
    if (y_size == 1) {
        diff_y = 0;
    }
    else {
        float uy = get_current_value(field, x, y), uyp1, uym1;
        uyp1 = get_current_value(field, x, y + 1);
        uym1 = get_current_value(field, x, y - 1);
        if (y == 0) {
            diff_y = (uyp1 - uy) /(dy*dy);
        }
        else if (y == (y_size - 1)) {
            diff_y = (uym1 - uy)/(dy*dy);
        }
        else {    
            diff_y = (uyp1 - 2*uy + uym1)/(dy*dy);
        }
    }
    
    return diff_x + diff_y;
}

float advection(float v, Field1D* field, int i, int j){
    float coeff = v*(1/field->x_axis.delta);
    float advec = 0;
    
    if (v > 0) {
        if (i > 0)
            advec = coeff *(get_current_value(field, i, j) - get_current_value(field, i - 1, j));
    }
    else if (v < 0) {
        if (i < field->x_axis.size - 1)
            advec = coeff *(get_current_value(field, i + 1, j) - get_current_value(field, i, j));
    }
    return advec;
}

void save_current_field(Field1D *field, float t) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = y_size - 1; j >= 0; j--) {            
        for (int x = 0; x < field->x_axis.size; x++)
            fprintf(fp, "%.5f,0,%.5f\n", field->x_axis.values[x], get_current_value(field, x, j));
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

void save_field(Field1D *field, float t) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = y_size - 1; j >= 0; j--) {            
        for (int x = 0; x < field->x_axis.size; x++)
            fprintf(fp, "%.5f,0,%.5f\n", field->x_axis.values[x], get_next_value(field, x, j));
        fprintf(fp,"\n");
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

Field1D* create_field(enum CellType cell_type, char *name, float x_min, float x_max, float dx){
    Field1D *field = (Field1D*) calloc(1, sizeof(Field1D));
    strcpy(field->name, name);
    field->cell_type = cell_type;
    Axis x_axis;
    x_axis.min = x_min;
    x_axis.max = x_max;
    x_axis.dim = (x_max - x_min);    
    x_axis.delta = dx;    
    x_axis.size = (x_axis.dim/dx) + 1;
    printf("size: %d\n", x_axis.size);
    field->cur_values = (float*) calloc(x_axis.size, sizeof(float));
    field->next_values = (float*) calloc(x_axis.size, sizeof(float));
    x_axis.values = (float*) calloc(x_axis.size, sizeof(float));
    x_axis.indexes = (int*) calloc(x_axis.size, sizeof(int));

    int i = 0;
    for (float x = 0; x < x_axis.dim; x += x_axis.delta) {
        x_axis.values[i] = new_precision(x, 2.0);
        x_axis.indexes[i] = i;
        i++;
    }    
    field->x_axis = x_axis;

    return field;
}

Field1D** create_fields(enum CellType cell_type[], char *name[], float x_min, float x_max, float dx){
    Field1D** fields = (Field1D**) calloc(N_FIELDS, sizeof(Field1D*));
    for (int i = 0; i < N_FIELDS; i++){
        Field1D *field = (Field1D*) calloc(1, sizeof(Field1D));
        strcpy(field->name, name[i]);
        fields[i] = create_field(cell_type[i], name[i], x_min, x_max, dx);
    }
    return fields;
}

int binary_search(float x, float *array, int low, int high) {    
    while (low <= high) {
        int mid = (low + high) / 2;
        //printf("mid = %d\n", mid);
        //printf("array[mid] = %f\n", array[mid]);

        if (fabsf(array[mid] - x) < powf(10.0,-3.0)) {
            return mid;
        }

        if (x > array[mid])
            low = mid + 1;

        else
            high = mid - 1;
    }
    return -1;
}

float derivation_x(Field1D* field, int i, int j){
    float central_diff = 0;
    if (i == 0){
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i, j))/ (field->x_axis.delta);
    }
    else if (i == field->x_axis.size - 1){
        central_diff = (get_current_value(field, i, j) - get_current_value(field, i - 1, j))/ (field->x_axis.delta);
    }
    else {
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i - 1, j))/ (2*field->x_axis.delta);
    }
    return central_diff;
}

float get_fields_sum_at_pos(Field1D** fields, int i, int j, int indexes[], int len){
    float sum  = 0.0;
    for (int i = 0; i < len; i++){
        sum += get_current_value(fields[indexes[i]], i, j);
    }
    return sum;
}


int get_index_coord(float x, float *array, int len){ 
    return binary_search(x, array, 0, len - 1); 
}

/*float Pmin(Field1D* field, float x[2], int i, int j){
    if (x[0] == -1) 
        return 0.01;
    
    if ((get_index_coord(x[0], field->x_axis.values, field->x_axis.size) == i) 
        || (get_index_coord(x[1], field->x_axis.values, field->x_axis.size) == i))   
        return 0.01;
    else
        return 0.0;
}*/

float Pmin(Field1D* field, float x[2], int i, int j){
    if (x[0] == -1) 
        return pmin;
    
    if ( (x[0]/field->x_axis.delta) == i 
        || (x[1]/field->x_axis.delta) == i)   
        return pmin;
    else
        return pmin;
}

float Pmax(Field1D* field, float x[2], int i, int j){
    if (x[0] == -1) 
        return pmax;
    //printf("%d\n", get_index_coord(x[1], field->x_axis.values, field->x_axis.size));

    if ( (x[0]/field->x_axis.delta) == i 
        || (x[1]/field->x_axis.delta) == i) 
        return pmax;
    else 
        return pmax;
}

int main(){    
    int nsteps = (int)((T_FINAL+dt)/dt);
    printf("steps: %d\n", nsteps);   
    int interval = nsteps/10;
    
    FILE *tFile = fopen("results/t.csv", "w");
    for (float t = 0; t < T_FINAL+dt; t += dt)
        fprintf(tFile, "%.6lf\n", t);
    fclose(tFile);    

    enum CellType field_types[N_FIELDS] = {B, N, Ch};
    char *field_names[N_FIELDS] = {"results/B", "results/N", "results/Ch"};

    float x_min = 0.0, x_max = 150.0, dx = 0.1;
    Field1D** fields = create_fields(field_types, field_names, x_min, x_max, dx);
    Field1D* b_field = fields[0];
    Field1D* n_field = fields[1];
    Field1D* ch_field = fields[2];
    
    FILE *xFile = fopen("results/x.csv", "w"); 
    for (float x = 0; x < b_field->x_axis.dim ; x +=  b_field->x_axis.delta) {
        fprintf(xFile, "%.2lf\n", x);        
    }
    fclose(xFile);    
    int x_size = b_field->x_axis.size;

    //condicao inicial 
    for (int j = 0; j < y_size; j++) {
        #pragma omp parallel for num_threads(6) 
        for (int i = 0; i < x_size; i++) {
            if ((i >= (x_size/2.0 - 5)) && (i <= (x_size/2.0 + 5)))
                set_current_value(b_field, i, j, 1.0);
            else 
                set_current_value(b_field, i, j, 0.0);
            set_current_value(ch_field, i, j, 0.0);
            set_current_value(n_field, i, j, 0.0);
        }
    }
   
    for (int i =0; i < N_FIELDS; i++)
        save_current_field(fields[i], 0);    

    float t = 0;
    for(int step = 1; step <= nsteps; ++step) {
                    
        for (int j = 0; j < y_size; j++) {

            #pragma omp parallel for num_threads(6) 
            for (int i = 0; i < x_size; i++) {

                float b = get_current_value(b_field, i, j);
                float ch = get_current_value(ch_field, i, j);
                float n = get_current_value(n_field, i, j);
                                
                float b_next = r*b/(1 + 0.1*b) - lambda_nb*n*b + D_b*diffusion(b_field, i, j); 

                float n_next = (Pmin(n_field,source_points,i,j) 
                    + Pmax(n_field,source_points,i,j)*ch/(keq + ch))* NBlood 
                    - m_n*n + D_n*diffusion(n_field, i, j) - advection(v_n, n_field, i, j);

                float ch_next = beta_ch*b*n - mu_ch*ch + D_ch*diffusion(ch_field, i, j);
                
                set_next_value(b_field, i, j,  b + (b_next * dt));                
                set_next_value(ch_field, i, j, ch + (ch_next * dt));
                set_next_value(n_field, i, j, n + (n_next * dt));
                
            }
        }

        if (step % interval == 0) {
            for (int i =0; i < N_FIELDS; i++)
                save_field(fields[i], t);            
        }

        for (int i =0; i < N_FIELDS; i++) {
            float *aux = fields[i]->cur_values;
            fields[i]->cur_values = fields[i]->next_values;
            fields[i]->next_values = aux;
        }

        t += dt;  
    }

    for (int i = 0; i < N_FIELDS; i++) {
        free(fields[i]->cur_values);
        free(fields[i]->next_values);
        free(fields[i]);
    }
    free(fields);
    return 0; 
}
