#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#define N_FIELDS 3
#define T_FINAL 50.0
#define dt 0.0001

int x_size;
int y_size = 1;

float initial_damage = 0.05, 
    beta_td = 0.5, 
    beta_ch = 0.4,
    mu_ch = 0.6,
    s_m1 = 1.,
    keq_m1 = 1.0,
    mu_m1 = 0.,
    D_ch = 0.8,
    D_m1 = 0.1;

/*enum CellType {
    Cells = 0,
    Fluid = 1,
    Tissue = 2,
    Nutrient = 3,
    Velocity,
    Pressure,
    Other
};*/

enum CellType {
    TissueDamage,
    Ch,
    M1
};

typedef struct TField1D {
    enum CellType cell_type;
    float *cur_values;
    float *next_values;
    float x_min, x_max;
    float y_min, y_max;
    float x_dim, y_dim;
    float dx, dy;
    char name[100];
} Field1D;

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
    float dx = p_field->dx, dy = p_field->dy;

    if (x > 0){
        float pm1 = get_current_value(p_field, x - 1, y);
        float chm1 = get_current_value(ch_field, x - 1, y);
        if(ch -  chm1 > 0)        
            flux_left = -((ch - chm1) * pm1)/dx;
        else        
            flux_left = -((ch - chm1) * p)/dx;
    }
    if(x < x_size - 1){
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
    float dx = field->dx;
    float dy = field->dy;

    if (x == 0) {
        uxp1 = get_current_value(field, x + 1, y);
        diff_x = (uxp1 - ux) /(dx*dx);
    }
    else if (x == (x_size - 1)) {        
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
    float coeff = v*(dt/field->dx);
    float advec = 0;
    
    if (v > 0) {
        if (i > 0)
            advec = coeff *(get_current_value(field, i, j) - get_current_value(field, i - 1, j));
    }
    else if (v < 0) {
        if (i < x_size - 1)
            advec = coeff *(get_current_value(field, i + 1, j) - get_current_value(field, i, j));
    }
    return advec;
}

void save_current_field(Field1D *field, float t, float *x_values) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = y_size - 1; j >= 0; j--) {            
        for (int x = 0; x < x_size; x++)
            fprintf(fp, "%.5f,0,%.5f\n", x_values[x], get_current_value(field, x, j));
        //fprintf(fp,"\n");
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

void save_field(Field1D *field, float t, float *x_values) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = y_size - 1; j >= 0; j--) {            
        for (int x = 0; x < x_size; x++)
            fprintf(fp, "%.5f,0,%.5f\n", x_values[x], get_next_value(field, x, j));
        fprintf(fp,"\n");
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

Field1D* create_field(enum CellType cell_type, char *name, float x_min, float x_max, float dx, float dy){
    Field1D *field = (Field1D*) calloc(1, sizeof(Field1D));
    strcpy(field->name, name);
    field->cell_type = cell_type;
    field->x_min = x_min;
    field->x_max = x_max;
    field->x_dim = (x_max - x_min);    
    field->dx = dx;
    field->dy = dy;
    field->cur_values = (float*) calloc((field->x_dim /dx) + 1, sizeof(float));
    field->next_values = (float*) calloc((field->x_dim /dx) + 1, sizeof(float));
    return field;
}

Field1D** create_fields(enum CellType cell_type[], char *name[], float x_min, float x_max, float dx, float dy){
    Field1D** fields = (Field1D**) calloc(N_FIELDS, sizeof(Field1D*));
    for (int i = 0; i < N_FIELDS; i++){
        Field1D *field = (Field1D*) calloc(1, sizeof(Field1D));
        strcpy(field->name, name[i]);
        field->cell_type = cell_type[i];
        field->x_min = x_min;
        field->x_max = x_max;
        field->x_dim = (x_max - x_min);    
        field->dx = dx;
        field->dy = dy;
        field->cur_values = (float*) calloc((field->x_dim /dx) + 1, sizeof(float));
        field->next_values = (float*) calloc((field->x_dim /dx) + 1, sizeof(float));
        fields[i] = field;
    }
    
    return fields;
}

int binary_search(float *array, float x, int low, int high) {
    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (fabsf(array[mid] - x) < pow(10,-5))
            return mid;

        if (array[mid] < x)
            low = mid + 1;

        else
            high = mid - 1;
    }
    return -1;
}

float derivation_x(Field1D* field, int i, int j){
    float central_diff = 0;
    if (i == 0){
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i, j))/ (field->dx);
    }
    else if (i == x_size - 1){
        central_diff = (get_current_value(field, i, j) - get_current_value(field, i - 1, j))/ (field->dx);
    }
    else {
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i - 1, j))/ (2*field->dx);
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

float get_fields_sum(Field1D** fields, int indexes[], int len){
    float sum = 0.0;
    for (int j = 0; j < y_size; j++) {
        for (int i = 0; i < x_size; i++) {
            sum += get_fields_sum_at_pos(fields, i, j, indexes, len);
        }
    }
    return sum;
}

/* 
float g_right_border_medium(int x, float** w, int t){
    return g((w[t][x+1] + w[t][x])/2);
}

float g_left_border_medium(int x, float** w, int t){
    return g((w[t][x] + w[t][x-1])/2);
}

float u_right_border(float** u, int x, int t){
    return u[t][x+1] - u[t][x];
}

float u_left_border(float** u, int x, int t){
    return u[t][x] - u[t][x-1];
}

float u_right_border_medium(float** u, int x, int t){
    return (u[t][x+1] + u[t][x])/2;
}

float u_left_border_medium(float** u, int x, int t){
    return (u[t][x] + u[t][x-1])/2;
}

float non_linear_diffusion(float** u, int x, float** w, int t) {
    if(x == 0)
        return (gRightBorderMedium(x, w, t)*(uRightBorder(u, x, t))) /(deltaX*deltaX);
    else if(x == size_-1)
        return - (gLeftBorderMedium(x, w, t)*(uLeftBorder(u, x, t))) /(deltaX*deltaX);
    else
        return (gRightBorderMedium(x, w, t)*uRightBorder(u, x, t) -  gLeftBorderMedium(x, w, t)*uLeftBorder(u, x, t) ) /(deltaX*deltaX);
}
*/

//TO DO: methods for inserting random values in the field

//x, y -> i, j
int get_index_coord(float w, float *array, int len){ 
    return binary_search(array, w, 0, len - 1); 
}

int main(){    
    int nsteps = (int)((T_FINAL+dt)/dt);
    printf("steps: %d\n", nsteps);   
    int interval = nsteps/10;
    
    FILE *tFile = fopen("results/t.csv", "w");
    for (float t = 0; t < T_FINAL+dt; t += dt)
        fprintf(tFile, "%.6lf\n", t);
    fclose(tFile);    

    enum CellType field_types[N_FIELDS] = {TissueDamage, Ch, M1};
    char *field_names[N_FIELDS] = {"results/TissueDamage", "results/Ch", "results/M1"};

    float x_min = 0.0, x_max = 10.0, dx = 0.05, dy = 0.05;
    Field1D** fields = create_fields(field_types, field_names, x_min, x_max, dx, dy);
    Field1D* td_field = fields[0];
    Field1D* ch_field = fields[1];
    Field1D* m1_field = fields[2];
    
    x_size = (td_field->x_dim/td_field->dx) + 1;
    y_size = 1;
    
    float* x_values = (float*) calloc(x_size, sizeof(float));
    int i =0;
    FILE *xFile = fopen("results/x.csv", "w"); 
    for (float x = 0; x < td_field->x_dim + td_field->dx; x += td_field->dx) {
        fprintf(xFile, "%.3lf\n", x);
        x_values[i++] = x;
    }
    fclose(xFile);    

    /*
    //tests 
    printf("index: %d\n", get_index_coord(0.0, x_values, x_size));
    printf("index: %d\n", get_index_coord(0.2, x_values, x_size));
    printf("index: %d\n", get_index_coord(1.0, x_values, x_size));*/
    
    //set_current_value(c_field, 0.0, 0.0, 0.8);
    
    for (int j = 0; j < y_size; j++) {
        #pragma omp parallel for num_threads(6) 
        for (int i = 0; i < x_size; i++) {
            if ((i >= x_size/2.0 - 5) && (i <= x_size/2.0 + 5))
                set_current_value(td_field, i, j, initial_damage);
            else 
                set_current_value(td_field, i, j, 0.0);
            set_current_value(ch_field, i, j, 0.0);
            set_current_value(m1_field, i, j, 0.0);
        }
    }
    //set_current_value(c_field, x_size/2, y_size/2, 5*c_0);
   
    for (int i =0; i < N_FIELDS; i++)
        save_current_field(fields[i], 0, x_values);

    float t = 0;
    for(int step = 1; step <= nsteps; ++step) {
                    
        for (int j = 0; j < y_size; j++) {
            #pragma omp parallel for num_threads(6) 
            for (int i = 0; i < x_size; i++) {
                float td = get_current_value(td_field, i, j);
                float ch = get_current_value(ch_field, i, j);
                float m1 = get_current_value(m1_field, i, j);
                
                float td_next = beta_td*m1;
                float ch_next = beta_ch*(td + m1) - mu_ch*ch + D_ch*diffusion(ch_field, i, j);
                float m1_next = (s_m1*ch)/(keq_m1 + ch) -mu_m1*m1 -beta_td*m1 + D_m1*diffusion(m1_field, i, j);
                
                set_next_value(td_field, i, j,  td + (td_next * dt));                
                set_next_value(ch_field, i, j, ch + (ch_next * dt));
                set_next_value(m1_field, i, j, m1 + (m1_next * dt));
                
            }
        }

        if (step % interval == 0) {
            for (int i =0; i < N_FIELDS; i++)
                save_field(fields[i], t, x_values);            
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
    free(x_values);
    return 0; 
}
