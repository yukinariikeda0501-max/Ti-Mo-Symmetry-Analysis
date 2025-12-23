#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define NUM_OPE 48
#define A_X 2.93573192
#define A_Y 2.93573192
#define A_Z 4.62578688 
#define EPS 1e-2

typedef struct {
    int    ope_type;
    double n;
    double axis[3];
    double translation[3];
} SymmetryOperation;

SymmetryOperation seitz[NUM_OPE] = {
/*1*/    {0, 1, {1,1,1}, {0,0,0}}, 
/*2*/	{0, 2, {0,0,1}, {0,0,0}}, 
/*3*/	{0, 2, {0,1,0}, {0,0,0}}, 
/*4*/	{0, 2, {1,0,0}, {0,0,0}},
/*5*/    {0, 3, {1,1,1}, {0,0,0}}, 
/*6*/	{0, 3, {-1,1,-1}, {0,0,0}}, 
/*7*/	{0, 3, {1,-1,-1}, {0,0,0}}, 
/*8*/	{0, 3, {-1,-1,1}, {0,0,0}},
/*9*/    {0, -3, {1,1,1}, {0,0,0}}, 
/*10*/	{0, -3, {1,-1,-1}, {0,0,0}}, 
/*11*/	{0, -3, {-1,-1,1}, {0,0,0}}, 
/*12*/	{0, -3, {-1,1,-1}, {0,0,0}},
/*13*/    {0, 2, {1,1,0}, {0,0,0}}, 
/*14*/	{0, 2, {1,-1,0}, {0,0,0}}, 
/*15*/	{0, -4, {0,0,1}, {0,0,0}}, 
/*16*/	{0, 4, {0,0,1}, {0,0,0}},
/*17*/    {0, -4, {1,0,0}, {0,0,0}}, 
/*18*/	{0, 2, {0,1,1}, {0,0,0}}, 
/*19*/	{0, 2, {0,1,-1}, {0,0,0}}, 
/*20*/	{0, 4, {1,0,0}, {0,0,0}},
/*21*/    {0, 4, {0,1,0}, {0,0,0}}, 
/*22*/	{0, 2, {1,0,1}, {0,0,0}}, 
/*23*/	{0, -4, {0,1,0}, {0,0,0}}, 
/*24*/	{0, 2, {-1,0,1}, {0,0,0}},
/*25*/    {1, 1, {1,1,1}, {0,0,0}}, 
/*26*/	{1, 2, {0,0,1}, {0,0,0}}, 
/*27*/	{1, 2, {0,1,0}, {0,0,0}}, 
/*28*/	{1, 2, {1,0,0}, {0,0,0}},
/*29*/    {1, 3, {1,1,1}, {0,0,0}}, 
/*30*/	{1, 3, {-1,1,-1}, {0,0,0}}, 
/*31*/	{1, 3, {1,-1,-1}, {0,0,0}}, 
/*32*/	{1, 3, {-1,-1,1}, {0,0,0}},
/*33*/    {1, -3, {1,1,1}, {0,0,0}}, 
/*34*/	{1, -3, {1,-1,-1}, {0,0,0}}, 
/*35*/	{1, -3, {-1,-1,1}, {0,0,0}}, 
/*36*/	{1, -3, {-1,1,-1}, {0,0,0}},
/*37*/    {1, 2, {1,1,0}, {0,0,0}}, 
/*38*/	{1, 2, {1,-1,0}, {0,0,0}}, 
/*39*/	{1, -4, {0,0,1}, {0,0,0}}, 
/*40*/	{1, 4, {0,0,1}, {0,0,0}},
/*41*/    {1, -4, {1,0,0}, {0 ,0,0}}, 
/*42*/	{1, 2, {0,1,1}, {0,0,0}}, 
/*43*/	{1, 2, {0,1,-1}, {0,0,0}}, 
/*44*/	{1, 4, {1,0,0}, {0,0,0}},
/*45*/    {1, 4, {0,1,0}, {0,0,0}}, 
/*46*/	{1, 2, {1,0,1}, {0,0,0}}, 
/*47*/	{1, -4, {0,1,0}, {0,0,0}}, 
/*48*/	{1, 2, {-1,0,1}, {0,0,0}}
};

double original_position[9][3] = {
/*1*/	{0,0,0},
/*2*/	{1,0,0},
/*3*/	{1,1,0},
/*4*/	{0,1,0},
/*5*/	{0.5,0.5,0.5},
/*6*/	{0,0,1},
/*7*/	{1,0,1},
/*8*/	{1,1,1},
/*9*/	{0,1,1},
};

void compute_center(double center[3], double position[][3], int ATOM) {
    for(int j=0;j<3;j++) center[j] = 0;
    for(int i=0;i<ATOM;i++)
        for(int j=0;j<3;j++)
            center[j] += position[i][j];
    for(int j=0;j<3;j++) center[j] /= ATOM;
}

void affine_transform(double R[4][4], double axis[3], double n,
                      int ope_type, double translation[3]) {
    double norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (norm < EPS) {
        for (int i=0;i<4;i++) for(int j=0;j<4;j++) R[i][j] = (i==j)?1:0;
    } else {
        double x = axis[0]/norm, y = axis[1]/norm, z = axis[2]/norm;
        double theta = 2*PI/fabs(n);
        double c = cos(theta);
        double s = (n>0?1:-1)*sin(theta);
        double v = 1-c;

        translation[0] *= A_X;
        translation[1] *= A_Y;
        translation[2] *= A_Z;

        R[0][0]=x*x*v+c;     R[0][1]=x*y*v-z*s; R[0][2]=z*x*v+y*s; R[0][3]=translation[0];
        R[1][0]=x*y*v+z*s;   R[1][1]=y*y*v+c;   R[1][2]=y*z*v-x*s; R[1][3]=translation[1];
        R[2][0]=z*x*v-y*s;   R[2][1]=y*z*v+x*s; R[2][2]=z*z*v+c;   R[2][3]=translation[2];
        R[3][0]=0; R[3][1]=0; R[3][2]=0; R[3][3]=1;
    }
}

void apply_operation(double transformed[][3], double R[4][4], double position[][3],
                     int ATOM, double center[3], int ope_type) {
    for(int i=0;i<ATOM;i++) {
        double shifted[3] = {
            position[i][0] - center[0],
            position[i][1] - center[1],
            position[i][2] - center[2]
        };
        for(int j=0;j<3;j++) {
            double sum = 0;
            for(int m=0;m<3;m++) sum += R[j][m] * shifted[m];
            if(ope_type == 1) sum = -sum;
            transformed[i][j] = sum + center[j] + R[j][3];
        }
    }
}

bool check_match(double position[][3], double transformed[][3], int ATOM) {
    bool matched[ATOM];
    for(int i=0;i<ATOM;i++) matched[i] = false;

    for(int i=0;i<ATOM;i++) {
        bool found = false;
        for(int j=0;j<ATOM;j++) {
            if(matched[j]) continue;
            double dx = position[i][0] - transformed[j][0];
            double dy = position[i][1] - transformed[j][1];
            double dz = position[i][2] - transformed[j][2];
            if(dx*dx + dy*dy + dz*dz < EPS*EPS) {
                matched[j] = true;
                found = true;
				printf("  Match found: position[%d] <--> transformed[%d]\n", i, j);
                break;
            }
        }
        if(!found) return false;
    }
    return true;
}

int run_all_operations(int points[], double position[][3], int ATOM, int valid[], FILE *fp) {
    double R[4][4], transformed[ATOM][3];
    int valid_count = 0;

    for(int k=0;k<NUM_OPE;k++) {
        SymmetryOperation op = seitz[k];
        double center[3];
        compute_center(center, position, ATOM);

        affine_transform(R, op.axis, op.n, op.ope_type, op.translation);
        apply_operation(transformed, R, position, ATOM, center, op.ope_type);

        if(check_match(position, transformed, ATOM)) {
            valid[valid_count++] = k;
            fprintf(fp, "=> Operation No.%d is valid.\n", k+1);
        }
    }
    return valid_count;
}

void recheck_valid_operations(int points[], double position[][3], int ATOM, int valid[], int valid_count, FILE *fp) {
    int sub_ATOM;
    printf("Enter the number of ATOMS: ");
    scanf("%d", &sub_ATOM);

    int sub_points[sub_ATOM];
    double sub_position[sub_ATOM][3];
    for (int i=0;i<sub_ATOM;i++) {
        printf("Enter the point number for atom %d: ", i+1);
        scanf("%d", &sub_points[i]);
        for(int j=0;j<3;j++)
            sub_position[i][j] = original_position[sub_points[i]-1][j];
    }

    double R[4][4], transformed[ATOM][3];
    double sub_center[3];
    compute_center(sub_center, sub_position, sub_ATOM);

    int local_valid = 0;
    fprintf(fp, "\n=== Rechecking around Atom ");
    for (int i=0;i<sub_ATOM;i++) fprintf(fp, "%d ", sub_points[i]);
    fprintf(fp, "===\n");

    for(int i=0;i<valid_count;i++) {
        int k = valid[i];
        SymmetryOperation op = seitz[k];

        affine_transform(R, op.axis, op.n, op.ope_type, op.translation);
        apply_operation(transformed, R, position, ATOM, sub_center, op.ope_type);

        if(check_match(position, transformed, ATOM)) {
            fprintf(fp, "✓ Still valid (Op. No.%d)\n", k+1);
            local_valid++;
        } else {
            fprintf(fp, "✗ Invalid now (Op. No.%d)\n", k+1);
        }
    }

    fprintf(fp, "→ Valid operations: %d / %d\n", valid_count, local_valid);
}

int main(void) {
    char filename[256];
    printf("出力ファイル名を入力してください（例: result.txt）：");
    scanf("%s", filename);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("ファイルオープンに失敗しました");
        return 1;
    }

    int ATOM;
    printf("Enter the number of ATOMS: ");
    scanf("%d", &ATOM);

    int points[ATOM];
    double position[ATOM][3];
    for (int i=0;i<ATOM;i++) {
        printf("Enter the point number for atom %d: ", i+1);
        scanf("%d", &points[i]);
        for(int j=0;j<3;j++)
            position[i][j] = original_position[points[i]-1][j];
    }

    int valid[NUM_OPE], valid_count;
    valid_count = run_all_operations(points, position, ATOM, valid, fp);

    fprintf(fp, "\nTotal valid operations: %d\n", valid_count);
    for(int i=0;i<valid_count;i++) 
        fprintf(fp, "Valid No.%d\n", valid[i]+1);

    recheck_valid_operations(points, position, ATOM, valid, valid_count, fp);
    fclose(fp);
    return 0;
}
