#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define NUM_OPE 24
#define A_X 2.93573192
#define A_Y 2.93573192
#define A_Z 4.62578688 
#define EPS 0.3

typedef struct {
    int    ope_type;
    double n;
    double axis[3];
    double translation[3];
} SymmetryOperation;

SymmetryOperation seitz[NUM_OPE] = {
/*1*/    {0,  1,  {1, 1, 1},    {0,   0,    0}},
/*2*/    {0,  3,  {0, 0, 1},    {0,   0,    0}},
/*3*/    {0, -3,  {0, 0, 1},    {0,   0,    0}},
/*4*/    {0,  2,  {0, 0, 1},    {0,   0,  0}},
/*5*/    {0, -6,  {0, 0, 1},    {0,   0,  0}},
/*6*/    {0,  6,  {0, 0, 1},    {0,   0,  0}},
/*7*/    {0,  2,  {0.5, sqrt(3)/2, 0},    {0,   0,    0}},
/*8*/    {0,  2,  {1, 0, 0},    {0,   0,    0}},
/*9*/    {0,  2,  {-0.5, sqrt(3)/2, 0},    {0,   0,    0}},
/*10*/    {0,  2,  {1.5, -sqrt(3)/2, 0},    {0,   0, 0}},
/*11*/    {0,  2,  {0, 1, 0},    {0,   0, 0}},
/*12*/    {0,  2,  {1.5, sqrt(3)/2, 0},    {0,   0, 0}},
/*13*/    {1,  1,  {1, 1, 1},    {0,   0,    0}},
/*14*/    {1,  3,  {0, 0, 1},    {0,   0,    0}},
/*15*/    {1, -3,  {0, 0, 1},    {0,   0,    0}},
/*16*/    {1,  2,  {0, 0, 1},    {0,   0,    0}},
/*17*/    {1, -6,  {0, 0, 1},    {0,   0,    0}},
/*18*/    {1,  6,  {0, 0, 1},    {0,   0,    0}},
/*19*/    {1,  2,  {0.5, sqrt(3)/2, 0},    {0, 0,    0}},
/*20*/    {1,  2,  {1, 0, 0},    {0, 0,    0}},
/*21*/    {1,  2,  {-0.5, sqrt(3)/2, 0},     {0, 0,    0}},
/*22*/    {1,  2,  {1.5, -sqrt(3)/2, 0},     {0,   0,    0}},
/*23*/    {1,  2,  {0, 1, 0},   {0,   0,    0}},
/*24*/    {1,  2,  {sqrt(3)/2,0.5, 0},    {0,   0,    0}}
};

double original_position[34][3] = {
{0, 0, 0},
{-2.27819362, -3.945947, 0},
{2.27819362, -3.945947, 0},
{4.55638737, 0, 0},
{2.27819368, 3.94594713, 0},
{-2.27819369, 3.94594713, 0},
{-4.55638724, 0, 0},
{2.25E-06, -2.63063263, 1.40600179},
{2.27819136, -1.31531436, 1.40600179},
{2.27819601, 1.31531437, 1.40600179},
{-2.33E-06, 2.63063276, 1.40600179},
{-2.27819137, 1.31531437, 1.40600179},
{-2.27819588, -1.31531436, 1.40600179},
{0, 0, 2.81200358},
{-2.27819362, -3.945947, 2.81200358},
{2.27819362, -3.945947, 2.81200358},
{4.55638737, 0, 2.81200358},
{2.27819368, 3.94594713, 2.81200358},
{-2.27819369, 3.94594713, 2.81200358},
{-4.55638724, 0, 2.81200358},
{2.25E-06, -2.63063263, 4.21800537},
{2.27819136, -1.31531436, 4.21800537},
{2.27819601, 1.31531437, 4.21800537},
{-2.33E-06, 2.63063276, 4.21800537},
{-2.27819137, 1.31531437, 4.21800537},
{-2.27819588, -1.31531436, 4.21800537},
{0, 0, 5.62400716},
{-2.27819362, -3.945947, 5.62400716},
{2.27819362, -3.945947, 5.62400716},
{4.55638737, 0, 5.62400716},
{2.27819368, 3.94594713, 5.62400716},
{-2.27819369, 3.94594713, 5.62400716},
{-4.55638724, 0, 5.62400716},
{-4.55638724, -2.63063263, 1.40600179},
};

int original_id[34]=
{
0,0,0,0,0,0,0,
1,1,1,1,1,1,
0,0,0,0,0,0,0,
1,1,1,1,1,1,
0,0,0,0,0,0,0,
1
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
                     int ATOM, double center[3], double center2[3], int ope_type) {
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
              transformed[i][j] = sum + center2[j] + R[j][3];        	
        }
    }
}

bool check_match(double position[][3], double transformed[][3], int ATOM, int p_id0[], int p_id1[]) {
    bool matched[ATOM];
    for(int i=0;i<ATOM;i++) matched[i] = false;

    for(int i=0;i<ATOM;i++) {
        bool found = false;
        for(int j=0;j<ATOM;j++) {
            if(matched[j]) continue;
            double dx = position[i][0] - transformed[j][0];
            double dy = position[i][1] - transformed[j][1];
            double dz = position[i][2] - transformed[j][2];
            if(dx*dx + dy*dy + dz*dz < EPS*EPS && p_id0[j] == p_id1[i]) {
           // if(dx*dx + dy*dy + dz*dz < EPS*EPS) {
                matched[j] = true;
                found = true;
			//	printf("  Match found: position[%d] <--> transformed[%d]\n", i, j);
                break;
            }
        }
        if(!found) return false;
    }
    return true;
}

int run_all_operations(int points[], double position[][3], double position2[][3], int ATOM, int valid[], int p_id0[], int p_id1[]) {
    double R[4][4], transformed[ATOM][3];
    int valid_count = 0;
    double center2[3];
	
    compute_center(center2, position2, ATOM);
    for(int k=0;k<NUM_OPE;k++) {
        SymmetryOperation op = seitz[k];
        double center[3];
        compute_center(center, position, ATOM);

        affine_transform(R, op.axis, op.n, op.ope_type, op.translation);
        apply_operation(transformed, R, position, ATOM, center, center2, op.ope_type);

        if(check_match(position2, transformed, ATOM, p_id0, p_id1)) {
            valid[valid_count++] = k;
           // printf("Operation No.%d is valid.\n", k+1);
        }
    }
    return valid_count;
}

int main(int argc, char *argv[]) {
    char filename[256];
	char s[256],s2[256];
	double temp[3];
    	double DIV=atof(argv[1]);

	FILE *fp1,*fp; 
	fp1=fopen("CONTCAR_omega", "r");
	fp=fopen("res.txt", "r");
    if (!fp1 || !fp) {
        perror("ファイルオープンに失敗しました");
        return 1;
    }


	int num1,num2;
	double boundary[9];
	fgets(s, BUFSIZ, fp1);
	fgets(s, BUFSIZ, fp1);
	fgets(s, BUFSIZ, fp1);
	sscanf(s, "%lf %lf %lf", &boundary[0], &boundary[1],&boundary[2]); 
	fgets(s, BUFSIZ, fp1);
	sscanf(s, "%lf %lf %lf", &boundary[3], &boundary[4],&boundary[5]); 
	fgets(s, BUFSIZ, fp1);
	sscanf(s, "%lf %lf %lf", &boundary[6], &boundary[7],&boundary[8]); 
	for(int i=0;i<9;i++) boundary[i]/=DIV;
	for(int i=0;i<9;i++) printf("%E\n",boundary[i]);
	
	fgets(s2, BUFSIZ, fp1);
	printf(s2);
	
	fgets(s, BUFSIZ, fp1);
	sscanf(s, "%d %d ", &num1,&num2);
	printf("%d %d \n", num1,num2);
	
	fgets(s, BUFSIZ, fp1);
	
	double mr[3*(num1+num2)];
	int p_id[(num1+num2)];
	int bd[4*75*(num1+num2)];
	double mr_abs[3*75*(num1+num2)];
	for(int i=0;i<num1+num2;i++){
		
		fgets(s, BUFSIZ, fp1);
		sscanf(s, "%lf%lf%lf%d", &mr[3*i],&mr[3*i+1],&mr[3*i+2],&p_id[i]);
		printf("%E %E %E %d\n", mr[3*i],mr[3*i+1],mr[3*i+2],p_id[i]);
		
		bd[4*i]=i;
		bd[4*i+1]=0;
		bd[4*i+2]=0;
		bd[4*i+3]=0;

		temp[0]=mr[3*i];		
		temp[1]=mr[3*i+1];		
		temp[2]=mr[3*i+2];		
		mr_abs[3*i]=temp[0]*boundary[0]+temp[1]*boundary[3]+temp[2]*boundary[6];
		mr_abs[3*i+1]=temp[0]*boundary[1]+temp[1]*boundary[4]+temp[2]*boundary[7];
		mr_abs[3*i+2]=temp[0]*boundary[2]+temp[1]*boundary[5]+temp[2]*boundary[8];
	}
	fclose(fp1);

	int ii=(num1+num2);
	for (int x2=-2; x2<3; x2++) {
        for (int y2=-2; y2<3; y2++) {
        for (int z2=-1; z2<2; z2++) {

		if(x2==0 && y2==0 && z2==0) continue;
		for(int i=0;i<num1+num2;i++){

	                bd[4*ii]=i;
	                bd[4*ii+1]=x2;
        	        bd[4*ii+2]=y2;
                	bd[4*ii+3]=z2;

                	temp[0]=mr[3*i]+x2;
                	temp[1]=mr[3*i+1]+y2;
                	temp[2]=mr[3*i+2]+z2;

                	mr_abs[3*ii]=temp[0]*boundary[0]+temp[1]*boundary[3]+temp[2]*boundary[6];
                	mr_abs[3*ii+1]=temp[0]*boundary[1]+temp[1]*boundary[4]+temp[2]*boundary[7];
                	mr_abs[3*ii+2]=temp[0]*boundary[2]+temp[1]*boundary[5]+temp[2]*boundary[8];
		ii++;
		}
	}
	}
	}
	
	
    int ATOM=4;
//    printf("Enter the number of ATOMS: ");
//   scanf("%d", &ATOM);
    int p_id0[ATOM];
	int p_id1[ATOM];

    int points[ATOM];
    double position[ATOM][3];
	double position2[ATOM][3];
	int bd_id[ATOM][3];
	for (int i=0;i<ATOM;i++) {
		points[i]=atoi(argv[i+2]);
		p_id0[i]=original_id[points[i]-1];
        for(int j=0;j<3;j++)
            position[i][j] = original_position[points[i]-1][j];
    }
	
	int valid[NUM_OPE], valid_count;
	double counter=0.0;
	for (int i1=0;i1<num1+num2;i1++) {
		
		position2[0][0] = mr_abs[3*i1];
		position2[0][1] = mr_abs[3*i1+1];
		position2[0][2] = mr_abs[3*i1+2];
		bd_id[0][0]=0;
		bd_id[0][1]=0;
		bd_id[0][2]=0;
		p_id1[0]=p_id[i1];

		for (int i2=i1+1;i2<75*(num1+num2);i2++) {
			
			
			position2[1][0] = mr_abs[3*i2];
			position2[1][1] = mr_abs[3*i2+1];
			position2[1][2] = mr_abs[3*i2+2];
			p_id1[1]=p_id[bd[4*i2]];
		
			bd_id[1][0]=bd[4*i2+1];
			bd_id[1][1]=bd[4*i2+2];
			bd_id[1][2]=bd[4*i2+3];
					
					
					
					
				
				
				for (int i3=i2+1;i3<75*(num1+num2);i3++) {

					position2[2][0] = mr_abs[3*i3];
					position2[2][1] = mr_abs[3*i3+1];
					position2[2][2] = mr_abs[3*i3+2];
					p_id1[2]=p_id[bd[4*i3]];
		
					bd_id[2][0]=bd[4*i3+1];
					bd_id[2][1]=bd[4*i3+2];
					bd_id[2][2]=bd[4*i3+3];
					
					
					for (int i4=i3+1;i4<75*(num1+num2);i4++) {

						position2[3][0] = mr_abs[3*i4];
						position2[3][1] = mr_abs[3*i4+1];
						position2[3][2] = mr_abs[3*i4+2];
						p_id1[3]=p_id[bd[4*i4]];

						bd_id[3][0]=bd[4*i4+1];
						bd_id[3][1]=bd[4*i4+2];
						bd_id[3][2]=bd[4*i4+3];
						
						valid_count = run_all_operations(points, position, position2, ATOM, valid,p_id0,p_id1);
						
						if(valid_count!=0){
						
							printf("%d %d %d %d ( %d, %d(%d,%d,%d), %d(%d,%d,%d), %d(%d,%d,%d) ) \n",i1+1,bd[4*i2]+1,bd[4*i3]+1,bd[4*i4]+1,i1+1,bd[4*i2]+1,bd[4*i2+1],bd[4*i2+2],bd[4*i2+3],bd[4*i3]+1,bd[4*i3+1],bd[4*i3+2],bd[4*i3+3],bd[4*i4]+1,bd[4*i4+1],bd[4*i4+2],bd[4*i4+3]);
	
							printf("%E %E %E\n", position2[0][0],position2[0][1],position2[0][2]);
		        			printf("%E %E %E\n", position2[1][0],position2[1][1],position2[1][2]);
		        			printf("%E %E %E\n", position2[2][0],position2[2][1],position2[2][2]);
		        			printf("%E %E %E\n", position2[3][0],position2[3][1],position2[3][2]);

							
							int bd_check[5][5][5];
							int bd_ct=0;	

							for(int i=-2;i<3;i++) 
								for(int j=-2;j<3;j++) 
									for(int k=-1;k<2;k++) bd_check[i+1][j+1][k+1]=0;
							

							for(int i_bd= 0;i_bd < ATOM;i_bd++){
								for(int i=-2;i<3;i++) 
									for(int j=-2;j<3;j++) 
										for(int k=-1;k<2;k++) 
											if(bd_id[i_bd][0]==i && bd_id[i_bd][1]==j && bd_id[i_bd][2]==k) bd_check[i+1][j+1][k+1]++;
							}
							
							printf("This cluster belongs to cell");
							for(int i=-2;i<3;i++){ 
								for(int j=-2;j<3;j++){ 
									for(int k=-1;k<2;k++){
										if(bd_check[i+1][j+1][k+1]!=0){
											bd_ct++;
											printf(" (%d, %d, %d), ",i,j,k);
										}
									}
								}
							}
							printf("\n%d cells \n\n",bd_ct);
							
							counter+=1.0/bd_ct;
							

					
						}
						
						valid_count=0;
					}

				}
				
		}
		
	}
	
	printf("counter=%f\n",counter);
    fclose(fp);
    return 0;
}
