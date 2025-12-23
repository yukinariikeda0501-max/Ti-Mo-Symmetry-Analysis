#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* 定数定義 */
#define PI 3.14159265358979323846
#define NUM_OPE 24        /* HCPの対称操作数 */
#define A_X 2.93573192
#define A_Y 2.93573192
#define A_Z 4.62578688 
#define EPS 1e-2          /* 判定精度 */

/* 構造体定義 */
typedef struct {
    int   ope_type;
    double n;
    double axis[3];
    double translation[3];
} SymmetryOperation;

/* グローバル変数（HCP座標データ） */
double original_position[18][3] = {
/*1*/   {1.461490465, 0.843791834, -1.156446725},
/*2*/   {1.461490386, -0.843791883, 1.156446725},
/*3*/   {-2.85482E-08, 1.687583764, 1.156446725},
/*4*/   {2.922980801, 1.68758377, 1.156446725},
/*5*/   {1.461490465, 0.843791834, 3.469340155},
/*6*/   {1.44518E-08, -1.687583741, -1.156446725},
/*7*/   {2.922980847, -1.687583749, -1.156446725},
/*8*/   {4.384471222, 0.843791835, -1.156446725},
/*9*/   {2.922980799, 3.375167488, -1.156446725},
/*10*/  {5.64518E-08, 3.37516749, -1.156446725},
/*11*/  {-1.461490359, 0.843791837, -1.156446725},
/*12*/  {2.922980847, -1.687583749, 3.469340155},
/*13*/  {4.384471222, 0.843791835, 3.469340155},
/*14*/  {2.922980799, 3.375167488, 3.469340155},
/*15*/  {5.64518E-08, 3.37516749, 3.469340155},
/*16*/  {-1.461490359, 0.843791837, 3.469340155},
/*17*/  {1.44518E-08, -1.687583741, 3.469340155},
/*18*/  {-1.461490358, -0.843791881, 1.156446725}
};

/* Seitz Operator Database (HCP用) */
SymmetryOperation seitz[NUM_OPE] = {
/*1*/   {0,  1,  {1, 1, 1},    {0,   0,    0}},
/*2*/   {0,  3,  {0, 0, 1},    {0,   0,    0}},
/*3*/   {0, -3,  {0, 0, 1},    {0,   0,    0}},
/*4*/   {0,  2,  {0, 0, 1},    {0,   0,  0.5}},
/*5*/   {0, -6,  {0, 0, 1},    {0,   0,  0.5}},
/*6*/   {0,  6,  {0, 0, 1},    {0,   0,  0.5}},
/*7*/   {0,  2,  {0.5, sqrt(3)/2, 0},    {0,   0,    0}},
/*8*/   {0,  2,  {1, 0, 0},    {0,   0,    0}},
/*9*/   {0,  2,  {-0.5, sqrt(3)/2, 0},    {0,   0,    0}},
/*10*/  {0,  2,  {1.5, -sqrt(3)/2, 0},    {0,   0, 0}},
/*11*/  {0,  2,  {0, 1, 0},    {0,   0, 0}},
/*12*/  {0,  2,  {1.5, sqrt(3)/2, 0},    {0,   0, 0}},
/*13*/  {1,  1,  {1, 1, 1},    {0,   0,    0}},
/*14*/  {1,  3,  {0, 0, 1},    {0,   0,    0}},
/*15*/  {1, -3,  {0, 0, 1},    {0,   0,    0}},
/*16*/  {1,  2,  {0, 0, 1},    {0,   0,    0}},
/*17*/  {1, -6,  {0, 0, 1},    {0,   0,    0}},
/*18*/  {1,  6,  {0, 0, 1},    {0,   0,    0}},
/*19*/  {1,  2,  {0.5, sqrt(3)/2, 0},    {0, 0,    0}},
/*20*/  {1,  2,  {1, 0, 0},    {0, 0,    0}},
/*21*/  {1,  2,  {-0.5, sqrt(3)/2, 0},     {0, 0,    0}},
/*22*/  {1,  2,  {1.5, -sqrt(3)/2, 0},     {0,   0,    0.5}},
/*23*/  {1,  2,  {0, 1, 0},    {0,   0,    0.5}},
/*24*/  {1,  2,  {sqrt(3)/2,0.5, 0},    {0,   0,    0.5}}
};

/* 重心計算 */
void compute_center(double center[3], double position[][3], int ATOM) {
    for(int j=0;j<3;j++) center[j] = 0;
    for(int i=0;i<ATOM;i++)
        for(int j=0;j<3;j++)
            center[j] += position[i][j];
    for(int j=0;j<3;j++) center[j] /= ATOM;
}

/* アフィン変換行列の計算 (HCP版は並進を含むため affine_transform を使用) */
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

        /* ローカル変数にコピーしてスケーリング */
        double t[3];
        t[0] = translation[0] * A_X;
        t[1] = translation[1] * A_Y;
        t[2] = translation[2] * A_Z;

        R[0][0]=x*x*v+c;     R[0][1]=x*y*v-z*s; R[0][2]=z*x*v+y*s; R[0][3]=t[0];
        R[1][0]=x*y*v+z*s;   R[1][1]=y*y*v+c;   R[1][2]=y*z*v-x*s; R[1][3]=t[1];
        R[2][0]=z*x*v-y*s;   R[2][1]=y*z*v+x*s; R[2][2]=z*z*v+c;   R[2][3]=t[2];
        R[3][0]=0; R[3][1]=0; R[3][2]=0; R[3][3]=1;
    }
}

/* 座標変換適用 */
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

/* 補助関数: パターンlからTi/Moの文字列を生成する */
void get_config_strings(long l, int ATOM, int points[], char* ti_str, char* mo_str) {
    ti_str[0] = '\0';
    mo_str[0] = '\0';
    
    for(int b=0; b<ATOM; b++) {
        char temp[16];
        int atom_num = points[b]; 
        
        if ( ((l >> b) & 1) == 0 ) {
            /* Bit 0 -> Ti */
            if(strlen(ti_str) > 0) strcat(ti_str, ",");
            sprintf(temp, "%d", atom_num);
            strcat(ti_str, temp);
        } else {
            /* Bit 1 -> Mo */
            if(strlen(mo_str) > 0) strcat(mo_str, ",");
            sprintf(temp, "%d", atom_num);
            strcat(mo_str, temp);
        }
    }
}

/* alpha_l の計算 */
void calculate_alpha_l(double position[][3], int points[], int ATOM, FILE *fp) {
    fprintf(fp, "\n=== Calculation of Multiplicity (alpha_l) for HCP ===\n");
    printf("\n=== Calculation of Multiplicity (alpha_l) for HCP ===\n");

    /* 1. 対称操作（Valid Operations）と置換マップの作成 */
    int valid_ops[NUM_OPE];
    int permutation_map[NUM_OPE][ATOM];
    int valid_count = 0;
    
    double center[3];
    compute_center(center, position, ATOM);
    double transformed[ATOM][3];
    double R[4][4];

    for(int k=0; k<NUM_OPE; k++) {
        SymmetryOperation op = seitz[k];
        affine_transform(R, op.axis, op.n, op.ope_type, op.translation);
        apply_operation(transformed, R, position, ATOM, center, op.ope_type);

        bool is_valid_op = true;
        int current_map[ATOM];

        for(int i=0; i<ATOM; i++) {
            bool found_dest = false;
            for(int j=0; j<ATOM; j++) {
                double d2 = 0;
                for(int d=0; d<3; d++) d2 += pow(transformed[i][d] - position[j][d], 2);
                if(d2 < EPS*EPS) {
                    current_map[i] = j;
                    found_dest = true;
                    break;
                }
            }
            if(!found_dest) {
                is_valid_op = false;
                break;
            }
        }

        if(is_valid_op) {
            valid_ops[valid_count] = k;
            for(int i=0; i<ATOM; i++) permutation_map[valid_count][i] = current_map[i];
            valid_count++;
        }
    }

    fprintf(fp, "Symmetry of Cluster N_A(G): Order %d\n", valid_count);
    
    /* 2. 全組み合わせの生成と alpha_l の計算 */
    long num_combinations = 1L << ATOM; // 2^ATOM
    bool *processed = (bool*)calloc(num_combinations, sizeof(bool));
    if(!processed) {
        printf("Memory allocation failed.\n");
        return;
    }

    /* ヘッダー出力 */
    fprintf(fp, "\n%-4s | %-25s | %-25s | %s\n", "ID", "Ti Atoms", "Mo Atoms", "Alpha_l");
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    printf("%-4s | %-25s | %-25s | %s\n", "ID", "Ti Atoms", "Mo Atoms", "Alpha_l");

    int orbit_id_counter = 0;

    for(long l = 0; l < num_combinations; l++) {
        if(processed[l]) continue;

        orbit_id_counter++;
        long *orbit_patterns = (long*)malloc(sizeof(long) * valid_count);
        int orbit_size = 0;

        /* 軌道探索 */
        for(int v=0; v<valid_count; v++) {
            long l_prime = 0;
            for(int i=0; i<ATOM; i++) {
                if((l >> i) & 1) {
                    int target_pos = permutation_map[v][i];
                    l_prime |= (1L << target_pos);
                }
            }

            bool seen_in_orbit = false;
            for(int o=0; o<orbit_size; o++) {
                if(orbit_patterns[o] == l_prime) {
                    seen_in_orbit = true;
                    break;
                }
            }
            if(!seen_in_orbit) {
                orbit_patterns[orbit_size++] = l_prime;
                processed[l_prime] = true;
            }
        }

        /* 代表パターンの出力 */
        char ti_str[256], mo_str[256];
        get_config_strings(l, ATOM, points, ti_str, mo_str);

        fprintf(fp, "%-4d | %-25s | %-25s | %d\n", orbit_id_counter, ti_str, mo_str, orbit_size);
        printf("%-4d | %-25s | %-25s | %d\n", orbit_id_counter, ti_str, mo_str, orbit_size);

        /* 消去された等価パターン（軌道内の他のパターン）を出力（矢印なし） */
        for(int k=0; k<orbit_size; k++) {
            long l_prime = orbit_patterns[k];
            if(l_prime == l) continue; /* 自分自身はスキップ */

            get_config_strings(l_prime, ATOM, points, ti_str, mo_str);
            
            /* IDやAlpha_lは空欄 */
            fprintf(fp, "     | %-25s | %-25s |\n", ti_str, mo_str);
            printf("     | %-25s | %-25s |\n", ti_str, mo_str);
        }
        
        fprintf(fp, "     |                           |                           |\n");

        free(orbit_patterns);
    }

    free(processed);
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "Total distinct configurations (orbits): %d\n", orbit_id_counter);
}

int main(void) {
    char filename[256];
    printf("出力ファイル名を入力してください（例: hcp_alpha_result.txt）：");
    scanf("%255s", filename);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("ファイルオープン失敗");
        return 1;
    }

    int ATOM;
    printf("原子の数を入力してください: ");
    scanf("%d", &ATOM);
    if(ATOM > 18) {
        printf("原子数が多すぎます（Max 18）。処理を中断します。\n");
        fclose(fp);
        return 1;
    }

    int points[ATOM];
    double position[ATOM][3];

    printf("各原子の番号を入力してください（HCPは 1-18 まで）。\n");
    
    for (int i=0; i<ATOM; i++) {
        printf("Atom %d No: ", i+1);
        scanf("%d", &points[i]);
        if (points[i] < 1 || points[i] > 18) {
            printf("Error: Invalid atom number.\n");
            fclose(fp);
            return 1;
        }
        for(int j=0; j<3; j++)
            position[i][j] = original_position[points[i]-1][j];
    }

    calculate_alpha_l(position, points, ATOM, fp);

    fclose(fp);
    printf("処理が完了しました。%s を確認してください。\n", filename);
    return 0;
}