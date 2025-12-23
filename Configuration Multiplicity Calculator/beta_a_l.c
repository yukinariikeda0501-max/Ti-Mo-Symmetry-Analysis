#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* 定数定義 */
#define PI 3.14159265358979323846
#define NUM_OPE 48
#define EPS 1e-3

/* 構造体定義 */
typedef struct {
    int   ope_type;
    double n;
    double axis[3];
    double translation[3];
} SymmetryOperation;

/* グローバル変数（座標データ） */
double original_position[9][3] = {
/*1*/   {0,0,0},
/*2*/   {1,0,0},
/*3*/   {1,1,0},
/*4*/   {0,1,0},
/*5*/   {0.5,0.5,0.5},
/*6*/   {0,0,1},
/*7*/   {1,0,1},
/*8*/   {1,1,1},
/*9*/   {0,1,1},
};

/* Seitz Operator Database (一部抜粋) */
SymmetryOperation seitz[NUM_OPE] = {
    {0, 1, {1,1,1}, {0,0,0}}, {0, 2, {0,0,1}, {0,0,0}}, {0, 2, {0,1,0}, {0,0,0}}, {0, 2, {1,0,0}, {0,0,0}},
    {0, 3, {1,1,1}, {0,0,0}}, {0, 3, {-1,1,-1}, {0,0,0}}, {0, 3, {1,-1,-1}, {0,0,0}}, {0, 3, {-1,-1,1}, {0,0,0}},
    {0, -3, {1,1,1}, {0,0,0}}, {0, -3, {1,-1,-1}, {0,0,0}}, {0, -3, {-1,-1,1}, {0,0,0}}, {0, -3, {-1,1,-1}, {0,0,0}},
    {0, 2, {1,1,0}, {0,0,0}}, {0, 2, {1,-1,0}, {0,0,0}}, {0, -4, {0,0,1}, {0,0,0}}, {0, 4, {0,0,1}, {0,0,0}},
    {0, -4, {1,0,0}, {0,0,0}}, {0, 2, {0,1,1}, {0,0,0}}, {0, 2, {0,1,-1}, {0,0,0}}, {0, 4, {1,0,0}, {0,0,0}},
    {0, 4, {0,1,0}, {0,0,0}}, {0, 2, {1,0,1}, {0,0,0}}, {0, -4, {0,1,0}, {0,0,0}}, {0, 2, {-1,0,1}, {0,0,0}},
    {1, 1, {1,1,1}, {0,0,0}}, {1, 2, {0,0,1}, {0,0,0}}, {1, 2, {0,1,0}, {0,0,0}}, {1, 2, {1,0,0}, {0,0,0}},
    {1, 3, {1,1,1}, {0,0,0}}, {1, 3, {-1,1,-1}, {0,0,0}}, {1, 3, {1,-1,-1}, {0,0,0}}, {1, 3, {-1,-1,1}, {0,0,0}},
    {1, -3, {1,1,1}, {0,0,0}}, {1, -3, {1,-1,-1}, {0,0,0}}, {1, -3, {-1,-1,1}, {0,0,0}}, {1, -3, {-1,1,-1}, {0,0,0}},
    {1, 2, {1,1,0}, {0,0,0}}, {1, 2, {1,-1,0}, {0,0,0}}, {1, -4, {0,0,1}, {0,0,0}}, {1, 4, {0,0,1}, {0,0,0}},
    {1, -4, {1,0,0}, {0,0,0}}, {1, 2, {0,1,1}, {0,0,0}}, {1, 2, {0,1,-1}, {0,0,0}}, {1, 4, {1,0,0}, {0,0,0}},
    {1, 4, {0,1,0}, {0,0,0}}, {1, 2, {1,0,1}, {0,0,0}}, {1, 4, {-1,0,1}, {0,0,0}} 
};

/* 重心計算 */
void compute_center(double center[3], double position[][3], int ATOM) {
    for(int j=0; j<3; j++) center[j] = 0;
    for(int i=0; i<ATOM; i++)
        for(int j=0; j<3; j++)
            center[j] += position[i][j];
    for(int j=0; j<3; j++) center[j] /= ATOM;
}

/* 変換行列計算 */
void compute_transform_matrix(double R[4][4], SymmetryOperation op) {
    double norm = sqrt(op.axis[0]*op.axis[0] + op.axis[1]*op.axis[1] + op.axis[2]*op.axis[2]);
    if (norm < EPS) {
        for (int i=0; i<4; i++) for(int j=0; j<4; j++) R[i][j] = (i==j) ? 1 : 0;
    } else {
        double x = op.axis[0]/norm, y = op.axis[1]/norm, z = op.axis[2]/norm;
        double theta = 2*PI/fabs(op.n);
        double c = cos(theta);
        double s = (op.n > 0 ? 1 : -1) * sin(theta);
        double v = 1 - c;
        R[0][0]=x*x*v+c;   R[0][1]=x*y*v-z*s; R[0][2]=z*x*v+y*s; R[0][3]=0;
        R[1][0]=x*y*v+z*s; R[1][1]=y*y*v+c;   R[1][2]=y*z*v-x*s; R[1][3]=0;
        R[2][0]=z*x*v-y*s; R[2][1]=y*z*v+x*s; R[2][2]=z*z*v+c;   R[2][3]=0;
        R[3][0]=0;         R[3][1]=0;         R[3][2]=0;         R[3][3]=1;
    }
}

/* 座標変換適用 */
void apply_operation(double transformed[][3], double R[4][4], double position[][3],
                     int ATOM, double center[3], int ope_type) {
    for(int i=0; i<ATOM; i++) {
        double shifted[3] = {
            position[i][0] - center[0],
            position[i][1] - center[1],
            position[i][2] - center[2]
        };
        for(int j=0; j<3; j++) {
            double sum = 0;
            for(int m=0; m<3; m++) sum += R[j][m] * shifted[m];
            if(ope_type == 1) sum = -sum; /* Inversion */
            transformed[i][j] = sum + center[j];
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
    fprintf(fp, "\n=== Calculation of Multiplicity (alpha_l) ===\n");
    printf("\n=== Calculation of Multiplicity (alpha_l) ===\n");

    /* 1. 対称操作（Valid Operations）と置換マップの作成 */
    int valid_ops[NUM_OPE];
    int permutation_map[NUM_OPE][ATOM];
    int valid_count = 0;
    
    double center[3];
    compute_center(center, position, ATOM);
    double transformed[ATOM][3];
    double R[4][4];

    for(int k=0; k<NUM_OPE; k++) {
        compute_transform_matrix(R, seitz[k]);
        apply_operation(transformed, R, position, ATOM, center, seitz[k].ope_type);

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

        /* 消去された等価パターン（軌道内の他のパターン）を出力（矢印削除） */
        for(int k=0; k<orbit_size; k++) {
            long l_prime = orbit_patterns[k];
            if(l_prime == l) continue; /* 自分自身はスキップ */

            get_config_strings(l_prime, ATOM, points, ti_str, mo_str);
            
            /* IDやAlpha_lは空欄にし、位置合わせのみ行う */
            fprintf(fp, "     | %-25s | %-25s |\n", ti_str, mo_str);
            printf("     | %-25s | %-25s |\n", ti_str, mo_str);
        }
        
        /* 視認性のための区切り線（オプション） */
        fprintf(fp, "     |                           |                           |\n");

        free(orbit_patterns);
    }

    free(processed);
    fprintf(fp, "--------------------------------------------------------------------------------\n");
    fprintf(fp, "Total distinct configurations (orbits): %d\n", orbit_id_counter);
}

int main(void) {
    char filename[256];
    printf("出力ファイル名を入力してください（例: alpha_result.txt）：");
    scanf("%255s", filename);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("ファイルオープン失敗");
        return 1;
    }

    int ATOM;
    printf("原子の数を入力してください: ");
    scanf("%d", &ATOM);
    if(ATOM > 20) {
        printf("原子数が多すぎます。処理を中断します。\n");
        fclose(fp);
        return 1;
    }

    int points[ATOM];
    double position[ATOM][3];

    printf("各原子の番号を入力してください（例: 5, 6, 7...）。\n");
    
    for (int i=0; i<ATOM; i++) {
        printf("Atom %d No: ", i+1);
        scanf("%d", &points[i]);
        if (points[i] < 1 || points[i] > 9) {
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
