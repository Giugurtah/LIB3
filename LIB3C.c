//cc -fPIC -shared -o LIB3C.so LIB3C.c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h> 

#define EPSILON 1e-5

//*TREE FUNCTIONS (SLBA is in a separate section of this file)
//Helper functions
void mat_print_double(int I, int J, double mat[I][J]){ //Function used to print a double matrix
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }   
}

void mat_rand(int I, int J, double mat[I][J]){ //Function to randomly generate a double matrix
    for (int i=0; i<I; i++){
        for (int j=0; j<J; j++){
            mat[i][j] = ((double)rand()/RAND_MAX)*5;
        }
    }
}

void mat_tras(int I, int J, double mat[I][J], double mat_res[J][I]){ //Function that given a matrix returns its transpose
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            mat_res[j][i] = mat[i][j];
        }
    }
}

void mat_tras_ptr(int I, int J, double** mat, double mat_res[J][I]){ //Function that given a matrix returns its transpose
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            mat_res[j][i] = mat[i][j];
        }
    }
}

void prod_matrix(int I, int J, int Z, double mat1[I][J], double mat2[J][Z], double mat_res[I][Z]){ //Given two matrices this function returns their product
    double sum;
    for(int i=0; i<I; i++){
        for(int z=0; z<Z; z++){
            sum = 0.0;
            for(int j=0; j<J; j++){
                sum += ((double)mat1[i][j] * (double)mat2[j][z]);
            }
            mat_res[i][z] = sum;
        }
    }
}

void prod_mat_vect(int I, int J, double mat[I][J], double arr[J], double arr_res[I]){
    double sum;
    for(int i=0; i<I; i++){
        sum = 0;
        for(int j=0; j<J; j++){
            sum += mat[i][j]*arr[j];
        }
        arr_res[i] = sum;
    }
}

void prod_mat_ptr_vect(int I, int J, double **mat, double arr[J], double arr_res[I]){
    double sum;
    for(int i=0; i<I; i++){
        sum = 0;
        for(int j=0; j<J; j++){
            sum += mat[i][j]*arr[j];
        }
        arr_res[i] = sum;
    }
}

void prod_matrix_ptr(int I, int J, int Z, double** mat1, double mat2[J][Z], double mat_res[I][Z]){ //Given two matrices this function returns their product
    double sum;
    for(int i=0; i<I; i++){
        for(int z=0; z<Z; z++){
            sum = 0.0;
            for(int j=0; j<J; j++){
                sum += ((double)mat1[i][j] * (double)mat2[j][z]);
            }
            mat_res[i][z] = sum;
        }
    }
}

double calc_tr(int I, double mat[I][I]){ //Function that given a matrix returns its trace
    double sum = 0;
    for(int i=0; i<I; i++){
        sum += mat[i][i];
    }
    return sum;
}

void conv_UA(int I, int J, double mat1[I][J], double mat2[I][J]){ //Function used, given U, to get A 
    double sum;
    for(int i=0; i<I; i++){
        sum = 0.0;
        for(int j=0; j<J; j++){
            sum += exp(mat1[i][j]);
        }
        for(int j=0; j<J; j++){
            mat2[i][j] = exp(mat1[i][j])/sum;
        }
    }
}

void conv_ZB(int I, int J, double mat1[I][J], double mat2[I][J]){ //Function used, given Z, to get B
    double sum;
    for(int j=0; j<J; j++){
        sum = 0.0;
        for(int i=0; i<I; i++){
            sum += exp(mat1[i][j]);
        }
        for(int i=0; i<I; i++){
            mat2[i][j] = exp(mat1[i][j])/sum;
        }
    }
}

void DX_calc(int I, int J, double mat1[I][J], double mat2[J][J]){ //Given the transpose of F this function returns DX
    double sum = 0.0;
    for(int j=0; j<J; j++){
        for(int i=0; i<J; i++){
            mat2[i][j] = 0.0;
        }
    }
    for(int j=0; j<J; j++){
        sum = 0.0;
        for(int i=0; i<I; i++){
            sum += mat1[i][j];
        }
        mat2[j][j] = sum;
    }
}

void conv_df_dU(int Jb, int K, double mat1[Jb][K], double mat2[Jb][K], double mat_res[Jb][K]){ //Given df_dA this function returns df_dU
    double sum = 0;
    double check;
    for (int i=0; i<Jb; i++){
        for (int j=0; j<K; j++){
            sum = 0;
            for (int k=0; k<K; k++){
                check = 0.0;
                if(j == k){
                    check = 1.0;
                }
                sum += mat1[i][k]*(check*mat2[i][k] - mat2[i][k]*mat2[i][j]);

            }
            mat_res[i][j] = sum;
        }
    }
}

void conv_df_dZ(int Jy, int K, double mat1[Jy][K], double mat2[Jy][K], double mat_res[Jy][K]){ //Given df_dB this function returns df_dZ
    double sum = 0;
    double check;
    for (int i=0; i<Jy; i++){
        for (int j=0; j<K; j++){
            sum = 0;
            for (int k=0; k<Jy; k++){
                check = 0.0;
                if(i == k){
                    check = 1.0;
                }
                sum += mat1[k][j]*(check*mat2[k][j] - mat2[k][j]*mat2[i][j]);
            }
            mat_res[i][j] = sum;
        }
    }
}

void update_UZ(int J, int K, double rate, double mat1[J][K], double mat2[J][K]){ //Given U(or Z) and df_dU(or df_dZ) this function updates U(or Z) 
    for(int i=0; i<J; i++){
        for(int j=0; j<K; j++){
            mat1[i][j] = mat1[i][j] - rate*mat2[i][j];
        }
    }
}

int check_stop(int I, int J, int K, double mat1[I][K], double mat2[J][K]){ //Given df_dU and df_dZ this function returns one if the derivatives are sufficiently close to zero, zero otherwise
    int stop = 1;
    for(int k=0; k<K; k++){
        for(int j=0; j<J; j++){
            if(mat2[j][k] > 0.001 || mat2[j][k] < -0.001){
                stop = 0;
            }
        }
        for(int i=0; i<I; i++){
            if(mat1[i][k] > 0.001 || mat1[i][k] < -0.001){
                stop = 0;
            }
        }
    }
    return stop;
}

void incr_S(int I, int S[I], int S_N[I]){ // Given S and S_N as boolean arrays, this function increases S by one
    for(int i=0; i<I; i++){
        if(S[i] == 0){
            S[i] = 1;
            S_N[i] = 0;
            return;
        } else {
            S[i] = 0;
            S_N[i] = 1;
        }
    }
}

//Impurity and predicatbality functions 
double gini_impurity(int I, int J, int S[I], double** F){ //Function used to evaluate the gini impurity
    double sumI, sumJ = 0.0;
    double N = 0.0;
    for(int j=0; j<J; j++){
        for(int i=0; i<I; i++){
            N += (double)S[i]*F[i][j];
        }
    }
    for(int j=0; j<J; j++){
        sumI = 0.0;
        for(int i=0; i<I; i++){
            sumI += (double)S[i]*F[i][j]/N;
        }
        sumJ += pow(sumI, 2);
    }
    return N*(1-sumJ);
}

double gini_impurity_twoing(int I, int J, int Sx[I], int Sy[J], int Sy_N[I], double** F){ //Function to evaluate the gini impurity in the twoing
    double sumI_1 = 0.0, sumI_0 = 0.0;
    double N = 0.0;

    for(int j=0; j<J; j++){
        for(int i=0; i<I; i++){
            N += (double)Sx[i]*F[i][j];
        }
    }
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            sumI_1 += (double)Sy[j]*Sx[i]*F[i][j]/N;
            sumI_0 += (double)Sy_N[j]*Sx[i]*F[i][j]/N;
        }
    }

    return N*(1-pow(sumI_1,2)-pow(sumI_0,2));
}

double gpi_c(int I, int J, double** F){ //Given F this function returns the gpi
    int S[I];
    double aux = 0.0;
    double Y_impurity = 0.0, x_impurity = 0.0, gpi = 0.0;

    for(int i=0; i<I; i++){
        S[i] = 1;
    }
    Y_impurity = gini_impurity(I, J, S, F);

    for(int i=0; i<I; i++){
        S[i] = 0;
    }
    for(int i=0; i<I; i++){
        S[i] = 1;
        x_impurity += gini_impurity(I, J, S, F);
        S[i] = 0;
    }
    gpi = (Y_impurity - x_impurity)/Y_impurity;
    return gpi;
}

//Function models
void twoStage_c(int I, int J, double** F, double* S_best, double** alpha, double** beta){ //This function, given F a contingency matrix, executes a two stage algorithm
    // Algo variables
    int S_N[I], S[I];
    int count = 1;
    double parent_impurity = 0.0, left_impurity = 0.0, right_impurity = 0.0, ppi;

    // S and S_Negative inizialization
    for(int i=0; i<I; i++){
        S[i] = 0;
        S_N[i] = 1;
    }
    // Parent impurity
    parent_impurity = gini_impurity(I, J, S_N, F);

    while(count< pow(2,I-1)){
        incr_S(I, S, S_N);
        // Children impurity
        left_impurity = gini_impurity(I,J,S,F);
        right_impurity = gini_impurity(I,J,S_N,F);
        // PPI improvement
        ppi = (parent_impurity - left_impurity - right_impurity)/(parent_impurity);
        if(ppi > S_best[0]){
            S_best[0] = ppi;
            for(int i=1; i<I+1; i++){
                S_best[i] = (double)S[i-1]*ppi;
            }
        }        
        // Iteration increase
        count += 1;
    }
}

void twoing_c(int I, int J, double** F, double* S_best, double** alpha, double** beta){ //This function, given F a contingency matrix, executes a twoing algorithm
    // Algo variables
    int Sx[I], Sx_N[I];
    int Sy[J], Sy_N[J];
    int count_x = 1, count_y = 1;
    double parent_impurity = 0.0, left_impurity = 0.0, right_impurity = 0.0, ppi;

    // Sx, Sy, Sx_Negative and Sy_Negative inizialization
    for(int i=0; i<I; i++){
        Sx[i] = 0;
        Sx_N[i] = 1;
    }
    for(int j=0; j<J; j++){
        Sy[j] = 0;
        Sy_N[j] = 1;
    }

    // Parent impurity
    parent_impurity = gini_impurity(I, J, Sx_N, F);

    while(count_y < pow(2, J-1)){
        incr_S(J, Sy, Sy_N);
        while(count_x < pow(2, I-1)){
            incr_S(I, Sx, Sx_N);
            // Children impurity
            left_impurity = gini_impurity_twoing(I,J,Sx,Sy,Sy_N,F);
            right_impurity = gini_impurity_twoing(I,J,Sx_N,Sy,Sy_N,F);
            // PPI improvement
            ppi = (parent_impurity - left_impurity - right_impurity)/(parent_impurity);
            if(ppi > S_best[0]){
                left_impurity = gini_impurity(I,J,Sx,F);
                right_impurity = gini_impurity(I,J,Sx_N,F);
                ppi = (parent_impurity - left_impurity - right_impurity)/(parent_impurity);
                
                S_best[0] = ppi;
                for(int i=1; i<I+1; i++){
                    S_best[i] = (double)Sx[i-1]*ppi;
                }
            }        
            // Iteration increase
            count_x += 1;
        }
        count_x = 1;
        // Reset Sx and Sx_Negative
        for(int i=0; i<I; i++){
            Sx[i] = 0;
            Sx_N[i] = 1;
        }
        count_y += 1;
    }
}

void lba_execute(int Jy, int Jb, double** F, double** alpha, double** beta){ //This function is called by lba_c to perform the LBA 
    int K=2, count=0, stop=0;
    double r = 0.1, d_stop = 0.0000001, f_stop=0.005;
    double f=0, f_p, d;
    double F_T[Jy][Jb], U[Jb][K], Z[Jy][K], DX[Jb][Jb], A[Jb][K], B[Jy][K];
    double df_dA[Jb][K], df_dB[Jy][K], df_dU[Jb][K], df_dZ[Jy][K];
    
    double B_T[K][Jy];
    double A_T[K][Jb];
    double DX_A[Jb][K], DX_A_BT_B[Jb][K], F_B[Jb][K];
    double FT_A[Jy][K], B_AT_DX_A[Jy][K];
    double BT_B[K][K];
    double B_AT[Jy][Jb], B_AT_DX[Jy][Jb];
    double A_BT[Jb][Jy]; 
    double FT_A_BT[Jy][Jy], B_AT_DX_A_BT[Jy][Jy];
    
    // Randomly generated matrices U and Z
    mat_rand(Jb, K, U);
    mat_rand(Jy, K, Z);
    
    for(int i=0; i<Jb; i++){
        for(int j=0; j<Jy; j++){
            F_T[j][i] = F[i][j];
        }
    }

    while(count < 100000 && stop==0){
    // while(count<10000 && stop==0){
        f_p = f;
        f = 1;

        //A and B are evaluated using U and Z
        conv_UA(Jb, K, U, A);
        conv_ZB(Jy, K, Z, B);

        //2tr[(F_t)(A)(B_t)]
        mat_tras(Jy, K, B, B_T);
        prod_matrix(Jb, K, Jy, A, B_T, A_BT); //Product between A and transpose of B
        prod_matrix(Jy, Jb, Jy, F_T, A_BT, FT_A_BT); //Product between transpose of F and A_BT 
        f -= 2 * calc_tr(Jy, FT_A_BT);

        //tr[(B)(A_T)(Dx)(A)(B_T)]
        mat_tras(Jb, K, A, A_T);
        prod_matrix(Jy, K, Jb, B, A_T, B_AT); //Product between B and transpose of A
        DX_calc(Jy, Jb, F_T, DX); //DX as F*1
        prod_matrix(Jy, Jb, Jb, B_AT, DX, B_AT_DX); 
        prod_matrix(Jy, Jb, Jy, B_AT_DX, A_BT, B_AT_DX_A_BT); //Product between B, transpose of A, Dx, A and transpose of B
        f += calc_tr(Jy, B_AT_DX_A_BT);

        //df_dA
        prod_matrix(Jb, Jb, K, DX, A, DX_A); //Product between DX and A
        prod_matrix(K, Jy, K, B_T, B, BT_B); //Product between transpose of B and B
        prod_matrix(Jb, K, K, DX_A, BT_B, DX_A_BT_B); //Product between DX, A, transpose of B and B
        prod_matrix_ptr(Jb, Jy, K, F, B, F_B); //Product between F and B 
        for(int i=0; i<Jb; i++){
            for(int j=0; j<K; j++){
                df_dA[i][j] = 2.0*DX_A_BT_B[i][j] - 2.0*F_B[i][j]; //df_dA
            }
        }

        //df_dB
        prod_matrix(Jy, Jb, K, F_T, A, FT_A); //Product between transpose of F and A 
        prod_matrix(Jy, Jb, K, B_AT, DX_A, B_AT_DX_A); //Product between B, transpose of A, DX and A
        for(int i=0; i<Jy; i++){
            for(int j=0; j<K; j++){
                df_dB[i][j] = 2.0*B_AT_DX_A[i][j] - 2.0*FT_A[i][j]; //df_dB
            }
        }

        //df_dU and df_dZ
        conv_df_dU(Jb, K, df_dA, A, df_dU);
        conv_df_dZ(Jy, K, df_dB, B, df_dZ);

        //U and Z update
        update_UZ(Jb, K, r ,U, df_dU);
        update_UZ(Jy, K, r, Z, df_dZ);

        count++;
        d = f_p - f;

        stop = check_stop(Jb, Jy, K, df_dU, df_dZ);
        if((d < d_stop && d > (-1)*d_stop)||f<f_stop){
            stop = 1;
        }
    }

    //TODO Cancellare
    printf("Il valore finale di f è: %f. L'ultimo decr è stato di: %f. Il numero totale di iterazioni è stato: %d\n", f, d, count);
    printf("La matrice A è:\n");
    mat_print_double(Jb, K, A);
    printf("La matrice B è:\n");
    mat_print_double(Jy, K, B);

    for(int i=0; i<Jb; i++){
        for(int k=0; k<K; k++){
            alpha[i][k] = A[i][k];
        }
    }
    for(int j=0; j<Jy; j++){
        for(int k=0; k<K; k++){
            beta[j][k] = B[j][k];
        }
    }
}
void lba_c(int I, int J, double** F, double* S_best, double** alpha, double** beta){ //This function, given F a contingency matrix, executes an LBA algorithm
    // Algo variables
    double sum;
    double parent_impurity;
    double left_impurity, right_impurity;
    double ppi;
    int S[I], S_N[I];
    int sumS=0, sumS_N=0; 

    // S and S_Negative inizialization
    for(int i=0; i<I; i++){
        S[i] = 1;
        S_N[i] = 0;
    }
    parent_impurity = gini_impurity(I,J,S,F);

    // LBA is executed
    lba_execute(J, I, F, alpha, beta);

    // Eval S and S_Negative based on the LBA results
    for(int i=0; i<I; i++){
        if(alpha[i][0] < alpha[i][1]){
            S[i] = 0;
            S_N[i] = 1;
            sumS += S[i];
            sumS_N += S_N[i];
        } 
    }

    // If the resulting split is such that all modalities go to just one child node the ppi is equal to 0
    if(sumS == I || sumS_N == I){
        ppi = 0;
    } else {
        left_impurity = gini_impurity(I,J,S,F);
        right_impurity = gini_impurity(I,J,S_N,F);

        // PPI eval
        ppi = (parent_impurity - left_impurity - right_impurity)/(parent_impurity);
    }

    S_best[0] = ppi;
    for(int i=1; i<I+1; i++){
        S_best[i] = S[i-1]*ppi;
    }
}

//*CATEGORIZATION
void kplusplus_init(int K, double clusters[K], int I, double X[I]){ //This function is used to initialize the centroids for a KMeans algorithm 
    double sum, dist, dist_min;
    double distances[I];

    srand(time(NULL));

    //First centroid init
    clusters[0] = X[rand() % I];

    //K-1 remaining centroids
    for(int k=1; k<K; k++){
        sum = 0.0;

        //DIstances evaluation
        for(int j=0; j<I; j++){
            //Find distance from closest centroid for each point of X 
            dist_min = DBL_MAX;
            for(int c=0; c<k; c++){
                dist = pow(X[j] - clusters[c],2);
                if(dist < dist_min){
                    dist_min = dist;
                }
            }
            distances[j] = dist_min*dist_min;
            sum += distances[j];
        }

        //Distances-proportional new centroid initialization
        double r = ((double)rand() / RAND_MAX) * sum;
        double cumulative = 0.0;
        for (int j = 0; j<I; j++) {
            cumulative += distances[j];
            if (cumulative >= r) {
                clusters[k] = X[j];
                break;
            }
        }
    }
    return;
}

int check_N(int I, int K, int minN, int labels[I]){ //This function checks whether all clusters have a lower size than minN or not
    int counts[K];
    int check = 1;

    for(int c=0; c<K; c++){
        counts[c] = 0; 
    }
    for(int i=0; i<I; i++){
        counts[labels[i]] += 1;
    }
    for(int c=0; c<K; c++){
        if(counts[c] > minN){
            check = 0;
        }
    }
    return check;
}

double kmeans_c(int I, int K, double X[I], int labels[I], double var){ //This function executes the Kmeans algorithm
    double centers[K], new_centers[K];
    double min_dist, dist;
    int count[K];
    int changed = 1;
    int counter = 0, MAX_ITER = 1000;

    //Centroids initialization
    kplusplus_init(K, centers, I, X);
    for(int c=0; c<K; c++){
        new_centers[c] = 0.0;
        count[c] = 0;
    }

    //Itearative execution
    while(counter<MAX_ITER && changed == 1){
        //Assign centroids and update centroids
        printf("Counter exe: %d\n", counter);
        for(int i=0; i<I; i++){
            min_dist = DBL_MAX;
            for(int c=0; c<K; c++){
                dist = fabs(centers[c] - X[i]);
                if(dist < min_dist){
                    min_dist = dist;
                    labels[i] = c;
                }
            }
            new_centers[labels[i]] += X[i];
            count[labels[i]] += 1;
        }

        //Check convergence and reset cycle
        changed = 0;
        for(int c=0; c<K; c++){
            new_centers[c] = (double) new_centers[c]/count[c];
            if(new_centers[c] != centers[c]){
                changed = 1;
                centers[c] = new_centers[c];
            }
            new_centers[c] = 0.0;
            count[c] = 0;
        }
        counter += 1;
    }

    // Elbow method
    if(var == 1.0){
        for(int i=0; i<I; i++){
            var += pow(X[i]-centers[labels[i]],2);
        }
        return var;
    }

    // Silhouette method
    if(var == 2.0){
        double dist[K];
        double min_b, tot_score = 0.0;
        int count[K];

        for(int i=0; i<I; i++){
            min_b = DBL_MAX;
            for(int c=0; c<K; c++){
                dist[c] = 0.0;
                count[c] = 0;
            }
            for(int j=0; j<I; j++){
                dist[labels[j]] += fabs(X[i] - X[j]);
                count[labels[j]] += 1;
            }
            count[labels[i]] -= 1;
            dist[labels[i]] = (double)dist[labels[i]]/count[labels[i]];

            for(int c=0; c<K; c++){
                if(c!=labels[i]){
                    dist[c] = (double) dist[c]/count[c];
                    if(dist[c] < min_b){
                        min_b = dist[c];
                    }
                }
            }
            if(min_b > dist[labels[i]]){
                tot_score += (min_b - dist[labels[i]])/min_b;
            } else {
                tot_score += (min_b - dist[labels[i]])/dist[labels[i]];
            }
        }
        return (double) tot_score/I;
    }

    //Sequential execution
    if(var == 3.0){
        double center;
        if(centers[0] > centers[1]){
            center = centers[1] + (centers[0]-centers[1])/2;
        } else {
            center = centers[0] + (centers[1]-centers[0])/2;
        }   
        for(int i=0; i<I; i++){
            if(X[i] > center){
                return (double) i-1;
            }
        }
    }

    return 1.0;
}

int kmeans_elbow_c(int I, int Kmax, int Kmin, int minN, double X[I], int labels_1[I], int labels_2[I]){ //This function is called to perform the Kmeans algorithm selecting the optimal number k of cluestes using the elbow method
    int current_k = 1, alt = 0;
    int check;
    double ssq[Kmax];

    ssq[0] = kmeans_c(I, current_k, X, labels_1, 1.0);
    ssq[Kmax-1] = kmeans_c(I, Kmax, X, labels_2, 1.0);
    double m = (double) (ssq[Kmax-2] - ssq[0])/(Kmax-1);

    printf("Valore di m: %f\n", m);
    current_k += 1;
    while(current_k <= Kmax){
        printf("Iterazione su K=%d ", current_k);
        if(alt == 0){
            ssq[current_k-1] = kmeans_c(I, current_k, X, labels_2, 1.0);
            check = check_N(I, current_k, minN, labels_2);
            if(((ssq[current_k-1] - ssq[current_k-2]) > m && current_k > Kmin) || check == 1) {
                return current_k; 
            }
            alt = 1;
        } else {
            ssq[current_k-1] = kmeans_c(I, current_k, X, labels_1, 1.0);
            check = check_N(I, current_k, minN, labels_1);
            if(((ssq[current_k-1] - ssq[current_k-2]) > m && current_k > Kmin) || check == 1) {
                return current_k; 
            }
            alt = 0;
        }
        current_k += 1;
    }
    return Kmax;
}

double kmeans_silhouette_c(int I, int Kmax, int Kmin, int minN, double X[I], int labels_1[I], int labels_2[I]){ //This function is called to perform the Kmeans algorithm selecting the optimal number k of cluestes using the silhouette method
    int current_k = 2, alt = 0, check = 0;
    double shiluette_score, highest_siluette = 0.0;
    int optimal_k = 1;

    while(current_k <= Kmax && check == 0){
        printf("Iterazione su K=%d ", current_k);
        if(alt == 0){
            shiluette_score = kmeans_c(I, current_k, X, labels_2, 2.0);
            check = check_N(I, current_k, minN, labels_2);
            printf("slihouette score: %f\n", shiluette_score);
            if(shiluette_score > highest_siluette && current_k >= Kmin) {
                highest_siluette = shiluette_score;
                optimal_k = current_k;
                alt = 1;
            }
        } else {
            shiluette_score = kmeans_c(I, current_k, X, labels_1, 2.0);
            check = check_N(I, current_k, minN, labels_1);
            printf("slihouette score: %f\n", shiluette_score);
            if(shiluette_score > highest_siluette && current_k >= Kmin) {
                highest_siluette = shiluette_score;
                optimal_k = current_k;
                alt = 0;
            }
        }
        current_k += 1;
    }
    return (double) 0.5*alt + optimal_k;
}

//*SIMULTANEOUS LATENT BUDGET ANALYSIS
//Helper function
void print_double_bytes(double num) { //This function is used for debugging
    unsigned char *byte_ptr = (unsigned char*)&num; 
    for (int i = 0; i < sizeof(double); i++) {
        printf("%02X ", byte_ptr[i]);  
    }
}

double round_with_tolerance(double x) { //This function is used to prevent errors that might occur because of small approximatons
    double rounded = round(x);  

    if (fabs(x - rounded) < EPSILON) {
        return rounded;  
    }

    return x; 
}

void mat_inv_GJ(int J, double A[J][J], double Ai[J][J]){ //This function, given A, returns its inverse using Gauss-Johnson method
    double prop;
    //Inizialization of the inverted matrix
    for(int i=0; i<J; i++){
        for(int j=0; j<i; j++){
            Ai[i][j] = 0.0;
        }
        Ai[i][i] = 1.0;
        for(int j=i+1; j<J; j++){
            Ai[i][j] = 0.0;
        }
    }
    //Gauss-Jordan algorithm for the inverted matrix 
    for(int i=0; i<J; i++){
        for(int j=i+1; j<J; j++){
            prop = (double)A[j][i]/A[i][i];
            for(int k=0; k<J; k++){
                A[j][k] = (double)(A[j][k] - A[i][k]*prop);
                Ai[j][k] = (double)(Ai[j][k] - Ai[i][k]*prop);
            }
        }
    }
    for(int i=J-1; i>=0; i--){
        for(int j=i-1; j>=0; j--){
            prop = A[j][i]/A[i][i];
            for(int k=0; k<J; k++){
                A[j][k] = (double)(A[j][k] - A[i][k]*prop);
                Ai[j][k] = (double)(Ai[j][k] - Ai[i][k]*prop);
            }
        }
    }
    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            Ai[i][j] = Ai[i][j]/A[i][i];
        }
    }
    return;
}

void mat_inv_NS(int J, double A[J][J], double Ai[J][J]){ //This function, given A, returns its pseudo-inverse using Newton-Shultz method
    double X[J][J], X_temp[J][J], I_2[J][J], AX[J][J], I_2AX[J][J], norm = 0, err;
    int MAX_ITER = 200;

    //Trace and 2I calculations
    prod_matrix(J, J, J, A, A, X);
    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            I_2[i][j] = 0;
        }
        norm += X[i][i];
        I_2[i][i] = 2;
    }
    mat_tras(J, J, A, X);
    //X_0 initialization
    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            X[i][j] = X[i][j]/norm;
        }
    }

    // printf("La matrice X è:\n");
    // mat_print_double(J,J,X);

    //Iterative execution
    for(int k=0; k<MAX_ITER; k++){
        prod_matrix(J, J, J, A, X, AX);

        // printf("La matrice AX è:\n");
        // mat_print_double(J,J,AX);

        for(int i=0; i<J; i++){
            for(int j=0; j<J; j++){
                I_2AX[i][j] = I_2[i][j] - AX[i][j];
            }
        }

        // printf("La matrice I2_AX è:\n");
        // mat_print_double(J,J,I_2AX);

        prod_matrix(J, J, J, X, I_2AX, X_temp);
        //Check convergence
        err = 0;
        for(int i=0; i<J; i++){
            for(int j=0; j<J; j++){
                err += fabs(X_temp[i][j] - X[i][j]);
            }
        }
        printf("L'errore è: %f\n", err);

        if(err < 0.00000005){
            for(int i=0; i<J; i++){
                for(int j=0; j<J; j++){
                    Ai[i][j] = X_temp[i][j];
                }
            }
            return;
        }
  
        //Update X 
        for(int i=0; i<J; i++){
            for(int j=0; j<J; j++){
                X[i][j] = X_temp[i][j];
            }
        }
    }

    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            Ai[i][j] = X_temp[i][j];
        }
    }
    return;
}

int check_det(int J, double A[J][J]){ //This function, given A, returns its determinant
    double prop, det = 1.0;
    double A_aux[J][J];
    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            A_aux[i][j] = A[i][j];
        }
    }
    for(int i=0; i<J; i++){
        for(int j=i+1; j<J; j++){
            prop = A_aux[j][i]/A[i][i];
            for(int k=0; k<J; k++){
                A_aux[j][k] = (double)(A_aux[j][k] - A_aux[i][k]*prop);
            }
        }
    }
    for(int i=0; i<J; i++){
        det *= A_aux[i][i];
    }
    printf("determinante: %f\n", det);

    det = round_with_tolerance(det);

    if(det == 0.0 || det == -0.0){
        return 0;
    }
    return 1;
}

void get_AIj(int I, int J, int K, double A[I][K], double AIj[I*J][K*J]){
    for(int i=0; i<I*J; i++){
        for(int k=0; k<K*J; k++){
            AIj[i][k] = 0.0;
        }
    }
    for(int j=0; j<J; j++){
        for(int i=0; i<I; i++){
            for(int k=0; k<K; k++){
                AIj[I*j + i][K*j + k] = A[i][k];
            }
        }
    }
}

void vectorize_P(int I, int J, int t, double ***mat, double vect[I*J]){
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            vect[i*J + j] = mat[j][i][t];
        }
    }
}

//The following functions are used to manage C and d, i.e. add/remove rows/values to/from C and d or free the memory when C and d are no longer needed
void add_row_C(int K, int J, double ***C, double arr[J]){
    // Add row to C
    *C = (double **)realloc(*C, (K+1) * sizeof(double *));
    (*C)[K] = (double *)malloc(J * sizeof(double));

    for(int i=0; i<J; i++){
        (*C)[K][i] = arr[i];
    }
    return;
}

void add_val_d(int K, double **d, double val){
    // Add value to d
    *d = (double *)realloc(*d, (K + 1) * sizeof(double));
    (*d)[K] = val;
    return;
}

void remove_row_C(int K, int J, int rem, double ***C){
    for(int i=rem; i<K-1; i++){
        for(int j=0; j<J; j++){
            (*C)[i][j] = (*C)[i+1][j];
        }
    }

    //free(*C[K-1]);
    *C = (double **)realloc(*C, (K - 1) * sizeof(double *));
    return;
}

void remove_val_d(int K, int rem, double **d){
    for(int i=rem; i<K-1; i++){
        (*d)[i] = (*d)[i+1];
    }
    *d = (double *)realloc(*d, (K - 1) * sizeof(double));
    return;
}

void free_C(int K, double **C){
    //Free the memory
    for(int i=0; i<K; i++){
        free(C[i]);
    }
    free(C);
    return;
}

//Function used in the execution of the algorithm
void calc_x0(int I, int J, double Q[I][J], double QtQi[J][J], double r[I], double x_0[J]){ //This function evaluates X0
    double Qt[J][I], QtQ[J][J], QtQiQt[J][I];
    int det;

    mat_tras(I, J, Q, Qt); //Transpose of the matrix Q
    prod_matrix(J, I, J, Qt, Q, QtQ); //Product of Qt and Q

    printf("Matrice QtQ:\n");
    mat_print_double(J,J,QtQ);

    det = check_det(J, QtQ);
    if(det == 0.0){
        mat_inv_NS(J, QtQ, QtQi); //QtQi is calculated as the inverted matrix of QtQ
    } else {
        mat_inv_GJ(J, QtQ, QtQi); //QtQi is calculated as the inverted matrix of QtQ
    }
    
    printf("\nMatrice QtQi:\n");
    mat_print_double(J,J,QtQi);

    prod_matrix(J, J, I, QtQi, Qt, QtQiQt); //(Q'Q)^(-1)Q'

    printf("\nMatrice QtQiQt:\n");
    mat_print_double(J,I,QtQiQt);

    prod_mat_vect(J, I, QtQiQt, r, x_0); // x_0 = (Q'Q)^(-1)Q'r

    for(int i=0; i<J; i++){
        x_0[i] = round_with_tolerance(x_0[i]);
        if(x_0[i] == -0.0){
            x_0[i] = 0.0;
        }
    }

    printf("Array x_0: "); //!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE

    return;
}

void calc_x0_exp(int I, int J, int T, double Q[T][I][J], double QtQi[J][J], double r[T][I], double x_0[J]){
    double Qt[J][I], QtQ[J][J], Qtr[J], QtQ_temp[J][J], Qtr_temp[J];
    int det;

    //Sum[t from 1 to T]:Q'(t)Q(t) and Sum[t from 1 to T]:Q'(t)r(t) inzialization
    for(int i=0; i<J; i++){
        for(int j=0; j<J; j++){
            QtQ[i][j] = 0.0;
        }
        Qtr[i] = 0.0;
    }

    for(int t=0; t<T; t++){
        mat_tras(I, J, Q[t], Qt); //Transpose of the matrix Q(t)
        prod_matrix(J, I, J, Qt, Q[t], QtQ_temp); //Product of Q'(t) and Q(t)

        prod_mat_vect(J, I, Qt, r[t], Qtr_temp); //Product of Q'(t) and r(t)

        for(int i=0; i<J; i++){
            for(int j=0; j<J; j++){
                QtQ[i][j] += QtQ_temp[i][j]; //Sum[t from 1 to T]:Q'(t)Q(t)
            }
            Qtr[i] += Qtr_temp[i];//Sum[t from 1 to T]:Q'(t)r(t)
        }
    }

    printf("Matrice QtQ:\n");
    mat_print_double(J,J,QtQ);

    det = check_det(J, QtQ);
    if(det == 0.0){
        mat_inv_NS(J, QtQ, QtQi); //QtQi is calculated as the inverted matrix of QtQ
    } else {
        mat_inv_GJ(J, QtQ, QtQi); //QtQi is calculated as the inverted matrix of QtQ
    }

    printf("\nMatrice QtQi:\n");
    mat_print_double(J,J,QtQi);

    printf("\nVettore Qtr: ");
    for(int i=0; i<J; i++){
        printf("%f ", Qtr[i]);
    }

    prod_mat_vect(J, J, QtQi, Qtr, x_0); // x_0 = sum[Q(t)'Q(t)]^(-1)sum[Q(t)'r(t)]

    for(int i=0; i<J; i++){
        x_0[i] = round_with_tolerance(x_0[i]);
        if(x_0[i] == -0.0){
            x_0[i] = 0.0;
        }
    }

    printf("\nArray x_0: "); //!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE

    return;
}

void calc_lamb(int K, int J, double QtQi[J][J], double **C, double *d, double x_0[J], double lambda[K]){ //This function evaluates Lambda
    double Ct[J][K], CQtQi[K][J], CQtQiCt[K][K], CQtQiCti[K][K], Cx0[K], d_Cx0[K];
    int det; 

    printf("\nVETTORE X_0 All'interno della funzione lambda: "); //!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE
    
    mat_tras_ptr(K, J, C, Ct); //Transpose of the matrix C

    prod_matrix_ptr(K, J, J, C, QtQi, CQtQi); //Product of C and (Q'Q)^-1

    prod_matrix(K, J, K, CQtQi, Ct, CQtQiCt); //Product of CQtQiCt

    det = check_det(K, CQtQiCt);
    if(det == 0.0){
        mat_inv_NS(K, CQtQiCt, CQtQiCti); //QtQi is calculated as the inverted matrix of QtQ
    } else {
        mat_inv_GJ(K, CQtQiCt, CQtQiCti); //QtQi is calculated as the inverted matrix of QtQ
    }

    prod_mat_ptr_vect(K, J, C, x_0, Cx0); //Product of Cx_0

    for(int i=0; i<K; i++){
        d_Cx0[i] = d[i] - Cx0[i]; //Array of (d-Cx_0)
    }

    printf("\nVETTORE X_0 dopo calcolato d_Cx0: "); //!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE

    prod_mat_vect(K, K, CQtQiCti, d_Cx0, lambda); //lambda = (((C(Q'Q)^-1)C')^-1)(d-Cx_0)

    printf("\nVETTORE X_0 dopo calcolato lambda: "); //!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE

    for(int i=0; i<K; i++){
        lambda[i] = round_with_tolerance(lambda[i]);
        if(lambda[i] == -0.0){
            lambda[i] = 0.0;
        }
    }

    printf("Array lambda: ");  //!DA CANCELLARE
    for(int i=0; i<K; i++){  //!DA CANCELLARE
        printf("%f ", lambda[i]);  //!DA CANCELLARE
    }  //!DA CANCELLARE
    printf("\n");  //!DA CANCELLARE
    return;
}

void calc_x(int K, int J, double QtQi[J][J], double **C, double lambda[K], double x_0[J], double x[J]){ //This function evaluates X
    double Ct[J][K], QtQiCt[J][K], QtQiCtL[J];

    mat_tras_ptr(K, J, C, Ct);
    prod_matrix(J, J, K, QtQi, Ct, QtQiCt);

    printf("\nVETTORE X_0 All'interno della funzione calc x: ");//!DA CANCELLARE
    for(int i=0; i<J; i++){ //!DA CANCELLARE
        printf("%f ", x_0[i]); //!DA CANCELLARE
    } //!DA CANCELLARE
    printf("\n"); //!DA CANCELLARE

    prod_mat_vect(J, K, QtQiCt, lambda, QtQiCtL);
    
    for(int i=0; i<J; i++){
        x[i] = round_with_tolerance(x_0[i] + QtQiCtL[i]);
        if(x[i] == -0.0){
            x[i] = 0.0;
        }
    }

    printf("\nArray x: ");
    for(int i=0; i<J; i++){
        printf("%f ", x[i]);
    }
    printf("\n");
    return;
}

void gen_ssq(int I, int J, int K, double Q[I][J], double r[I], double x[J], double ***C, double **d){ //This function, given [Q,r,C,d], performs the SSQ method and returns X
    int stop = 0, feasible, i_f, counter, k_temp = K; 
    double add_array[J], QtQi[J][J], x_0[J], lambda[K], x_0_saved[J];

    printf("INIZIO SSQ METHOD:\n");
    while(stop == 0){
        printf("INIZIO ITERAZIONE N.%d\n", stop+1);
        printf("Constraints C and d:\n");
        printf("K: %d\n", k_temp);
        for(int i=0; i<k_temp; i++){
            for(int j=0; j<J; j++){
                printf("%f ", (*C)[i][j]);
            }
            printf(" | %f\n", (*d)[i]);
        }

        calc_x0(I, J, Q, QtQi, r, x_0);

        for(int i=0; i<J; i++){
            x_0_saved[i] = x_0[i];
        }

        calc_lamb(k_temp, J, QtQi, *C, *d, x_0, lambda);

        calc_x(k_temp, J, QtQi, *C, lambda, x_0_saved, x);

        counter = 0;
        feasible = 1;
        while(counter < J && feasible == 1){
            if(x[counter] < 0){
                i_f = counter;
                feasible = 0;
            }
            counter += 1;
        }

        if(feasible == 0){
            printf("X NON FEASIBLE\n");
            for(int i=0; i<i_f; i++){
                add_array[i] = 0.0;
            }
            add_array[i_f] = 1.0;
            for(int i=i_f+1; i<J; i++){
                add_array[i] = 0.0;
            }

            add_row_C(k_temp, J, C, add_array);
            add_val_d(k_temp, d, 0.0);
            k_temp += 1;
            printf("C and d updated\n\n");
        }
        if(feasible == 1){
            printf("X FEASIBLE\n");
            for(int i=J; i<K; i++){
                if(lambda[i] < 0){
                    i_f = i;
                    feasible = 0;
                }
            }
            if(feasible == 0){
                printf("LAMBDA NON FEASIBLE\n");
                remove_row_C(k_temp, J, i_f, C);
                remove_val_d(k_temp, i_f, d);
                k_temp -= 1;
            }
        }
        if(feasible == 1){
            printf("LAMBDA FEASIBLE\n");
            stop = 1;
        }
        getchar();
    }
    while(k_temp != K){
        printf("Err 1\n");
        remove_row_C(k_temp, K, J, C);
        printf("Err 2\n");
        remove_val_d(k_temp, K, d);
        k_temp -= 1;       
    }
}

void gen_ssq_exp(int I, int J, int K, int T, double Q[T][I][J], double r[T][I], double x[J], double ***C, double **d){
    int stop = 0, feasible, i_f, counter, k_temp = K; 
    double add_array[J], QtQi[J][J], x_0[J], lambda[K];

    printf("INIZIO SSQ METHOD:\n");
    while(stop == 0){
        printf("INIZIO ITERAZIONE N.%d\n", stop+1);
        printf("Constraints C and d:\n");
        printf("K: %d\n", k_temp);
        for(int i=0; i<k_temp; i++){
            for(int j=0; j<J; j++){
                printf("%f ", (*C)[i][j]);
            }
            printf(" | %f\n", (*d)[i]);
        }

        calc_x0_exp(I, J, T, Q, QtQi, r, x_0);
        calc_lamb(k_temp, J, QtQi, *C, *d, x_0, lambda);
        calc_x(k_temp, J, QtQi, *C, lambda, x_0, x);

        counter = 0;
        feasible = 1;
        while(counter < J && feasible == 1){
            if(x[counter] < 0){
                i_f = counter;
                feasible = 0;
            }
            counter += 1;
        }

        if(feasible == 0){
            printf("X NON FEASIBLE\n");
            for(int i=0; i<i_f; i++){
                add_array[i] = 0.0;
            }
            add_array[i_f] = 1.0;
            for(int i=i_f+1; i<J; i++){
                add_array[i] = 0.0;
            }

            add_row_C(k_temp, J, C, add_array);
            add_val_d(k_temp, d, 0.0);
            k_temp += 1;
            printf("C and d updated\n\n");
        }
        if(feasible == 1){
            printf("X FEASIBLE\n");
            for(int i=J; i<K; i++){
                if(lambda[i] < 0){
                    i_f = i;
                    feasible = 0;
                }
            }
            if(feasible == 0){
                printf("LAMBDA NON FEASIBLE\n");
                remove_row_C(k_temp, J, i_f, C);
                remove_val_d(k_temp, i_f, d);
                k_temp -= 1;
            }
        }
        if(feasible == 1){
            printf("LAMBDA FEASIBLE\n");
            stop = 1;
        }
        getchar();
    }
    while(k_temp != K){
        remove_row_C(k_temp, K, J, C);
        remove_val_d(k_temp, K, d);
        k_temp -= 1;       
    }
}

//Functions called by Python
void omo1_slba_c(int J, int I, int T, double ***mat){
    //First step: Minimizing with respect to the mixing parameters
    printf("Inzio esecuzione in C");
    int K = 3;
    int t = 0;
    int Kc = 1, Kc_l;
    double randmat[J][K], B[J][K], ai[K], p[J], A[I][K];
    mat_rand(J, K, randmat);
    conv_ZB(J, K, randmat, B);

    printf("No errore 1\n");
    double **C;
    double *d;

    C = (double **)malloc(Kc * sizeof(double *)); 
    for(int i=0; i<K; i++){ 
        C[i] = (double *)malloc(K * sizeof(double));     
    } 
    d = (double *)malloc(Kc * sizeof(double)); 
    C[0][0] = 1.0; C[0][1] = 1.0; C[0][2] = 1.0;
    d[0] = 1.0;

    printf("No errore 2\n");
    for(int i=0; i<I; i++){
        Kc_l = Kc;
        for(int j=0; j<J; j++){
            p[j] = mat[j][i][t];
        }
        printf("No errore 3\n");
        gen_ssq(J, K, Kc_l, B, p, ai, &C, &d);
        for(int k=0; k<K; k++){
            A[i][k] = ai[k];
        }
        while(Kc_l != Kc){
            remove_row_C(Kc_l, K, Kc, &C);
            remove_val_d(Kc_l, Kc, &d);
            Kc_l -= 1;
        }
        printf("\n");
    }

    printf("Matrice A ottimizzata: \n");
    mat_print_double(I, K, A);
    printf("\n");

    printf("No errore 4\n");
    free_C(Kc, C); //! NON CANCELLARE!!!!!!!!!
    free(d); //! NON CANCELLARE!!!!!!!!!
    return;
}

void slba_omogen_A(int J, int I, int T, double ***mat){
    //First step: Minimizing with respect to the mixing parameters
    srand(time(NULL));
    printf("Inzio esecuzione in C slba omogen\n");
    int K = 3, Ke = K*J;
    int Kc = 1, Kc_l;
    int iter = 0, MAX_ITER=1;
    double randmat[J][K], B[T][J][K], p[T][J], A[I][K], AIj[I*J][K*J], Bv[J*K], Pv[I*J], arr[Ke];

    for(int t=0; t<T; t++){
        mat_rand(J, K, randmat);
        conv_ZB(J, K, randmat, B[t]);
        printf("Matrice inizializzata B(%d):\n", t+1);
        mat_print_double(J, K, B[t]);
    }

    printf("No errore 1\n");

    //Inizialization of C and d used to optimize A
    double **C;
    double *d;

    C = (double **)malloc(Kc * sizeof(double *)); 
    for(int i=0; i<Kc; i++){ 
        C[i] = (double *)malloc(K * sizeof(double));     
    } 
    d = (double *)malloc(Kc * sizeof(double)); 
    C[0][0] = 1.0; C[0][1] = 1.0; C[0][2] = 1.0;
    d[0] = 1.0;

    //Inizialization of C_e and d_e used to optimize B
    double **C_e = NULL;
    double *d_e = NULL;

    C_e = (double **)malloc(K * sizeof(double *)); 
    for(int i=0; i<K; i++){ 
        C_e[i] = (double *)malloc(Ke * sizeof(double));     
    } 
    for(int i=0; i<K; i++){
        for(int j=0; j<J*K; j++){
            C_e[i][j] = 0.0;
        }
    }
    for(int i=0; i<K; i++){
        for(int j=0; j<J; j++){
            C_e[i][J*i + j]= 1.0;
        }
    }

    d_e = (double *)malloc(K * sizeof(double)); 
    for(int i=0; i<K; i++){
        d_e[i] = 1.0;
    }

    printf("No errore 2\n");
    while(iter < MAX_ITER){
        //Optimization with respect of A
        for(int i=0; i<I; i++){ 
            Kc_l = Kc;
            for(int t=0; t<T; t++){
                for(int j=0; j<J; j++){
                    p[t][j] = mat[j][i][t];
                }
            }

            printf("No errore 3\n");
            gen_ssq_exp(J, K, Kc_l, T, B, p, A[i], &C, &d);

            printf("\n");
        }

        printf("Matrice A ottimizzata: \n");
        mat_print_double(I, K, A);
        printf("\n");

        //Optimization with respect of B(t)
        get_AIj(I, J, K, A, AIj);

        for(int t=0; t<T; t++){
            vectorize_P(I, J, t, mat, Pv); 
            gen_ssq(I*J, K*J, K, AIj, Pv, Bv, &C_e, &d_e);

            for(int k=0; k<K; k++){
                for(int j=0; j<J; j++){
                    B[t][j][k] = Bv[k*J + j];
                }
            }
        }

        for(int t=0; t<T; t++){
            printf("Matrice B(%d) ottimizzata: \n", t+1);
            mat_print_double(J, K, B[t]);
            printf("\n");
        }

        iter += 1;
    }

    printf("No errore 4\n");
    free_C(Kc, C); //! NON CANCELLARE!!!!!!!!!
    free(d); //! NON CANCELLARE!!!!!!!!!
    free_C(K, C_e); //! NON CANCELLARE!!!!!!!!!
    free(d_e); //! NON CANCELLARE!!!!!!!!!

    return;
}

// int main() {
//     // int I = 4;
//     // int J = 3;
//     // int K = 3;
//     // double Q[I][J], r[I], x_0[J], QtQi[J][J], lambda[K], x[J];

//     // Q[0][0] = 11.0; Q[0][1] = 9.0; Q[0][2] = 3.0;
//     // Q[1][0] = 5.0; Q[1][1] = 0.0; Q[1][2] = 4.0;
//     // Q[2][0] = 7.0; Q[2][1] = 9.0; Q[2][2] = 9.0;
//     // Q[3][0] = 1.0; Q[3][1] = 2.0; Q[3][2] = 3.0;
//     // r[0] = 1.0; r[1] = 3.0; r[2] = 2.0; r[3] = 1.0; 

//     // // Initialization of C and d
//     // double **C;
//     // double *d;
//     // C = (double **)malloc(K * sizeof(double *)); 
//     // for(int i=0; i<K; i++){ 
//     //     C[i] = (double *)malloc(J * sizeof(double));     
//     // } 
//     // d = (double *)malloc(K * sizeof(double)); 
//     // C[0][0] = 1.0; C[0][1] = 2.0; C[0][2] = 1.0; 
//     // C[1][0] = 3.0; C[1][1] = 2.0; C[1][2] = 1.0;
//     // C[2][0] = 1.0; C[2][1] = 4.0; C[2][2] = 7.0;
//     // d[0] = 2; d[1] = 4; d[2] = 6; 

//     int I = 2, J = 2, K = 1;
//     double Q[I][J], r[I], x_0[J], QtQi[J][J], lambda[K], x[J];

//     Q[0][0] = -1.0; Q[0][1] = 2.0; Q[1][0] = -3.0; Q[1][1] = 4.0; 
//     r[0] = 3.0, r[1] = 5.0;

//     double **C;
//     double *d;

//     C = (double **)malloc(K * sizeof(double *)); 
//     for(int i=0; i<K; i++){ 
//         C[i] = (double *)malloc(J * sizeof(double));     
//     } 
//     d = (double *)malloc(K * sizeof(double)); 
//     C[0][0] = -2.0; C[0][1] = 2.0;
//     d[0] = 3.0;

//     gen_ssq(I, J, K, Q, r, x, C, d);

//     //Free the memory
//     free_C(K, C); //! NON CANCELLARE!!!!!!!!!
//     free(d); //! NON CANCELLARE!!!!!!!!!
//     return 0;
// }

