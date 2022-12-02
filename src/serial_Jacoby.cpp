#include "utils.h"

void JMstep(Matrix &A, int i, int j) {
    float xi = (A[j][j] - A[i][i])/(2*A[i][j]);
    float t = 1./(abs(xi)+sqrt(1.+xi*xi));
    if(xi<0) {
        t = -t;
    };
    float cos = 1./sqrt(1.+t*t);
    float sin = t*cos;
    float oldA_ii = A[i][i];
    float oldA_jj = A[j][j];
    float oldA_ij = A[i][j];
    for(int q=0; q<A.size(); q++) {
        float oldA_iq = A[i][q];
        // float oldA_qi = A[q][i];
        float oldA_jq = A[j][q];
        // float oldA_qj = A[q][j];
        A[i][q] = cos * oldA_iq - sin * oldA_jq;
        A[q][i] = cos * oldA_iq - sin * oldA_jq;

        A[j][q] = sin * oldA_iq + cos * oldA_jq;
        A[q][j] = sin * oldA_iq + cos * oldA_jq;
    };
    A[i][i] = oldA_ii - t*oldA_ij;
    A[j][j] = oldA_jj + t*oldA_ij;
    A[i][j] = 0;
    A[j][i] = 0;
}





int main(int argc, char **argv) {
    int N = 50;
    if(argc>1) {
        N = atoi(argv[1]);
    };
    int debug = 0;
    if(argc>2) {
        debug = atoi(argv[2]);
    };
    

    Matrix A = genRandomMatrix(N);
    saveMatrix(A, "A.txt");
    double init_sum = offDiagSum(A);
    printf("sum(A) = %f\n", offDiagSum(A));

    std::vector< std::array<int, 2> > ijs = genIJs(N);
    int step = ijs.size()/10;
    std::ofstream progress_file("../results/sjResults/progress.txt");
    std::string debug_fileName = "../results/sjResults/matrices/";
    int it_num = 0;
    for(int it = 0; it<3; ++it) {
        printf("== %d ",it);
        for(int k = 0; k<ijs.size(); ++k) {
            it_num += 1;
            if( k % step == 0 || it_num < 50) {
                std::cout<<"."<<std::flush;
                saveMatrix(A, debug_fileName+std::to_string(it_num) + ".txt");
            };
            int i = ijs[k][0];
            int j = ijs[k][1];
            JMstep(A, i, j);
            double offSum = offDiagSum(A);
            progress_file << offSum << "\n";
        };
        std::cout << std::endl;
    }
    progress_file.close();
    double final_sum = offDiagSum(A);
    printf("sum(A)_final  = %f\n", final_sum);
    printf(" ratio = %f\n", final_sum/init_sum);
    saveMatrix(A, "JA_20");


    return 0;
}