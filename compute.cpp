#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include <vector>
#include <tuple>
#include <gmp.h>
#include <chrono>
#include <stdint.h>
#include <omp.h>




void poly_mult_update(mpz_t* p1, mpz_t* p2, uint32_t max_degree) {
    // Computes p1 = p1 * p2
    mpz_t* product = new mpz_t[1+max_degree];
    
    mpz_t term;
    mpz_init(term);

    mpz_t tempProd;
    mpz_init(tempProd);
    
    for (size_t i = 0; i <= max_degree; i++) {
        mpz_init(product[i]); // initialize product[i] to 0
        for (size_t j = 0; j <= i; j++) {
            mpz_set (tempProd, product[i]);

            mpz_mul(term, p1[j], p2[i-j]); // multiply coefficients
            mpz_add(product[i], tempProd, term); // add to product[i]
        }
    }

    mpz_clear(term);
    mpz_clear(tempProd);
    
    for (size_t i = 0; i <= max_degree; i++) {
        mpz_set (p1[i], product[i]);
        mpz_clear(product[i]);
    }

    delete product;
}

void poly_mult(mpz_t* p1, mpz_t* p2, mpz_t* product, uint32_t max_degree) {
    
    /*
        Expects product to be additive identity
        Computes product = p1 * p2
    */
    mpz_t term;
    mpz_init(term);

    mpz_t tempProd;
    mpz_init(tempProd);

    for (size_t i = 0; i <= max_degree; i++) {
        for (size_t j = 0; j <= i; j++) {
            mpz_set (tempProd, product[i]);

            mpz_mul(term, p1[j], p2[i-j]); // multiply coefficients
            mpz_add(product[i], tempProd, term); // add to product[i]
            
        }
    }

    mpz_clear(term);
    mpz_clear(tempProd);
}




int main(int argc, char *argv[]){

    int option;

    int32_t precision = 128;
    int32_t genomeSize = 1000;
    int32_t readCountHap1 = -1;
    int32_t readCountHap2 = -1;
    int32_t hetLocus = 200;
    int32_t threadCount = 32;
    
    std::string readDistFileHap1;
    std::string readDistFileHap2;

    std::cout << "inputs: Genome Size, Read count, Het Locus, Read Distribution File, Precision" << std::endl;

    while ((option = getopt(argc, argv, "g:R:r:h:D:d:p:t:")) != -1) {
        switch (option){
            case 'g':
                genomeSize = atoi(optarg);
                std::cout << "Genome size " << genomeSize << std::endl;
                break;
            case 'R':
                readCountHap1 = atoi(optarg);
                std::cout << "Read Count Haplotype 1 " << readCountHap1 << std::endl;
                break;
            case 'r':
                readCountHap2 = atoi(optarg);
                std::cout << "Read Count Haplotype 2 " << readCountHap2 << std::endl;
                break;
            case 'h':
                hetLocus = atoi(optarg);
                std::cout << "Heterozygous Locus " << hetLocus << std::endl;
                break;
            case 'D':
                readDistFileHap1 = optarg;
                std::cout << "Hap 1 Distribution in " << readDistFileHap1 << std::endl;
                break;
            case 'd':
                readDistFileHap2 = optarg;
                std::cout << "Hap 2 Distribution in " << readDistFileHap2 << std::endl;
                break;
            case 'p':
                precision = atoi(optarg);
                std::cout << "Precision " << precision << std::endl;
                break;
            case 't':
                threadCount = atoi(optarg);
                std::cout << "Number of threads " << threadCount << std::endl;
                break;
            default:
                std::cout << "One of the inputs is invalid." << std::endl;
                std::cerr << "Usage: " << argv[0] << "-g genomeSize -R readCountHap1 -r readCountHap2 -h heterozygousLocus -D readDistributionFileHap1 -d readDistributionFileHap2 -p precision -t numThreads" << std::endl;
                return 1;
        }
    }

    if (precision <= 0 
        || genomeSize <= 0
        || readCountHap1 <= 0
        || readCountHap2 <= 0
        || hetLocus <= 0
        || hetLocus >= genomeSize
        || readDistFileHap1.empty()
        || readDistFileHap2.empty()
        || threadCount <= 0) {
            std::cout << "One of the inputs is invalid." << std::endl;
            std::cerr << "Error: Missing required options." << std::endl;
            std::cerr << "Usage: " << argv[0] << "-g genomeSize -R readCountHap1 -r readCountHap2 -h heterozygousLocus -D readDistributionFileHap1 -d readDistributionFileHap2 -p precision -t threadCount" << std::endl;
            return 1;
    }

    // Setting OpenMP parameters
    omp_set_dynamic (0);
    omp_set_num_threads (threadCount);
    
    // Setting precision
    mpf_set_default_prec(precision);
    std::cout << "Set precision to " << precision << " bits" << std::endl;

    /*
        readDistFile contains lines of the form:
        readLength readCount
        It is assumed that each read length appears exactly once. 
        Read lengths appear in increasing order.
    */

    typedef std::tuple<int32_t, int32_t> distElement;

    std::vector<distElement> readLengthDistHap1;
    std::vector<distElement> readLengthDistHap2;

    std::ifstream readLengthFileHap1(readDistFileHap1);
    if (!readLengthFileHap1) {
        std::cout << "Error concerning read length distribution for Haplotype 1." << std::endl;
        std::cerr << "Can not open read distribution file for Haplotype 1" << std::endl;
        return 1;
    }
    std::ifstream readLengthFileHap2(readDistFileHap2);
    if (!readLengthFileHap2) {
        std::cout << "Error concerning read length distribution for Haplotype 2." << std::endl;
        std::cerr << "Can not open read distribution file for Haplotype 2" << std::endl;
        return 1;
    }

    int32_t tempReadLength = -1;
    int32_t tempReadCount = -1;
    int32_t maxReadLengthHap1 = -1;
    int32_t maxReadLengthHap2 = -1;

    std::string line;
    
    while (std::getline(readLengthFileHap1, line)) {
        std::istringstream iss(line);
        if (iss >> tempReadLength >> tempReadCount){
            if (tempReadLength <= 0) {
                std::cout << "Read lengths can not be negative. Please vet the input file for haplotype 1." << std::endl;
                std::cerr << "Invalid read length " << tempReadLength << std::endl;
                return 1;
            }
            if (tempReadCount <= 0) {
                std::cout << "Read counts can not be negative. Please vet the input file for haplotype 1." << std::endl;
                std::cerr << "Invalid read count " << tempReadCount << std::endl;
                return 1;
            }

            std::cout << tempReadLength << " " << tempReadCount << std::endl;

            if (tempReadLength > maxReadLengthHap1)
                maxReadLengthHap1 = tempReadLength;
            
            std::cout << "Max Read Length Haplotype 1 " << maxReadLengthHap1 << std::endl;
            readLengthDistHap1.emplace_back (tempReadLength, tempReadCount);
        }
    }
    readLengthFileHap1.close();

    while (std::getline(readLengthFileHap2, line)) {
        std::istringstream iss(line);
        if (iss >> tempReadLength >> tempReadCount){
            if (tempReadLength <= 0) {
                std::cout << "Read lengths can not be negative. Please vet the input file for haplotype 1." << std::endl;
                std::cerr << "Invalid read length " << tempReadLength << std::endl;
                return 1;
            }
            if (tempReadCount <= 0) {
                std::cout << "Read counts can not be negative. Please vet the input file for haplotype 1." << std::endl;
                std::cerr << "Invalid read count " << tempReadCount << std::endl;
                return 1;
            }

            std::cout << tempReadLength << " " << tempReadCount << std::endl;

            if (tempReadLength > maxReadLengthHap2)
                maxReadLengthHap2 = tempReadLength;
            
            std::cout << "Max Read Length Haplotype 2 " << maxReadLengthHap2 << std::endl;
            readLengthDistHap2.emplace_back (tempReadLength, tempReadCount);
        }
    }
    readLengthFileHap2.close();

    int32_t readDistHap1[1+maxReadLengthHap1] = {};
    for (const auto& tuple1 : readLengthDistHap1){
        int32_t getReadLength = std::get<0>(tuple1);
        int32_t getReadCount = std::get<1>(tuple1);
        readDistHap1 [getReadLength] = getReadCount;
    }
    int32_t readDistHap2[1+maxReadLengthHap2] = {};
    for (const auto& tuple2 : readLengthDistHap2) {
        int32_t getReadLength = std::get<0>(tuple2);
        int32_t getReadCount = std::get<1>(tuple2);
        readDistHap2 [getReadLength] = getReadCount;
    }


    int32_t distinctReadLengthsHap1 = 0;
    for (int32_t i = 0; i <= maxReadLengthHap1; i++) {
        if (readDistHap1[i] > 0)
            distinctReadLengthsHap1 += 1;
    }
    int32_t distinctReadLengthsHap2 = 0;
    for (int32_t i = 0; i <= maxReadLengthHap2; i++) {
        if (readDistHap2[i] > 0)
            distinctReadLengthsHap2 += 1;
    }

    std::cout << "Distinct read lengths on haplotype 1 " << distinctReadLengthsHap1 << std::endl;
    std::cout << "Distinct read lengths on haplotype 2 " << distinctReadLengthsHap2 << std::endl;

    int32_t validReadLengthsHap1[distinctReadLengthsHap1] = {};
    for (int32_t i = 0, j = 0; i <= maxReadLengthHap1; i++) {
        if (readDistHap1[i] == 0)
            continue;
        
        validReadLengthsHap1[j] = i;
        j += 1;
    }
    int32_t validReadLengthsHap2[distinctReadLengthsHap2] = {};
    for (int32_t i = 0, j = 0; i <= maxReadLengthHap2; i++) {
        if (readDistHap2[i] == 0)
            continue;
        
        validReadLengthsHap2[j] = i;
        j += 1;
    }

    // DEBUG
    std::cout << "Valid Read Lengths Hap 1";
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++)
        std::cout << validReadLengthsHap1[i] << " ";
    std::cout << std::endl;
    std::cout << "Valid Read Lengths Hap 2";
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++)
        std::cout << validReadLengthsHap2[i] << " ";
    std::cout << std::endl;
    // END DEBUG

    
    /*
    We compute total number of read sequencing outputs in the following manner:
    On hap1 reads covering hetLocus have to have generating functions 1...
    On hap2 reads covering hetLocus have to have generating functions 1...
    */
    mpz_t cVar;
    mpz_init (cVar);

    // These permit 0 or more on hap1 covering hetLocus
    mpz_t* genFuncHap1[distinctReadLengthsHap1][1+genomeSize];
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        genFuncHap1[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
            genFuncHap1[i][j] = new mpz_t[size];
            for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                mpz_init (genFuncHap1[i][j][k]);
                mpz_set_ui (genFuncHap1[i][j][k], 1);
            }
        }
    }
    // These permit 0 or more on hap2 covering hetLocus
    mpz_t* genFuncHap2[distinctReadLengthsHap2][1+genomeSize];
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        genFuncHap2[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
            genFuncHap2[i][j] = new mpz_t[size];
            for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                mpz_init (genFuncHap2[i][j][k]);
                mpz_set_ui (genFuncHap2[i][j][k], 1);
            }
        }
    }
    std::cout << "Initialised generating functions" << std::endl;


    /*
    generatingFunctions contains all the polynomials for reads of length 
    validReadLengths[i] stopping at position j. 
    For a fixed i, all polynomials will be multiplied. 
    The result is stored in an array of polynomials "products" indexed by i.
    */

    mpz_t* prodHap1[distinctReadLengthsHap1];
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
        prodHap1[i] = new mpz_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpz_init (prodHap1[i][j]);
        }
        mpz_set_ui (prodHap1[i][0], 1);
    }
    mpz_t* prodHap2[distinctReadLengthsHap2];
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
        prodHap2[i] = new mpz_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpz_init (prodHap2[i][j]);
        }
        mpz_set_ui (prodHap2[i][0], 1);
    }
    std::cout << "Initialised products" << std::endl;

    /*
    Multiplication occurs below
    PARALLELIZED
    */
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            poly_mult_update (
                prodHap1[i], 
                genFuncHap1[i][j], 
                readDistHap1[validReadLengthsHap1[i]]
            );
        }
    }
    #pragma omp parallel for
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            poly_mult_update (
                prodHap2[i], 
                genFuncHap2[i][j], 
                readDistHap2[validReadLengthsHap2[i]]
            );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
        end 
        - start
    );
    std::cout << "All products computed in " << duration.count() << " ms" << std::endl;

    // Extract total number of permutations for a haplotype
    mpz_t totalReadSeqOutputsHap1;
    mpz_init (totalReadSeqOutputsHap1);
    mpz_set_ui (totalReadSeqOutputsHap1, 1);

    mpz_t totalReadSeqOutputsHap2;
    mpz_init (totalReadSeqOutputsHap2);
    mpz_set_ui (totalReadSeqOutputsHap2, 1);

    mpz_set_ui (cVar, 0);
    
    start = std::chrono::high_resolution_clock::now();
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        mpz_set (cVar, totalReadSeqOutputsHap1);
        mpz_mul (
            totalReadSeqOutputsHap1, 
            cVar, 
            prodHap1[i][readDistHap1[validReadLengthsHap1[i]]]
        );
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        mpz_set (cVar, totalReadSeqOutputsHap2);
        mpz_mul (
            totalReadSeqOutputsHap2, 
            cVar, 
            prodHap2[i][readDistHap2[validReadLengthsHap2[i]]]
        );
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (
        end
        - start
    );
    std::cout << "Total count computed in " << duration.count() << " ms" << std::endl;

for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpz_clear (prodHap1[i][j]);
        }
        delete prodHap1[i];
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpz_clear (prodHap2[i][j]);
        }
        delete prodHap2[i];
    }

    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                mpz_clear (genFuncHap1[i][j][k]);
            }
            delete genFuncHap1[i][j];
        }
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                mpz_clear (genFuncHap2[i][j][k]);
            }
            delete genFuncHap2[i][j];
        }
    }

    // These are for zero on hap1
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        genFuncHap1[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
            genFuncHap1[i][j] = new mpz_t[size];
            if ( (hetLocus <= j) && (j <= ((hetLocus + maxReadLengthHap1) - 1)) ) {
                if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                    mpz_init (genFuncHap1[i][j][0]);
                    mpz_set_ui (genFuncHap1[i][j][0], 1);
                    for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                        mpz_init (genFuncHap1[i][j][k]);
                        mpz_set_ui (genFuncHap1[i][j][k], 0);
                    }
                }
                else{
                    for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                        mpz_init (genFuncHap1[i][j][k]);
                        mpz_set_ui (genFuncHap1[i][j][k], 1);
                    }
                }
            }
            else {
                for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                    mpz_init (genFuncHap1[i][j][k]);
                    mpz_set_ui (genFuncHap1[i][j][k], 1);
                }
            }
        }
    }
    // These are for 0 on hap2
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        genFuncHap2[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
            genFuncHap2[i][j] = new mpz_t[size];
            if ( (hetLocus <= j) && (j <= ((hetLocus + maxReadLengthHap2) - 1)) ) {
                if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                    mpz_init (genFuncHap2[i][j][0]);
                    mpz_set_ui (genFuncHap2[i][j][0], 1);
                    for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                        mpz_init (genFuncHap2[i][j][k]);
                        mpz_set_ui (genFuncHap2[i][j][k], 0);
                    }
                }
                else {
                    for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                        mpz_init (genFuncHap2[i][j][k]);
                        mpz_set_ui (genFuncHap2[i][j][k], 1);
                    }
                }
            }
            else {
                for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                    mpz_init (genFuncHap2[i][j][k]);
                    mpz_set_ui (genFuncHap2[i][j][k], 1);
                }
            }
        }
    }
    std::cout << "Initialised generating functions" << std::endl;


    /*
    generatingFunctions contains all the polynomials for reads of length 
    validReadLengths[i] stopping at position j. 
    For a fixed i, all polynomials will be multiplied. 
    The result is stored in an array of polynomials "products" indexed by i.
    */

    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
        prodHap1[i] = new mpz_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpz_init (prodHap1[i][j]);
        }
        mpz_set_ui (prodHap1[i][0], 1);
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
        prodHap2[i] = new mpz_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpz_init (prodHap2[i][j]);
        }
        mpz_set_ui (prodHap2[i][0], 1);
    }
    std::cout << "Initialised products" << std::endl;

    /*
    Multiplication occurs below
    PARALLELIZED
    */
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            poly_mult_update (
                prodHap1[i], 
                genFuncHap1[i][j], 
                readDistHap1[validReadLengthsHap1[i]]
            );
        }
    }
    #pragma omp parallel for
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            poly_mult_update (
                prodHap2[i], 
                genFuncHap2[i][j], 
                readDistHap2[validReadLengthsHap2[i]]
            );
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (
        end 
        - start
    );
    std::cout << "All products computed in " << duration.count() << " ms" << std::endl;

    // Extract total number of permutations for a haplotype
    mpz_t totalReadSeqOutputsZeroHap1;
    mpz_init (totalReadSeqOutputsZeroHap1);
    mpz_set_ui (totalReadSeqOutputsZeroHap1, 1);

    mpz_t totalReadSeqOutputsZeroHap2;
    mpz_init (totalReadSeqOutputsZeroHap2);
    mpz_set_ui (totalReadSeqOutputsZeroHap2, 1);

    mpz_set_ui (cVar, 0);
    
    start = std::chrono::high_resolution_clock::now();
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        mpz_set (cVar, totalReadSeqOutputsZeroHap1);
        mpz_mul (
            totalReadSeqOutputsZeroHap1, 
            cVar, 
            prodHap1[i][readDistHap1[validReadLengthsHap1[i]]]
        );
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        mpz_set (cVar, totalReadSeqOutputsZeroHap2);
        mpz_mul (
            totalReadSeqOutputsZeroHap2, 
            cVar, 
            prodHap2[i][readDistHap2[validReadLengthsHap2[i]]]
        );
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (
        end
        - start
    );
    std::cout << "Total count computed in " << duration.count() << " ms" << std::endl;

    // char strP[256];
    // mp_exp_t expP;
    // mpz_get_str(strP, &expP, 10, 0, totalReadSeqOutputs);
    mpz_t T11;
    mpz_init (T11);
    mpz_t T10;
    mpz_init (T10);
    mpz_t T01;
    mpz_init (T01);
    mpz_t T00;
    mpz_init (T00);

    mpz_mul (T11, totalReadSeqOutputsHap1, totalReadSeqOutputsHap2);
    mpz_mul (T10, totalReadSeqOutputsHap1, totalReadSeqOutputsZeroHap2);
    mpz_mul (T01, totalReadSeqOutputsZeroHap1, totalReadSeqOutputsHap2);
    mpz_mul (T00, totalReadSeqOutputsZeroHap1, totalReadSeqOutputsZeroHap2);

    std::cout << "T11 " << T11 << std::endl;
    std::cout << "T10 " << T10 << std::endl;
    std::cout << "T01 " << T01 << std::endl;
    std::cout << "T00 " << T00 << std::endl;

    // Deleting genFuncHap1, genFuncHap2, prodHap1, and prodHap2
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpz_clear (prodHap1[i][j]);
        }
        delete prodHap1[i];
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpz_clear (prodHap2[i][j]);
        }
        delete prodHap2[i];
    }

    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                mpz_clear (genFuncHap1[i][j][k]);
            }
            delete genFuncHap1[i][j];
        }
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                mpz_clear (genFuncHap2[i][j][k]);
            }
            delete genFuncHap2[i][j];
        }
    }

    /*
    Computing number of permutations for a fixed (x1, x2). Multiple steps here:
    1. Count permutations on haplotype 1 where we only condition on classes 2 
        and 4 : N1
    2. Count permutations on haplotype 1 where we condition on classes 2 and 
        4, and force class 1 to have no reads : N11
    3. Count permutations on haplotype 1 where we condition on classes 2 and 
        4, and force class 3 to have no reads : N13
    4. Count permutations on haplotype 1 where we condition on classes 2 and 
        4, and force classes 1 and 3 to have no reads : N113
    5. Count permutations on haplotype 2 where we only condition on classes 2 
        and 4 : N2
    6. Count permutations on haplotype 2 where we condition on classes 2 and 
        4, and force class 1 to have no reads : N21
    7. Count permutations on haplotype 2 where we condition on classes 2 and 
        4, and force class 3 to have no reads : N23
    8. Count permutations on haplotype 2 where we condition on classes 2 and 
        4, and force classes 1 and 3 to have no reads : N213

    From this point, the code is parllelized again for loop x1.
    Below is a list of variables defined inside the loop x1.
    */

    mpz_t* aggregateCount = new mpz_t[threadCount];
    for (int32_t i = 0; i < threadCount; i++) {
        mpz_init (aggregateCount[i]);
    }

    mpz_t* aggregateError = new mpz_t[threadCount];
    for (int32_t i = 0; i < threadCount; i++) {
        mpz_init (aggregateError[i]);
    }

    /*
    For Haplotype 1:
    *   x1 lies in [hetLocus + 1, (hetLocus + maxReadLengthHap1) - 1]
    1.  Class 1 Reads stop at x1 and cover hetLocus.
        Thus, their lengths are at least (x1 - hetLocus) + 1.
        Test condition: validReadLengths[i] >= (x1 - hetLocus) + 1.
        At least one of these reads must exist.
    2.  Class 2 Reads stop in [hetLocus, x1 - 1].
        Fixing the stop position as j, their lengths are at least (j - hetLocus) + 1.
        Test condition: validReadLengths[i] >= (j - hetLocus) + 1.
    3.  Class 3 Reads stop in [hetLocus+1, x1 - 1] and also start in [hetLocus+1, x1 - 1].
        Of these, the reads which start and stop in [hetLocus+1, x2] don't cause coverage gaps.
        Reads which start and stop in [x2+1, x1-1] may cause coverage gaps, but doesn't fit our definition.
        Reads which start in [hetLocus+1, x2] and stop in [x2+1, x1-1] cause assembly gaps.
        At least one of the reads which cause assembly gaps must exist on either haplotype.
    4.  Class 4 Reads stop in [1, hetLocus - 1] or [x1+1, genomeSize]. 
        They must start in either [1, hetLocus - 1] or [x2+1, genomeSize].
    For Haplotype 2:
    *   x2 lies in [hetLocus, (hetLocus + maxReadLengthHap2) - 1].
    *   If the right limit is x1 - 1, nothing changes in the code.
    *   If the right limit is hetLocus + maxReadLengthHap2 - 1, then every occurrence of x1 - 1 must be replaced.
    1.  Class 1 Reads stop at x2 and cover hetLocus.
        Thus, their lengths are at least (x2 - hetLocus) + 1
        Test condition: validReadLengths[i] >= (x2 - hetLocus) + 1.
        At least one of these reads must exist.
    2.  Class 2 Reads stop in [hetLocus, x2 - 1].
        Fixing the stop position as j, their lengths are at least (j - hetLocus) + 1.
        Test condition: validReadLengths[i] >= (j - hetLocus) + 1.
    3.  Class 3 Reads stop in [hetLocus+1, x1 - 1] and also start in [hetLocus+1, x1 - 1].
        Of these, the reads which start and stop in [hetLocus+1, x2] don't cause coverage gaps.
        Reads which start and stop in [x2+1, x1-1] cause coverage gaps but don't fit our definition.
        Reads which start in [hetLocus+1, x2] and stop in [x2+1, x1-1] case assembly gaps.
        At least one of the reads which cause assembly gaps must exist on either haplotype.
    4.  Class 4 Reads stop in [1, hetLocus - 1] or [x1+1, genomeSize].
        They must start in either [1, hetLocus - 1] or [x2+1, genomeSize].
    If x2 == x1:
        Do nothing. Continue.
    If x2 > x1:
        Reverse the behaviours of haplotype 1 and haplotype 2. 
        Initialise polynomials for haplotype 1 with the conditions of haplotype 2 and vice versa.
    */

    #pragma omp parallel for
    for (int32_t x1 = hetLocus; x1 <= ((hetLocus + maxReadLengthHap1) - 1); x1++) {
        
        // Define and initialise generatingFunctions
        mpz_t* generatingFunctionsHap1[distinctReadLengthsHap1][1+genomeSize];
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            generatingFunctionsHap1[i][0] = NULL;
            for (int32_t j = 1; j <= genomeSize; j++) {
                int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
                generatingFunctionsHap1[i][j] = new mpz_t[size];
                for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                    mpz_init (generatingFunctionsHap1[i][j][k]);
                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                }
            }
        }
        mpz_t* generatingFunctionsHap2[distinctReadLengthsHap2][1+genomeSize];
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            generatingFunctionsHap2[i][0] = NULL;
            for (int32_t j = 1; j <= genomeSize; j++) {
                int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
                generatingFunctionsHap2[i][j] = new mpz_t[size];
                for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                    mpz_init (generatingFunctionsHap2[i][j][k]);
                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                }
            }
        }

        // Define and initialise calcVar
        mpz_t calcVar;
        mpz_init (calcVar);

        // Define and initialise products
        mpz_t* productsHap1[distinctReadLengthsHap1];
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
            productsHap1[i] = new mpz_t[size];
            for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                mpz_init (productsHap1[i][j]);
            }
            mpz_set_ui (productsHap1[i][0], 1);
        }
        mpz_t* productsHap2[distinctReadLengthsHap2];
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
            productsHap2[i] = new mpz_t[size];
            for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                mpz_init (productsHap2[i][j]);
            }
            mpz_set_ui (productsHap2[i][0], 1);
        }

        // Inner for loop x2
        // int32_t limitHap2 = std::min((x1 - 1), ((hetLocus + maxReadLengthHap2) - 1));
        for (int32_t x2 = hetLocus; x2 <= ((hetLocus + maxReadLengthHap2) - 1); x2++) {
            
            if (x1 == x2){
                continue;
            }

            mpz_t N1;
            mpz_init (N1);
            mpz_t N11;
            mpz_init (N11);
            mpz_t N13;
            mpz_init (N13);
            mpz_t N113;
            mpz_init (N113);
            mpz_t N2;
            mpz_init (N2);
            mpz_t N21;
            mpz_init (N21);
            mpz_t N23;
            mpz_init (N23);
            mpz_t N213;
            mpz_init (N213);

            if (x2 < x1) {
                /* 
                    Compute N1
                    j == x1 means the reads are in class 1, or class 3b which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }

                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N1
                mpz_set_ui (N1, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N1);
                    mpz_mul (
                        N1, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N11
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (invalid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }
                }
                
                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N11
                mpz_set_ui (N11, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N11);
                    mpz_mul (
                        N11, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N13
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N13
                mpz_set_ui (N13, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N13);
                    mpz_mul (
                        N13, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N113
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N113
                mpz_set_ui (N113, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N113);
                    mpz_mul (
                        N113, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }



                /* 
                    Compute N2
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid) 
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (validreadLengthsHap2[i] >= ((j - hetLocus) + 1))
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N2
                mpz_set_ui (N2, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N2);
                    mpz_mul (
                        N2, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N21
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3a (valid ) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid ) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid ) if (i < j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid ) 
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid ) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] <= (j - x2))) {
                                // Group 3c (valid ) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N21
                mpz_set_ui (N21, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N21);
                    mpz_mul (
                        N21, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N23
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid) 
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }

                }
                
                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N23
                mpz_set_ui (N23, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N23);
                    mpz_mul (
                        N23, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N213
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (invalid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }

                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N213
                mpz_set_ui (N213, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N213);
                    mpz_mul (
                        N213, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }
            }
            else {
                /* 
                    Compute N2
                    j == x2 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x1) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x1 + 1) && (j <= (x2-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < j))
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }

                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N2
                mpz_set_ui (N2, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N2);
                    mpz_mul (
                        N2, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N21
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3c (invalid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < j))
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }
                }
                
                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N21
                mpz_set_ui (N21, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N21);
                    mpz_mul (
                        N21, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N23
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N23
                mpz_set_ui (N23, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N23);
                    mpz_mul (
                        N23, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }


                /* 
                    Compute N213
                    j == x1 means the reads are in class 1, or class 3 which upon 
                    deletion create a coverage gap.
                    (j >= hetLocus) && (j <= x2) accounts for some of the reads in 
                    class 2, and some reads in class 3 which upon deletion do not 
                    create a coverage gap.
                    (j >= x2 + 1) && (j <= (x1-1)) accounts for the remaining reads
                    in class 2, and those reads in class 3 which upon deletion 
                    cause a coverage gap.
                    The rest are class 4.
                */
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x1) {
                            // Group 1, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((x1 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else if (((x1 - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((x1 - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (x1 - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= (x1 - 1))) {
                            // Group 2, 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else if (((j - x2) < validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3c (valid) if 1 <= i <= (j - x2)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) < j) && (j <= x2)) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid)
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if ((1 <= validReadLengthsHap2[i]) && (validReadLengthsHap2[i] < ((j - hetLocus) + 1)))
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid)
                            // Account for 1 <= i <= maxReadLengthHap2
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4  j > x1
                            // Account for 1 <= i <= maxReadLengthHap2
                            if (validReadLengthsHap2[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap2[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set_ui (productsHap2[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                        mpz_set_ui (productsHap2[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap2[i], 
                            generatingFunctionsHap2[i][j], 
                            readDistHap2[validReadLengthsHap2[i]]
                        );
                    }
                }

                // Extract N213
                mpz_set_ui (N213, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                    mpz_set (calcVar, N213);
                    mpz_mul (
                        N213, 
                        calcVar, 
                        productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                    );
                }



                /* 
                    Compute N1
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (valid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid) 
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N1
                mpz_set_ui (N1, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N1);
                    mpz_mul (
                        N1, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N11
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (invalid)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3a (valid ) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid ) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid ) if (i < j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid ) 
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (valid ) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] <= (j - x2))) {
                                // Group 3c (valid ) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }
                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N11
                mpz_set_ui (N11, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N11);
                    mpz_mul (
                        N11, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N13
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2 (valid) 
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }

                }
                
                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N13
                mpz_set_ui (N13, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N13);
                    mpz_mul (
                        N13, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }


                /* 
                    Compute N113
                */
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        if (j == x2) {
                            // Group 1, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((x2 - hetLocus) + 1)) {
                                // Group 1 (invalid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < x2)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (((hetLocus + 1) <= j) && (j <= (x2 - 1))) {
                            // Group 2, 3a
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                                // Group 2 (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Group 3a (valid) if (i < j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                        }
                        else if (j == hetLocus) {
                            // Group 2
                            // Account for 1 <= i <= maxReadLengthHap1
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                            }
                        }
                        else if (((x2 + 1) <= j) && (j <= x1)) {
                            // Group 3b, 3c
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (((j - x2) < validReadLengthsHap1[i]) && (validReadLengthsHap1[i] < ((j - hetLocus) + 1))) {
                                // Group 3b (invalid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                            else if ((1 <= validReadLengthsHap1[i]) && (validReadLengthsHap1[i] <= (j - x2))) {
                                // Group 3c (valid) 
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                // Invalid if (i >= j)
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else if (j < hetLocus) {
                            // Group 4
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                        else {
                            // Group 4 if j > x1
                            // Account for 1 <= i <= maxReadLengthHap1
                            if (validReadLengthsHap1[i] <= (j - x2)) {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 1);
                                }
                            }
                            else {
                                for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                    mpz_set_ui (generatingFunctionsHap1[i][j][k], 0);
                                }
                            }
                        }
                    }
                }

                // Product polynomials
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set_ui (productsHap1[i][0], 1);
                    for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                        mpz_set_ui (productsHap1[i][j], 0);
                    }

                }

                // Multiplication
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    for (int32_t j = 1; j <= genomeSize; j++) {
                        poly_mult_update (
                            productsHap1[i], 
                            generatingFunctionsHap1[i][j], 
                            readDistHap1[validReadLengthsHap1[i]]
                        );
                    }
                }

                // Extract N113
                mpz_set_ui (N113, 1);
                mpz_set_ui (calcVar, 0);
                
                for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                    mpz_set (calcVar, N113);
                    mpz_mul (
                        N113, 
                        calcVar, 
                        productsHap1[i][readDistHap1[validReadLengthsHap1[i]]]
                    );
                }
            }



            // Computing (aggregate) the numerator
            mpz_t tempVar1;
            mpz_init (tempVar1);
            mpz_t tempVar2;
            mpz_init (tempVar2);
            mpz_t tempVar3;
            mpz_init (tempVar3);
            mpz_t tempVar4;
            mpz_init (tempVar4);

            mpz_t prodVar1;
            mpz_init (prodVar1);
            mpz_t prodVar2;
            mpz_init (prodVar2);
            mpz_t coverageGapPermutations;
            mpz_init (coverageGapPermutations);

            mpz_sub (tempVar1, N1, N11);
            mpz_sub (tempVar2, N2, N21);
            mpz_mul (prodVar1, tempVar1, tempVar2);

            mpz_sub (tempVar3, N13, N113);
            mpz_sub (tempVar4, N23, N213);
            mpz_mul (prodVar2, tempVar3, tempVar4);

            mpz_sub (calcVar, prodVar1, prodVar2);
            mpz_abs (coverageGapPermutations, calcVar);

            mpz_set (calcVar, aggregateCount[omp_get_thread_num()]);
            mpz_add (aggregateCount[omp_get_thread_num()], calcVar, coverageGapPermutations);

            mpz_clear (tempVar1);
            mpz_clear (tempVar2);
            mpz_clear (tempVar3);
            mpz_clear (tempVar4);
            mpz_clear (prodVar1);
            mpz_clear (prodVar2);
            mpz_clear (coverageGapPermutations);
            mpz_clear (N1);
            mpz_clear (N11);
            mpz_clear (N13);
            mpz_clear (N113);
            mpz_clear (N2);
            mpz_clear (N21);
            mpz_clear (N23);
            mpz_clear (N213);
        }

        mpz_clear(calcVar);

        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            for (int32_t j = 1; j <= genomeSize; j++) {
                for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                    mpz_clear (generatingFunctionsHap1[i][j][k]);
                }
                delete generatingFunctionsHap1[i][j];
            }
        }
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            for (int32_t j = 1; j <= genomeSize; j++) {
                for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                    mpz_clear (generatingFunctionsHap2[i][j][k]);
                }
                delete generatingFunctionsHap2[i][j];
            }
        }
    
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                mpz_clear (productsHap1[i][j]);
            }
            delete productsHap1[i];
        }
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                mpz_clear (productsHap2[i][j]);
            }
            delete productsHap2[i];
        }


    }

    // Computing the probability
    mpz_t numerator;
    mpz_init (numerator);

    mpz_t denominator;
    mpz_init (denominator);

    mpq_t ratio;
    mpq_init (ratio);

    for (int32_t i = 0; i < threadCount; i++) {
        mpz_set (cVar, numerator);
        mpz_add (numerator, cVar, aggregateCount[i]);
    }

    mpf_t probability;
    mpf_init (probability);

    // (denominator) cVar = T11 - T10 - T01 + T00
    mpz_sub (cVar, T11, T10);
    mpz_set (denominator, cVar);
    mpz_sub (cVar, denominator, T01);
    mpz_set (denominator, cVar);
    mpz_add (cVar, denominator, T00);
    mpz_set (denominator, cVar);
    
    mpq_set_num (ratio, numerator);
    mpq_set_den (ratio, denominator);
    mpq_canonicalize (ratio);

    mpf_set_q (probability, ratio);


    std::cout << "Total number of read sequencing outputs: " << denominator << std::endl;
    std::cout << "No of outputs with coverage gap event: " << numerator << std::endl;
    std::cout << "probability = " << probability << std::endl;
    
    // Free initialised variables
    for (int32_t i = 0; i < threadCount; i++) {
        mpz_clear (aggregateCount[i]);
        // mpz_clear (aggregateError[i]);
    }
    delete aggregateCount;
    // delete aggregateError;
    
    mpf_clear (probability);
    // mpz_clear (totalError);
    // mpz_clear (relativeError);
    mpz_clear (cVar);
    mpz_clear (numerator);
    mpz_clear (denominator);
    mpq_clear (ratio);
    mpz_clear (totalReadSeqOutputsHap1);
    mpz_clear (totalReadSeqOutputsHap2);
    mpz_clear (totalReadSeqOutputsZeroHap1);
    mpz_clear (totalReadSeqOutputsZeroHap2);
    mpz_clear (T11);
    mpz_clear (T10);
    mpz_clear (T01);
    mpz_clear (T00);

return 0;
}
