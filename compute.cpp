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




void poly_mult_update(mpf_t* p1, mpf_t* p2, uint32_t max_degree) {
    // Computes p1 = p1 * p2
    mpf_t* product = new mpf_t[1+max_degree];
    
    mpf_t term;
    mpf_init(term);

    mpf_t tempProd;
    mpf_init(tempProd);
    
    for (size_t i = 0; i <= max_degree; i++) {
        mpf_init(product[i]); // initialize product[i] to 0
        for (size_t j = 0; j <= i; j++) {
            mpf_set (tempProd, product[i]);

            mpf_mul(term, p1[j], p2[i-j]); // multiply coefficients
            mpf_add(product[i], tempProd, term); // add to product[i]
        }
    }

    mpf_clear(term);
    mpf_clear(tempProd);
    
    for (size_t i = 0; i <= max_degree; i++) {
        mpf_set (p1[i], product[i]);
        mpf_clear(product[i]);
    }

    delete product;
}

void poly_mult(mpf_t* p1, mpf_t* p2, mpf_t* product, uint32_t max_degree) {
    
    /*
        Expects product to be additive identity
        Computes product = p1 * p2
    */
    mpf_t term;
    mpf_init(term);

    mpf_t tempProd;
    mpf_init(tempProd);

    for (size_t i = 0; i <= max_degree; i++) {
        for (size_t j = 0; j <= i; j++) {
            mpf_set (tempProd, product[i]);

            mpf_mul(term, p1[j], p2[i-j]); // multiply coefficients
            mpf_add(product[i], tempProd, term); // add to product[i]
            
        }
    }

    mpf_clear(term);
    mpf_clear(tempProd);
}




int main(int argc, char *argv[]){

    int option;

    int32_t precision = -1;
    int32_t genomeSize = -1;
    int32_t readCountHap1 = -1;
    int32_t readCountHap2 = -1;
    int32_t hetLocus = -1;
    int32_t threadCount = -1;
    
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

    
    mpf_t cVar;
    mpf_init (cVar);

    mpf_t* genFuncHap1[distinctReadLengthsHap1][1+genomeSize];
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        genFuncHap1[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
            genFuncHap1[i][j] = new mpf_t[size];
            for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                mpf_init (genFuncHap1[i][j][k]);
                mpf_set_d (genFuncHap1[i][j][k], 1.0);
            }
        }
    }
    mpf_t* genFuncHap2[distinctReadLengthsHap2][1+genomeSize];
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        genFuncHap2[i][0] = NULL;
        for (int32_t j = 1; j <= genomeSize; j++) {
            int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
            genFuncHap2[i][j] = new mpf_t[size];
            for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                mpf_init (genFuncHap2[i][j][k]);
                mpf_set_d (genFuncHap2[i][j][k], 1.0);
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

    mpf_t* prodHap1[distinctReadLengthsHap1];
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
        prodHap1[i] = new mpf_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpf_init (prodHap1[i][j]);
        }
        mpf_set_d (prodHap1[i][0], 1.0);
    }
    mpf_t* prodHap2[distinctReadLengthsHap2];
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
        prodHap2[i] = new mpf_t[size];
        // j goes from 0 upto n_i
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpf_init (prodHap2[i][j]);
        }
        mpf_set_d (prodHap2[i][0], 1.0);
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
    mpf_t totalReadSeqOutputsHap1;
    mpf_init (totalReadSeqOutputsHap1);
    mpf_set_d (totalReadSeqOutputsHap1, 1.0);

    mpf_t totalReadSeqOutputsHap2;
    mpf_init (totalReadSeqOutputsHap2);
    mpf_set_d (totalReadSeqOutputsHap2, 1.0);

    mpf_set_d (cVar, 0.0);
    
    start = std::chrono::high_resolution_clock::now();
    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        mpf_set (cVar, totalReadSeqOutputsHap1);
        mpf_mul (
            totalReadSeqOutputsHap1, 
            cVar, 
            prodHap1[i][readDistHap1[validReadLengthsHap1[i]]]
        );
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        mpf_set (cVar, totalReadSeqOutputsHap2);
        mpf_mul (
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

    // char strP[256];
    // mp_exp_t expP;
    // mpf_get_str(strP, &expP, 10, 0, totalReadSeqOutputs);
    mpf_mul (cVar, totalReadSeqOutputsHap1, totalReadSeqOutputsHap2);
    std::cout << "Number of read sequencing outputs " << cVar << std::endl;


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

    mpf_t* aggregateCount = new mpf_t[threadCount];
    for (int32_t i = 0; i < threadCount; i++) {
        mpf_init (aggregateCount[i]);
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
        Of these, the reads which stop in [hetLocus+1, x2] don't cause coverage gaps.
        Reads which stop in [x2+1, x1-1] cause coverage gaps.
        At least one of the reads which cause coverage gaps must exist on either haplotype.
    4.  Class 4 Reads stop in [1, hetLocus - 1] or [x1+1, genomeSize]. 
        They must start in either [1, hetLocus - 1] or [x2+1, genomeSize].
    For Haplotype 2:
    *   x2 lies in [hetLocus, min (hetLocus + maxReadLengthHap2 - 1, x1 - 1)].
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
        Of these, the reads which stop in [hetLocus+1, x2] don't cause coverage gaps.
        Reads which stop in [x2+1, x1-1] cause coverage gaps.
        At least one of the reads which cause coverage gaps must exist on either haplotype.
    4.  Class 4 Reads stop in [1, hetLocus - 1] or [x1+1, genomeSize].
        They must start in either [1, hetLocus - 1] or [x2+1, genomeSize].
    */

    #pragma omp parallel for
    for (int32_t x1 = hetLocus+1; x1 <= ((hetLocus + maxReadLengthHap1) - 1); x1++) {
        
        // Define and initialise generatingFunctions
        mpf_t* generatingFunctionsHap1[distinctReadLengthsHap1][1+genomeSize];
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            generatingFunctionsHap1[i][0] = NULL;
            for (int32_t j = 1; j <= genomeSize; j++) {
                int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
                generatingFunctionsHap1[i][j] = new mpf_t[size];
                for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                    mpf_init (generatingFunctionsHap1[i][j][k]);
                    mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                }
            }
        }
        mpf_t* generatingFunctionsHap2[distinctReadLengthsHap2][1+genomeSize];
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            generatingFunctionsHap2[i][0] = NULL;
            for (int32_t j = 1; j <= genomeSize; j++) {
                int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
                generatingFunctionsHap2[i][j] = new mpf_t[size];
                for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                    mpf_init (generatingFunctionsHap2[i][j][k]);
                    mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                }
            }
        }

        // Define and initialise calcVar
        mpf_t calcVar;
        mpf_init (calcVar);

        // Define and initialise products
        mpf_t* productsHap1[distinctReadLengthsHap1];
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            int32_t size = 1 + readDistHap1[validReadLengthsHap1[i]];
            productsHap1[i] = new mpf_t[size];
            for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                mpf_init (productsHap1[i][j]);
            }
            mpf_set_d (productsHap1[i][0], 1.0);
        }
        mpf_t* productsHap2[distinctReadLengthsHap2];
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            int32_t size = 1 + readDistHap2[validReadLengthsHap2[i]];
            productsHap2[i] = new mpf_t[size];
            for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                mpf_init (productsHap2[i][j]);
            }
            mpf_set_d (productsHap2[i][0], 1.0);
        }

        // Inner for loop x2
        int32_t limitHap2 = std::min((x1 - 1), ((hetLocus + maxReadLengthHap2) - 1));
        for (int32_t x2 = hetLocus; x2 <= limitHap2; x2++) {
            /* 
                Compute N1
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x1) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                        }
                    }
                    else if ((j <= (x1-1)) && (j >= hetLocus)) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap1[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set_d (productsHap1[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                    mpf_set_d (productsHap1[i][j], 0.0);
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
            mpf_t N1;
            mpf_init (N1);
            mpf_set_d (N1, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set (calcVar, N1);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x1) {
                        if (validReadLengthsHap1[i] <= (x1 - hetLocus)){
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else if ((j <= (x1-1)) && (j >= hetLocus)) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap1[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set_d (productsHap1[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                    mpf_set_d (productsHap1[i][j], 0.0);
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
            mpf_t N11;
            mpf_init (N11);
            mpf_set_d (N11, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set (calcVar, N11);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x1) {
                        if (validReadLengthsHap1[i] >= ((x1 - hetLocus) + 1)){
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else if ((j >= hetLocus) && (j <= x2)) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= x2 + 1) && (j <= (x1-1))) {
                        if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap1[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set_d (productsHap1[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                    mpf_set_d (productsHap1[i][j], 0.0);
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
            mpf_t N13;
            mpf_init (N13);
            mpf_set_d (N13, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set (calcVar, N13);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x1) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                        }
                    }
                    else if ((j >= hetLocus) && (j <= x2)) {
                        for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= x2 + 1) && (j <= (x1-1))) {
                        if (validReadLengthsHap1[i] >= ((j - hetLocus) + 1)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap1[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap1[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap1[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set_d (productsHap1[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                    mpf_set_d (productsHap1[i][j], 0.0);
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
            mpf_t N113;
            mpf_init (N113);
            mpf_set_d (N113, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
                mpf_set (calcVar, N113);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x2) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= hetLocus) && (j <= (x2-1))) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= (x2+1)) && (j <= x1)) {
                        if (validReadLengthsHap2[i] <= (j - hetLocus)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap2[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set_d (productsHap2[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                    mpf_set_d (productsHap2[i][j], 0.0);
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
            mpf_t N2;
            mpf_init (N2);
            mpf_set_d (N2, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set (calcVar, N2);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x2) {
                        if (validReadLengthsHap2[i] <= (x2 - hetLocus)){
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else if ((j >= hetLocus) && (j <= (x2-1))) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= (x2+1)) && (j <= x1)) {
                        if (validReadLengthsHap2[i] <= (j - hetLocus)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap2[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set_d (productsHap2[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                    mpf_set_d (productsHap2[i][j], 0.0);
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
            mpf_t N21;
            mpf_init (N21);
            mpf_set_d (N21, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set (calcVar, N21);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x2) {
                        if (validReadLengthsHap2[i] >= ((x2 - hetLocus) + 1)){
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                    }
                    else if ((j >= hetLocus) && (j <= (x2-1))) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= (x2+1)) && (j <= x1)) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                        }
                    }
                    else if (j < hetLocus) {
                        if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap2[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set_d (productsHap2[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                    mpf_set_d (productsHap2[i][j], 0.0);
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
            mpf_t N23;
            mpf_init (N23);
            mpf_set_d (N23, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set (calcVar, N23);
                mpf_mul (
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
                    mpf_set_d (calcVar, 0.0);
                    if (j == x2) {
                        if (validReadLengthsHap2[i] <= x2 - hetLocus) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else if ((j >= hetLocus) && (j <= (x2-1))) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                        }
                    }
                    else if ((j >= (x2+1)) && (j <= x1)) {
                        for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                            mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                        }
                    }   
                    else if (j < hetLocus) {
                        if (validReadLengthsHap2[i] <= ((genomeSize + j) - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                    else {
                        if (validReadLengthsHap2[i] <= (j - x2)) {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 1.0);
                            }
                        }
                        else {
                            for (int32_t k = 1; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                                mpf_set_d (generatingFunctionsHap2[i][j][k], 0.0);
                            }
                        }
                    }
                }
            }

            // Product polynomials
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set_d (productsHap2[i][0], 1.0);
                for (int32_t j = 1; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                    mpf_set_d (productsHap2[i][j], 0.0);
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
            mpf_t N213;
            mpf_init (N213);
            mpf_set_d (N213, 1.0);
            mpf_set_d (calcVar, 0.0);
            
            for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
                mpf_set (calcVar, N213);
                mpf_mul (
                    N213, 
                    calcVar, 
                    productsHap2[i][readDistHap2[validReadLengthsHap2[i]]]
                );
            }



            // Computing (aggregate) the numerator
            mpf_t tempVar1;
            mpf_init (tempVar1);
            mpf_t tempVar2;
            mpf_init (tempVar2);
            mpf_t tempVar3;
            mpf_init (tempVar3);
            mpf_t tempVar4;
            mpf_init (tempVar4);
            mpf_t prodVar1;
            mpf_init (prodVar1);
            mpf_t prodVar2;
            mpf_init (prodVar2);
            mpf_t coverageGapPermutations;
            mpf_init (coverageGapPermutations);

            mpf_sub (tempVar1, N1, N11);
            mpf_sub (tempVar2, N2, N21);
            mpf_mul (prodVar1, tempVar1, tempVar2);

            mpf_sub (tempVar3, N13, N113);
            mpf_sub (tempVar4, N23, N213);
            mpf_mul (prodVar2, tempVar3, tempVar4);

            mpf_sub (coverageGapPermutations, prodVar1, prodVar2);

            mpf_set (calcVar, aggregateCount[omp_get_thread_num()]);
            mpf_add (aggregateCount[omp_get_thread_num()], calcVar, coverageGapPermutations);

            mpf_clear (tempVar1);
            mpf_clear (tempVar2);
            mpf_clear (tempVar3);
            mpf_clear (tempVar4);
            mpf_clear (prodVar1);
            mpf_clear (prodVar2);
            mpf_clear (coverageGapPermutations);
            mpf_clear (N1);
            mpf_clear (N11);
            mpf_clear (N13);
            mpf_clear (N113);
            mpf_clear (N2);
            mpf_clear (N21);
            mpf_clear (N23);
            mpf_clear (N213);
        }

        mpf_clear(calcVar);

        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            for (int32_t j = 1; j <= genomeSize; j++) {
                for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                    mpf_clear (generatingFunctionsHap1[i][j][k]);
                }
                delete generatingFunctionsHap1[i][j];
            }
        }
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            for (int32_t j = 1; j <= genomeSize; j++) {
                for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                    mpf_clear (generatingFunctionsHap2[i][j][k]);
                }
                delete generatingFunctionsHap2[i][j];
            }
        }
    
        for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
            for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
                mpf_clear (productsHap1[i][j]);
            }
            delete productsHap1[i];
        }
        for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
            for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
                mpf_clear (productsHap2[i][j]);
            }
            delete productsHap2[i];
        }


    }

    // Computing the probability
    mpf_t numerator;
    mpf_init (numerator);

    mpf_t denominator;
    mpf_init (denominator);

    for (int32_t i = 0; i < threadCount; i++) {
        mpf_set (cVar, numerator);
        mpf_add (numerator, cVar, aggregateCount[i]);
    }

    mpf_t probability;
    mpf_init (probability);

    mpf_mul (denominator, totalReadSeqOutputsHap1, totalReadSeqOutputsHap2);
    mpf_div (probability, numerator, denominator);

    std::cout << "Total number of read sequencing outputs: " << denominator << std::endl; 

    std::cout << "probability = " << probability << std::endl;
    
    // Free initialised variables
    for (int32_t i = 0; i < threadCount; i++) {
        mpf_clear (aggregateCount[i]);
    }
    delete aggregateCount;
    
    mpf_clear (probability);
    mpf_clear (cVar);
    mpf_clear (numerator);
    mpf_clear (denominator);
    mpf_clear (totalReadSeqOutputsHap1);
    mpf_clear (totalReadSeqOutputsHap2);

    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 0; j <= readDistHap1[validReadLengthsHap1[i]]; j++) {
            mpf_clear (prodHap1[i][j]);
        }
        delete prodHap1[i];
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 0; j <= readDistHap2[validReadLengthsHap2[i]]; j++) {
            mpf_clear (prodHap2[i][j]);
        }
        delete prodHap2[i];
    }

    for (int32_t i = 0; i < distinctReadLengthsHap1; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap1[validReadLengthsHap1[i]]; k++) {
                mpf_clear (genFuncHap1[i][j][k]);
            }
            delete genFuncHap1[i][j];
        }
    }
    for (int32_t i = 0; i < distinctReadLengthsHap2; i++) {
        for (int32_t j = 1; j <= genomeSize; j++) {
            for (int32_t k = 0; k <= readDistHap2[validReadLengthsHap2[i]]; k++) {
                mpf_clear (genFuncHap2[i][j][k]);
            }
            delete genFuncHap2[i][j];
        }
    }


return 0;
}