#include <iostream>
#include <vector>
#include <fstream>

#include "suk/SolidKmers.hpp"
#include "kmc_file.h"
#include "kmer_api.h"
#include <zlib.h>

std::string reverseKmer(std::string kmer) {
    std::string rev;
    for (int i = 16; i >= 0; i--) {
        char c;
        switch (kmer[i]) {
        case 'A': 
            c = 'T';
            break;
        case 'T': 
            c = 'A';
            break;
        case 'C': 
            c = 'G';
            break;
        case 'G':
                c = 'C';
            break;
        }
        rev.push_back(c);
    }
    return rev;
}

int plotKmcKmerHist(CKMCFile kmc_database, std::string kmc_output_path, 
                    std::string output_path, size_t histSize) {
    if(!kmc_database.OpenForListing(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return false;
    }
	CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    // size_t histSize = 1000;
    
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
	
    std::vector<size_t> histArray(histSize + 1);
    uint64 counter;
    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter <= histSize) histArray[counter]++;
    }
    std::ofstream histogram(output_path);
    for (int i : histArray)
        histogram << i << std::endl;
    histogram.close();
    return 0;
}

int plotRefKmerHist(CKMCFile kmc_database, std::string kmc_output_path, 
                    std::string input_path, std::string output_path, size_t histSize) {
    if (!kmc_database.OpenForRA(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return -1;
    }
    CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    // size_t histSize = 1000;

    std::vector<size_t> histArray(histSize, 0);
    CKmerAPI* kmer_object = new CKmerAPI(kmc_info.kmer_length);
    std::ifstream input(input_path);
    std::string line;
    while(std::getline(input, line)) {    
        // if (line.empty()) break;
            // std::cout << line << std::endl;
        CKmerAPI* kmer_object = new CKmerAPI(kmc_info.kmer_length);
        kmer_object->from_string(line);
        uint64 counter = 0;
        kmc_database.CheckKmer(*kmer_object, counter);
        if (counter != 0) {
            histArray[counter]++;
            continue;
	    }
        // std::cout << counter << std::endl;
        std::string rev = reverseKmer(line);
        // std::cout << rev << std::endl;
        kmer_object->from_string(rev);
            counter = 0;
            kmc_database.CheckKmer(*kmer_object, counter);
            histArray[counter]++;
    }
    input.close();
    std::ofstream output(output_path);
    for (auto i : histArray) {
        output << i << std::endl;
    }
    output.close();
    return 0;
}

int main(int argc, char* argv[]) {
	
    if (argc < 2) {
        std::cout << "use: suk_cnt k [path]" << std::endl;
    }
    CKMCFile kmc_database;
    std::string path = argc > 2 ? argv[2] : "./";
    std::string k = argv[1];
    std::string kmc_output_path = path + k + "mer" + "/kmc_result.res";

    size_t histSize = 1000;
    
    /* Plot histogram */
    std::string file = path + "ref_" + k + "mer.txt";
    std::string output = path + "hist/" + k + "mer_hist.txt";
    plotRefKmerHist(kmc_database, kmc_output_path, file, output, histSize);
    
    file = path + "ref_" + k + "mer_rep.txt";
    output = path + "hist/" + k + "mer_rep_hist.txt";
    plotRefKmerHist(kmc_database, kmc_output_path, file, output, histSize);
    
    output = path + "hist/" + k + "mer_kmc_hist.txt";
    plotKmcKmerHist(kmc_database, kmc_output_path, output, histSize);

    return 0;
}