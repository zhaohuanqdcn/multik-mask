#include <iostream>
#include <vector>
#include <fstream>

#include "suk/SolidKmers.hpp"
#include "kmc_file.h"
#include "kmer_api.h"
#include <zlib.h>

std::string reverseKmer(size_t k, std::string kmer) {
    std::string rev;
    for (int i = k - 1; i >= 0; i--) {
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

int plotRefKmerHist(size_t k, CKMCFile kmc_database, std::string kmc_output_path, 
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
        std::string rev = reverseKmer(k, line);
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
        return 0;
    }
    CKMCFile kmc_database;
    std::string path = argc > 2 ? argv[2] : "./";
    std::string k_str = argv[1];
    std::string kmc_output_path = path + k_str + "mer" + "/kmc_result.res";

    size_t k = std::stoi(k_str);
    size_t histSize = 1000;
    
    /* Plot histogram */
    // std::string file_ref = path + "ref_" + k_str + "mer.txt";
    // std::string output_ref = path + "hist/" + k_str + "mer_hist.txt";
    // plotRefKmerHist(k, kmc_database, kmc_output_path, file_ref, output_ref, histSize);
    
    // std::string file_rep = path + "ref_" + k_str + "mer_rep.txt";
    // std::string output_rep = path + "hist/" + k_str + "mer_rep_hist.txt";
    // plotRefKmerHist(k, kmc_database, kmc_output_path, file_rep, output_rep, histSize);
    
    std::string output_path = path + "hist/" + k_str + "mer_kmc_hist.txt";
    std::vector<size_t> histArray(histSize + 1);
    if(!kmc_database.OpenForListing(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return false;
    }
	CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
	
    uint64 counter;
    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter <= histSize) histArray[counter]++;
    }
    std::ofstream histogram(output_path);
    for (int i : histArray)
        histogram << i << std::endl;
    histogram.close();
    printf("histogram of %ld-mer saved to: ", k);
    std::cout << output_path << std::endl;
    
    // idenfity cutoffs
    int sum = 0, prefix = 0; 
    float min = 1;
    int l = 0, r = 0; 
    for (int i = 1; i < int(histSize / 5); i++)
        if (histArray[i - 1] != 0) {
            float ratio = 1 - (histArray[i] / histArray[i - 1]);
            if (ratio < min) {
                min = ratio;
                l = i;
            }
        }
    for (int i = l; i < histSize; i++) 
        sum += histArray[i];
    for (int i = l; i < histSize; i++) {
        prefix += histArray[i];
        if (prefix >= float(sum) * 0.95) {
            r = i;
            break;
        }
    }
    
    printf("identified cutoff of %ld-mer: (%d, %d)\n", k, l, r);
    std::string cutoff_path = path + "hist/" + k_str + "mer_cutoff.txt";
    std::ofstream cutoff(cutoff_path);
    cutoff << l << std::endl;
    cutoff << r << std::endl;
    cutoff.close();
    printf("cutoff of %ld-mer saved to: ", k);
    std::cout << cutoff_path << std::endl;
    
    return 0;
}