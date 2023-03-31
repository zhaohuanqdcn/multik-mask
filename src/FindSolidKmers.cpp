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

int writeKmcKmersInRange(CKMCFile kmc_database, std::string kmc_output_path, 
                         std::string output_path, size_t lower, size_t upper) {
    
    if(!kmc_database.OpenForListing(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return false;
    }
	CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
    std::ofstream output(output_path);
    uint64 counter;

    int kcount = 0;
    std::string s;
    s.reserve(30 * 1e4);
    
    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter < upper && counter >= lower) {
            std::string kmer = kmer_object->to_string();
            kcount += 1;
            s.append(">\n" + kmer + "\n");
            if (kcount >= 1e4) {
                kcount = 0;
                output << s;
                s = "";
            }
        }
    }
    output << s;
    output.close();
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "use: suk_kmer k lower upper [path]" << std::endl;
        return 0;
    }
    CKMCFile kmc_database;
    std::string path = argc > 4 ? argv[4] : "./";
    std::string k = argv[1];
    std::string kmc_output_path = path + k + "mer/kmc_result.res";
    // std::cout << kmc_output_path << std::endl;
    
    /* Take KMC kmers in range */
    std::string output = path + k + "mers_solid.fa";
    size_t lower = std::stoi(argv[2]);
    size_t upper = std::stoi(argv[3]);

    std::cout << "writing solid " << k << "-mers to: " << output << std::endl;
    writeKmcKmersInRange(kmc_database, kmc_output_path, output, lower, upper);
    
    return 0;
}