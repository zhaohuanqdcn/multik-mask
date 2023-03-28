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

int writeRefKmersInRange(CKMCFile kmc_database, std::string kmc_output_path, 
                        std::string input_paths[], std::string output_path, 
                        size_t lower, size_t upper, size_t size) {
    if (!kmc_database.OpenForRA(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return -1;
    }
    CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    int kcount = 0;

    std::string s;
    s.reserve(40 * 1e4);
    CKmerAPI* kmer_object = new CKmerAPI(kmc_info.kmer_length);
    std::ofstream output(output_path);
    std::string line;

    for (int i = 0; i < size; i++) {
        std::ifstream input(input_paths[i]); 
        while(std::getline(input, line)) { 
            // if (line.empty()) break;
                // std::cout << line << std::endl;
            CKmerAPI* kmer_object = new CKmerAPI(kmc_info.kmer_length);
            kmer_object->from_string(line);
            uint64 counter0 = 0;
            kmc_database.CheckKmer(*kmer_object, counter0);
            // std::cout << counter << std::endl;
            std::string rev = reverseKmer(line);
            kmer_object->from_string(rev);
            uint64 counter1 = 0;
            kmc_database.CheckKmer(*kmer_object, counter1);
            // std::cout << counter1 << std::endl;
            
            // batch flushing
            if ((counter0 < upper && counter0 >= lower) || 
                (counter1 < upper && counter1 >= lower)) {
                kcount += 1;
                s.append(line + "\n" + rev + "\n"); // write in both seq and its rev
                if (kcount >= 1e4) {
                    kcount = 0;
                    output << s;
                    s = "";
                }
            }
        }
        kcount = 0;
        output << s;
        s = "";
        input.close();
    }
    output.close();
    return 0;
}

int writeErrorKmersInRange(CKMCFile kmc_database, std::string kmc_output_path, 
                         std::string input_path, std::string output_path,
                         size_t lower, size_t upper) {

    std::ifstream input(input_path);
    std::string line;
    std::set<std::string> ref_kmers;

    while(std::getline(input, line)) { 
        ref_kmers.insert(line);
    }
    input.close();
    
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
    s.reserve(25 * 1e4);

    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter < upper && counter >= lower) {
            std::string kmer = kmer_object->to_string();
            if (ref_kmers.find(kmer) == ref_kmers.end()) { // not in refs
                kcount += 1;
                s.append(">\n" + kmer + "\n");
                if (kcount >= 1e4) {
                    kcount = 0;
                    output << s;
                    s = "";
                }
            }
        }
    }
    output << s;
    output.close();
    return 0;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
        std::cout << "use: suk_err k lower upper" << std::endl;
    }
    CKMCFile kmc_database;
   std::string k = argv[1];
    std::string kmc_output_path = k + "mer/kmc_result.res";

    /* Write error kmers at peak */
    std::string file1 = "ref_" + k + "mer.txt";
    std::string file2 = "ref_" + k + "mer_rep.txt";
    std::string file3 = "peak_" + k + "mer_ref.txt";
    std::string file4 = "peak_" + k + "mer_err.fa";
    std::string inputs[] {file1, file2};
    size_t lower = std::stoi(argv[2]);
    size_t upper = std::stoi(argv[3]);
    writeRefKmersInRange(kmc_database, kmc_output_path, inputs , file3, lower, upper, 2);
    writeErrorKmersInRange(kmc_database, kmc_output_path, file3, file4, lower, upper);

    return 0;
}