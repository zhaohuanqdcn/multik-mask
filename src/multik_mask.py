import os
import sys
import argparse
from subprocess import PIPE, Popen
from threading import Lock
from threading import Thread

def find_solid_kmers(filtered, k, output_dir):
    cmd_kmc = ["../build/bin/suk", "-i", filtered, "-k", str(k), "-f",  output_dir]
    p = Popen(cmd_kmc, stdout=PIPE)
    p.wait()
    cmd_cnt = ["../build/bin/suk_cnt", str(k), output_dir]
    p = Popen(cmd_cnt, stdout=PIPE)
    p.wait()
    cutoff = open(f"{output_dir}/hist/{k}mer_cutoff.txt", "rt")
    l = eval(cutoff.readline())
    r = eval(cutoff.readline())
    cmd_kmer = ["../build/bin/suk_kmer", str(k), str(l), str(r), output_dir]
    p = Popen(cmd_kmer, stdout=PIPE)
    p.wait()
    return


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", type=int, nargs="+", dest="K", help="A list of kmer sizes", required=True)
    parser.add_argument("-i", type=str, dest="input_path", help="Input file of reads", required=True)
    parser.add_argument("-f", type=str, dest="output_folder", help="Output folder name", required=True)
    parser.add_argument("-o", type=str, dest="prefix_name", default="assembly", help="Prefix of output files (default: assembly)", required=False)
    parser.add_argument("--no-pipeline", default=False, action="store_true", help="Mask reads without assembling", required=False)
    args = parser.parse_args()
    
    K = args.K
    input_path = args.input_path
    output_dir = args.output_folder
    prefix = args.prefix_name
    # print(K, input_path, output_dir, prefix)
    os.mkdir(output_dir)
    os.mkdir(output_dir + "/hist")
    filtered = output_dir + "/filtered.fa"
    masked = output_dir + "/masked.fa"

    # Phred filter
    print("phred filtering...")
    cmd_filter = ["python", "phred_filter.py", input_path, filtered]
    p = Popen(cmd_filter, stdout=PIPE)
    p.wait()
    print(f"filtered reads are saved to: {filtered}")

    # solid kmers
    print("getting kmers...")
    threads = [Thread(
                    target=find_solid_kmers, \
                    args=(filtered, k, output_dir)) \
                for k in K]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    print(f"solid kmers are saved to: {output_dir}/kmers_solid.fa")
    
    # masking
    cmd_mask = ["python", "mask.py", filtered, masked]
    for k in K:
        cmd_mask.append(str(k))
        cmd_mask.append(f"{output_dir}/{k}mers_solid.fa")
    p = Popen(cmd_mask, stdout=PIPE)
    p.wait()
    print(f"masked reads are saved to: {masked}")
    if args.no_pipeline:
        exit
    print(f"Pipeline starting...")
    cmd_minia = ["../gatb-minia-pipeline/gatb", "-s", masked, "--kmer-sizes", "21", "-o", prefix, "--no-error-correction"]
    p = Popen(cmd_mask, stdout=PIPE)
    p.wait()
    print(f"done. assembly is saved: {output_dir}/{prefix}_final.contigs.fa")

