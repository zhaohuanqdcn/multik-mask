# Multi-k Error Masking
This repository presents a multi-k error masking tool for high-coverage long reads.

The framework is updated from [SUK](https://github.com/Ritu-Kundu/suk).

## Dependencies
- Linux 64 bits

- gcc 4.8+ or clang 3.4+

- python 3.7+ 

- python >=2.7 and <3 (for BESST)

- cmake 3.10+

- [BESST](https://github.com/ksahlin/BESST)

  ```
  python2 -m pip install --user BESST
  ```

- [KMC3](https://github.com/refresh-bio/KMC)

  ```
  conda install -c bioconda kmc
  ```

- [GATB-minia](https://github.com/GATB/gatb-minia-pipeline)

  ```
  git clone --recursive https://github.com/GATB/gatb-minia-pipeline
  cd gatb-minia-pipeline; make test

- (KMC, GATB-minia and MultiKMasking should be in path (`$PATH`))

## Building
CmakeLists is provided in the project root folder.

Run the following commands to install and build:
```console
  git clone --recursive https://github.com/zhaohuanqdcn/multik-mask
  cd multik-mask
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make
```
**Notes:** 

* If `--recursive` was omitted from `git clone`, please run `git submodule init` and `git submodule update` before compiling.
* If the target machine is different from the one on which Hypo is being compiled, exclude the flag `-Doptimise_for_native=ON`.

## Usage of the Tool

To run the tool from `multik-mask/src` folder:

```
multik_mask.py [-h] -k K+ -i INPUT_PATH -f OUTPUT_FOLDER [-o PREFIX_NAME] [--no-pipeline]

optional arguments:
  -h, --help        show this help message and exit
  -k K [K ...]      A list of kmer sizes
  -i INPUT_PATH     Input file of reads
  -f OUTPUT_FOLDER  Output folder name
  -o PREFIX_NAME    Prefix of output files (default: assembly)
  --no-pipeline     Mask reads without assembling
```

Subroutines are not recommended since they are not intended to be standalone components. However, in case you want to run a component separately:

 ```bash
  python phred_filter.py src_path dest_path
  bin/suk_cnt k [path (default ./)]
  bin/suk_kmer k lower upper [path (default ./)]
  python mask.py filter_path dest_path (k kmer_path)+
 ```
**Notes:** 

- The input file should be in `.fastq` or `.fq` format. Compressed variant `.fq.gz` is also accepted, but the performance may be compromised as `gzip` format does not support random access.
- `suk` has been modified to integrate into the pipeline. Please refer to [SUK](https://github.com/Ritu-Kundu/suk) instead if you want to use `suk` alone.

### Output Files

All output files will reside in `<OUTPUT_FOLDER>`. Masked reads are stored in `/masked.fa` and the assembly is stored in `/<PREFIX_NAME>_final.contigs.fa` (if the pipeline is run). The KMC database of `k`mers can be found at `/<k>mer` and its frequency histogram will be stored in `/hist/<k>mer_kmc.hist.txt`.

### Example  

`python multik_mask.py -k 13 17 21 -i F2426_reads.fq -f output -o f2426`

The above command runs MultiKMasking with `k=[13,17,21]` on reads `F2426_reads.fq` stored in the same folder (`multik-mask/src`), and directs the output to folder `output`. The assembly will be available at `./output/f2426_final.contigs.fa`.

## External Libraries

 * [sdsl](https://github.com/simongog/sdsl-lite) library has been used for bit-vector implementation.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.



