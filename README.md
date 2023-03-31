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

- (KMC and GATB-minia should be in the path (`$PATH`))

## Building
CmakeLists is provided in the project root folder.
### Building Library
Run the following commands to build a library `suk.a` in `build/lib` :
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
placeholder
```

The usage of subroutines is not recommended since they are not intended to be standalone components. However, in case you want to do so:

 ```bash
  python phred_filter.py src_path dest_path
  bin/suk_cnt k [path (default ./)]
  bin/suk_kmer k lower upper [path (default ./)]
  python mask.py k kmer_path filtered_path dest_path
 ```
**Notes:** 

- The input file should be in `.fastq` or `.fq` format. Compressed variant `.fq.gz` is also accepted, but the performance may be compromised as `gzip` format does not support random access.
- `suk` has been modified to integrate into the pipeline. Please refer to [SUK](https://github.com/Ritu-Kundu/suk) instead if you want to use `suk` alone.

### Output Files

placeholder

### Example  

placeholder

## External Libraries

 * [sdsl](https://github.com/simongog/sdsl-lite) library has been used for bit-vector implementation.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.



