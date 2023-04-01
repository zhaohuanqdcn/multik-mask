import os
import sys
from queue import Queue
from threading import Lock
from threading import Thread

def file_writer(filepath, queue):
    buffer_size = 3*1024*1024 # 3M
    with open(filepath, 'at', buffering=buffer_size) as file:
        while True:
            line = queue.get()
            file.write(line)
            file.flush()
            queue.task_done()

encode = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
decode = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
compl = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def complement(kmer):
    res = []
    for c in kmer[::-1]:
        res.append(compl[c])
    return ''.join(res)

def to_int(kmer):
    res = 0
    for i in kmer:
        res = res << 2
        res += encode[i]
    return res

def from_int(k, kmer):
    res = []
    while kmer > 0:
        res.append(decode[kmer & 3])
        kmer = (kmer >> 2)
    while len(res) < k:
        res.append('A')
    return ''.join(res[::-1])

def mask_line(bps):
    n = len(bps) - 1
    res = []
    next = 0
    seq = [0 for k in K]
    for i in range(n):
        for j in range(len(K)):
            k = K[j]
            seq[j] = seq[j] << 2
            seq[j] += encode[bps[i]]
            seq[j] = seq[j] & ((1 << (2 * k)) - 1)
            if i >= k - 1 and seq[j] in kmers:
                res.append(bps[max(next, i - k + 1): i + 1])
                next = i + 1
                # print(from_int(k, seq[j]))
    return "".join(res)


def mask_file(thread_no, src, start, end):
    batch = 1e3
    res = []
    fin = open(src, "rt")
    fin.seek(start)
    while True:
        if fin.tell() > end:
            break
        header = fin.readline()
        if not header:
            break
        while header[0] != '>':
           header = fin.readline()
        # print(str(thread_no), fin.tell())
        bps = fin.readline()
        masked = mask_line(bps)
        if len(masked) > 0:
            res.append(">\n")
            res.append(masked)
            res.append("\n")
        if len(res) > batch:
            queue.put(''.join(res))
            res = []
            if thread_no == 0: # showing progress on thread 0
                print(f"estimated progress: {int(100 * fin.tell() / end)}%", end="\r")
            # break
    fin.close()
    queue.put(''.join(res))
    res = []
    return
    

lock = Lock()
queue = Queue()
kmers = set()
K = list()

if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) % 2 != 1:
        print("use: python mask.py filter_path dest_path (k kmer_path)+")
    # for i, arg in enumerate(sys.argv):
    #   print(f"Argument {i}: {arg}")

    src = sys.argv[1]
    dst = sys.argv[2]

    for i in range(3, len(sys.argv), 2):
        k = eval(sys.argv[i])
        K.append(k)
        kmer_path = sys.argv[i + 1]
        print(f"reading in solid {k}-mers...")
        count = 0
        with open(kmer_path, "rt") as file:
            for kmer in file.read().splitlines()[1::2]:
                count += 1
                kmers.add(to_int(kmer))
                kmers.add(to_int(complement(kmer)))
        print(f"total # of solid {k}-mers: {count}")
    
    number = 48
    size = os.path.getsize(src)
    size_per_thread = int(size/number)
    print(f"no. of threads: {number}")
    print(f"initing threads...")

    writer = Thread(target=file_writer, args=(dst, queue), daemon=True)
    writer.start()
    threads = [Thread(
                    target=mask_file, \
                    args=(i, src, size_per_thread * i, size_per_thread * (i+1))) \
                for i in range(number)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    queue.join()
    print(f"done. masked reads are saved to: {dst}")
    