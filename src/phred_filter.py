import os
import sys
import gzip
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

def filter_line(bps, qs, q0=13):
    res = []
    start = 0
    end = -2
    for i, q in enumerate(qs[:-1]):
        if ord(q) - ord('!') >= q0:
           end = i
        elif end == i - 1:
            if i - start > 17:
                res.append(">\n")
                res.append(bps[start: i])
                res.append("\n")
            start = i + 1
        else:
            start = i + 1
    if start < end and len(qs) - start > 17:
        res.append(">\n")
        res.append(bps[start: -1])
        res.append("\n")
    return "".join(res)


def filter_file(thread_no, src, start, end, q0=13):
    batch = 1e3
    res = []
    codec = os.path.splitext(src)[1]
    if codec == ".gz":
        fin = gzip.open(src, "rt")
    else:
        fin = open(src, "rt")
    fin.seek(start)
    while True:
        if fin.tell() > end:
            break
        header = fin.readline()
        if not header:
            break
        while header[0] != '@':
           header = fin.readline()
        # print(str(thread_no), fin.tell())
        bps = fin.readline()
        fin.readline()
        qs = fin.readline()
        filtered = filter_line(bps, qs)
        if len(filtered) >= 20:
            res.append(filtered)
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

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("use: python phred_filter.py source_path dest_path")
    # for i, arg in enumerate(sys.argv):
    #   print(f"Argument {i}: {arg}")
    src = sys.argv[1]
    dst = sys.argv[2]
    number = 24
    size = os.path.getsize(src)
    size_per_thread = int(size/number)
    print(f"size of {src}: {size}")
    print(f"no. of threads: {number}")
    print(f"initing threads...")
    writer = Thread(target=file_writer, args=(dst, queue), daemon=True)
    writer.start()
    
    threads = [Thread(
                    target=filter_file, \
                    args=(i, src, size_per_thread * i, size_per_thread * (i+1))) \
                for i in range(number)]
    for thread in threads:
        thread.start()
    print(f"threads running...")
    
    for thread in threads:
        thread.join()
    queue.join()
    print(f"done. filtered reads are saved to: {dst}")
    
