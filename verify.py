#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import os
import subprocess
import pandas as pd
from datetime import datetime
from pathlib import Path
import libverify

LIM_SECOND= 20 * 60

re_expr = re.compile(r'Valid dat file. (?P<dist>\d+\.\d+)')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("exec", help="executable file name")
    parser.add_argument("tsps", nargs='+', help="tsp files")
    parser.add_argument(
        "-o", "--output-directory", required=False, default=None,
        help="output directory name"
    )
    args = parser.parse_args()

    exe_path = Path(args.exec).resolve()
    if not exe_path.exists() or exe_path.is_dir():
        sys.stderr.write("ERROR: {0} does not exist.\n".format(exe_path))
        sys.exit(1)

    current_dir = Path.cwd().resolve()
    if args.output_directory is None:
        output_dir_parent = current_dir / (exe_path.stem + "_d")
    else:
        output_dir_parent = Path(args.output_directory).resolve()

    stats = []

    for tsp in args.tsps:
        tsp_path = Path(tsp).resolve()
        if not tsp_path.exists():
            continue

        output_dir = output_dir_parent / tsp_path.stem
        output_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(str(output_dir))
        print("\n**************************************************")
        print("Executing: {0} {1}".format(args.exec, tsp))

        dats = execute(exe_path, tsp_path, LIM_SECOND);
        if dats:
            dat_path, t_delta = dats[-1]
            dat_path = Path(dat_path).resolve()
            process = libverify.check(tsp_path, dat_path)
            stdout, stderr = process.communicate()

            stdout = stdout.decode('utf-8').strip()
            stderr = stderr.decode('utf-8').strip()
            print(stdout, stderr)

            if (process.returncode == 0):
                match = re_expr.match(stdout)
                if match:
                    dist = match.group('dist')
                    print("{0} is valid. Distance: {1}, Time: {2}".format(
                        dat_path.name, dist, t_delta)
                    )
                    stats.append((tsp_path.name, dist, t_delta))
                else:
                    print("{0} is not valid.".format(dat_path.name))
        else:
            print("No tsp file generated over {0}".format(str(tsp_path.name)))
            stats.append((tsp_path.name, None, None))

        os.chdir(str(current_dir))

    stats_path = output_dir_parent / (exe_path.stem + ".csv")
    save_stats(stats, stats_path)

def execute(exe_path, tsp_path, timeout):
    t1 = datetime.now()
    process = subprocess.Popen([str(exe_path), str(tsp_path)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    dats = []
    while True:
        t_delta = (datetime.now() - t1).total_seconds()
        output = process.stdout.readline().decode('utf-8')
        if output:
            if output.startswith("tour"):
                e = output.find("dat") + 3
                datname = output[:e]
                dats.append((datname, t_delta))
            print(output.strip())
        elif process.poll() is not None:
            print("[Finished '{0}']".format(tsp_path.stem), end=" ", flush=True)
            break

        if t_delta >= timeout:
            process.kill()
            print("[Terminated '{0}']".format(tsp_path.stem), end=" ", flush=True)
            break
    return dats

def save_stats(stats, stats_path):
    df = pd.DataFrame([e[1:] for e in stats],
                      columns=["Distance", "Time"],
                      index=[e[0] for e in stats])
    df.index.name = str(datetime.now())
    with stats_path.open('a') as f:
        f.write(df.to_csv())

if __name__ == '__main__':
    main()
