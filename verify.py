#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import subprocess
from datetime import datetime
from pathlib import Path

LIM_SECOND= 20 * 60

parser = argparse.ArgumentParser()
parser.add_argument("exec", help="executable file name")
parser.add_argument("tsps", nargs='+', help="tsp files")
parser.add_argument(
    "-o", "--output-directory", required=False, default=None,
    help="output directory name"
)
args = parser.parse_args()

path_exec = Path(args.exec).resolve()
if not path_exec.exists() or path_exec.is_dir():
    sys.stderr.write("ERROR: {0} does not exist.\n".format(path_exec))
    sys.exit(1)

current_dir = Path.cwd().resolve()
if args.output_directory is None:
    output_dir_parent = current_dir / Path(path_exec.stem)
else:
    output_dir_parent = Path(args.output_directory).resolve()

for tsp in args.tsps:
    tspfile = Path(tsp).resolve()
    if not tspfile.exists():
        continue

    output_dir = output_dir_parent / Path(tspfile.stem)
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(str(output_dir))
    print("Executing: {0} {1}".format(args.exec, tsp))

    t1 = datetime.now()
    process = subprocess.Popen([str(path_exec), str(tspfile)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    datfile = ""
    while True:
        t_delta = datetime.now() - t1
        output = process.stdout.readline().decode('utf-8')
        if output:
            if output.startswith("tour"):
                e = output.find("dat") + 3
                datfile = output[:e]
            print(output.strip())
        elif process.poll() is not None:
            print("Finished in {0} sec:".format(t_delta), end=" ")
            break

        if t_delta.seconds > LIM_SECOND:
            process.kill()
            print("Terminated:", end=" ")
            break
    print("The last dat file is", datfile)
    os.chdir(str(current_dir))
