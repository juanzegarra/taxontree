import os
import subprocess
import time
import sys
import re
import multiprocessing
from pathlib import Path
current_dir = os.getcwd()
directories = [d for d in os.listdir(current_dir) if os.path.isdir(d)]

print(directories)

def trim_alignment(alignment_file, out_dir):
    input_name = Path(alignment_file).stem
    trimmed_output = out_dir / f"{input_name}_trimmed.fasta"

    start = time.perf_counter()

    command = [
        "trimal",
        "-in", str(alignment_file),
        "-out", str(trimmed_output),
        "-strictplus",
    ]

def monitor_resources(pid, program, iteration):
    pidstat_cmd = ["pidstat", "-T", "ALL", "-H", "-u", "-r", "-p", str(pid), "1"]
    log = f"resource_log_{program}_i{iteration}.txt"
    with open(log,"w") as log_file:
        process = subprocess.Popen(pidstat_cmd, stdout=log_file, stderr=subprocess.STDOUT)
    return process, log

def run_command(program, input_file, i):
    print(f"running {program} on file {input_file}")
    if program == "fasttree" or program == "veryfasttree":
        output = input_file.replace("_aln.fasta", f"_{program}.nwk")
    else:
        output = input_file.replace(".fasta", f"_{program}_aln.fasta")

    # Use micromamba run -n <env> to execute within the correct environment
    famsa2_command = ["micromamba", "run", "-n", "famsa_env", "famsa", input_file, output, "-t", "8", "-refine_mode", "auto"]
    muscle5_command = ["micromamba", "run", "-n", "famsa_env", "muscle", "-super5", input_file, "-out", output, "-t", "8"]
    muscle5_1core_command = ["micromamba", "run", "-n", "famsa_env", "muscle", "-super5", input_file, "-out", output, "-t", "1"]
    muscle3_command = ["micromamba", "run", "-n", "muscle_env", "muscle", "-in", input_file, "-out", output, "-maxiters", "8"]
    fasttree_command = ["fasttree", str(input_file)]
    veryfast_command = ["veryfasttree", str(input_file)]

    start = time.perf_counter()

    process = None
    try:
        if program == "fasttree":
            with open(output, "w", encoding="utf-8") as fout:
                process = subprocess.Popen(
                    fasttree_command,
                    stdout=fout,
                    stderr=subprocess.PIPE
                )
        elif program == "veryfasttree":
            with open(output, "w", encoding="utf-8") as fout:
                process = subprocess.Popen(
                    veryfast_command,
                    stdout=fout,
                    stderr=subprocess.PIPE
                )
        elif program == "muscle5":
                process = subprocess.Popen(muscle5_command)

        elif program == "muscle5_1core":
                process = subprocess.Popen(muscle5_1core_command)

        elif program == "muscle3":
               process = subprocess.Popen(muscle3_command)

        elif program == "famsa":
               process = subprocess.Popen(famsa2_command)
        else:
            print(f"{program} program not in list of executable programs, exiting...")
            sys.exit()

        pid, log_file = monitor_resources(process.pid, program,i)

    except subprocess.CalledProcessError as e:
        print(f"Error performing {program} on input {input_file}: {e}")
        output = None


    elapsed = time.perf_counter() - start if output else "Nan"
    return output, elapsed

def multiprocess_worker(program, seq_file, aligned_out, directory, i):
    output, exec_time = run_command(program, seq_file, i)
    aligned_out.add((directory, output))
    result = f"{directory}\t{program}\t{i}\t{exec_time}\n"
    return result

def main():
    iters = 10
    aln_programs = ["muscle3", "famsa", "muscle5", "muscle5_1core"]
    aligned_out = set()

    # Open files outside the loop so they aren't overwritten for each directory
    with open("summary_times_programs_msa.tsv", "w") as f_msa:
        f_msa.write("Directory\tProgram\tIteration\tExecution\n")

        for directory in directories:
            prefix = directory
            seq_file = os.path.join(directory, f"{prefix}_blast_output_sequences.fasta")

            tasks = []
            # Check if path actually exists
            if os.path.isfile(seq_file):
                for program in aln_programs:
                    for i in range(iters):
                        # Create tasks to be executed by worker
                        tasks.append((program,seq_file, aligned_out, directory, i))
            with multiprocessing.Pool(processes=2) as pool:
                results = pool.starmap(multiprocess_worker, tasks)

            for line in results:
                f_msa.write(line)

    tree_programs = ["fasttree","veryfasttree"]
    # Run FastTree on the generated alignments

    with open("summary_times_trees.tsv", "a") as f_trees:
        f_trees.write("Directory\tAlignment\tIteration\tExecution\n")
        for directory, aln in aligned_out:
            # Fixed regex to capture 1 or more alphanumeric characters
            aln_program_match = re.search(r"_([A-Za-z0-9_]+)_aln\.fasta", aln)

            if aln_program_match:
                aln_program = aln_program_match.group(1)
                for program in tree_programs:
                    for i in range(iters):
                        output, exec_time = run_command(program, aln)
                        f_trees.write(f"{directory}\t{program}\t{i}\t{exec_time}\n")


if __name__ == "__main__":
    main()
