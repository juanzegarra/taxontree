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
    # Ensure out_dir is a Path object so the '/' operator works
    trimmed_output = Path(out_dir) / f"{input_name}_trimmed.fasta"

    command = [
        "trimal",
        "-in", str(alignment_file),
        "-out", str(trimmed_output),
        "-strictplus",
    ]

    print(f"Trimming alignment: {alignment_file}")
    # Actually execute the trimal command
    subprocess.run(command, check=True)

    # Return as string because run_command uses .replace() on the input_file
    return str(trimmed_output)

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
    muscle5_command = ["micromamba", "run", "-n", "famsa_env", "muscle", "-super5", input_file, "-output", output, "-t", "8"]
    muscle5_1core_command = ["micromamba", "run", "-n", "famsa_env", "muscle", "-super5", input_file, "-output", output, "-t", "1"]
    muscle3_command = ["micromamba", "run", "-n", "muscle_env", "muscle", "-in", input_file, "-out", output, "-maxiters", "8"]
    fasttree_command = ["fasttree", str(input_file)]
    veryfast_command = ["VeryFastTree", str(input_file)]

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

        pid, log_file = monitor_resources(process.pid, program, i)

        # CRITICAL: Wait for the alignment/tree program to actually finish!
        process.wait()

    except subprocess.CalledProcessError as e:
        print(f"Error performing {program} on input {input_file}: {e}")
        output = None

    elapsed = time.perf_counter() - start if output else "Nan"
    return output, elapsed

def multiprocess_worker(args):
    program, seq_file, directory, i = args
    output, exec_time = run_command(program, seq_file, i)
    result = f"{directory}\t{program}\t{i}\t{exec_time}\n"
    # Return the result string AND the alignment output data
    return result, directory, output

def main():
    iters = 10
    aln_programs = ["muscle3", "famsa", "muscle5", "muscle5_1core"]
    aligned_out = set()

    with open("summary_times_programs_msa.tsv", "w") as f_msa:
        f_msa.write("Directory\tProgram\tIteration\tExecution\n")
        f_msa.flush() # Force write header immediately

        for directory in directories:
            prefix = directory
            seq_file = os.path.join(directory, f"{prefix}_blast_output_sequences.fasta")

            tasks = []
            if os.path.isfile(seq_file):
                for program in aln_programs:
                    for i in range(iters):
                        # Pack arguments into a single tuple; remove aligned_out
                        tasks.append((program, seq_file, directory, i))

            with multiprocessing.Pool(processes=2) as pool:
                # imap_unordered yields results the moment a worker finishes
                for line, dir_out, out_file in pool.imap_unordered(multiprocess_worker, tasks):
                    f_msa.write(line)
                    f_msa.flush() # Force write loop output immediately

                    if out_file: # Add to set in the main process
                        aligned_out.add((dir_out, out_file))

    tree_programs = ["fasttree","veryfasttree"]

    with open("summary_times_trees.tsv", "a") as f_trees:
        f_trees.write("Directory\tAlignment\tIteration\tExecution\n")
        f_trees.flush() # Force write header immediately

        for directory, aln in aligned_out:
            aln_program_match = re.search(r"_([A-Za-z0-9_]+)_aln\.fasta", aln)

            if aln_program_match:
                # aln is already the full path to the alignment file.
                trim_in = aln

                try:
                    # Run trimming ONCE per alignment file
                    trim_out = trim_alignment(trim_in, directory)

                    for program in tree_programs:
                        for i in range(iters):
                            # run_command now uses the trimmed alignment
                            output, exec_time = run_command(program, trim_out, i)
                            f_trees.write(f"{directory}\t{program}\t{i}\t{exec_time}\n")
                            f_trees.flush() # Force write loop output immediately

                except subprocess.CalledProcessError as e:
                    print(f"Trimming failed for {trim_in}, skipping trees: {e}")

if __name__ == "__main__":
    main()
