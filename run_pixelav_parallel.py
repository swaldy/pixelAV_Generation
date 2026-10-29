import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
 
# ============================================================
# Configuration
# ============================================================
 
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
 
DATA_DIR = "./pixelav_48x192_10ps/gzip_files_48_192_center"
#DATA_DIR = "/eos/user/s/swaldych/smart_pix/pixelav_48x192_10ps/gzip_files_48_192_center"
 
BINARY = os.path.join(
    BASE_DIR,
    "ppixelav2_list_trkpy_n_2f_10ps_48x192"
)
 
# 200,000 tracks / 5,000 tracks per block = 40 blocks
NUM_BLOCKS = 40
RUNSIZE = 5000
 
# filebase from ppixel2.init
FILEBASE = 16400
 
# Number of PixelAV jobs running simultaneously.
# Start conservatively.
MAX_PARALLEL = 5
 
 
def run_block(frun):
    """
    Run one PixelAV block.
 
    frun=1  -> tracks 1-5000
    frun=2  -> tracks 5001-10000
    ...
    frun=40 -> tracks 195001-200000
    """
 
    file_index = FILEBASE + frun
 
    output_file = os.path.join(
        DATA_DIR,
        f"pixel_clusters_d{file_index:05d}.out"
    )
 
    gzip_file = output_file + ".gz"
 
    log_file = os.path.join(
        DATA_DIR,
        f"pixelav_block_{frun:02d}.log"
    )
 
    # Don't accidentally regenerate a completed block
    if os.path.exists(gzip_file):
        print(f"[Block {frun:02d}] Already finished. Skipping.")
        return frun, True
 
    command = [
        BINARY,
        str(frun),
        str(RUNSIZE)
    ]
 
    print(
        f"[Block {frun:02d}] Starting "
        f"tracks {(frun - 1) * RUNSIZE + 1}"
        f"-{frun * RUNSIZE}"
    )
 
    with open(log_file, "w") as log:
        result = subprocess.run(
            command,
            cwd=DATA_DIR,
            stdout=log,
            stderr=subprocess.STDOUT
        )
 
    if result.returncode != 0:
        print(
            f"[Block {frun:02d}] FAILED "
            f"(return code {result.returncode})"
        )
        return frun, False
 
    if not os.path.exists(output_file):
        print(
            f"[Block {frun:02d}] FAILED: "
            f"{output_file} was not created."
        )
        return frun, False
 
    # PixelAV does not stop when an input table is missing (e.g. it prints
    # "cant find SIRUTH.SPR" and keeps going with a placeholder), so check the log.
    with open(log_file, errors="replace") as log:
        if "cant find" in log.read():
            print(
                f"[Block {frun:02d}] FAILED: PixelAV could not find an input file "
                f"(see {log_file}); output left uncompressed, do not use it."
            )
            return frun, False
 
    # gzip only after PixelAV successfully finishes
    gzip_result = subprocess.run(
        ["gzip", "-f", output_file]
    )
 
    if gzip_result.returncode != 0:
        print(f"[Block {frun:02d}] gzip FAILED.")
        return frun, False
 
    print(f"[Block {frun:02d}] Finished successfully.")
 
    return frun, True
 
 
def main():
 
    # --------------------------------------------------------
    # Sanity checks before launching anything
    # --------------------------------------------------------
 
    if not os.path.isfile(BINARY):
        raise FileNotFoundError(
            f"PixelAV binary not found:\n{BINARY}"
        )
 
    tracklist = os.path.join(DATA_DIR, "track_list_200k_cotBeta_lt2.txt")
    ppixel_init = os.path.join(DATA_DIR, "ppixel2.init")
    wgt_init = os.path.join(DATA_DIR, "wgt_pot.init")
    # Energy-loss cross sections read by dsdlne(); without it PixelAV runs anyway
    # but deposits almost no charge. It lives in the repo root: copy it to DATA_DIR.
    siruth = os.path.join(DATA_DIR, "SIRUTH.SPR")
 
    for required_file in [tracklist, ppixel_init, wgt_init, siruth]:
        if not os.path.exists(required_file):
            raise FileNotFoundError(
                f"Required input file not found:\n{required_file}"
            )
 
    # Leftovers from an earlier run make PixelAV resume instead of starting fresh:
    # a seedfile stores the last event done (a finished block's says "all done", so
    # nothing new is simulated) and the .out file is opened in append mode.
    stale = []
    for frun in range(1, NUM_BLOCKS + 1):
        file_index = FILEBASE + frun
        out = os.path.join(DATA_DIR, f"pixel_clusters_d{file_index:05d}.out")
        if os.path.exists(out + ".gz"):
            continue                      # finished block: run_block skips it
        for leftover in (out, os.path.join(DATA_DIR, f"seedfile{file_index:05d}")):
            if os.path.exists(leftover):
                stale.append(leftover)
    if stale:
        raise RuntimeError(
            "Leftover files from a previous run would make PixelAV resume or append "
            "instead of starting fresh. Move them out of the way first:\n  "
            + "\n  ".join(stale)
        )
 
    print("=" * 60)
    print("PixelAV 48x192 / 10 ps production")
    print("=" * 60)
    print(f"Binary:             {BINARY}")
    print(f"Dataset directory:  {DATA_DIR}")
    print(f"Number of blocks:   {NUM_BLOCKS}")
    print(f"Tracks per block:   {RUNSIZE}")
    print(f"Total tracks:       {NUM_BLOCKS * RUNSIZE}")
    print(f"Parallel processes: {MAX_PARALLEL}")
    print("=" * 60)
 
    failed_blocks = []
 
    # --------------------------------------------------------
    # Launch jobs
    # --------------------------------------------------------
 
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as executor:
 
        futures = {
            executor.submit(run_block, frun): frun
            for frun in range(1, NUM_BLOCKS + 1)
        }
 
        for future in as_completed(futures):
 
            frun = futures[future]
 
            try:
                block, success = future.result()
 
                if not success:
                    failed_blocks.append(block)
 
            except Exception as exc:
                print(
                    f"[Block {frun:02d}] Exception: {exc}"
                )
                failed_blocks.append(frun)
 
    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------
 
    print()
    print("=" * 60)
    print("PixelAV production finished")
    print("=" * 60)
 
    if failed_blocks:
        print("Failed blocks:")
        print(sorted(failed_blocks))
    else:
        print("All 40 blocks completed successfully.")
        print("Total tracks requested: 200,000")
 
    print("=" * 60)
 
 
if __name__ == "__main__":
    main()
