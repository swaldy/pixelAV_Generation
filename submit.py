#!/usr/bin/env python3
#

# Submission script for parallel processing of pixel_clusters .out.gz files on lxplus

import os, sys
from datetime import datetime

def main():

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = f"logs/run_{timestamp}"
    os.makedirs(logdir, exist_ok=True)
    print(f"Logs will be written to: {logdir}")

    #for i in range(16401,16440):
    for i in [16401]:
        subfile = f"{logdir}/d{i}.sub"
        f = open(subfile,"w")

        f.write("universe                = vanilla \n")
        f.write("executable              = datagen.sh \n")
        f.write("arguments               = "+str(i)+" \n")
        f.write("request_memory          = 128 GB \n")
        f.write("transfer_input_files    = datagen_10ps_48x192.py, datagen.sh \n")
        f.write("transfer_output_files   = \"\" \n")
        f.write(f"output                  = {logdir}/d{i}.out \n")
        f.write(f"error                   = {logdir}/d{i}.err \n")
        f.write(f"log                     = {logdir}/d{i}.log \n")

        # Job flavour determines job wall time
        # https://batchdocs.web.cern.ch/local/submit.html#job-flavours
        f.write("+JobFlavour             = \"longlunch\" \n")
        f.write("queue \n")

        f.close()

        # submit the job
        os.system("condor_submit " + subfile)

    return

if __name__ == "__main__":
    main()
