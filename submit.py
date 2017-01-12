import os, shutil
import subprocess

SLURMJobTemplate = """#!/bin/bash -l 
##execution shell environment 

## name of your job
#SBATCH -J %job
## system error message output file
#SBATCH -e %outDir/stderr.txt
## system message output file
#SBATCH -o %outDir/stdout.txt
## a per-process (soft) memory limit
## limit is specified in MB
## example: 1 GB is 1000
#SBATCH --mem-per-cpu=%memory
## how long a job takes, wallclock time hh:mm:ss
#SBATCH -t %wallTime
## number of processes
#SBATCH -n %cores

mkdir -p %outDir

%command"""

def submit(command, outDir, job, memory=4000, cores=1, wallTime="48:00:00", dummy=False, clear=False):
    if not dummy:
        if os.path.exists(outDir):
            if clear:
                print "Removing output directory", outDir
                shutil.rmtree(outDir)
            else:
                print "Output directory", outDir, "already exists"
                raise Exception()
        print "Making output directory", outDir
        os.makedirs(outDir)
    
    command = command.replace("%outDir")
    template = SLURMJobTemplate
    for param, value in [("%command", command), ("%outDir", outDir), ("%job", job), ("%memory", memory), ("%cores", cores), ("%wallTime", wallTime)]:
        template = template.replace(param, value)
    print "==========", "Template", "=========="
    print template
    print "===================================="
    if not dummy:
        with open(os.path.join(outDir, "template.txt"), "wt") as f:
            f.write(template)
        print "Submitting job", job
        p = subprocess.Popen("sbatch", stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.communicate(input=template)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option('-c','--command', help='')
    optparser.add_option("-o", "--outDir", default=None, help="")
    optparser.add_option('-j','--job', default=None, help='')
    optparser.add_option('-m','--memory', default=4000, type=int, help='')
    optparser.add_option('-r','--cores', default=1, type=int, help='')
    optparser.add_option('-t','--time', default="48:00:00", help='')
    optparser.add_option("--dummy", default=False, action="store_true", help="")
    optparser.add_option("--clear", default=False, action="store_true", help="")
    (options, args) = optparser.parse_args()
    
    submit(options.command, options.outDir, options.job, options.memory, options.cores, options.time, options.dummy, options.clear)