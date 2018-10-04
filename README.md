# CAFA3
University of Turku CAFA3 project

Files are in the new machine in address: /home/sukaew/CAFA3

CNN experiment can be run with python train.py
You'll need to copy the data folder from /home/kahaka/CAFA3/

Running preprocessing and sequence analysis
-------------------------------------------
All preprocessing steps and sequence analyses can be run within the directory 'sequence_features' using the following command line. 

`python3 target_process.py -o [out_folder] -s [ori_seq]`

The program needs two inputs, the out_folder and the ori_seq. The out_folder should be absolute path ended with '/' for the input ori_seq file/directory and the output features directory. The input ori_seq can be in these four formats (folder of non-compressed fasta files, tar.gz, gz and zip). The sequence analysese include Blast Protein, DeltaBlast, Interproscan5, NetAcet, predGPI, nucPred and Taxonomy hierarchy. All analysis results are in folder called `feature`

Running the Feature Generation, Classification and Analysis
-----------------------------------------------------------
All experiments can be run using the program `run.py`. The experimental code uses a three-step system. One or more of these actions can be performed using the command line option `--action` or `--a`. By default, all three actions (`build`, `classify` and `statistics`) are performed.

The run.py program can be called like this:

`python run.py -e [TASK] -o [OUTPUT] --targets external`

The `[TASK]` value can be one of `cafa3`, `cafa3hpo` or `cafapi`. Depending on task, different input files are used. The `--targets` option defines how to handle CAFA targets.

Cross-validation
----------------
By default, the scikit-learn classification will use the train/devel/test split for the learning data. To use n-fold cross-validation instead, use the `--fold` option of `run.py`. To do 10-fold cross-validation, the program can be run 10 times using a script like this:

`for FOLD in 0 1 2 3 4 5 6 7 8 9; do python run.py -o /tmp/CAFA10fold/fold$FOLD --fold $FOLD; done`

Ensemble
--------
The program `ensemble.py` can be used to combine predictions from different systems and the BLAST fallback baseline. To run the ensemble, use a command like:

`python ensemble.py -a [PRED1_DIR]  -b [PRED2_DIR] -o [OUTPUT] --baseline 4 --simple --terms 1000000 --write --cafa --clear`
