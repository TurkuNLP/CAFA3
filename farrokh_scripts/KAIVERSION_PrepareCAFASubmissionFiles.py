# example input parameters:
#
#-i /home/farmeh/Desktop/PROJECTS/GIT/CAFA3/data/testpred/LASTKAI_cafa_targets.tsv.gz -o /home/farmeh/Desktop/PROJECTS/GIT/CAFA3/data/testpred/ExampleSubmissionFolder/ -m 2 -l /home/farmeh/Desktop/PROJECTS/GIT/CAFA3/data/testpred/ExampleSubmissionFolder/

TARGET_ORGANISMS = ["moon"] + [str(i) for i in [10090 , 10116 , 160488 , 170187 , 208963 , 223283 , 224308 , 237561 , 243232 , 243273, 273057 ,284812, 321314, 3702 , 44689 , 559292 , 7227 , 7955 , 83333 , 8355 , 85962 , 9606 , 99287]] 
PARAM_KEYWORDS   = "sequence alignment, sequence-profile alignment, profile-profile alignment, predicted properties, protein structure, machine learning, other functional information" + "." 

import shutil , datetime , argparse , gzip , csv , zipfile
import MySQLdb as mdb 

ErrLogFileHandler = None 

def DATETIME_GetNowStr():
    return datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S") ;  

def Errlog(text):
    global LogFileHandler 
    text = "[" + DATETIME_GetNowStr() + "]\t" + text 
    print text 
    ErrLogFileHandler.write (text + "\n") 

def func_is_number(s):
    try:
        if s == None: 
            return False ;
        float(s)
        return True
    except ValueError:
        return False

def func_os_MakeDir (FolderAddress):
    if not shutil.os.path.exists(FolderAddress):
        shutil.os.makedirs(FolderAddress)

def func_os_RemoveDirectoryWithContent(FolderAddress):
    if shutil.os.path.exists(FolderAddress):
        shutil.rmtree(FolderAddress)

def func_os_IsDirectory (FolderAddress):
    return shutil.os.path.isdir (FolderAddress)

def func_os_FileExists (FileAddress):
    return shutil.os.path.exists (FileAddress)
    
def func_db_connectToEVEXDB (_DB_ADDR, _UserName, _Password, _DBName):
    try:
        EVEXDBcon = mdb.connect(_DB_ADDR,_UserName, _Password, _DBName);
        return EVEXDBcon
    except mdb.Error, e:
        Errlog ("Error :" + str(e.args[0]) + " , " + str (e.args[1])) 
        shutil.sys.exit(1)

def func_EXIT (ExitMSG , code=-1):
    print ExitMSG ; 
    shutil.sys.exit (code); 

def func_db_GetUniprot_Info (UniprotID , EVEXDBcon):
    SQLSTR = 'select CAFA_id , ncbitax_id from map_cafa2protein_organism where symbol = "' + UniprotID + '";' 
    CUR = EVEXDBcon.cursor()
    CUR.execute(SQLSTR)
    DATA = CUR.fetchall()
    CUR.close ()
    return DATA 

def func_openInputFile (InputFileAddress):
    try:    
        FileExtention = InputFileAddress.rsplit(".", 1)[-1].lower()
        if FileExtention in ("tsv" , "txt"):
            return open(InputFileAddress , "rt")
        elif FileExtention == "gz":
            return gzip.open (InputFileAddress , "rt")
        else:
            raise Exception ("Invalid file format:" + FileExtention) 
    except Exception as E:
        Errlog ("Error in opening input file.\nError:" + E.message)
        Errlog ("Abnormal program exit!!!")
        shutil.sys.exit(-1)

def func_FloatNumberWithNDecimalPoints (NUM, N):
    if N <= 0: return NUM ; #Safe    
    from decimal import Decimal ; 
    try:
        L = float(Decimal(NUM).quantize (Decimal(10) ** -N)) ;
        return L;
    except:
        L = str(NUM);
        return [L[0:L.find(".")+N+1]];
        
if __name__== "__main__": 
    CurrentExecutionFolder = shutil.os.path.dirname(shutil.os.path.realpath(__file__)) 
    if CurrentExecutionFolder[-1] <> "/":
        CurrentExecutionFolder+= "/" 

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile"  , help= "absolute address for the input prediction result file.") 
    parser.add_argument("-o", "--outputDir"  , help= "output folder absolute address") 
    parser.add_argument("-m", "--modelNumber", help= "Should be 1, 2, or 3")
    parser.add_argument("-l", "--logDir"     , help= "logfile folder absolute address"   , default=CurrentExecutionFolder)
    parser.add_argument("-t", "--teamName"   , help= "name of the team"                  , default="TurkuBioNLP1") 
    args = parser.parse_args()
    
    if (args.inputFile == None) or (not(func_os_FileExists (args.inputFile))):
        func_EXIT ("Invalid input file:" + str(args.inputFile)) 
    
    if (args.outputDir == None) or (not(func_os_IsDirectory (args.outputDir))):
        func_EXIT ("Invalid output folder argument:" + str(args.outputDir)) 

    if (args.logDir == None) or (not(func_os_IsDirectory (args.logDir))):
        func_EXIT ("Invalid log folder argument:" + str(args.logDir)) 

    if (not func_is_number(args.modelNumber)) or (not int(args.modelNumber) in (1,2,3)):
        func_EXIT ("Model number argument should be either 1, 2 or 3. Now it is:" + str(args.modelNumber)) 
    
    if args.outputDir[-1] <> "/": args.outputDir+="/" ;
    if args.logDir[-1]    <> "/": args.logDir   +="/" ;

    print "\n\n\n" 
    print "1) Connecting to EVEXDB ..." 
    EVEXDBcon = func_db_connectToEVEXDB  ('evex.utu.fi','evex_write' , 'EVEXopensesame156', 'CAFA3');
    print "Done.\n" 
    
    ErrLogFileAddress = args.logDir + "MODEL" + args.modelNumber + "_" + DATETIME_GetNowStr() + ".log" 
    print "2) Creating error log file: " , ErrLogFileAddress
    ErrLogFileHandler = open (ErrLogFileAddress , "wt") 
    print "Done.\n" 
    
    print "3) Opening Input file:" , args.inputFile
    INPUT_FileHandler = func_openInputFile (args.inputFile) 
    print "Done.\n" 
    
    OUTPUT_FOLDER = args.outputDir + "MODEL" + args.modelNumber + "_" + DATETIME_GetNowStr() + "/" 
    print "4) Deleting and Creating Output Directory:\n" , OUTPUT_FOLDER
    func_os_RemoveDirectoryWithContent(OUTPUT_FOLDER) 
    func_os_MakeDir(OUTPUT_FOLDER)
    print "Done.\n" 
        
    print "5) Creating output files ..." 
    ALL_FILES = {}
    for ORG_ID in TARGET_ORGANISMS:
        filename = args.teamName + "_" + args.modelNumber + "_" + ORG_ID + ".txt" 
        ALL_FILES[ORG_ID] = open(OUTPUT_FOLDER + filename , "wt") 
        HEADER = "AUTHOR "   + args.teamName    + "\n" + \
                 "MODEL "    + args.modelNumber + "\n" + \
                 "KEYWORDS " + PARAM_KEYWORDS   + "\n" 
        ALL_FILES[ORG_ID].write (HEADER)
    print "Done.\n" 
    
    print "6) Writing prediction outputs ..." ; 
        
    cnt = 0 
    csv_reader = csv.DictReader(INPUT_FileHandler, delimiter='\t')
    Last = {"UniProtID": None , "Info": None }
    
    for row in csv_reader:
        cnt+= 1 
        if cnt%10000 == 0:
            print "row processed:" , cnt 

        UniProtID         = row["id"]
        PredicatedGoTerm  = row["label"]
        ConfidenceScore   = float(row["confidence"])
        IsPositivePred    = int  (row["predicted"]) == 1
        Internal_cafa_ids = None if len(row["cafa_ids"])==0 else set(row["cafa_ids"].split(","))
        
        if not IsPositivePred:
            continue ; 
        
        if (ConfidenceScore <= 0) or (ConfidenceScore > 1):
            Errlog ("Invalid Confidence Score. Should be in the interval (0.00, 1.00]. Info: " + str(row))
            continue

        if UniProtID == Last["UniProtID"]:
            Info = Last["Info"]
        else:
            Info = func_db_GetUniprot_Info (UniProtID , EVEXDBcon)
            Last["UniProtID"] = UniProtID
            Last["Info"] = Info 
            
        if len (Info) == 0: #information is not found in DB. 
            if Internal_cafa_ids == None: #No info in DB, no info in file. That's not CAFA Target, skip
                continue 
            else:#No info in DB, but CAFA id in File 
                Errlog ("Unknown CAFA TargetID (Defined in file, not in MySQL): " + row["cafa_ids"]) 
                continue 
        else: #info is in the DB, let's check if CAFA id in DB and file agree together... 
            DB_cafa_ids = set(i[0] for i in Info)
            
            for c_id in Internal_cafa_ids:
                if not c_id in DB_cafa_ids:
                    Errlog ("cafa_id in FILE, but not in DB :" + str(Internal_cafa_ids - DB_cafa_ids) + str(row))
                    continue 
            """
            if len(DB_cafa_ids - Internal_cafa_ids)>0:
                Errlog ("cafa_id in DB, but not in FILE :" + str(Internal_cafa_ids - DB_cafa_ids) + str(row))
                continue 
            """
                    
        
        for cafa_id in Internal_cafa_ids:
            if cafa_id.startswith("M"):
                ORG_TYPE = "moon" 
            else:
                ORG_TYPE = str(Info[0][1])# ncbitax_id) #change to STR ... by default it is Long int

            if not ORG_TYPE in ALL_FILES:
                Errlog ("Prediction not among CAFA targets. Info: " + ORG_TYPE + "     " + str(row))
                continue
                
            ALL_FILES[ORG_TYPE].write (cafa_id + "\t" + PredicatedGoTerm + "\t" + '%.2f' % ConfidenceScore + "\n") 
            #print CAFA_id , ORG_TYPE , PredicatedGoTerm , ConfidenceScore 
    
    print "row processed:" , cnt      
    print "-------------------------------------------------" 
    print "CHECK LOG FILE TO SEE IF THERE IS ANY ERRORS!" ; 
    print "log:" + ErrLogFileAddress
    print "-------------------------------------------------" 

    #CLOSE INPUT/OUTPUT FILES ...     
    INPUT_FileHandler.close ()
    for FileHandler in ALL_FILES:
        ALL_FILES[FileHandler].write ("END")
        ALL_FILES[FileHandler].close() 

    """
    Zip each file into a separate zip file ...
    # find . -type f -execdir zip '{}.zip' '{}' \; -exec rm '{}' \;
    shutil.os.chdir (OUTPUT_FOLDER) 
    shutil.os.system ("find . -type f -execdir zip '{}.zip' '{}' \; -exec rm '{}' \;") 
    
    print "COMPRESSING FILES FOR CAFA SUBMISSION ..." 
    ProteinCentricPredFiles = [] 
    TermCentricAndMoonlightPredFiles = []
    for f in ALL_FILES:
        name = ALL_FILES[f].name 
        if not "moon" in name:
            ProteinCentricPredFiles.append (name)
        else:
            TermCentricAndMoonlightPredFiles.append (name)

    #ProteinCentricPredictions
    ZipFileName = args.teamName + "_" + args.modelNumber + "_ProteinCentric" 
    Z1 = zipfile.ZipFile(OUTPUT_FOLDER + ZipFileName + ".zip" , "w" )
    for f in ProteinCentricPredFiles:
        Z1.write(f, ZipFileName + "/" + f.split("/")[-1], compress_type=zipfile.ZIP_DEFLATED)
    Z1.close()
    
    #TermCentricandMoonlightingPredictions
    ZipFileName = args.teamName + "_" + args.modelNumber + "_TermCentricAndMoonlighting"
    Z2 = zipfile.ZipFile(OUTPUT_FOLDER + ZipFileName + ".zip" , "w")
    for f in TermCentricAndMoonlightPredFiles:
        Z2.write(f, ZipFileName + "/" + f.split("/")[-1], compress_type=zipfile.ZIP_DEFLATED)
    Z2.close()

    #delte txt files ... 
    shutil.os.chdir (OUTPUT_FOLDER) 
    shutil.os.system ("rm -rf *.txt") 
    """
    
    print "-------------------------------------------------" 
    print "EXITING PROGRAM..." 
    EVEXDBcon.close ()
    ErrLogFileHandler.close ()
    print "END." 


    """
    VALID keywords: 
    --------------------------------------
    [x] sequence alignment
    [x] sequence-profile alignment
    [x] profile-profile alignment
    [ ] phylogeny
    [?] sequence properties
    [ ] physicochemical properties
    [x] predicted properties
    [ ] protein interactions
    [ ] gene expression
    [ ] mass spectrometry
    [ ] genetic interactions
    [x] protein structure
    [ ] literature
    [?] genomic context
    [ ] synteny
    [ ] structure alignment
    [ ] comparative model
    [?] predicted protein structure
    [ ] de novo prediction
    [x] machine learning
    [ ] genome environment
    [ ] operon
    [ ] ortholog
    [ ] paralog
    [ ] homolog
    [ ] hidden Markov model
    [ ] clinical data
    [ ] genetic data
    [ ] natural language processing
    [x] other functional information
    """