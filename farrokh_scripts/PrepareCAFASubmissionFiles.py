# example input parameters:
#
# -i /home/farmeh/Desktop/MY_ALL_PROJECTS/CAFA/CAFA3/PrepareSubmissionFiles/Results_Jari_Example/devel-predictions.tsv -o /home/farmeh/Desktop/MY_ALL_PROJECTS/CAFA/CAFA3/PrepareSubmissionFiles/ExampleSubmissionFolder/ -m 1 -l /home/farmeh/Desktop/MY_ALL_PROJECTS/CAFA/CAFA3/PrepareSubmissionFiles/ExampleSubmissionFolder/ -t TurkuBioNLP


TARGET_ORGANISMS = ["moon"] + [str(i) for i in [10090 , 10116 , 160488 , 170187 , 208963 , 223283 , 224308 , 237561 , 243232 , 243273, 273057 ,284812, 321314, 3702 , 44689 , 559292 , 7227 , 7955 , 83333 , 8355 , 85962 , 9606 , 99287]] 
PARAM_KEYWORDS   = "sequence alignment, sequence-profile alignment, profile-profile alignment" + "." 

import shutil , datetime , argparse , gzip
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
        FileExtention = InputFileAddress.split("/")[-1].split(".")[1].lower()
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
    parser.add_argument("-t", "--teamName"   , help= "name of the team"                  , default="EVEX") 
    parser.add_argument("-s", "--skipFirstLine", help = "Should I Skip first line of input file." , default="False")
    args = parser.parse_args()

    
    if (args.inputFile == None) or (not(func_os_FileExists (args.inputFile))):
        func_EXIT ("Invalid input file:" + str(args.inputFile)) 
    
    if (args.outputDir == None) or (not(func_os_IsDirectory (args.outputDir))):
        func_EXIT ("Invalid output folder argument:" + str(args.outputDir)) 

    if (args.logDir == None) or (not(func_os_IsDirectory (args.logDir))):
        func_EXIT ("Invalid log folder argument:" + str(args.logDir)) 

    if (not func_is_number(args.modelNumber)) or (not int(args.modelNumber) in (1,2,3)):
        func_EXIT ("Model number argument should be either 1, 2 or 3. Now it is:" + str(args.modelNumber)) 
    
    #Jari writes column names into his files. That must be skipped ... 
    if not str(args.skipFirstLine).lower() in ("true" , "false"):
        func_EXIT ("skipFirstLine argument should be either true or false.")
    if args.skipFirstLine.lower() == "true":
        args.skipFirstLine = True
    else:
        args.skipFirstLine = False 
    
    if args.outputDir[-1] <> "/": args.outputDir+="/" ;
    if args.logDir[-1]    <> "/": args.logDir   +="/" ;

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
    if args.skipFirstLine:
        INPUT_FileHandler.readline() 
        
    cnt = 0 
    for line in INPUT_FileHandler:
        cnt+= 1 
        if cnt == 4000:
            break ;

        if args.skipFirstLine and (cnt==1):
            continue ; 
        
        origline = line 
        line = line.split("\n")[0].split("\t") 
        UniProtID , PredicatedGoTerm , ConfidenceScore = line[0] , line[1] , line[2]
        
        Info = func_db_GetUniprot_Info (UniProtID , EVEXDBcon)
        if len (Info) == 0:
            Errlog ("Unknown UniProtID :" + UniProtID)
            continue 

        #The score must be in the interval (0.00, 1.00] and contain two significant figures. A score of 0.00 is not allowed        
        if (float(ConfidenceScore) <= 0) or (float(ConfidenceScore) > 1):
            print ConfidenceScore 
            Errlog ("Invalid Confidence Score. Should be in the interval (0.00, 1.00]. Info: " + origline.split("\n")[0])
            continue
        
        print UniProtID , PredicatedGoTerm , ConfidenceScore  
        for record in Info:
            CAFA_id , ncbitax_id = record[0] , record[1] 
            if CAFA_id[0] == "M":
                ORG_TYPE = "moon" 
            else:
                ORG_TYPE = ncbitax_id 
            

            if not ORG_TYPE in ALL_FILES:
                Errlog ("Prediction not among CAFA targets. Info: " + origline.split("\n")[0])
                continue
                
            ALL_FILES[ORG_TYPE].write (CAFA_id + "\t" + PredicatedGoTerm + "\t" + '%.2f' % ConfidenceScore + "\n") 
            
    print "EXITING PROGRAM..." 
    INPUT_FileHandler.close ()
    for FileHandler in ALL_FILES:
        ALL_FILES[FileHandler].write ("END")
        ALL_FILES[FileHandler].close() 
    EVEXDBcon.close ()
    ErrLogFileHandler.close ()
    print "END." 
