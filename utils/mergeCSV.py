import sys, os
import glob
import gzip
import csv
import Stream

def readCSVs(pathPattern, idName, delimiter="\t"):
    rows = []
    rowById = {}
    for pathName in glob.glob(pathPattern):
        print "Loading input file", pathName
        #fileName = os.path.basename(pathName)
        with (gzip.open(pathName, "rt") if pathName.endswith(".gz") else open(pathName, "rt")) as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for line in reader:
                rowId = line[idName]
                if rowId not in rowById:
                    rowById[rowId] = {key:[] for key in line.keys()}
                    rows.append(rowById[rowId])
                row = rowById[rowId]
                for key in line:
                    if key != idName:
                        row[key].append(line[key])
                    else:
                        row[key] = idName
                rowById[rowId] = row
    return rows

def mergeRows(rows, columns):
    print "Merging the rows"
    keys = sorted(rows[0].keys())
    for row in rows:
        for key in keys:
            if key not in columns or columns[key]["action"] == "mean":
                row[key] = sum([columns[key]["type"](x) for x in row[key]]) / len(row[key])
            elif columns[key]["action"] == "id":
                pass
            elif columns[key]["action"] == "range":
                values = [columns[key]["type"](x) for x in row[key]]
                minVal = min(values)
                maxVal = max(values)
                if (minVal != maxVal):
                    row[key] = str(minVal) + "-" + str(maxVal)
                else:
                    row[key] = minVal
            else:
                raise Exception("Unknown column action '" + str(columns[key]) + "'")
    return rows

def noType(value):
    return value

def getColumns(columnPattern):
    columns = {}
    splits = [x.split(":") for x in columnPattern.split(",")]
    for split in splits:
        if len(split) == 2:
            columns[split[0]] = {"action":split[1], "type":noType}
        else:
            columns[split[0]] = {"action":split[1], "type":eval(split[2])}
    return columns

def validateColumns(columns):
    print "Validating the column rules:", columns
    idName = None
    for key in columns:
        if columns[key]["action"] == "id":
            assert idName == None
            idName = key
    assert idName != None
    return idName

def mergeCSVs(pathPattern, outPath, columns, delimiter="\t"):
    Stream.openLog(outPath + "-log.txt")
    columns = getColumns(columns)
    idName = validateColumns(columns)
    rows = readCSVs(pathPattern, idName, delimiter=delimiter)
    if len(rows) == 0:
        print "No input rows"
        return
    rows = mergeRows(rows, columns)
    print "Saving results to", outPath
    with (gzip.open(outPath, "wt") if outPath.endswith(".gz") else open(outPath, "wt")) as f:
        writer = csv.DictWriter(f, fieldnames=sorted(columns.keys()), delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-i", "--input", help="Glob pattern. Remember to quote so the shell doesn't expand it!")
    optparser.add_option("-o", "--output", help="Output file path")
    optparser.add_option("-c", "--columns", help="Comma-separated list of column:action[:type] duplets/triplets, e.g. 'mycolumn:mean:float'.")
    optparser.add_option("-d", "--delimiter", default="\t", help="CSV file delimiter")
    (options, args) = optparser.parse_args()
    
    mergeCSVs(options.input, options.output, options.columns, options.delimiter)