import sys, os
import glob
import gzip
import csv
import Stream

def readCSVs(pathPattern, columns, sortBy=None, delimiter="\t"):
    rows = []
    rowById = {}
    idName = getIdName(columns)
    for pathName in glob.glob(pathPattern):
        print "Loading input file", pathName
        #fileName = os.path.basename(pathName)
        with (gzip.open(pathName, "rt") if pathName.endswith(".gz") else open(pathName, "rt")) as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            fileRows = []
            for line in reader:
                rowId = line[idName]
                if rowId not in rowById:
                    rowById[rowId] = {key:[] for key in line.keys()}
                    if sortBy != None:
                        rowById[rowId]["rank"] = []
                    rows.append(rowById[rowId])
                row = rowById[rowId]
                for key in line:
                    if key != idName:
                        line[key] = columns[key]["type"](line[key])
                        row[key].append(line[key])
                    else:
                        row[key] = line[key]
                rowById[rowId] = row
                fileRows.append(line)
            if sortBy != None:
                print "Sorting by column", sortBy
                fileRows.sort(key=lambda x: x[sortBy], reverse=True)
                rank = 0
                for fileRow in fileRows:
                    fileRow["rank"] = rank
                    rowById[fileRow[idName]]["rank"].append(rank)
                    rank += 1
    #if sortBy != None:
    #    rows.sort(key=lambda x: min(x["rank"]))
    return rows

def mergeRows(rows, columns):
    print "Merging the rows"
    keys = sorted(rows[0].keys())
    for row in rows:
        for key in keys:
            if key not in columns or columns[key]["action"] == "mean":
                row[key] = sum(x for x in row[key]) / len(row[key])
            elif columns[key]["action"] == "max":
                row[key] = max(row[key])
            elif columns[key]["action"] == "min":
                row[key] = min(row[key])
            elif columns[key]["action"] == "id":
                pass
            elif columns[key]["action"] == "unique":
                values = set(row[key])
                assert len(values) == 1, row
                row[key] = list(values)[0]
            elif columns[key]["action"] == "range":
                minVal = min(row[key])
                maxVal = max(row[key])
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

def getIdName(columns):
    print "Validating the column rules:", columns
    idName = None
    for key in columns:
        if columns[key]["action"] == "id":
            assert idName == None
            idName = key
    assert idName != None
    return idName

def mergeCSVs(pathPattern, outPath, columns, sortBy=None, delimiter="\t"):
    Stream.openLog(outPath + "-log.txt")
    columns = getColumns(columns)
    rows = readCSVs(pathPattern, columns, sortBy=sortBy, delimiter=delimiter)
    if len(rows) == 0:
        print "No input rows"
        return
    rows = mergeRows(rows, columns)
    if sortBy != None:
        print "Sorting results by column", sortBy
        rows.sort(key=lambda x: x["importance"], reverse=True)
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
    optparser.add_option("-s", "--sortBy", default=None, help="Column to sort by (generates the 'rank' column)")
    optparser.add_option("-d", "--delimiter", default="\t", help="CSV file delimiter")
    (options, args) = optparser.parse_args()
    
    mergeCSVs(options.input, options.output, options.columns, options.sortBy, options.delimiter)