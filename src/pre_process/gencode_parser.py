import sys
import gzip


def readGTF(gtfFile):
    l = []
    r = open(gtfFile, "r")
    for line in r:
        if line[0] != "#":
            line = line.strip().split("\t")
            if line[2] == "gene":
                chrom = line[0]
                atts = line[8].split(";")
                del atts[-1]
                gene_id = atts[0].split(" ")[1].replace('"','').split(".")[0]
                gene_type = atts[1].split(" ")[2].replace('"','')
                #gene_status = atts[3].split(" ")[2].replace('"','')
                gene_status = "NULL"
                gene_name = atts[2].split(" ")[2].replace('"','')
                level = atts[3].split(" ")[2]
                l.append((chrom, gene_id, gene_name, gene_type, gene_status, level))
    r.close()
    return l


def writePositions(refGTF, outTXT):
    GTFlist = readGTF(refGTF)
    w = open(outTXT, "w")
    w.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("chrom", "gene_id", "gene_name", "gene_type", "gene_status", "level"))
    for gene in GTFlist:
        w.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene[0], gene[1], gene[2], gene[3], gene[4], gene[5]))
    w.close()


input_gtf_file = sys.argv[1]
output_txt_file = sys.argv[2]

writePositions(input_gtf_file, output_txt_file)
