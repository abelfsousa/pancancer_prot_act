#!/usr/bin/python

# import modules
import argparse
import sys

####################################################################
### notes: in this library, bed files for operations like merge  ###
### and substracting contains require at least four columns:     ###
### chr,start,end,strand.                                        ###
####################################################################

################################
######## function list  ########
################################
#1.	neighbor_merge: merge two ranges, say two exons
#2.	bedmerge: merge a list of ranges, say a few exons 
#3.	exon2intron: calculate intron bed from a list of exons.
#4.	masked_intron: identify introns which are covered by exons of other isoforms
#5.	get_UTR: get 5' and 3' UTR of isoforms.
#6.	get_gene_bed: generate bed files for a list of input genes given in a txt file
#7.	merge_exon: merge all exons of multiple isoforms from the same gene.
#8. 	get_intron: calculate intron bed for each isoform
#
#
#
#################################
############ functions:  ########
#################################


def neighbor_merge(range1,range2):
	if range2[1]<=range1[2]:
		merged =[(range1[0],range1[1],max(range1[2],range2[2]),range1[3])]
	else:
		merged=[range1,range2]
	return merged



def bedmerge(featureRange):
	# featureRange: a list of ranges in bed format, such as [(1,1000,2000,+),(1,2200,3000,-)]
	featureRange.sort(key = lambda x: (x[0],x[1]))	
	merged=[]
	nRange=len(featureRange)
	if nRange == 1:
		merged=featureRange
	elif nRange == 2:
		imerge=neighbor_merge(featureRange[0],featureRange[1])
		for each in imerge:
			merged.append(each)
	else:
		i = 2
		imerge=neighbor_merge(featureRange[0],featureRange[1])
		n_imerge=len(imerge)
		while n_imerge > 0:
			if n_imerge == 2:
				merged.append(imerge[0])
			
			imerge=neighbor_merge(imerge[n_imerge-1],featureRange[i])
			n_imerge=len(imerge)
			if i == nRange-1:
				for each in imerge:
					merged.append(each)
				n_imerge = -1
			i+=1

	return merged	

def exon2intron(featureRange):
	featureRange.sort(key = lambda x: (x[0],x[1]))	
	nRange=len(featureRange)
	intron=[]
	if nRange > 1:
		for i in range(1,nRange):
			exon0=featureRange[i-1]
			exon1=featureRange[i]
			thisintron=(exon0[0],exon0[2],exon1[1],exon1[3])
			intron.append(thisintron)
	return intron

def masked_intron(GTFfile,maskedfile="masked_intron.bed",chroms=map(str,range(1,23))+['X','Y']):
	merged = merge_exon(GTFfile,chroms=chroms)
	introns= get_intron(GTFfile,chroms=chroms)

	exon = merged['merged_exon']
	intron   = introns['intron']
	iso2gene = introns['iso2gene']
	fnew = open(maskedfile,'w')

	for isoform in intron:
		the_introns = intron[isoform]
		the_gene = iso2gene[isoform]
		the_exons = exon[the_gene]
		for i in range(len(the_introns)):
			for j in range(len(the_exons)):
				if the_introns[i][1] >= the_exons[j][1]  and the_introns[i][2] <= the_exons[j][2]:
					output = "\t".join(map(str,the_introns[i]))+'\t'+isoform+'\t'+the_gene+'\n'
					fnew.write(output)
	fnew.close()

def get_UTR(GTFfile, utr_file="",chroms=map(str,range(1,23))+['X','Y']):
	f=open(GTFfile)
	# record UTR, CDS and strand information for determining iUTR5 or iUTR3
	utr = {} 
	cds = {}
	strand = {} 
	for line in f:
		table = line.split('\t')
		if line[0] != '#' and table[0] in chroms:
			if table[2] == 'UTR' and 'transcript_id' in line:
				tcx  = line.split('transcript_id')[1].split('"')[1]
				gene = line.split('gene_id')[1].split('"')[1] 
				infor=[table[0],int(table[3]),int(table[4]),table[6],gene]
				if tcx in utr:
					utr[tcx].append(infor)
				else:
					utr[tcx]=[infor]
				strand[tcx] = table[6]
			if table[2] == 'CDS' and 'transcript_id' in line:
				tcx  = line.split('transcript_id')[1].split('"')[1]
				gene = line.split('gene_id')[1].split('"')[1] 
				cds[tcx] = int((int(table[3])+int(table[4]))/2)
				strand[tcx] = table[6]
		
	f.close()


	# determine 5UTR and 3UTR
	allUTR=[]
	for tcx in cds:
		if tcx in utr:
			istrand = strand[tcx]
			icds = cds[tcx]
			iutr = utr[tcx]
			iUTR5 = []
			iUTR3 = []
			for thisUTR in iutr:
				if istrand == '+':
					if thisUTR[1] > icds:
						iUTR3.append(thisUTR[1])
						iUTR3.append(thisUTR[2])
					else:
						iUTR5.append(thisUTR[1])
						iUTR5.append(thisUTR[2])
				elif istrand == '-':
					if thisUTR[1] > icds:
						iUTR5.append(thisUTR[1])
						iUTR5.append(thisUTR[2])
					else:
						iUTR3.append(thisUTR[1])
						iUTR3.append(thisUTR[2])
			iUTR5.sort()
			iUTR3.sort()
			iiutr=iutr[0]
			if len(iUTR5) > 1:
				allUTR.append([iiutr[0],iUTR5[0]-1,iUTR5[-1],iiutr[3],iiutr[4],'5UTR',tcx])
			if len(iUTR3) > 1:
				allUTR.append([iiutr[0],iUTR3[0]-1,iUTR3[-1],iiutr[3],iiutr[4],'3UTR',tcx])

	allUTR.sort(key = lambda x: (x[0],x[1]))				
	

	# print UTR to file if required
	if len(utr_file) > 1:
		f=open(utr_file,'w')
		for iUTR in allUTR:
			table = iUTR
			iUTR[1] = str(iUTR[1])
			iUTR[2] = str(iUTR[2])
			out = "\t".join(iUTR)+'\n'
			f.write(out)
		f.close()



	return(allUTR)


def get_tss_region(GTFfile,w=2000,tss_bed_file='',chroms=map(str,range(1,23))+['X','Y']):
	# assume genelistfile contains one gene per line.
	# get bed format
	TSSbed=[]
	f = open(GTFfile)
	for line in f:
		table = line.split('\t')
		if table[0] in chroms and table[2] == 'transcript':
			chrom  = table[0]
			strand = table[6] 
			tcx = line.split('transcript_id')[1].split('"')[1]
			geneid = line.split('gene_id')[1].split('"')[1]
			genesymbol = line.split('gene_name')[1].split('"')[1]
			if strand == "+":
				if int(table[3])-1-w >0:
					iregion = [chrom,int(table[3])-1-w,int(table[3])-1+w,strand,tcx,geneid,genesymbol]
			elif strand == '-':
				if int(table[4])-w > 0:
					iregion = [chrom,int(table[4])-w,int(table[4])+w,strand,tcx,geneid,genesymbol]
			TSSbed.append(iregion)
	f.close()
	# write to file
	if len(tss_bed_file) > 1:
		f=open(tss_bed_file,'w')
		for item in TSSbed:
			out = '\t'.join([item[0],str(item[1]),str(item[2]),item[3],item[4],item[5],item[6]])+'\n'
			f.write(out)
		f.close()


def get_gene_bed(GTFfile,gene_bed_file='',chroms=map(str,range(1,23))+['X','Y']):
	# assume genelistfile contains one gene per line.
	# get bed format
	genebed={}
	f = open(GTFfile)
	for line in f:
		table = line.split('\t')
		if table[0] in chroms:
			if table[2] == "gene":
				ensid  = line.split('gene_id')[1].split('"')[1]			
				symbol = line.split('gene_name')[1].split('"')[1].upper()			
				record = (table[0],int(table[3])-1,int(table[4]),table[6],ensid,symbol)
				if table[0] in genebed:
					genebed[table[0]].append(record)
				else:
					genebed[table[0]] = [record]
	f.close()

	# print to file if required 
	if len(gene_bed_file) > 1:
		f=open(gene_bed_file,'w')
		for ichrom in genebed.keys():
			ichrom_gene = genebed[ichrom]
			for igene in ichrom_gene:
				f.write('\t'.join([igene[0],str(igene[1]),str(igene[2]),igene[3],igene[4],igene[5],])+'\n')
		f.close()


	return(genebed)



def get_isoform_bed(GTFfile,isoform_bed_file='',chroms=map(str,range(1,23))+['X','Y']):
	# assume genelistfile contains one gene per line.

	# get bed format
	isoformbed={}
	f = open(GTFfile)
	for line in f:
		table = line.split('\t')
		if table[0] in chroms:
			if table[2] == "transcript":
				isoformid  = line.split('transcript_id')[1].split('"')[1]			
				ensid  = line.split('gene_id')[1].split('"')[1]			
				symbol = line.split('gene_name')[1].split('"')[1].upper()			
				record = (table[0],int(table[3])-1,int(table[4]),table[6],isoformid,ensid,symbol)
				if table[0] in isoformbed:
					isoformbed[table[0]].append(record)
				else:
					isoformbed[table[0]] = [record]
	f.close()

	# print to file if required 
	if len(isoform_bed_file) > 1:
		f=open(isoform_bed_file,'w')
		for ichrom in isoformbed.keys():
			ichrom_isoform = isoformbed[ichrom]
			for iisoform in ichrom_isoform:
				f.write('\t'.join([iisoform[0],str(iisoform[1]),str(iisoform[2]),iisoform[3],iisoform[4],iisoform[5],iisoform[6]])+'\n')
		f.close()
	return(isoformbed)






def merge_exon(GTFfile,merged_exon_file='',chroms=map(str,range(1,23))+['X','Y']):
	# record exon coordination
	exon={}
	f=open(GTFfile)
	for line in f:
		table=line.split('\t')
		if line[0] != '#' and table[0] in chroms:  # skip comment line
			if table[2] == 'exon':
				gene=line.split('gene_id')[1].split('"')[1]
				iexon=(table[0],int(table[3])-1,int(table[4]),table[6])  # gtf to bed coordination
				if gene in exon:
					exon[gene].append(iexon)
				else:
					exon[gene]=[iexon]
	f.close()
	# merge all exons of each gene
	merged_exon={}
	gene_length={}
	for gene in exon:
		merged=bedmerge(exon[gene])
		merged_exon[gene]=merged
	
		# calculate merged gene length(sum of non-overlapping exons)
		length=0
		for each in merged:
			length+=each[2]-each[1]
		gene_length[gene]=length

	# assign merged exons and gene lengths
	merged = {}
	merged['merged_exon']=merged_exon
	merged['merged_gene_length']=gene_length


	# print merged exons if required
	if len(merged_exon_file) > 1:
		m=open(merged_exon_file,'w')
		for gene in merged_exon:
			imerged=merged_exon[gene]
			for each in imerged:
				joined="\t".join([each[0],str(each[1]),str(each[2]),gene,'0',each[3]])+"\n"
				m.write(joined)
		m.close()

	return merged


def get_independent_intron(GTFfile,independent_intron_file='iintron.tmp',chroms=map(str,range(1,23))+['X','Y']):
	# record exon coordination
	exon={}
	f=open(GTFfile)
	for line in f:
		table=line.split('\t')
		if line[0] != '#' and table[0] in chroms:  # skip comment line
			ichr = table[0]
			if table[2] == 'exon':
				iexon=(ichr,int(table[3])-1,int(table[4]),table[6])  # gtf to bed coordination
				if ichr in exon:
					exon[ichr].append(iexon)
				else:
					exon[ichr]=[iexon]
	f.close()
	
	# get gene coordinates in bed formats
	genebed=get_gene_bed(GTFfile,chroms=chroms)

	# calculate and write independent introns
	f = open(independent_intron_file,'w')
	record=[]
	for ichrom in chroms:
		if ichrom in genebed.keys():
			ichrom_gene = genebed[ichrom]
			ichrom_exon = exon[ichrom]
			sub = bed_subtract(ichrom_gene,ichrom_exon)
			for item in sub:
				if item[2] - item[1] >= 10:
					uniqued = unique_judge(item,ichrom_gene)	
					if len(uniqued)>0:
						record.append(uniqued)
	# write to file	
	f = open(independent_intron_file,'w')
	record = list(set(record))
	record.sort(key = lambda x: (x[0],x[1]))
	for item in record:
		ilength = item[2]-item[1]
		f.write('\t'.join([item[0],str(item[1]),str(item[2]),item[4],str(ilength)])+'\n')
	f.close()


def unique_judge(intron,genes):
	# intron: a single intron
	# genes: a list of genes
	hostgene=[]
	ngene=len(genes)
	i = 0
	count = 0
	flag = 0
	while i < ngene and flag == 0:
		igene = genes[i]
		if intron[1] > igene[1] and intron[1] < igene[2]:
			count = count + 1
			hostgene = igene
			if count > 1:
				flag = 1
		i = i+1

	if  count == 1:
		#return((intron[0],intron[1],intron[2],intron[3]))
		return((intron[0],intron[1],intron[2],intron[3],hostgene[4]))
	else:
		return([])




def bed_subtract(bedA,bedB):
	# assume that A and B are on the same chromosome.
	# calculate bedA-bedB
	bedA=bedmerge(bedA)
	bedB=bedmerge(bedB)
	AminusB=[]
	for ibed in bedA:
		ichr    = ibed[0]
		istart  = ibed[1]
		iend    = ibed[2]
		istrand = ibed[3]

		# find overlapping regions in bedB
		overlapped=[]
		for ibedB in bedB:
			ibedB = list(ibedB)
			judgeleft  = (ibedB[1] > istart and ibedB[1] < iend)
			judgeright = (ibedB[2] > istart and ibedB[2] < iend)
			if judgeleft == 1 or judgeright == 1:
				if ibedB[1] <= istart:
					ibedB[1] = istart
				if ibedB[2] >= iend:
					ibedB[2] = iend
				overlapped.append(ibedB)

		# calculate independent introns
		nexons = len(overlapped)
		if nexons == 0:
			AminusB.append((ichr,istart,iend,istrand))
		else:
			for i in range(nexons):
				iexon = overlapped[i]
				if i == 0:
					if iexon[1] > istart:
						AminusB.append((ichr,istart,iexon[1],istrand))
				if i >= 1:
					pexon = overlapped[i-1]
					AminusB.append((ichr,pexon[2],iexon[1],istrand))
				if i == nexons-1:
					if iexon[2] < iend:
						AminusB.append((ichr,iexon[2],iend,istrand))
	
	# return
	return(AminusB)


def subtract(fileA, fileB):
	A=[]
	f = open(fileA)
	for line in f:
		t = line.strip().split('\t')
		A.append((t[0],int(t[1]),int(t[2]),t[3]))
	f.close()

	B=[]
	f = open(fileB)
	for line in f:
		t = line.strip().split('\t')
		B.append((t[0],int(t[1]),int(t[2]),t[3]))
	f.close()

	A=bedmerge(A)
	for item in A:
		print('\t'.join([item[0],str(item[1]),str(item[2]),item[3]]))


	B=bedmerge(B)

	R=bed_subtract(A,B)
	length = len(R)
	print('length:'+str(length))
	for item in R:
		print('\t'.join([item[0],str(item[1]),str(item[2]),item[3]]))

def get_isoform_length(GTFfile,isoformlength_file='',chroms=map(str,range(1,23))+['X','Y']):
	# calculate length of isoforms
	f=open(GTFfile)
	isolength = {}
	gene2iso= {}
	iso2gene= {}
	for line in f:
		table=line.split('\t')
		if line[0] != '#' and table[0] in chroms:  # skip comment line
			if table[2] == 'exon':
				gene= line.split('gene_id')[1].split('"')[1]  # gene ID
				tcx = line.split('transcript_id')[1].split('"')[1]  # tcx ID
				exon_length = int(table[4])-int(table[3])+1

				# isoform length calculate
				if tcx in isolength:
					isolength[tcx] = isolength[tcx] + exon_length
				else:
					isolength[tcx] =  exon_length

				# record isoforms of genes
				if gene in gene2iso:
					gene2iso[gene].append(tcx)
				else:
					gene2iso[gene] = [tcx]
				iso2gene[tcx] = gene

	if len(isoformlength_file) > 1:
		f=open(isoformlength_file,'w')
		f.write('isoform\tgene\tlength\n')
		for thisiso in iso2gene:
			out=thisiso+'\t'+iso2gene[thisiso]+'\t'+str(isolength[thisiso])+'\n'
			f.write(out)
		f.close()	

	# return
	ret={}
	ret['gene2iso']  = gene2iso
	ret['isoform_length'] = isolength
	return(ret)

def get_gene_length(GTFfile,genelength_file='',chroms=map(str,range(1,23))+['X','Y']):
	# record exon coordination
	# calculate length of merged exons
	merged_gene_length = merge_exon(GTFfile,chroms=chroms)['merged_gene_length']

	# calculate length of each transcrpt isoform
	ret_isoform_length = get_isoform_length(GTFfile)
	gene2iso =ret_isoform_length['gene2iso']
	isolength=ret_isoform_length['isoform_length']

	# calculate gene length as mean, median or max of isoforms
	gene_length = {}
	for igene in gene2iso:
		isoforms = gene2iso[igene]
		length_set = []
		for iisoform in isoforms:
			length_set.append(isolength[iisoform])
		gene_length[igene] = [int(list_mean(length_set)),int(list_median(length_set)),max(length_set),merged_gene_length[igene]]

	# print
	if len(genelength_file) > 1:
		f = open(genelength_file,'w')
		f.write('gene\tmean\tmedian\tlongest_isoform\tmerged\n')
		for igene in gene_length:
			tmp =  gene_length[igene]
			tmplst = [igene,str(tmp[0]),str(tmp[1]),str(tmp[2]),str(tmp[3])]
			f.write('\t'.join(tmplst)+'\n')
		f.close()



def get_exon(GTFfile,exon_file='',chroms=map(str,range(1,23))+['X','Y']):
	# record exon coordination
	f=open(GTFfile)
	exon = {}
	iso2gene= {}
	for line in f:
		table=line.split('\t')
		if line[0] != '#' and table[0] in chroms:  # skip comment line
			if table[2] == 'exon':
				gene= line.split('gene_id')[1].split('"')[1]  # gene ID
				tcx = line.split('transcript_id')[1].split('"')[1]  # tcx ID
				iso2gene[tcx] = gene   # map isoforms to genes
				iexon=(table[0],int(table[3])-1,int(table[4]),table[6])  # gtf to bed coordination
				if tcx in exon:
					exon[tcx].append(iexon)
				else:
					exon[tcx]=[iexon]
	f.close()
	ret={}
	ret['exon']=exon
	ret['iso2gene']=iso2gene


	if len(exon_file) > 1:
		f=open(exon_file,'w')
		alltcx = exon.keys()
		for itcx in alltcx:
			exonlist = exon[itcx]
			for item in exonlist:
				out = map(str,item)+[itcx, iso2gene[itcx]]
				out = '\t'.join(out)+'\n'
				f.write(out)
	#return
	return(ret)




def cal_exon_skip_tags(exonlist):
	nexon = len(exonlist)
	nexon_middle = nexon-2
	#if nexon == 3:
	#	read_left  = 
	#	read_right =
	#	read_skip  = 

	#if nexon >=4:
	#for i in range(1,nexon_middle):
	#	for j in range(i+1,n



#def get_exon_skip_db(GTFfile,skipdb_file='',chroms=map(str,range(1,23))+['X','Y']):
#	# record exon coordination
#	exon_raw = get_exon(GTFfile,exon_file='',chroms=chroms)
#	exon = exon_raw['exon']
#	iso2gene= exon_raw['iso2gene']
#	alltcx = exon.keys()
#	for itcx in alltcx:
#		exonlist = exon[itcx]
#		nexon = len(exonlist)
#		if nexon >= 3:
#			xxx=cal_exon_skip_tags(exonlist)
		

def get_intron(GTFfile,intron_file='',chroms=map(str,range(1,23))+['X','Y']):
	# record exon coordination
	f=open(GTFfile)
	exon = {}
	iso2gene= {}
	for line in f:
		table=line.split('\t')
		if line[0] != '#' and table[0] in chroms:  # skip comment line
			if table[2] == 'exon':
				gene= line.split('gene_id')[1].split('"')[1]  # gene ID
				tcx = line.split('transcript_id')[1].split('"')[1]  # tcx ID
				iso2gene[tcx] = gene   # map isoforms to genes
				iexon=(table[0],int(table[3])-1,int(table[4]),table[6])  # gtf to bed coordination
				if tcx in exon:
					exon[tcx].append(iexon)
				else:
					exon[tcx]=[iexon]
	f.close()

	# calculate intron location based on exons for each tcx
	intron={}
	for tcx in exon:
		iexon = exon[tcx]
		iintron=exon2intron(iexon)
		if len(iintron) > 0:
			intron[tcx]=iintron
	to_return={}
	to_return['intron']=intron
	to_return['iso2gene']=iso2gene


	# print intron to file if required
	if len(intron_file) > 0:
		f=open(intron_file,'w')
		for tcx in intron:
			iintron=intron[tcx]
			for each in iintron:
				joined="\t".join([each[0],str(each[1]),str(each[2]),tcx,iso2gene[tcx],each[3]])+"\n"
				f.write(joined)
		f.close()


	# return
	return to_return

def chroms_interpreter(chroms_input):
	t = chroms_input.strip().split(',')
	chroms = []
	for it in t:
		if '-' in it:
			start  = int(it.split('-')[0].strip())
			end    = int(it.split('-')[1].strip())+1
			irange = map(str,range(start,end))
			chroms = chroms + irange
		else:
			chroms = chroms + [it]
	return(chroms)	



# convert GENCODE to ENSEMBL format
def gencode2ensembl(gtf1,gtf2):
#gtf1: ensembl GTF
#gtf2: gencode GTF

	f  = open(gtf1)
	fo = open(gtf2,'w')

	for line in f:
		if line[0:3] == 'chr':
			t = line.strip().split('\t')
			# update
			cname = t[0].split('chr')[1]
			if cname == 'M':
				cname = 'MT'

			# replace
			t[0] =  cname
			out = '\t'.join(t)+'\n'
			fo.write(out)
		else:
			fo.write(line)

	f.close()
	fo.close()


# check format of GTF: ensembl or gencode
def gtf_format_check(gtf):
	f=open(gtf)
	f.seek(10000)
	f.readline()
	chrom=f.readline().split('\t')[0]
	if len(chrom) >= 4:
		ftype = 'GENCODE'
	else:
		ftype = 'ENSEMBL'
	return(ftype)

############################################
##############   END OF FUNCTIONS  #########
############################################



############################################
##############   START: SIMPLE MATH ########
############################################
def list_sum(x):
	# assume x is lists of numbers
	s = 0
	for i in x:
		s = s + i
	return(s)

def list_mean(x):
	m = list_sum(x)/float(len(x))
	return m

def list_median(x):
	x.sort()
	length = len(x)
	if length % 2 == 0:
		half = length/2
		return((x[half-1]+x[half])/2)
	else:
		return(x[(length-1)/2])

############################################
##############   END: SIMPLE MATH  #########
############################################



############################################
#############  ENTRY  ######################
############################################


# argument parser
parser = argparse.ArgumentParser()  # create a parser
parser.add_argument('GTFfile',help="GTF file: only ENSEMBL or GENCODE GTF file accepted")
parser.add_argument('-c','--chroms',metavar='chromosomes',help="chromosome list to analyze. Chromosomes can be separated by comma(,) or dash(-). For example, '-c 1-5,X,Y' means chromosomes 1 to 5 plus X and Y. Default is: 1-22,X,Y")
parser.add_argument('-m','--merged_exon',metavar='merged_exon',help="file name for outputing merged exons from all isoforms of a gene in bed format")
parser.add_argument('-e','--exon',metavar='exon',help="file name for outputing exon coordination of splice isoforms in bed format")
parser.add_argument('-i','--intron',metavar='intron',help="file name for outputing intron coordination of splice isoforms in bed format")
parser.add_argument('-d','--independent_intron',metavar='independent_intron',help="file name for outputing independent intron coordination of genes. Independent introns refer to those introns that do not overlap with any exon of isoforms. It is calcualted by merging all exons of a chromosome followed by substracting them from gene regions.")
parser.add_argument('-l','--gene_length',metavar='gene_length',help="file name for gene length file. Gene length is calculated as the mean, median and max of isoforms. The length by merging exons of all isoforms is also calculated.")
parser.add_argument('-r','--isoform_length',metavar='isoform_length',help="file name for isoform length file. Isoform length is calculated as the summed length of its exons.mean, median and max of isoforms.")
parser.add_argument('-k','--masked_intron',metavar='masked_intron',help="file name for intron completely overlapping exons")
parser.add_argument('-u','--UTR',metavar='UTR',help="file name for UTR data")
parser.add_argument('-s','--isoform',metavar='isoform',help="file name for isoform bed data.")
parser.add_argument('-g','--gene',metavar='gene',help="file name for gene bed data.")
parser.add_argument('-t','--TSS',metavar='TSS',help="file name for a region centering at transcription start site (TSS). It is calculated as (TSS-w,TSS+w) where w is a user-specified distance, say 2000bp.")
parser.add_argument('-w','--window',metavar='window_size',help="the value of w in calculating TSS regions as described in '-t'. Default: 2000")
parser.add_argument('-v','--version',action='version',version="GTFtools version:0.6.5")
#parser.add_argument('-a','--file_a',help="bed file a.")
#parser.add_argument('-b','--file_b',help="bed file b.")
#parser.add_argument('-p','--merged_file',help="merge two bed files.")
#parser.add_argument('-q','--subtracted_file',help="subtract b from A, i.e. a-b.")


args = parser.parse_args()   # parse command-line arguments


######################################
#################  main   ############
######################################
GTFfile=args.GTFfile
ftype = gtf_format_check(GTFfile)
print('GTF format is '+ftype)
if ftype == 'GENCODE':
	print('Converting GTF to ENSEMBL format')
	GTFfile_ensembl=GTFfile+'.ensembl'
	gencode2ensembl(GTFfile,GTFfile_ensembl)
	GTFfile = GTFfile_ensembl


#print '------> Analyzing '+GTFfile

if args.chroms:
	chroms = chroms_interpreter(args.chroms)
	print 'Analyzing chromosomes:'
	print chroms
else:
	chroms = map(str,range(1,23))+['X','Y'] 

# print merged exon in bed format
if args.merged_exon:
	merge_exon(GTFfile,merged_exon_file=args.merged_exon,chroms=chroms)

# print gene length of merged exons
if args.gene_length:
	get_gene_length(GTFfile,genelength_file=args.gene_length,chroms=chroms)

# print isoform length of summed exons
if args.isoform_length:
	get_isoform_length(GTFfile,isoformlength_file=args.isoform_length,chroms=chroms)


# print introns of splice isoforms in bed format
if args.exon:
	get_exon(GTFfile,exon_file=args.exon,chroms=chroms)


# print introns of splice isoforms in bed format
if args.intron:
	get_intron(GTFfile,intron_file=args.intron,chroms=chroms)


# print introns of splice isoforms in bed format
if args.independent_intron:
	get_independent_intron(GTFfile,independent_intron_file=args.independent_intron,chroms=chroms)


# print masked introns in bed format
if args.masked_intron:
	masked_intron(GTFfile,maskedfile=args.masked_intron)


# print UTR in bed format
if args.UTR:
	get_UTR(GTFfile,utr_file=args.UTR,chroms=chroms)


# print gene bed
if args.gene:
	get_gene_bed(GTFfile,gene_bed_file=args.gene,chroms=chroms)


# print isoform bed
if args.isoform:
	get_isoform_bed(GTFfile,isoform_bed_file=args.isoform,chroms=chroms)


# print subtracted bed
#if args.subtracted_file and args.file_a and args.file_b:
#	subtract(args.file_a,args.file_b)


# print TSS neighborhood bed
if args.window:
	w = int(args.window)
else:
	w = 2000
if args.TSS:
	get_tss_region(GTFfile,w=w,tss_bed_file=args.TSS,chroms=chroms)

