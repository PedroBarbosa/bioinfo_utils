import sys
table=sys.argv[1]
with open(table, 'r') as infile:
    for line in infile:
        if line.startswith("#List of"):
            print("Please remove the first two lines of the table")
            exit(1)
        elif line.startswith("#"):
            header=line.rstrip().split("\t")
            index_snv=header.index("SpliceAI")        
            new_header=header[:index_snv] + header[index_snv+1 :] 
            if "SpliceAI_ind" in header:
                index_ind=header.index("SpliceAI_ind")
                index_ind_new_header = new_header.index("SpliceAI_ind")
                new_header=new_header[:index_ind_new_header] + new_header[index_ind_new_header+1 :]
            print('\t'.join(new_header) + '\t' + '\t'.join(["gene_spliceAI","DS_AG","DS_AL","DS_DG","DS_DL"]))
        else:
            values = line.rstrip().split('\t')
            spliceai=line.rstrip().split('\t')[index_snv]
            values = values[:index_snv] + values[index_snv+1 :]  
            if index_ind:     
                values = values[:index_ind_new_header] + values[index_ind_new_header+1 :]
                if spliceai == ".":
                    spliceai=line.rstrip().split('\t')[index_ind]

            split_spliceai=spliceai.split("|")[1:6]
            #print(split_spliceai)
            print('\t'.join(values) +  '\t' + '\t'.join(split_spliceai))
infile.close()
