from collections import Counter
import sys

counter = Counter({'a','t','c','g'})
allgenome=0
i=0
with open(sys.argv[1], "r") as f:  # Always open files like this
    for line in f:  # no need for readlines
        line = line.strip("\n")		
        if not line.startswith(">"):
            counter.update({'a':line.count('a'),'t':line.count('t'),'c':line.count('c'),'g':line.count('g')})
            allgenome+=len(line)
        #else:
        #    i+=1
        #    print("%i contigs processed." % i)
    f.close()
    print("Total size of genome:\t%i" % allgenome)
    maskedgenome=sum(counter.values())
    print("Total number of bases masked:\t%i" % maskedgenome)
    print("Percentage genome masked:\t%f" % round(maskedgenome/allgenome*100,2))
#    print("Frequency of each masked base:\ta=%f,t=%f,c=%f,g=%f" % (round(counter['a']/maskedgenome*100,2),round(counter['t']/maskedgenome*100,2),round(counter['c']/maskedgenome*100,2),round(counter['g']/maskedgenome*100,2))                
