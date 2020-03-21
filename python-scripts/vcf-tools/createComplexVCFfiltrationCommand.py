import argparse
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def processGroupsGATK(samples,genotypes):
    cmd=""
    for group in samples:
        genotype = genotypes[samples.index(group)]

        with open(group,'r') as file:
            if genotype == "homRef":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"vc.getGenotype('" + sample +"').isHomRef()"
                    else:
                        cmd = cmd + " && vc.getGenotype('" + sample +"').isHomRef()"

            elif genotype == "homAlt":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"vc.getGenotype('" + sample +"').isHomVar()"
                    else:
                        cmd = cmd + " && vc.getGenotype('" + sample +"').isHomVar()"

            elif genotype == "het":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"vc.getGenotype('" + sample +"').isHet()"
                    else:
                        cmd = cmd + " && vc.getGenotype('" + sample +"').isHet()"
    cmd = cmd + "\""
    return cmd


def processGroupsSNPSift(samples,genotypes):

    cmd=""
    for group in samples:
        genotype = genotypes[samples.index(group)]

        with open(group,'r') as file:
            if genotype == "homRef":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"isHom[" + sample + "] & isRef[" + sample + "]"
                    else:
                        cmd = cmd + " & isHom[" + sample + "] & isRef[" + sample + "]"

            elif genotype == "homAlt":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"isHom[" + sample + "] & isVariant[" + sample + "]"
                    else:
                        cmd = cmd + " & isHom[" + sample + "] & isVariant[" + sample + "]"

            elif genotype == "het":
                for line in file:
                    sample = line.rstrip()
                    if not cmd:
                        cmd = "\"isHet[" + sample + "]"
                    else:
                        cmd = cmd + " & isHet[" + sample + "]"
    cmd = cmd + "\""
    return cmd



def main():

    parser = argparse.ArgumentParser(description='Script to generate complex queries to filter variants based on the genotype of groups of samples.')
    parser.add_argument(dest='VCF_file', metavar='vcf_file', nargs=1, help='VCF file to filter.')
    parser.add_argument(dest='FASTA_file', metavar='fasta_file', nargs=1, help='Fasta file of the reference used to call variants.')
    parser.add_argument('-s','--samples', required=True, nargs='+', help='List of files, each one representing the samples belonging to each group. One sample per line.')
    parser.add_argument('-g','--genotypes', required=True, nargs='+',choices=['homRef','homAlt','het'],
                        help='Genotypes desired for each group. List should be the same length as "-s" list.')
    parser.add_argument('-p','--program', required=True, choices=['snpsift','gatk'], help='Program to base the generation of the command. Check script to change the path of executables.')

    args = parser.parse_args()

    if len(args.samples) != len(args.genotypes):
        logging.error("Number of samples is different from genotypes provided. Please set 1 genotype per group. Exiting..")
        exit(1)



    if args.program == "snpsift":
        exec_snpsift = "cat " + args.VCF_file[0] + " | java -jar /opt/tools/snpEff-v4.2/SnpSift.jar filter "
        cmd = processGroupsSNPSift(args.samples, args.genotypes)
        logging.info("Please copy the command you should run:\n")
        print("%s%s" % (exec_snpsift,cmd))

    elif args.program == "gatk":
        exec_gatk = "java -jar -Xmx100G /opt/tools/GATK-v3.6/GenomeAnalysisTK.jar -T SelectVariants -R " + args.FASTA_file[0] + " -V " + args.VCF_file[0] +  " -select "
        cmd = processGroupsGATK(args.samples, args.genotypes)
        logging.info("Please copy the command you should run:\n")
        print("%s%s" % (exec_gatk,cmd))

if __name__ == "__main__":
    main()

