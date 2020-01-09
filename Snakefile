from Bio import SeqIO
import pandas as pd
def calcRSCU(cds_file):

    codonDict = {'Ala': {'GCT': {}, 'GCC': {}, 'GCA': {}, 'GCG': {}},
                 'Arg': {'CGT': {}, 'CGC': {}, 'CGA': {}, 'CGG': {}, 'AGA': {}, 'AGG': {}},
                 'Asn': {'AAT': {}, 'AAC': {}},
                 'Asp': {'GAT': {}, 'GAC': {}},
                 'Cys': {'TGT': {}, 'TGC': {}},
                 'Gln': {'CAA': {}, 'CAG': {}},
                 'Glu': {'GAA': {}, 'GAG': {}},
                 'Gly': {'GGT': {}, 'GGC': {}, 'GGA': {}, 'GGG': {}},
                 'His': {'CAT': {}, 'CAC': {}},
                 'Ile': {'ATT': {}, 'ATC': {}, 'ATA': {}},
                 'Leu': {'TTA': {}, 'TTG': {}, 'CTT': {}, 'CTC': {}, 'CTA': {}, 'CTG': {}},
                 'Lys': {'AAA': {}, 'AAG': {}},
                 'Phe': {'TTT': {}, 'TTC': {}},
                 'Pro': {'CCT': {}, 'CCC': {}, 'CCA': {}, 'CCG': {}},
                 'Ser': {'TCA': {}, 'TCT': {}, 'TCC': {}, 'TCG': {}, 'AGT': {}, 'AGC': {}},
                 'Thr': {'ACT': {}, 'ACC': {}, 'ACA': {}, 'ACG': {}},
                 'Stop': {'TAA': {}, 'TAG': {}, 'TGA': {}},
                 'Val': {'GTA': {}, 'GTT': {}, 'GTC':{}, 'GTG': {}},
                 'Tyr': {'TAT': {}, 'TAC': {}}
                }

    isoform_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        header, seq = record.description,str(record.seq)
        if "complete" in header:
            gene = header.split("_i")[0]
            if gene not in isoform_dict:
                isoform_dict[gene] = (header,seq)
            elif len(seq) > len(isoform_dict[gene][1]):
                isoform_dict[gene] = (header,seq)
    for gene in isoform_dict:
        # print(isoform_dict[gene][0])
        header = gene
        seq = isoform_dict[gene][1]
            # filter by longest isoform, somehow
        n = 3
        codons = [seq[i:i+n] for i in range(0, len(seq), n)]

        # print(header)
        # print(codons)
        for aa in codonDict:
            # print(aa)
            aa_codonCount = len(codonDict[aa])
            sum_redundantCodons = sum([codons.count(codon) for codon in codonDict[aa] ])
            for codon in codonDict[aa]:
                observed_codonCount = codons.count(codon)

                if ((1/aa_codonCount)*sum_redundantCodons) != 0:
                    rscu = observed_codonCount / ((1/aa_codonCount)*sum_redundantCodons)
                else:
                    rscu = 0

                codonDict[aa][codon][header] = (rscu,seq)

    return codonDict

def fop_func(optimal_codon_file,rscu_file):
    optimal_codon_dict = {}

    with open(optimal_codon_file) as f:
        line1 = 1
        for line in f:
            if line1:
                line1 = 0
                continue
            row = line.strip().split(",")
            optimal_codon_dict[row[1]] = True # if weighted average, put deltaRSCU_norm row here
    # quant_file = input.quant

    # quant_df  = pd.read_csv(quant_file, sep='\t', header=0)
    # quant_over2TPM = quant_df[quant_df["TPM"] > 2]
    # num_seqs = quant_over2TPM.sort_values("TPM").shape[0]
    # five_percent = round(quant_over2TPM.sort_values("TPM").shape[0] * 0.05)
    # top5 = quant_over2TPM.sort_values("TPM")[num_seqs-five_percent:]

    cds_dict = {}
    with open(rscu_file) as f:
        line1 = 1


        for line in f:
            nop = 0
            if line1:
                line1 = 0
                continue
            row = line.strip().split(",")
            header = row[0]
            seq = row[4]
            if header not in cds_dict:
                n = 3
                seq_codons = [seq[i:i+n] for i in range(0, len(seq), n)]
                for codon in optimal_codon_dict:
                    nop += seq_codons.count(codon) # for weight, this would be multiplied by optimal_codon_dict[codon] (deltaRSCU_norm)
                fop = nop / len(seq_codons)
                cds_dict[header] = fop
    return pd.DataFrame.from_dict(list(cds_dict,columns=['header', 'fop']))
    # cds seq  will be in file rscu output csv



SAMPLES_venom, = glob_wildcards("{sample}_venom_1.fq")  # read in file list
SAMPLES_body, = glob_wildcards("{sample}_body_1.fq")  # read in file list
SAMPLES, = glob_wildcards("{sample}_venom_1.fq")


print(SAMPLES)
rule final:
    input:
        expand("{sample}.fop.tmp", sample = SAMPLES)
        # expand("{sample}_body.optimalCodon.csv", sample = SAMPLES)
        # expand("{sample}_body.rscu.csv", sample = SAMPLES)
        # expand("{sample}_trinity.Trinity.fasta.transdecoder.cds", sample = SAMPLES)
        # expand("{sample}_supertranscript.fasta", sample = SAMPLES)
        # expand("{sample}_trinity.Trinity.fasta", sample = SAMPLES)


rule fastp:
    input:
        raw1 = "{sample}_{tissue}_1.fq",
        raw2 = "{sample}_{tissue}_2.fq"
    output:
        p1 = "{sample}_{tissue}_1.processed.fq",
        p2 = "{sample}_{tissue}_2.processed.fq"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -c -r -M 3 -i {input.raw1} -I {input.raw2} -o {output.p1} -O {output.p2}"

rule banana:
    input:
        p1 = "{sample}_{tissue}_1.processed.fq",
        p2 = "{sample}_{tissue}_2.processed.fq"
    output:
        banana1 = "{sample}_{tissue}_1.processed_banana.fq",
        banana2 = "{sample}_{tissue}_2.processed_banana.fq"
    run:
        with open(output.banana1,"w") as out:
            number = 0
            with open(input.p1) as f:
                for line in f:

                    if number%4:
                        out.write(line) #out.write("@banana"+str(number)+"/1\n")
                    else:
                        out.write("@banana"+str(number)+"/1\n") #out.write(line)
                    number+=1

        with open(output.banana2,"w") as out:
            number = 0
            with open(input.p2) as f:
                for line in f:
                    if number%4: #line[0] == "@":
                        out.write(line) #out.write("@banana"+str(number)+"/2\n")
                    else:
                        out.write("@banana"+str(number)+"/2\n") #out.write(line)
                    number+=1

# venom / body needed
rule trinity:
    input:
        body_banana1 = "{sample}_body_1.processed_banana.fq",
        body_banana2 = "{sample}_body_2.processed_banana.fq",
        venom_banana1 = "{sample}_venom_1.processed_banana.fq",
        venom_banana2 = "{sample}_venom_2.processed_banana.fq"
    output:
         "{sample}_trinity.Trinity.fasta",

    conda:
        "envs/trinity.yaml"
    threads: 64
    shell:
        "out={output};Trinity --seqType fq --max_memory 150G  --left {input.venom_banana1},{input.body_banana1} --right {input.venom_banana2},{input.body_banana2} --CPU {threads} --full_cleanup --output  ${{out%.Trinity.fasta}}"
        # "echo {input.body_banana1} {input.body_banana2} {input.venom_banana1} {input.venom_banana2} > {output}"
        # "echo {input.body_banana1} {input.body_banana2} {input.venom_banana1} {input.venom_banana2} > {output}"
# this will also be venom and body combined
rule transdecoder:
    input:
        "{sample}_trinity.Trinity.fasta"
    output:
        "{sample}_trinity.Trinity.fasta.transdecoder.cds"
    conda:
        "envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.LongOrfs -t {input}  -m 30
        TransDecoder.Predict -t {input} --single_best_only
        """

# this will also be venom and body combined
rule supertranscript:
    input:
        "{sample}_trinity.Trinity.fasta"
    output:
        "{sample}_supertranscript.fasta"
    conda:
        "envs/trinity.yaml"
    shell:
        "out={output};Trinity_gene_splice_modeler.py --trinity_fasta {input} --out_prefix ${{out%.fasta}}"

# body / venom doesn't matter
# map venom reads to combined, and map body reads to combined
rule salmon_index:
    input:
        "{sample}_supertranscript.fasta"
    output:
        directory("{sample}_index")
    conda:
        "envs/trinity.yaml"
    shell:
        "salmon index -t {input} -i {output}"
rule salmon_venom_quant:
    input:
        index = "{sample}_index",
        banana1 = "{sample}_venom_1.processed_banana.fq",
        banana2 = "{sample}_venom_2.processed_banana.fq"
    output:
        "{sample}_venom_quant/quant.sf"
    conda:
        "envs/trinity.yaml"
    shell:
        "out={output};salmon quant -i {input.index} -l A -1 {input.banana1} -2 {input.banana2} -o ${{out%quant.sf}}"

rule salmon_body_quant:
    input:
        index = "{sample}_index",
        banana1 = "{sample}_body_1.processed_banana.fq",
        banana2 = "{sample}_body_2.processed_banana.fq"
    output:
        "{sample}_body_quant/quant.sf"
    conda:
        "envs/trinity.yaml"
    shell:
        "out={output};salmon quant -i {input.index} -l A -1 {input.banana1} -2 {input.banana2} -o ${{out%quant.sf}}"


rule rscu:
    input:
        quant = "{sample}_body_quant/quant.sf",
        cds = "{sample}_trinity.Trinity.fasta.transdecoder.cds" # just body
    output:
        "{sample}_body.rscu.csv"
    run:
        quant_file = input.quant
        quant_df  = pd.read_csv(quant_file, sep='\t', header=0)
        quant_over2TPM = quant_df[quant_df["TPM"] > 2]
        num_seqs = quant_over2TPM.sort_values("TPM").shape[0]
        five_percent = round(quant_over2TPM.sort_values("TPM").shape[0] * 0.05)
        bottom5 = quant_over2TPM.sort_values("TPM")[:five_percent]
        top5 = quant_over2TPM.sort_values("TPM")[num_seqs-five_percent:]

        rscu_dict = calcRSCU(input.cds)
        with open(output[0],"w") as out:
            out.write("header,aa,codon,rscu,seq,class\n")
            for aa in rscu_dict:
                for codon in rscu_dict[aa]:
                    for header in rscu_dict[aa][codon]:
                        if header in list(bottom5["Name"]):
                            out.write(header + "," + aa + "," + codon + "," + str(rscu_dict[aa][codon][header][0]) + "," + rscu_dict[aa][codon][header][1] + ",bottom_5percent\n")
                        elif header in list(top5["Name"]):
                            out.write(header + "," + aa + "," + codon + "," + str(rscu_dict[aa][codon][header][0]) + "," + rscu_dict[aa][codon][header][1] + ",top_5percent\n")

rule optimal_codon:
    input:
        script = "src/rscu.R",
        rscu = "{sample}_body.rscu.csv"
    output:
        deltaRSCU = "{sample}_body.deltaRSCU.csv",
        optimalCodon = "{sample}_body.optimalCodon.csv",
        optimalCodonFig = "{sample}_body.optimalCodon.png"
    conda:
        "envs/r.yaml"
    shell:
        "Rscript {input.script} -d {input.rscu} -r {output.deltaRSCU} -o {output.optimalCodon} -f {output.optimalCodonFig}"


# fop
# optimial_codon
# cds
# for seq in cds
    # count optimal codons

# quant
# subset cds fop that is only in the top 5%
rule body_fop:
    input:
        optimal_codon = "{sample}_body.optimalCodon.csv",
        cds = "{sample}_trinity.Trinity.fasta.transdecoder.cds"
    output:
        "{sample}_body.fop.csv"
    run:
        fop_pd = fop_func(optimal_codon,cds)
        fop_pd.to_csv(output[0],index=False)

rule venom_fop:
    input:
        optimal_codon = "{sample}_body.optimalCodon.csv",
        cds = "{sample}_trinity.Trinity.fasta.transdecoder.cds"
    output:
        "{sample}_venom.fop.csv"
    run:
        fop_pd = fop_func(optimal_codon,cds)
        fop_pd.to_csv(output[0],index=False)
rule fop:
    input:
        venom = "{sample}_venom.fop.csv",
        body = "{sample}_body.fop.csv"
    output:
        "{sample}.fop.tmp"
    shell:
        "touch {output}"



        # rscu_panda = pd.DataFrame.from_dict({(i,j): rscu_dict[i][j]
        #                    for i in rscu_dict.keys()
        #                    for j in rscu_dict[i].keys()})
        # rscu_panda.to_csv(output[0])
    #     with open([output],"w") as out:
    #         out.write("temp")

        # get quant file, split into top and bottom 5 percent. Write to file

        # columns will be (header,aa,codon,rscu, high/low)
        # Sort tpm values from highest to lowest,
        # then count total number of tpm values,
        # set n = total number of tpm values * 0.05,
        # Top5% = first n values, low5% = last n values
