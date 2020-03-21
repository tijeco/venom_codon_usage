from Bio import SeqIO
import pandas as pd
import json

def calcAminoUsage(cds_file):

    aaDict = {'GCT':['Ala',4.76],'GCC':['Ala',4.76],'GCA':['Ala', 4.76],'GCG':['Ala',4.76],
                'CGT':['Arg', 56.34],'CGC':['Arg', 56.34],'CGA':['Arg', 56.34],'CGG':['Arg', 56.34],'AGA':['Arg', 56.34],'AGG':['Arg', 56.34],
                'AAT':['Asn', 33.72],'AAC':['Asn', 33.72],
                'GAT':['Asp', 32.72],'GAC':['Asp', 32.72],
                'TGT':['Cys', 57.16],'TGC':['Cys', 57.16],
                'CAA':['Gln', 37.48],'CAG':['Gln', 37.48],
                'GAA':['Glu', 36.48],'GAG':['Glu', 36.48],
                'GGT':['Gly', 1.00],'GGC':['Gly', 1.00],'GGA':['Gly', 1.00],'GGG':['Gly', 1.00],
                'CAT':['His', 58.70],'CAC':['His', 58.70],
                'ATT':['Ile', 16.04],'ATC':['Ile', 16.04],'ATA':['Ile', 16.04],
                'TTA':['Leu', 16.04],'TTG':['Leu', 16.04],'CTT':['Leu', 16.04],'CTC':['Leu', 16.04],'CTA':['Leu', 16.04],'CTG':['Leu', 16.04],
                'AAA':['Lys', 30.14],'AAG':['Lys', 30.14],
                'ATG':['Met', 64.68],
                'TTT':['Phe', 44.0],'TTC':['Phe', 44.0],
                'CCT':['Pro', 31.8],'CCC':['Pro', 31.8],'CCA':['Pro', 31.8],'CCG':['Pro', 31.8],
                'TCA':['Ser', 17.86],'TCT':['Ser', 17.86],'TCC':['Ser', 17.86],'TCG':['Ser', 17.86],'AGT':['Ser', 17.86],'AGC':['Ser', 17.86],
                'ACT':['Thr', 21.62],'ACC':['Thr', 21.62],'ACA':['Thr', 21.62],'ACG':['Thr', 21.62],
                'TAA':['Stop', 0],'TAG':['Stop', 0],'TGA':['Stop', 0],
                'TGG':['Trp', 73.0],
                'TAT':['Tyr', 57.00],'TAC':['Tyr', 57.00],
                'GTA':['Val', 12.28],'GTT':['Val', 12.28],'GTC':['Val', 12.28],'GTG':['Val', 12.28],
                }

    sc_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        header = record.description
        seq = str(record.seq)
        n = 3
        codons = [seq[i:i+n] for i in range(0, len(seq), n)] #an array including the values for the sequence ie "ATG","AAA",etc.
        for codon, aa in aaDict.items():    # for every codon in aaDict access both values for the codon
            aa_count = len(codons) - 1 #the total number of aa in a sequence which is the same as the number of codons in the codons array
        total_sc = 0.0 #have to give the variable a value before the loop
        for codon in codons: #for each codon in the codons array
            total_sc += aaDict[codon][1] # the total sc score of a sequence is calculated by searching the dictionary for the values of the codon and adding it to the running total
        mean_sc = total_sc / aa_count # The average sc formula
        sc_dict[header] = [mean_sc,total_sc] # creating a dictionary using header as the key and meansc, and total sc as values

    return sc_dict

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
    for record in SeqIO.parse(cds_file, "fasta"):
        header = record.description
        seq = str(record.seq)
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

def fop_func(optimal_codon_file,cds_file):
    optimal_codon_dict = {}


    with open(optimal_codon_file) as f:
        line1 = 1
        for line in f:
            if line1:
                line1 = 0
                continue
            row = line.strip().split(",")
            optimal_codon_dict[row[1]] = True # if weighted average, put deltaRSCU_norm row here

    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        header = record.description
        seq = str(record.seq)
        nop = 0
        if header not in cds_dict:
            n = 3
            seq_codons = [seq[i:i+n] for i in range(0, len(seq), n)]
            for codon in optimal_codon_dict:
                nop += seq_codons.count(codon) # for weight, this would be multiplied by optimal_codon_dict[codon] (deltaRSCU_norm)
            fop = nop / len(seq_codons)
            cds_dict[header] = fop
                    # continue
                # continue

    return pd.DataFrame(list(cds_dict.items()),columns=['Name', 'fop'])
    # cds seq  will be in file rscu output csv



# SAMPLES_venom, = glob_wildcards("{sample}_venom_1.fq")  # read in file list
# SAMPLES_body, = glob_wildcards("{sample}_body_1.fq")  # read in file list
SAMPLES, = glob_wildcards("{sample}_venom_1.fq")


print(SAMPLES)
rule final:
    input:
        expand("{sample}_aminoAcidUsage.csv", sample = SAMPLES) 
        # expand("{sample}.combined_5percent_fop.csv", sample = SAMPLES)
        # expand("{sample}_merged_quant.csv", sample = SAMPLES)
        # expand("{sample}.fop.csv", sample = SAMPLES)
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
        "fastp -i {input.raw1} -I {input.raw2} -o {output.p1} -O {output.p2}"
        # "fastp -c -r -M 3 -i {input.raw1} -I {input.raw2} -o {output.p1} -O {output.p2}"

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

rule longest_isoform:
    input:
        "{sample}_trinity.Trinity.fasta.transdecoder.cds"
    output:
        "{sample}_complete_longest_isoform.cds"
    run:
        isoform_dict = {}
        for record in SeqIO.parse(input[0], "fasta"):
            header = record.description
            seq = str(record.seq)
            # print("header",header)
            if "complete" in header:
                # print("Found complete gene!")
                gene = header.split("_i")[0]
                if gene not in isoform_dict:
                    isoform_dict[gene] = (header,seq)
                elif len(seq) > len(isoform_dict[gene][1]):
                    isoform_dict[gene] = (header,seq)
        print("Genes processed",len(isoform_dict))
        with open(output[0],"w") as out:
            for gene in isoform_dict:
                header = gene
                seq = isoform_dict[gene][1] # filter by longest isoform, somehow
                out.write(">" + header + "\n")
                out.write(seq + "\n")

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
        cds = "{sample}_complete_longest_isoform.cds"
    output:
        rscu = "{sample}_body.rscu.csv",
        json = "{sample}_body.rscu.json"
    run:
        quant_file = input.quant
        quant_df  = pd.read_csv(quant_file, sep='\t', header=0)
        quant_over2TPM = quant_df[quant_df["TPM"] > 2]
        num_seqs = quant_over2TPM.sort_values("TPM").shape[0]
        five_percent = round(quant_over2TPM.sort_values("TPM").shape[0] * 0.05)
        bottom5 = quant_over2TPM.sort_values("TPM")[:five_percent]
        top5 = quant_over2TPM.sort_values("TPM")[num_seqs-five_percent:]

        rscu_dict = calcRSCU(input.cds)
        with open(output.json, 'w') as json_file:
            json.dump(rscu_dict, json_file)

        with open(output.rscu,"w") as out:
            out.write("header,aa,codon,rscu,class\n")
            for aa in rscu_dict:
                for codon in rscu_dict[aa]:
                    for header in rscu_dict[aa][codon]:
                        if header in list(bottom5["Name"]):
                            out.write(header + "," + aa + "," + codon + "," + str(rscu_dict[aa][codon][header][0]) + ",bottom_5percent\n")
                        elif header in list(top5["Name"]):
                            out.write(header + "," + aa + "," + codon + "," + str(rscu_dict[aa][codon][header][0]) + ",top_5percent\n")

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
rule fop:
    input:
        optimal_codon = "{sample}_body.optimalCodon.csv",
        cds = "{sample}_complete_longest_isoform.cds"
    output:
        "{sample}.fop.csv"
    run:
        fop_pd = fop_func(input.optimal_codon,input.cds)
        fop_pd.to_csv(output[0],index=False)

rule merge_quant:
    input:
        script = "src/merge_quant.R",
        body_quant = "{sample}_body_quant/quant.sf",
        venom_quant = "{sample}_venom_quant/quant.sf"
    output:
        "{sample}_merged_quant.csv"
    conda:
        "envs/r.yaml"
    shell:
        "Rscript {input.script} -b {input.body_quant} -v {input.venom_quant} -o {output}"

rule merge_fop:
    input:
        script = "src/fop.R",
        quant = "{sample}_merged_quant.csv",
        fop = "{sample}.fop.csv"
    output:
        combined_5percent = "{sample}.combined_5percent_fop.csv",
        test = "{sample}.t_test.csv",
        violin = "{sample}.combined_5percent_fop.png"

    conda:
        "envs/r.yaml"
    shell:
        "Rscript {input.script} -q {input.quant} -f {input.fop} -o {output.combined_5percent} -v {output.violin} -s {output.test}"
rule aa_usage:
    input:
        cds = "{sample}_complete_longest_isoform.cds"
    output:
        aa_usage = "{sample}_aminoAcidUsage.csv"
    run:
        sc_dict = calcAminoUsage(input.cds)
        with open(output.sc,"w") as out:
            out.write("header,mean_sc,total_sc\n")
            for header in sc_dict:
                out.write(header + "," + sc_dict[header][0] + "," + sc_dict[header][1] + "\n")

##### Will update this soon #########
# rule merge_aa_usage
#     input:
#         script =
#     output:
#         test =
#     conda:
#         "envs/r.yaml"
###########################
