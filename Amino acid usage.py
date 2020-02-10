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
#aa_count = sum([codons.count(aa[0]) for aa in codonDict[codon]) # Sum up all amino acids in aaDict
#sc_score = sum((aa[1])for codon, aa in aaDict.items()) # Sum of all sc scores in aaDict
#for codon, aa in aaDict.items():
    #    total_sc += aa[1] # add the sc score for the amino acid to the running total sc score of sequence
    #    print(total_sc) # just a test to see if this is properly adding the sc scores 