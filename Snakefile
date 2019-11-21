SAMPLES, = glob_wildcards("{sample}_1.fq")  # read in file list
#TISSUE =  stuff

print(SAMPLES)
rule final:
    input:
        expand("{sample}_trinity", sample = SAMPLES)

rule fastp:
    input:
        raw1 = "{sample}_1.fq",
        raw2 = "{sample}_2.fq"
    output:
        p1 = "{sample}_1.processed.fq",
        p2 = "{sample}_2.processed.fq"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -c -r -M 3 -i {input.raw1} -I {input.raw2} -o {output.p1} -O {output.p2}"

rule banana:
    input:
        p1 = "{sample}_1.processed.fq",
        p2 = "{sample}_2.processed.fq"
    output:
        banana1 = "{sample}_1.processed_banana.fq",
        banana2 = "{sample}_2.processed_banana.fq"
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

rule trinity:
    input:
        banana1 = "{sample}_1.processed_banana.fq",
        banana2 = "{sample}_2.processed_banana.fq"
    output:
        trinity_dir = "{sample}_trinity",
        trinity_fasta = "{sample}_trinity/Trinity.fasta"
    conda:
        "envs/trinity.yaml"
    shell:
        "Trinity --seqType fq --max_memory 150G  --left {input.banana1} --right {input.banana2} --CPU 20 --full_cleanup --output {output.trinity_dir}"


rule transdecoder:
    input:
        "{sample}_trinity/Trinity.fasta"
    output:
        "{sample}_trinity/Trinity.fasta.TransDecoder_dir"
    conda:
        "envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.LongOrfs -t {input}  -m 30
        TransDecoder.Predict -t {input} --single_best_orf
        """
rule supertranscript:
    input:
        "{sample}_trinity/Trinity.fasta.TransDecoder_dir"
    output:
        "{sample}_supertranscript.fasta"
    conda:
        "envs/trinity.yaml"
    shell:
        "supertranscript -i {input} "
rule salmon:
    input:
        supertranscript = "{sample}_supertranscript.fasta"
        banana1 = "{sample}_1.processed_banana.fq"
        banana2 = "{sample}_2.processed_banana.fq"
    output:
        "{sample}_quant.sf"
    conda:
        "envs/trinity.yaml"
    shell:
        """
        salmon index {input.supertranscript}
        salmon quant {input.banana1}{input.banana2}
        """
