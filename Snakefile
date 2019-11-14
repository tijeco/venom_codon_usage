SAMPLES, = glob_wildcards("{sample}*_1.fq")  # read in file list

rule final:
    input:
        raw1 = "{sample}_1.fq",
        raw2 = "{sample}_2.fq"
    output:
        "{sample}_trinity/Trinity.fasta"

rule fastp:
    input:
        raw1 = "{sample}_1.fq",
        raw2 = "{sample}_2.fq"
    output:
        p1 = "{sample}_1.processed.fq",
        p2 = "{sample}_2.processed.fq"
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
        "{sample}_trinity"
    shell:
        "Trinity --seqType fq --max_memory 150G  --left {input.banana1} --right {input.banana2} --CPU 20 --full_cleanup --output {output}"
