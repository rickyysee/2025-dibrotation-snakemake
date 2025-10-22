SAMPLES = ["SRR1976948"]
RUNS = [1,2]
EMAIL = "rcantua@ucdavis.edu"

onsuccess:
    shell("sendmail {EMAIL} < {log}")

onerror:
    shell("sendmail {EMAIL} < {log}")

rule all:
    input:
        fastqczip = expand("fastqc/{sample}_{run}_fastqc.zip", sample=SAMPLES, run=RUNS),
        fastqchtml = expand("fastqc/{sample}_{run}_fastqc.html", sample=SAMPLES, run=RUNS),
        abund_gather = expand("sourmash/{sample}-abundtrim-gather.csv", sample=SAMPLES),
        raw_gather = expand("sourmash/{sample}-gather.csv", sample=SAMPLES),
        desulfo_sgc = "sgc/sgc-done.txt",

rule get_reads:
    output:
        out1 = expand("raw-data/{sample}_1.fastq.gz", sample=SAMPLES),
        out2 = expand("raw-data/{sample}_2.fastq.gz", sample=SAMPLES),
    params:
        id = expand("{sample}", sample=SAMPLES),
    shell:
        """
        mkdir -p raw-data
        wget -P raw-data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/{params.id}/{params.id}_1.fastq.gz
        wget -P raw-data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/{params.id}/{params.id}_2.fastq.gz
        """

rule fastqc:
    input:
        reads = expand("raw-data/{sample}_{run}.fastq.gz", sample=SAMPLES, run=RUNS),
    output:
        zip = expand("fastqc/{sample}_{run}_fastqc.zip", sample=SAMPLES, run=RUNS),
        html = expand("fastqc/{sample}_{run}_fastqc.html", sample=SAMPLES, run=RUNS),
    params:
        outdir="fastqc",
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p fastqc
        fastqc {input.reads} --outdir {params.outdir}
        """

rule trim_reads:
    input:
        in1 = expand("raw-data/{sample}_1.fastq.gz", sample=SAMPLES),
        in2 = expand("raw-data/{sample}_2.fastq.gz", sample=SAMPLES),
    output:
        out1 = expand("trim/{sample}_1.trim.fastq.gz", sample=SAMPLES),
        out2 = expand("trim/{sample}_2.trim.fastq.gz", sample=SAMPLES),
        json = expand("trim/{sample}.fastp.json", sample=SAMPLES),
        html = expand("trim/{sample}.fastp.html", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p trim
        fastp --in1 {input.in1}  --in2 {input.in2}  \
        --out1 {output.out1} --out2 {output.out2}  \
        --detect_adapter_for_pe  --qualified_quality_phred 4 \
        --length_required 31 --correction \
        --json {output.json} --html {output.html}
        """

rule khmer:
    input:
        in1 = expand("trim/{sample}_1.trim.fastq.gz", sample=SAMPLES),
        in2 = expand("trim/{sample}_2.trim.fastq.gz", sample=SAMPLES),
    output:
        out = expand("khmer/{sample}.abundtrim.fq.gz", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p khmer
        interleave-reads.py {input.in1} {input.in2} | trim-low-abund.py --gzip -C 3 -Z 18 -M 20e9 -V - -o {output}
        """

rule raw_sketch:
    input:
        in1 = expand("raw-data/{sample}_1.fastq.gz", sample=SAMPLES),
        in2 = expand("raw-data/{sample}_2.fastq.gz", sample=SAMPLES),
    output:
        sketch = expand("sourmash/{sample}.sig", sample=SAMPLES),
    params:
        id = expand("{sample}", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p sourmash
        sourmash sketch dna -o {output} \
        --name {params.id} \
        -p scaled=2000,k=21,k=31,k=51,abund \
        {input}
        """

rule get_database:
    output:
        "databases/genbank-k31.lca.json",
    shell:
        """
        mkdir -p databases
        curl -L https://osf.io/4f8n3/download -o databases/genbank-k31.lca.json.gz
        gunzip databases/genbank-k31.lca.json.gz
        """

rule raw_gather:
    input:
        sketch = expand("sourmash/{sample}.sig", sample=SAMPLES),
        database = expand("databases/genbank-k31.lca.json", sample=SAMPLES),
    output:
        gather = expand("sourmash/{sample}-gather.csv", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        sourmash gather -o {output} \
        {input.sketch} {input.database}
        """

rule abundtrim_sketch:
    input:
        in1 = expand("khmer/{sample}.abundtrim.fq.gz", sample=SAMPLES),
    output:
        sketch = expand("sourmash/{sample}-abundtrim.sig", sample=SAMPLES),
    params:
        id = expand("{sample}", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p sourmash
        sourmash sketch dna -o {output} \
        --name {params.id} \
        -p scaled=2000,k=21,k=31,k=51,abund \
        {input}
        """

rule abundtrim_gather:
    input:
        sketch = expand("sourmash/{sample}-abundtrim.sig", sample=SAMPLES),
        database = expand("databases/genbank-k31.lca.json", sample=SAMPLES),
    output:
        gather = expand("sourmash/{sample}-abundtrim-gather.csv", sample=SAMPLES),
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        sourmash gather -o {output} \
        {input.sketch} {input.database}
        """

rule get_sgc:
    output:
        out1 = "spacegraphcats/setup.py",
        out2 = "envs/sgc.yaml",
    shell:
        """
        git clone https://github.com/spacegraphcats/spacegraphcats/
        cp spacegraphcats/environment.yml envs/sgc.yaml
        mkdir -p sgc
        """

rule initialize_sgc:
    input:
        in1 = "spacegraphcats/setup.py",
        in2 = "envs/sgc.yaml",
    output:
        touch("sgc/sgc-initialized.txt"),
    conda: "envs/sgc.yaml"
    shell:
        """
        pip install spacegraphcats/
        """

rule desulfo_genomic:
    output:
        "sgc/GCA_001508995.1_ASM150899v1_genomic.fna.gz",
    conda: "envs/dib-rotation.yaml"
    shell:
        """
        mkdir -p sgc
        wget -P sgc/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/508/995/GCA_001508995.1_ASM150899v1/GCA_001508995.1_ASM150899v1_genomic.fna.gz
        """

rule desulfo_config_sgc:
    output:
        "sgc/desulfo-config.yml",
    params:
        id = expand("{sample}", sample=SAMPLES)
    shell:
        """
        mkdir -p sgc
        echo -e "catlas_base: '{params.id}'\ninput_sequences:\n- ../khmer/{params.id}.abundtrim.fq.gz\nksize: 31\nradius: 1\nsearch:\n- GCA_001508995.1_ASM150899v1_genomic.fna.gz\nsearchquick: GCA_001508995.1_ASM150899v1_genomic.fna.gz" > sgc/desulfo-config.yml
        """

rule desulfo_cDBG:
    input:
        in1 = "sgc/desulfo-config.yml",
        in2 = expand("khmer/{sample}.abundtrim.fq.gz", sample=SAMPLES),
        in3 = "sgc/GCA_001508995.1_ASM150899v1_genomic.fna.gz",
        in4 = "sgc/sgc-initialized.txt",
    output:
        out1 = directory(expand("sgc/{sample}{suffix}", sample=SAMPLES, suffix=["","_k31_r1","_k31_r1_search_oh0"])),
        out2 = touch("sgc/sgc-done.txt"),
    conda: "envs/sgc.yaml"
    shell:
        """
        cd sgc
        python -m spacegraphcats \
        run desulfo-config.yml extract_contigs \
        extract_reads --nolock 
        cd ..
        """