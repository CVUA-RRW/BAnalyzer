
shell.executable("bash")

# Settings ------------------------------------------------------------------------------------------------------------------

workdir: config["workdir"]

# Input rule ----------------------------------------------------------------------------------------------------------------
 
rule all:
    input: 
        "pairwise_alignment/table.tsv",
        "pairwise_alignment/alignments.txt",
        "reports/distances.tsv",
        "reports/sequence_lengths.txt",
        "reports/dereplication.tsv",
        "reports/report.html"
        
# Workflow ------------------------------------------------------------------------------------------------------------------

# Add derep report to html

rule export_dbinfo:
    output:
        temp("db_info.txt")
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Exporting database information"
    conda:
        "envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        blastdbcmd -db {params.blast_DB} -entry all -outfmt '%a\t%T\t%S' > {output}
        """
        
rule export_sequences:
    input:
        "db_info.txt"
    output:
        temp(dynamic("fastadump/{taxid}.fa"))
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Exporting sequences"
    conda:
        "envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        for tax in $(cut -d$'\t' -f2 {input} | sort -u); do
            # retrieving sequences per taxa
            blastdbcmd -db {params.blast_DB} -taxids $tax -outfmt '%f' \
                > fastadump/$tax.fa
        done
        """
        
rule dereplicate:
    input: 
        dynamic("fastadump/{taxid}.fa")
    output:
        report = temp("derepdump/dereplication.tsv"),
        fasta = temp("fasta/sequences_derep.fa"),
        tmpfa  = temp("derepdump/derep.fa"),
        tmptab = temp("derepdump/derep.txt")
    message: "Dereplicating sequences"
    conda:
        "envs/vsearch.yaml"
    shell:
        """
        for file in {input}; do
            if [ $(grep -c '^>' $file) -eq 1 ]; then
                cat $file >> {output.fasta}
            
            else
                vsearch --cluster_fast $file \
                        --id 1 \
                        --iddef 1 \
                        --centroids {output.tmpfa} \
                        --uc {output.tmptab} \
                        --quiet
                
                cat {output.tmpfa} >> {output.fasta}
                
                grep -E '^[S|H]' {output.tmptab} \
                     | cut -d$'\t' -f1,9,10 \
                     >> {output.report}
            fi
        done
        """

rule derep_stats:
    input:
        derep = "derepdump/dereplication.tsv",
        table = "db_info.txt"
    output:
        "reports/dereplication.tsv"
    message: "Collecting dereplication stats"
    shell:
        """
        join --nocheck-order -1 2 -2 1 -t $'\t' \
             <(sort -k2 {input.derep}) \
             <(sort -k1 {input.table}) \
             | sed -e 's/\tH\t/\thit\t/' -e 's/\tS\t/\tcentroid\t/' \
             | sort -k4n \
             > {output}
        
        sed -i '1 i\seqid\ttype\tcentroid\ttaxid\tname' {output}
    """

rule filter_seq:
    input:
        "fasta/sequences_derep.fa"
    output:
        temp("fasta/sequences.fa")
    message: "Filtering ambiguous sequences"
    params:
        max_n = config["max_n"]
    conda:
        "envs/cutadapt.yaml"
    shell:
        """
        cutadapt --max-n {params.max_n} {input} > {output}
        """


rule trim_primers:
    input:
        "fasta/sequences.fa"
    output:
        fasta = temp("fasta/sequences_trim.fa"),
        report = "reports/trimming_report.txt",
        primers_rc = temp("primers_rc.fa")
    params:
        primers = config["primers"]
    message: "Trimming primers"
    conda: 
        "envs/cutadapt.yaml"
    shell:
        """
        seqtk seq -r {params.primers} > {output.primers_rc}
        
        cutadapt -g file:{params.primers} \
                 {input} 2> {output.report} \
            | cutadapt -a file:{output.primers_rc} - \
                       > {output.fasta} 2>> {output.report}
        """

rule pariwise_alignement:
    input:
        "fasta/sequences_trim.fa" if config["trim_primers"] == True else "fasta/sequences.fa"
    output:
        table = "pairwise_alignment/table.tsv",
        align = "pairwise_alignment/alignments.txt"
    params:
        id = config["min_identity"]
    threads: workflow.cores
    message: "Performing pairwise global alignment"
    conda:
        "envs/vsearch.yaml"
    shell:
        """
        vsearch --allpairs_global {input} \
                --iddef 1 \
                --gapext 2I/2E \
                --gapopen 20I/20E \
                --threads {threads} \
                --id {params.id} \
                --userout {output.table} \
                --userfields query+target+id+mism+gaps+alnlen+qlo+qhi+ql+tlo+thi+tl \
                --alnout {output.align}
                
        sed -i '1 i\query\ttarget\tid\tmismatch\tgaps\taln_length\tqstart\tqend\tqlength\ttstart\ttend\ttlength' {output.align}
        """

rule collect_descriptors:
    input:
        aln = "pairwise_alignment/table.tsv",
        info = "db_info.txt"
    output:
        report = "reports/distances.tsv",
        queries = temp("queries.txt"),
        targets = temp("targets.txt"),
        dist = temp("distances.txt")
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Collecting sequence informations"
    conda: 
        "envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        # Separate query and target lists
        cat {input.aln} | tail -n+2 | cut -d$'\t' -f1 > {output.queries}
        cat {input.aln} | tail -n+2 | cut -d$'\t' -f2 > {output.targets}

        # calculate difference
        awk -F$'\t' '{{print $4+$5}}' {input.aln} | tail -n+2 > {output.dist}

        # Merge files
        paste <(blastdbcmd -entry_batch {output.queries} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(blastdbcmd -entry_batch {output.targets} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(cat {output.dist}) \
              > {output.report}
        # Add reverse orientation - Output from Vsearch is non-redundant TABLE TOO BIG
        # paste <(blastdbcmd -entry_batch {output.targets} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              # <(blastdbcmd -entry_batch {output.queries} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              # <(cat {output.dist}) \
              # >> {output.report}

        sed -i '1 i\query\tquery_taxid\tquery_name\ttarget\ttarget_taxid\ttarget_name\tdistance' {output.report}
        """

rule seq_sizes_raw:
    # FIXME: Export %a%T%S once and use join to collect infos
    input:
        raw = "fasta/sequences.fa",
        trimmed = "fasta/sequences_trim.fa" if config["trim_primers"] == True else "fasta/sequences.fa"
    output:
        raw = temp("reports/fasta_length_raw.tsv"),
        ids = temp("reports/seqids.txt"),
        length_raw = temp("reports/length_raw.txt"),
        length_trim = temp("reports/length_trim.txt"),
        table = "reports/sequence_lengths.txt"
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"],
        trim = config["trim_primers"]
    message: "Getting sequence length distribution"
    conda:
        "envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        # Get seqid length table
        cat {input.raw} \
            | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
            | sort -k1 \
            > {output.raw}
        
        # split seqids for blast input
        cat {output.raw} | cut -d$'\t' -f1 | cut -d'.' -f1 > {output.ids}
        cat {output.raw} | cut -d$'\t' -f2 > {output.length_raw}
        
        if [ {params.trim} = "True" ]; then
            # get trimmed lengths
            cat {input.trimmed} \
                | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
                | sort -k1 \
                | cut -d$'\t' -f2 \
                > {output.length_trim}
            
            # merge infos
            paste <(blastdbcmd -entry_batch {output.ids} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
                  <(cat {output.length_raw}) \
                  <(cat {output.length_trim}) \
                  > {output.table}
            # Add header
            sed -i '1 i\seqid\ttaxid\tname\tdb_length\ttrim_length' {output.table}
        
        else
            # if not trimming, just merge infos
            touch {output.length_trim} # otherwise snakemake complains about missing output
            paste <(blastdbcmd -entry_batch {output.ids} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
                  <(cat {output.length_raw}) \
                  > {output.table}
            # Add header
            sed -i '1 i\seqid\ttaxid\tname\tlength' {output.table}
        fi
        """

rule db_stats:
    input:
        unfiltered = "fasta/sequences_derep.fa",
        seq = "fasta/sequences.fa",
        table = "db_info.txt",
    output:
        taxids = temp("reports/taxids_number.txt"),
        nseq = temp("reports/seq_number.txt"),
        derep = temp("reports/derep_number.txt"),
        highN = temp("reports/high_N.txt")
    shell:
        """
        cut -d$'\t' -f2 {input.table} | sort -u | wc -l | cut -d$' ' -f1 > {output.taxids}
        wc -l {input.table} | cut -d$' ' -f1 > {output.nseq}
        grep -c "^>" {input.seq} > {output.derep}
        
        echo $(( $(grep -c "^>" {input.unfiltered}) - $(grep -c "^>" {input.seq}) )) > {output.highN}
        """

rule write_report:
    input:
        seq = "reports/seq_number.txt",
        taxids = "reports/taxids_number.txt",
        dist = "reports/distances.tsv",
        sizedist = "reports/sequence_lengths.txt",
        derep = "reports/dereplication.tsv",
        nderep = "reports/derep_number.txt",
        nNfilt = "reports/high_N.txt"
    output:
        "reports/report.html"
    params:
        database = config["blast_db"],
        workdir = config["workdir"],
        trimming = config["trim_primers"],
        id = config["min_identity"],
        primers = config["primers"] if config["trim_primers"] == True else config["trim_primers"]
    message: "Writting report"
    conda:
        "envs/rmarkdown.yaml"
    script:
        "scripts/write_report.Rmd"