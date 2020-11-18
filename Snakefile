
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
        "reports/report.html"
        
# Workflow ------------------------------------------------------------------------------------------------------------------

# Insert Selective retrieval of FASTA by taxid groups and dereplication
# Output a dereplication table looking like:
#     taxid sciname centroid_acc clustered_acc
# Could work in a loop to avoid generating 1000s of files.
# e.g. blasdbcmd -taxid -oufmt %f | vsearch merge_exact >> fasta

rule export_fasta:
    # FIXME: Export %a%T%S once and use join to collect infos
    output:
        seq = temp("fasta/sequences.fa"),
        taxids = temp("fasta/taxids.txt")
    params:
        blast_DB = config["blast_db"],
        taxdb = config["taxdb"]
    message: "Exporting fasta file"
    conda:
        "envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_DB} -entry all -outfmt '%f' > {output.seq}
        blastdbcmd -db {params.blast_DB} -entry all -outfmt '%T' | sort -u > {output.taxids}
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
    # FIXME: Export %a%T%S once and use join to collect infos
    input:
        "pairwise_alignment/table.tsv"
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
        cat {input} | tail -n+2 | cut -d$'\t' -f1 | cut -d'.' -f1 > {output.queries}
        cat {input} | tail -n+2 | cut -d$'\t' -f2 > {output.targets}

        # calculate difference
        awk -F$'\t' '{{print $4+$5}}' {input} | tail -n+2 > {output.dist}

        # Merge files
        paste <(blastdbcmd -entry_batch {output.queries} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(blastdbcmd -entry_batch {output.targets} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(cat {output.dist}) \
              > {output.report}
        # Add reverse orientation - Output from Vsearch is non-redundant
        paste <(blastdbcmd -entry_batch {output.targets} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(blastdbcmd -entry_batch {output.queries} -db {params.blast_DB} -outfmt '%a\t%T\t%S') \
              <(cat {output.dist}) \
              >> {output.report}

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
        seq = "fasta/sequences.fa",
        taxids = "fasta/taxids.txt"
    output:
        taxids = temp("reports/taxids_number.txt"),
        seq = temp("reports/seq_number.txt")
    shell:
        """
        wc -l {input.taxids} | cut -d$' ' -f1 > {output.taxids}
        grep -c "^>" {input.seq} > {output.seq}
        """

rule write_report:
    input:
        seq = "reports/seq_number.txt",
        taxids = "reports/taxids_number.txt",
        dist = "reports/distances.tsv",
        sizedist = "reports/sequence_lengths.txt"
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