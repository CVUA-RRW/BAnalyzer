import pandas as pd

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
        "reports/cluster_size.tsv",
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
    
rule cluster_size:
    input:
        "reports/dereplication.tsv"
    output: 
        "reports/cluster_size.tsv"
    message: "Collecting cluster sizes"
    run:
        df = pd.read_csv(input[0], sep= '\t')
        taxcount =  df.groupby('taxid').agg('count')['seqid'].rename('tax_size')
        centroids = df[df['type'] == 'centroid']
        centroidSize = df.groupby('centroid').nunique()['seqid']+1
        centroidSize = centroidSize[centroidSize.index != '*'].rename("cluster_size")
        dfout = centroids.set_index('seqid').join(centroidSize).fillna(1)
        dfout = dfout.join(taxcount, on = 'taxid').drop(["type", "centroid"], axis = 1)
        dfout = dfout.astype({'cluster_size' : 'int32', 'tax_size': 'int32'})
        dfout["size"] = dfout["cluster_size"].map(str) + "/" + dfout["tax_size"].map(str)
        dfout["rel_cluster_size"] = round(dfout["cluster_size"] / dfout["tax_size"] *100, 2)
        dfout.to_csv(output[0], sep = '\t')
        

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
        info = "db_info.txt",
        sizes = "reports/cluster_size.tsv"
    output:
        "reports/distances.tsv",
    message: "Collecting sequence informations"
    run:
        dfaln = pd.read_csv(input.aln, sep='\t', names= ["query", "target", "id", "mismatch", "gaps","aln_length", "qstart", "qend", "qlength", "tstart", "tend", "tlength"])
        dfinfo = pd.read_csv(input.info, sep = '\t', names = ["seqid", "taxid", "name"])
        dfsize = pd.read_csv(input.sizes, sep = '\t')

        dfout = pd.DataFrame({"query": dfaln["query"],
                              "target": dfaln["target"],
                              "distance": dfaln["mismatch"] + dfaln["gaps"]})

        dfout = dfout.join(dfinfo.set_index("seqid"), 
                            on = "query", how = "inner").rename(columns={'taxid' : 'query_taxid', 
                                                                         'name' : 'query_name'})
                                                                         
        dfout = dfout.join(dfinfo.set_index("seqid"), 
                            on = "target", how = "inner").rename(columns={'taxid' : 'target_taxid', 
                                                                         'name' : 'target_name'})
                                                                         
        dfout = dfout.join(dfsize.set_index("seqid")[['size', 'rel_cluster_size']], 
                            on = "query", how = "inner").rename(columns={'size' : 'query_size', 
                                                                         'rel_cluster_size' : 'query_relsize'})
                                                                         
        dfout = dfout.join(dfsize.set_index("seqid")[['size', 'rel_cluster_size']], 
                            on = "target", how = "inner").rename(columns={'size' : 'target_size', 
                                                                         'rel_cluster_size' : 'target_relsize'})
                                                                         
        dfout = dfout[['query', 'query_taxid', 'query_name', 'query_size', 'query_relsize',
                        'target', 'target_taxid', 'target_name', 'target_size', 'target_relsize',
                        'distance']]
                        
        dfout.to_csv(output[0], sep='\t', index=False)
        
rule get_seq_sizes:
    input:
        raw = "fasta/sequences.fa",
        trimmed = "fasta/sequences_trim.fa" if config["trim_primers"] == True else "fasta/sequences.fa"
    output:
        raw = temp("reports/fasta_length_raw.tsv"),
        length_trim = temp("reports/fasta_length_trim.txt"),
    params:
        trim = config["trim_primers"]
    message: "Getting sequence length distribution"
    shell:
        """
        # Get seqid length table
        cat {input.raw} \
            | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
            > {output.raw}
        
        if [ {params.trim} = "True" ]; then
            # get trimmed lengths
            cat {input.trimmed} \
                | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
                > {output.length_trim}
            
        else
            # if not trimming, just merge infos
            touch {output.length_trim} # otherwise snakemake complains about missing output
        
        fi
        """

rule seq_size_table:
    input:
        raw = "reports/fasta_length_raw.tsv",
        trim = "reports/fasta_length_trim.txt",
        info = "db_info.txt"
    output:
        "reports/sequence_lengths.txt"
    params:
        trim = config["trim_primers"]
    message: "Formatting sequence length table"
    run:
        dfinfo = pd.read_csv(input.info, sep = '\t', names = ["seqid", "taxid", "name"])
        dfraw = pd.read_csv(input.raw, sep = '\t', names = ["seqid", "length"])
        
        dfoutraw = dfraw.join(dfinfo.set_index("seqid"),
                              on = "seqid", how = "inner")
        
        if params.trim:
            dfout = dfoutraw.rename(columns={'length': 'db_length'})
            dftrim = pd.read_csv(input.trim, sep = '\t', names = ["seqid", "length"])
            dfout = dfout.join(dftrim.set_index('seqid'),
                               on = 'seqid', how = 'inner').rename(columns = {'length' : 'trim_length'})
            dfout = dfout[['seqid', 'taxid', 'name', 'db_length', 'trim_length']]
            dfout.to_csv(output[0], sep = '\t', index = False)
        
        else:
            dfoutraw = dfoutraw[['seqid', 'taxid', 'name', 'length']]
            dfoutraw.to_csv(output[0], sep = '\t', index = False)

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
        
rule get_consensus_level:
    input:
        distance_table = "reports/distances.tsv",
    output:
        cons = "reports/consensus.tsv"
    message:
        "Determining consensus ranks"
    params:
        lineage = config["taxonomy"]["rankedlineage_dmp"],
        nodes = config["taxonomy"]["nodes_dmp"],
    script:
        "../scripts/consensus_levels.py"

rule write_report:
    input:
        seq = "reports/seq_number.txt",
        taxids = "reports/taxids_number.txt",
        dist = "reports/distances.tsv",
        sizedist = "reports/sequence_lengths.txt",
        derep = "reports/dereplication.tsv",
        nderep = "reports/derep_number.txt",
        nNfilt = "reports/high_N.txt",
        clusterSize = "reports/cluster_size.tsv"
        consensus = "reports/consensus.tsv"
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