#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from taxidTools.taxidTools import Taxdump

def get_consensus_rank(txd, idlist):
	lca = txd.lowestCommonNode([str(t) for t in idlist])
	return(txd.getRank(lca))

def main(table, outfile, lineage, nodes):
	df = pd.read_csv(table, sep = '\t')
	
	df["entries"] = df["query"] + df["query_size"]
	
	dfout = pd.DataFrame()
	txd = Taxdump(lineage, nodes)
	
	for entry in set(df["entries"]):
		sub = df[df["entries"] == entry].reindex()
		
		new = sub.head(1)[["query_name", "query_taxid", "query_size", "query_relsize"]]
		
		for i in range(0,4):
			taxids = list(sub[sub["distance"] <= i]["target_taxid"])
			taxids.extend(list(sub.head(1)["query_taxid"]))
			rank = get_consensus_rank(txd, taxids)
			
			new[f"Consensus rank with {i} mismatches"] = rank
		
		dfout = dfout.append(new, ignore_index=True)
		
		dfout.to_csv(outfile, sep = '\t')

if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0], snakemake.params["lineage"], snakemake.params["nodes"])