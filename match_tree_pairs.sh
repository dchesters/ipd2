

 ####################################################################
 #
 # 
 # 20260906: complete replication of the tree processing and pairing
 # 
 #
 # Scripts: 
 #	process_catalogue_of_life.pl
 #	convert_scientific_notation_branchlengths.pl
 # 	process_newicks.pl	
 #	taxonomic_report.pl
 # 	relational_constraints.pl
 #	taxon_table.pl
 #	process_sourcetree_classes.pl
 # 	wrap_pair_terminals.pl
 # 	pair_terminals.pl
 # 	prune_tree_pair.pl
 #	prune_tree.pl
 #	integrate_PD_results.pl
 #
 # Files: 
 # 	current_list.202511B
 #	dataset-316115.txtree
 #	Newicks.202605.tar.xz
 #	newicks.sample_table.202511
 #	sourcetrees_table.202511B
 #
 # Software/libraries:
 # 	R-base
 #	TreeDist
 # 	ape
 # 
 # 
 ####################################################################








 #######################################################################################################
 # Catalog of Life is the single most comprehensive taxonomic database for the species level,
 # use this for correcting synonyms in phylogeny data,
 # first parse the database and make simple tabular format.
 # https://www.catalogueoflife.org/data/download
 # Base Release
CoL_database=dataset-316115.txtree
perl process_catalogue_of_life.pl $CoL_database insect_species_and_synonyms



 #######################################################################################################
 # Download the insect phylogeny database,
 # correct synonyms and remove branchlengths if present.
 # https://github.com/dchesters/insect_phylogeny_database
rm *.nwk
newick_list=current_list.202511B
synonyms=insect_species_and_synonyms
perl process_newicks.pl $newick_list $synonyms


 #######################################################################################################
 # Higher taxon names need to be inferred for species at terminals,
 # NCBI taxonomic database is more suitable for this task 
 # due to more detailed and well structured lineage information.
 # Download NCBI taxonomy (taxdump.tar.gz) from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
 # Unzip the object, there are several files in there, 2 are needed, 
 # copy names.dmp into the working directory
 # https://www.itis.gov/download.html
 # ITIS_database=/home/douglas/databases/ITIS/taxa_3842024.txt
 # -ITIS $ITIS_database 
NCBI_names=names.dmp
NCBI_nodes=nodes.dmp
taxon_of_interest=Hexapoda
rm Taxon_Table taxon_table.lineages
perl taxon_table.pl -node 6960 $taxon_of_interest -NCBI $NCBI_names $NCBI_nodes -outfile Taxon_Table.$taxon_of_interest


 ########################################################################################################
 # Assign higher taxon names to terminals for all newicks
 # This is one of the slowest steps of the pipeline, 1-2 days to cover all newicks.
 # includes a random selection of exemplars so PD results wont be exactly the same at later steps,
readarray -t newick_list < current_list.202511B
taxon_table=Taxon_Table.Hexapoda
for i in ${!newick_list[*]}
do
intree=${newick_list[$i]}
if test -f constraints.$intree; then
 echo "lineage_assignments.$intree file found."
else
 echo "lineage_assignments.$intree file not found."

# infer higher taxa and lineage for each terminal
# first, idiosyncrasy of a subsequent script, need a fasta file of taxon of interest
rm $intree.fas
perl taxonomic_report.pl -input $intree -output $intree.taxreport -newick
rm $intree.taxreport* taxonomic_report_verboseLOG
rm Backbone_Constraints_Newick.list_of_constraints lineage_assignments_to_barcodes lineage_assignments.$intree Backbone_Constraints_Newick.list_of_constraints constraints.$intree
perl relational_constraints.pl -seqfile $intree.fas -treefile $intree -outfile_prefix RelConsTEMP -backbone_terminal_format 0 -taxon_table $taxon_table
mv lineage_assignments_to_barcodes lineage_assignments.$intree;mv Backbone_Constraints_Newick.list_of_constraints constraints.$intree
rm RelConsTEMP.* tabulated_taxa all_taxon_assignments_forBAMM.txt backbone_constraints_newick_LOG constraint_phylogeny_pruned2.*
rm $intree.fas list_constrained_IDs MRMS_clade_members RelCons_VerboseLog.txt SRC_fasta SRC_tabulated

fi
done


 ########################################################################################################
 # Find candidate tree pairs, those that have taxonomic overlap
 # $job = 6;
sourcetrees_table=sourcetrees_table.202511B
rm match_newick_pairs_LOG match_newick_pairs_LOG2
perl process_sourcetree_classes.pl $sourcetrees_table


 ########################################################################################################
 # Candidate newick pair is modified as to try and increase number of terminal matches, by assignment of higher taxa.
 # Wrapper for pair_terminals.pl, which tests pair of phylogenies with overlapping taxa. Arguments are treeA treeB pairID
 # pair_terminals.pl invoked only if lineage_assignments file is found for both of newick pair.
 # Outputs are 2 newicks, modified with some higher taxon names where inexact matches have been made between tree pair,
rm *.nwk.pairID*
perl wrap_pair_terminals.pl match_newick_pairs_LOG2 NULL
wc -l treepair_match_LOG.txt


 ########################################################################################################
 # Read the log file, has one tree pair per line, need to prune the non overlapping terminal
 # output is All_prune_commands
rm All_prune_commands
perl prune_tree_pair.pl treepair_match_LOG.txt 1

# these might have been deleted if running previous step, just in case:
rm *pairID*.pruned
chmod +x All_prune_commands
 # branchlengths removed here, if present. nb sci-not bls would crash later steps anyway
 # running time 1 min:
./All_prune_commands

 # Script job 2 read log file again, input each pruned pair into r for calculation of phylogenetic distance
 # takes couple of hours
 # output has columns newick1 newick2 phylo_dist
rm disttopo_RESULTS disttopo_RESULTS2
perl prune_tree_pair.pl treepair_match_LOG.txt 2


 ########################################################################################################
 # Make complete R table
rm integrated_result list_all_paired_newicks
perl integrate_PD_results.pl disttopo_RESULTS2 treepair_match_LOG.txt newicks.sample_table.202511
wc -l integrated_result




#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



