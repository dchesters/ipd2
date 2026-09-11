

# 
# study specific script used for work in submission, provided for review
# 



$sourcetrees_table = $ARGV[0];


# this not used under default settings:
$newicks_path = "/home/douglas/databases/Phylo/Processed_Newick/";


# [ more can be added here ] scrap that, only one can be used due to possibility of misalignment
@COI_files = 	( "COI_NCBIandBOLD"
	#	"/home/douglas/scripted_analyses/Phylogenomic_Contribution/branchlengths/alignerOUT_hq.Insecta202309"
	#	"/home/douglas/scripted_analyses/Phylogenomic_Contribution/branchlengths/alignerOUT.sp"
	#	"/home/douglas/scripted_analyses/BEF_China_and_TreeDi/Lepidoptera/barcodes/alignerOUT.20230714.hq"
		);


$copy_newicks = 0; # default 0


		# 1 == read sourcetrees info table, get phylo class for each, and copy relevent newicks
$job = 6; 	# 2 == find terminal species which dont have coi from bold
		# 3 == get sequences and infer branchlengths for each newick
		# 4 == get sum branchlength for each sourcetree
		# 5 == look for tree pairs with potential overlapping terminals.
		# 6 == different approach for finding tree pairs based only on exact matches

			# default 3, higher values very unlikely to match additional tree pairs.
$max_rank_steps = 3;	# in checking overlap between tree pair, rank of focal taxon will be varied by this number of steps,
			# larger number of steps means less incidence of overlap.

$overlap_cutoff = 5;

#######################################################################################################################



#######################################################################################################################

unless($sourcetrees_table =~ /\w/){die "\nerror 35, no path given in command\n"};

if($job =~ /[23]/)
	{
	print "\njob == $job\n";
	get_sequences_for_newicks();
	exit;
	}elsif($job =~ /[4]/)
	{
	get_sum_branchlengths();
	exit;
	};


#######################################################################################################################






#################################################################################
$line_count = 0;
open(IN, $sourcetrees_table) || die "\nerror 3 , cant open $sourcetrees_table\n";

if($job =~ /[123]/)
	{
	open(OUT, ">newick_filenames") || die "\nerror 10.\n";
	open(OUT2, ">newicks.sample_table") || die "\nerror 304.\n";
	};

while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	my @splitline = split /\t/, $line;
	if($line_count == 0)
		{		# HEADER line only
		for my $column_index(0 .. $#splitline)
			{
			my $column_header = $splitline[$column_index];print "column_header:$column_header\n";
			push @column_titles, $column_header;
			};

		}elsif($line =~ /\.nwk/) # not header, and has tree file
		{

		my $current_datatype; my $current_filename;my $current_Lineage; my $current_Focal_taxon;my $current_Focal_order = "NA";
		my $current_CharsMol = "NA";my $current_CharsMorph = "NA";my $current_ClassMol = "NA";my $current_ClassMorph = "NA";
		for my $column_index(0 .. $#splitline)
			{
			my $column_header = $column_titles[$column_index];
			my $cell_entry = $splitline[$column_index]; # print "line_count:$line_count column_index:$column_index column_header:$column_header cell_entry:$cell_entry\n";

			if($column_header eq "ID" && $cell_entry =~ /\S\.nwk/)
				{ # here is newick file name
				$sourcetree_path = "$newicks_path$cell_entry";$current_filename = $cell_entry;
				# if($copy_newicks == 1){ system("cp $sourcetree_path $cell_entry") };
				}elsif( $column_header eq "Datatype.L1" && $cell_entry =~ /\w/)
				{ # tree type
				$all_datatypes{$cell_entry}++;$current_datatype = $cell_entry;
				unless($cell_entry =~ /Mol|Morph|Synth/){print "warning unexpected entry $line\n"};
				}elsif( $column_header eq "Class.Mol" && $cell_entry =~ /\w/)
				{
				$current_ClassMol=$cell_entry;
				}elsif( $column_header eq "Class.Morph" && $cell_entry =~ /\w/)
				{
				$current_ClassMorph=$cell_entry;
				}elsif( $column_header eq "Datatype.L2")
				{
				print "\nwarning, cant parse Datatype.L2 from line:$line\n";
				};


			if($column_header eq "Lineage")				{$current_Lineage	=$cell_entry};
			if($column_header eq "Focal_taxon")			{$current_Focal_taxon	=$cell_entry};
			if($column_header eq "Focal_Order")			{$current_Focal_order	=$cell_entry};
			if($column_header eq "Loci" && $cell_entry =~ /\d/)	{$current_CharsMol	=$cell_entry}; 		# Chars.Mol is not recorded much
			if($column_header eq "Chars.Morph"&&$cell_entry=~/\d/)	{$current_CharsMorph	=$cell_entry};

			};

		# print "\nline:$line\n\tcurrent_CharsMol:$current_CharsMol\n";die "";
# 20250524: will now class at two levels, first as before SG, MG, Mito, Omics, Morph, Comb
#		also molecular class SG, MG, Mito, transcriptome, UCE, AHE, Omics(....)

		if($current_datatype =~ /\w/ && $current_filename =~ /nwk/)
			{
			my $class="NA";my $classMol="NA";my $classMorph; # print "$current_datatype\t";

			######################################################################################################			
			# broad class
			if($current_datatype =~ /Mol[\,\/]Morph/)# mol and morph; combined:C
				{
				$class = "C"
				}elsif($current_datatype =~ /^SYNTH/i)# synthesis:S
				{
				$class = "S"
				}elsif($current_datatype =~ /^MORPH/i)
				{
				$class = "P"
				}else{
				if($current_ClassMol =~ /^AHE|^OME|^RADseq|^RAD\-Seq|^TR|^UCE|^TE|^NGS|^BUSCO|^USCO|^WGS|^LCG/)	# omics:O
					{
					$class = "O"
					}elsif($current_ClassMol =~ /^PM|^MG/) 				# multigene:G
					{
					$class = "G"
					}elsif($current_ClassMol =~ /^MT/) 					# mitogenome:M
					{ 
					$class = "M"
					}elsif($current_ClassMol =~ /^SG|^BC/)				# barcode:B
					{
					$class = "B"
					}else{
					print "warning, no class for $line\n"	
					};
			#	if($current_ClassMorph =~ /^IM|^AD|^LA|^PU|^LV|^MORPH|^Morph/) 	# morphology:P
			#		{
			#		$class = "P"
			#		};
				};
			######################################################################################################			
			# molecular subclass [edit, not interested in diff between e.g. gene trees and omics trees]

				if($current_ClassMol =~ /^AHE/) 
					{
					$classMol = "AHE"
					}elsif($current_ClassMol =~ /^UCE/)
					{
					$classMol = "UCE"
					}elsif($current_ClassMol =~ /TR/)
					{
					$classMol = "TR"
					}elsif($current_ClassMol =~ /^OME|^RADseq|^TE|^NGS|^BUSCO|^USCO/)
					{
					$classMol = "NA"; # $classMol = "OME"
					}elsif($current_ClassMol =~ /^PM|^MG/) 				
					{
					$classMol = "NA"; # $classMol = "MG"
					}elsif($current_ClassMol =~ /^MT/) 			
					{ 
					$classMol = "NA"; # $classMol = "MITO"
					}elsif($current_ClassMol =~ /^SG|^BC/)	
					{
					$classMol = "NA"; # $classMol = "SG"
					};

			######################################################################################################			


			if( $class =~ /^[CSPOGMB]$/ && $current_filename =~ /nwk/)
				{
				$class_counts{$class}++;
				if($copy_newicks == 1){system("cp $sourcetree_path $current_filename")}; # $class.$current_filename 
				print OUT "$current_filename\n"; 		# newick_filenames
				print OUT2 "$current_filename\t$class\t$current_Focal_order\t$current_CharsMol\t$current_CharsMorph\t$classMol\n";# newicks.sample_table
				}else{
				print "WARNING, sourcetree overlooked:$current_filename\n";
				};


			###################################################################################

			$Lineages{$current_filename} = $current_Lineage;$Focal_taxas{$current_filename} = $current_Focal_taxon;

			###################################################################################

			}elsif($line =~ /processed|manual|Mol/)
			{
			print "\nwarning, line not parsed:$line\n";		
			};
			
		};

	$line_count++;
	};
close IN;

if($job =~ /[123]/)
	{
	print "\nwritten sourcetree classes to newicks.sample_table\n";
close OUT;
close OUT2;
	};



#################################################################################


my @datatypes = keys %all_datatypes;@datatypes = sort @datatypes;
print "\ndata types found:\n";
foreach my $datatype(@datatypes){print "\t$datatype\t$all_datatypes{$datatype}\n"};

my @classes = keys %class_counts;@classes = sort @classes;
print "\nclassed to:\n";
foreach my $datatype(@classes)
	{
	print "\t$datatype\t$class_counts{$datatype}\n";$newicks_classed+=$class_counts{$datatype};
	};
print "newicks_classed:$newicks_classed\n";


#################################################################################


if($job =~ /[6]/)
{
print "\nsearching for sourcetree pairs based on exact terminal matches alone.\n";
my @parsed_treeIDs = keys %Lineages;@parsed_treeIDs = sort @parsed_treeIDs;
print "count parsed_treeIDs:$#parsed_treeIDs\n";
my %all_species;my $count_unique_sourcetrees=0;
#####################################################
foreach my $treeID(@parsed_treeIDs)
	{
	my $sourcetree_nonunique=0;
	my %nwk1_terminals = ();
	open(NWK1, $treeID) || die "\nerror 203, cant find $treeID.\n";
	while(my $line = <NWK1>)
		{
		while($line =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/)
			{
			my $terminal=$2; # $nwk1_terminals{$terminal}=1
			if($nwk1_terminals{$terminal}==1){print "Warning, $treeID has duplicate terminal $terminal\n"};
			$nwk1_terminals{$terminal}=1;
			$all_species{$terminal}+=1;
		#	if($all_species{$terminal} >= 2){$sourcetree_nonunique=1};
			};
		};
	close NKW1;
#	if($sourcetree_nonunique==0){$count_unique_sourcetrees++};
	};
#####################################################
# for each species, have summed how many sourcetrees they are in,
# now go through again, and any sourcetrees containing all species that are in 1 sourcetree only,
# can be discounted from further consideration.
my @candidate_pairable=();
foreach my $treeID(@parsed_treeIDs)
	{
	my $sourcetree_nonunique=0;
#	my %nwk1_terminals = ();
	open(NWK1, $treeID) || die "\nerror 203, cant find $treeID.\n";
	my $current_terminals = "";
	while(my $line = <NWK1>)
		{
		while($line =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/)
			{
			my $terminal=$2;$current_terminals .="$terminal\t"; # $nwk1_terminals{$terminal}=1
		#	if($nwk1_terminals{$terminal}==1){print "Warning, $treeID has duplicate terminal $terminal\n"};
		#	$nwk1_terminals{$terminal}=1;
		#	$all_species{$terminal}+=1;
			if($all_species{$terminal} >= 2){$sourcetree_nonunique=1};# tree has as least one terminal found in another tree
			};
		};
	close NKW1;
	if($sourcetree_nonunique==0)
		{}else{
		push @candidate_pairable, $treeID;$stored_terminals{$treeID}=$current_terminals;
		};
	};
print "candidate_pairable:$#candidate_pairable\n";
#####################################################
# go through tree pairs and look for exact matches
my $tree_index=0;
@sourcetree_pairs;
foreach my $treeID(@candidate_pairable)
	{
	if($tree_index =~ /00$/){print "tree_index $tree_index of $#candidate_pairable, count_overlapping_sourcetree_pairs:$count_overlapping_sourcetree_pairs\n"};
	my %nwk1_terminals = ();
	open(NWK1, $treeID) || die "\nerror 203, cant find $treeID.\n";
	while(my $line = <NWK1>)
		{
		while($line =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/)
			{
			my $terminal=$2;  $nwk1_terminals{$terminal}=1
			};
		};
	close NKW1;

	foreach my $treeID2(@candidate_pairable)
		{
		my @pairID_array = ($treeID,$treeID2);@pairID_array=sort@pairID_array;
		my $pairID = join '	',@pairID_array;# print "pairID:$pairID\n";
		# if($pairs_stored{$pairID}==1){print "\talready stored $pairID\n"}else{};
		unless($treeID eq $treeID2 || $pairs_stored{$pairID}==1) # 
			{
			my $overlap_terminals=0;
		#	open(NWK2, $treeID2) || die "\nerror 203, cant find $treeID2.\n";
		#	while(my $line2 = <NWK2>)
		#		{
		#		while($line2 =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/)
		#			{
		#			my $terminal=$2;  if($nwk1_terminals{$terminal}==1){$overlap_terminals++};
		#			};
		#		};
		#	close NKW2;
			my $terminal_set = $stored_terminals{$treeID2};$terminal_set =~ s/\t+$//;my @terminals_split=split /\t+/, $terminal_set;
			foreach my $terminal2(@terminals_split)
				{if($nwk1_terminals{$terminal2}==1){$overlap_terminals++};
				};
			if($overlap_terminals >= $overlap_cutoff) # 
				{
				$count_overlapping_sourcetree_pairs++;push @sourcetree_pairs, $pairID;$pairs_stored{$pairID}=1;
				};
			};
		};

	$tree_index++;
	};
print "count_overlapping_sourcetree_pairs:$count_overlapping_sourcetree_pairs\n";
open(MATCH_NEWICK_PAIRS_LOG2, ">match_newick_pairs_LOG2") || die "\nerror 281, cant open outfile \n";
foreach my $k(@sourcetree_pairs)
	{
	print MATCH_NEWICK_PAIRS_LOG2 "$k\n";
	};
close MATCH_NEWICK_PAIRS_LOG2;

#####################################################



exit;
};

#########################################################################################################


if($job =~ /[5]/)
{
print "\nLooking for tree pairs with overlapping terminals.\n";
@parsed_treeIDs = keys %Lineages;@parsed_treeIDs = sort @parsed_treeIDs;
foreach my $treeID(@parsed_treeIDs)
	{
	my $current_Lineage = $Lineages{$treeID};my $current_Focal_taxon = $Focal_taxas{$treeID}; # these originate in IPD table
#	print "treeID:$treeID Focal_taxon:$current_Focal_taxon Lineage:$current_Lineage\n";

	foreach my $treeID2(@parsed_treeIDs)  	# goes through all trees for all trees,
		{				# thus each pair assessed in both directions (focaltax-lineage, lineage-focaltax).

		my $current_Lineage2 = $Lineages{$treeID2};my $current_Focal_taxon2 = $Focal_taxas{$treeID2};
		# print "treeID:$treeID Focal_taxon:$current_Focal_taxon Lineage:$current_Lineage\n";
		unless($treeID eq $treeID2)
			{
			for my $steps(0..$max_rank_steps)
				{
				if($current_Lineage2 =~ s/\,(\w+)$//)
					{
					my $tax3 = $1;
					# if tax from lineage of treeID equals focal taxon of treeID1, store match
					if($tax3 eq $current_Focal_taxon)
						{$tree_matches{$treeID} .= "$treeID2:$steps,"};
					}
				};
			};
		};
	};

open(MATCH_NEWICK_PAIRS_LOG, ">match_newick_pairs_LOG") || die "\nerror 200, cant open outfile \n";
foreach my $treeID(@parsed_treeIDs)
	{

	if($tree_matches{$treeID} =~ /\w/)
		{
		$count_trees_paired++;
		my $current_pairs =$tree_matches{$treeID};$current_pairs =~ s/\,$//;my @array_pairs = split /\,/, $current_pairs;
		$sum_pairs += scalar @array_pairs;
	#	print "tree $treeID paired with $current_pairs\n";
		# tree Zhou2020Dacini.nwk paired with Valerio2021Bactrocera.nwk:1,Virgilio2009Dacus.nwk:1
		my %nwk1_terminals = ();
		open(NWK1, $treeID) || die "\nerror 203, cant find $treeID.\n";
		while(my $line = <NWK1>){while($line =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/){my $terminal=$2;$nwk1_terminals{$terminal}=1}};
		close NKW1;
	#	print "terminals read:",scalar keys %nwk1_terminals, "\n";

		# next open the newicks a check for exact matches in terminal IDs
		for my $steps(0..$max_rank_steps)
			{
			if($current_pairs =~ /\:$steps/)
				{
				while($current_pairs =~ s/([^\,\:\.]+\.nwk)\:$steps\,//)
					{
					my $nwk = $1; # print "\tsteps:$steps nwk:$nwk\n";
					my $exactmatch_count=0;
					open(NWK2, $nwk) || die "\nerror 218, cant find $nwk.\n";
					while(my $line = <NWK2>)
						{
						while($line =~ s/([\(\)\,])([A-Za-z_]+)([\:\,\)])/$1$3/)
							{
							my $terminal=$2;if($nwk1_terminals{$terminal}==1){$exactmatch_count++}
							}
						};
					close NKW2;
				#	print "\t\texactmatch_count:$exactmatch_count\n";
					if($exactmatch_count >= 4){print MATCH_NEWICK_PAIRS_LOG "$treeID\t$nwk\t$steps\t$exactmatch_count\n"};
					};
				};
			};

		}else{
		$count_trees_not_paired++;
		
		};
	$index3++;
	};

print "
JOB 5, searching for taxonomic overlap between tree pairs.
	max_rank_steps:$max_rank_steps
	count_trees_paired:$count_trees_paired
	count_trees_not_paired:$count_trees_not_paired
	sum_pairs:$sum_pairs (for each matched tree, how many others is it matched with)

printed match_newick_pairs_LOG

";

close MATCH_NEWICK_PAIRS_LOG;

};# if($job =~ /[5]/)

#################################################################################


	# 	SUBS








#######################################################################################################################

sub get_sum_branchlengths
{
print "\nsub get_sum_branchlengths\n";
open(IN, "newick_filenames") || die "\ncant open newick_filenames\n";
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /(\S+)/){my $filename = $1;push @newick_files, $filename};
	};
close IN;
print "\n read file newick_filenames, there are $#newick_files files\n";

system("rm sum_branchlengths");

for my $i(0 .. $#newick_files)
	{
	my $file = $newick_files[$i]; print "newick file $i of $#newick_files is $file\n";my $coi_not_available_for_terminal_current =0;
	####################################################################################
	# FOR EACH NEWICK FILE
	my $current_newick = "RAxML_bestTree.$file";
	# this script reads branchlengths and sums
	my $command = "perl /home/douglas/usr_scripts/newick/remove_branchlengths_from_newick.pl $current_newick";
	print "command:$command\n"; 
	 system("$command");

	####################################################################################
	};

open(RESULTS, "sum_branchlengths") || die "";
open(RESULTS2, ">sum_branchlengths2") || die "";
while(my $line = <RESULTS>)
	{
#	RAxML_bestTree.G.Letsch2016Anisoptera.nwk	43.2187509372781
	if($line =~ s/RAxML_bestTree\.(\w)\.(\S+)\t(\S+)/$2\t$1\t$3/)
		{
		$line =~ s/\tO\t/\tOmics\t/;
		$line =~ s/\tG\t/\tMultigene\t/;
		$line =~ s/\tM\t/\tMitogenome\t/;
		$line =~ s/\tP\t/\tMorphology\t/;
		$line =~ s/\tC\t/\tCombined\t/;
		print RESULTS2 "$line";
		};
	};
close RESULTS;
close RESULTS2;

}; # sub get_sum_branchlengths

#######################################################################################################################












#######################################################################################################################

sub get_sequences_for_newicks
{
print "\nsub get_sequences_for_newicks\n";
#######################################################

for my $index(0 .. $#COI_files) # Fasta reader, for COIs
	{
	my $COI_file = $COI_files[$index];
	#######################################################################
	open(FASTA_IN, $COI_file) || die "error 44. Cant open seq file $database_file.\n";print "seq file $index named $COI_file\n";
	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)	{$fasta_entry =~ s/^>//;process_entry($fasta_entry)};$fasta_entry = $fasta_line;
			}else{
			$fasta_entry .= $fasta_line;
			}
		};
	close(FASTA_IN);
	unless(length($fasta_entry)<=2)	{$fasta_entry =~ s/^>//; process_entry($fasta_entry)}
	sub process_entry
		{
		my $line = shift;my $current_id = "";
		if ($line =~ /^(.+)\n/ )
			{
			$current_id = $1;$line =~ s/^.+\n//;
			}else{die "\nerror 5631.\n"};
		$line =~ s/[\-\?n]/N/g;
		$sequences{$current_id} = $line; $sequences_stored++;# print ">$current_id\n";
		# also store exemplar for each genus
		if($current_id =~ /^(\w+)_/){my $genus = $1 . "_sp"; $sequences{$genus} = $line};
		}
	print "sequences_stored:$sequences_stored\n";
	######################################################################
	};
#######################################################
# read list of newick file names
open(IN, "newick_filenames") || die "\ncant open newick_filenames\n";
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /(\S+)/){my $filename = $1;push @newick_files, $filename};
	};
close IN;
print "\n read file newick_filenames, there are $#newick_files files\n";
#######################################################

if($job == 2)	{
		open(FIND_COIS, ">missings_COIs") || die "\nerror 83.\n";
		};

for my $i(0 .. $#newick_files)
	{
	my $file = $newick_files[$i]; print "newick file $i of $#newick_files is $file\n";my $coi_not_available_for_terminal_current =0;
	####################################################################################
	# FOR EACH NEWICK FILE

	# get list of terminals for current newick file; remove_fasta_entries with newick switch just reads terminals and lists
	system("rm $file.terminal_list");
	system("rm current_COIs current_COIs.phy current_Newick1 current_Newick2");

	my $command = "perl ~/usr_scripts/remove_fasta_entries.pl -newick $file"; system($command);print "command:$command\n";
	my @current_terminals = ();my %current_genera = ();
	open(LIST, "$file.terminal_list") || die "\nerror cant open file $file.terminal_list\n";
	while(my $line = <LIST>){$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /(.+)/)
		{
		my $ID = $1;# print "ID:$ID\n";
		if($ID =~ /^(\w+)_/){$current_genera{$1}++}; # this needs to go before push command, since latter disappears the variable
		push @current_terminals, $ID};
		};
	close LIST;
	# terminals are now read into hash %current_terminals

	open(TREE, "$file") || die "\nerror cant open newick file $file\n";
	my $current_tree = "";
	while(my $line = <TREE>){$current_tree .= $line};
	close TREE;

	# make fasta file for the available COIs, for current terminal set
	my $current_fasta_seqs = ""; my %prune_list = ();
	foreach my $terminal(@current_terminals) # note first item in list is probably outgroup.
		{
		# print "terminal:$terminal\n";
		my $terminalCOPY = $terminal;
		if($terminal =~ s/_cf_|_nr_/_/)
			{
			$current_tree =~ s/$terminalCOPY/$terminal/;
			print "\tremoving taxonomic shorthand\n";
			};
		if($terminal =~ s/_\w*\d+\w*$/_sp/){print "\tconverting non canonical to species unknown\n"};
		my $sourcetree_species_found =0;
		my $sequence_retrieved = "";
		if($sequences{$terminal} =~ /\w/) # simple, found sequence for exact match of binomial
			{
			$sourcetree_species_found=1;$coi_available_for_terminal++;
			$sequence_retrieved = ">$terminal\n$sequences{$terminal}\n";
			}elsif($terminal =~ /^(\w+)_.*\d/) # non canonical species, not many of these, should fix in original newick
			{
			$prune_list{$terminal} = 1;
			my $genus = $1 . "_sp";
			if($sequences{$genus} =~ /\w/){die "\nto implement.\n"};# non canonical species, but genus found. can only use if single genue exemplar

			}elsif($terminal =~ /^(\w+)_sp$/) # species ambiguous
			{
			my $genus = $1; # print "\t\t$terminal is species ambigous, of genus $genus\n";
			if($current_genera{$genus} == 1 && $sequences{$terminal} =~ /\w/)
				{ # can use genus level sequence as long as there are no congeners in the sourcetree. keep terminal id as is.
				$sequence_retrieved = ">$terminal\n$sequences{$terminal}\n";
				}else{
				$coi_not_available_for_terminal_current++;
				$prune_list{$terminal} = 1;
				}
			}elsif($terminal =~ /^[A-Z][a-z]+$/) # wierd terminal id, needs fixing earlier stage
			{
			print "format error terminal $terminal in souecetree $file\n";
			$coi_not_available_for_terminal_current++;
			$prune_list{$terminal} = 1;
			}elsif($terminal =~ /^([A-Z][a-z]+)_/)  # no exact match to binomial. terminal is not species ambiguous,
			{					# though can use congener sequence is only representative of genus
			my $genus = $1;my $newID = $genus . "_sp";
		#	print "\tno sequence for exact match $terminal checking genus:$genus current newick ","genus count:$current_genera{$genus}\n";
			if($current_genera{$genus} == 1 && $sequences{$newID} =~ /\w/)
				{
				$sequence_retrieved = ">$newID\n$sequences{$newID}\n";
				if($current_tree =~ s/$terminal([\:\)\,])/$newID$1/){}else{die "\nerror couldnt switch terminal $terminal\n"};
				$terminal_convert{$terminal} = $newID; # in this case, remove species name. instead it is genus rep
			#	print "current newick species $terminal is only member of genus $genus, ","for which COI is available for congener\n"
				}else{
				$coi_not_available_for_terminal_current++;
				$prune_list{$terminal} = 1;
				};
			}else{
			$coi_not_available_for_terminal_current++;
			$prune_list{$terminal} = 1;
			print FIND_COIS "$terminal\n";
			$queries_for_NCBI{$terminal} = 1;
			};

		print "\tsourcetree_species_found:$sourcetree_species_found\n";
		$current_fasta_seqs .= $sequence_retrieved;

		}; # foreach my $terminal(@current_terminals)
	print "coi_available_for_terminal:$coi_available_for_terminal\ncoi_not_available_for_terminal_current:$coi_not_available_for_terminal_current\n\n";

	if($job == 3)
		{

		open(SEQS, ">current_COIs") || die "\nerror 98.\n";
		print SEQS "$current_fasta_seqs";
		close SEQS;
		$current_tree =~ s/\:[\d\.]+//g;
		open(CURRENT_TREE, ">current_Newick1") || die "\nerror 178\n";
		print CURRENT_TREE $current_tree;
		close CURRENT_TREE;

		# need to prune the tree for species with no COI available ; edit may do this step later, since might find some of the missing species
		my @prunelist = keys %prune_list;@prunelist = sort @prunelist;
		 my $command = "perl ~/usr_scripts/prune_tree.pl -treefile current_Newick1 -output current_Newick2 -start_prune_list @prunelist -end_prune_list";
		print "command:$command\n";system($command);
		 my $command = "perl ~/usr_scripts/format_conversion.pl current_COIs current_COIs.phy fasta phylip";
		print "command:$command\n";system($command);
		 my $command = "raxmlHPC-8.2.4 -s current_COIs.phy -n $file -m GTRCAT -c 4 -p 123 -g current_Newick2";
		print "command:$command\n";system($command);
		};

	####################################################################################
	}; # for my $i(0 .. $#newick_files)

close FIND_COIS;

@NCBI_search_species = keys %queries_for_NCBI; @NCBI_search_species = sort @NCBI_search_species;

if($job == 2)	{
		open(BASH, ">bash_commands") || die "\nerro 144.\n";
		};
foreach my $binomial(@NCBI_search_species)
	{
	#print "search ncbi for missing:$binomial\n";
	if($binomial =~ /(.+)_(.+)/)
		{
		my $genus = $1; my $species = $2;
		my $command = "esearch -db nuccore -query \'\"$genus $species\" [ORGN] AND (\"COI\" [TITLE] OR \"COX1\" [TITLE] OR \"cytochrome c oxidase subunit I\" [TITLE] OR \"cytochrome c oxidase subunit 1\" [TITLE])\' -mindate 2000 -maxdate 2023 | efetch -format fasta > $binomial.NCBIspecies";
		print BASH "$command\n";
#esearch -db nuccore -query '"Eucalyptus" [ORGN] AND "astringens" [ORGN] AND ("COI" [TITLE] OR "COX1" [TITLE] OR "cytochrome c oxidase subunit I" [TITLE] OR "cytochrome c oxidase subunit 1" [TITLE])' -mindate 2000 -maxdate 2023 | efetch -format fasta > Eucalyptus_astringens.NCBI.fas
		};

	};
close BASH;


}; # sub get_sequences_for_newicks
#######################################################################################################################












