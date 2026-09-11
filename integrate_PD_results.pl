
# 
# study specific script used for work in submission, provided for review
# 

# 
# 20251226:	Option to skip plotting pairs with zero tree distance
# 
# 
# 
# 
# 
# 
# 
# 
##############################################################################################################

$in1 = $ARGV[0]; # disttopo_RESULTS
$in2 = $ARGV[1]; # treepair_match_LOG.txt
$in3 = $ARGV[2]; # newicks.sample_table.202504

$supplementary_figures 			= 0;  # default 0
	$sourcetree_limit		= 20; # 40 gives 9.9 mb, 30 gives 9.5
	$sourcetree_size_limit 		= 8;
	$plot_only_omics_pairs		= 0;
	$dont_plot_identical_pair	= 1;

# plot_only_omics_pairs	== 1, sourcetree_limit	= 40, sourcetree_size_limit 	= 10;  final_list:104

###############################################################################################################
# PAIRWISE PHYLOGENETIC DISTANCE VALUES
open(IN1, $in1) || die "\ncant open file $in1\n";
print "opened $in1\n";
while (my $line = <IN1>)
	{
#	print $line;
# Mengual2023Syrphidae.nwk.pairID973.pruned	Moran2022Eristalinae.nwk.pairID973.pruned	10
# now with multiple PD indices:
# DjernaesKlass2012Blattodea.nwk.pairID996.pruned	DjernaesMurienne2022Blattoidea.nwk.pairID996.pruned	0.16	4.3	2	4.5	4
# here are the distance indices (ie last 5 columns):
# TreeDistance, NyeSimilarity, MatchingSplitDistance, KendallColijn, MASTSize

	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(\S+)\t(\S+)\t(\d+.+)$/)
		{
		my $treeA = $1; my $treeB= $2; my $PD = $3;  print "\ttreeA:$treeA treeB:$treeB PD:($PD)\n";
# treeA:Cruaud2024Chalcidoidea.nwk.pairID997.pruned treeB:Sharkey2011Hymenoptera.nwk.pairID997.pruned PD:(0.244686983423239	12.000980392156918	7.87400787401181	14)
		if($line =~ /NA/){print "no value:$line\n"};
		$pair_PD{ "$treeA"."____"."$treeB" } = $PD;$values_stored++;
		$figure_filename{ "$treeA"."____"."$treeB" } =  "$treeA.$treeB";
		if($PD =~ /^([0-9\.]+)\t/){my $tree_distance=$1;$pair_TreeDistance{ "$treeA"."____"."$treeB" } = $tree_distance};
		};
	};
close IN1;
print "\tvalues_stored:$values_stored\n"; # die "";
###############################################################################################################
#
#
# 	PAIRWISE DETAILS (treepair_match_LOG.txt)
#
#
open(IN2, $in2) || die "\ncant open file $in2\n";
print "opened $in2\n";
while (my $line = <IN2>)
	{
#	print $line;
# pairID1738	Zwick2011Bombycoidea.nwk	RegierCook2008Bombycoidea.nwk	42	4	superfamily	Bombycoidea	genusArtace(Artace_cribraria,Artace_cribaria);subfamilyCeratocampinae(Citheronia_sepulcralis,Eacles_imperialis);subfamilyCercophaninae(Janiodes_laverna,Janiodes_sp);subfamilySaturniinae(Saturnia_naessigi,Antheraea_polyphemus);
# pairID$pairID treeA  treeB      scalar@exact_matches,    scalar@equivelents   shared_taxon   print_string, 

	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\w+)\t(\S+)\t([\d\.]+)/)
		{
		my $pairID = $1; my $nwk1 = $2; my $nwk2 = $3; my $exact_matches = $4; my $equivelents = $5;my $rank = $6;my $taxon = $7; my $proportion_overlap = $8; 
		$sum_overlap{$pairID} = $exact_matches+$equivelents;
		$proporion_overlaps{$pairID} = $proportion_overlap;

		# print "line:$line\n\tproportion_overlap:$proportion_overlap\n";
		if($exact_matches >= 1)
			{
			$overlap_index{$pairID} = $exact_matches/($exact_matches+$equivelents);
			};
		$which_rank{$pairID} = $rank;$all_ranks{$rank}++;
		$store_taxon{$pairID} = "$rank $taxon";

		$match_counts_parsed++;

		if($rank eq "NA")	{$which_rank_value{$pairID}="NA"};

		if($rank eq "clade")	{$which_rank_value{$pairID}="NA";print "$line";print "\nerror unquantifyable rank, go fix pair match script\n"}; # 
		if($rank eq "cohort"){$which_rank_value{$pairID}="NA";print "$line";print "\nerror unquantifyable rank, go fix pair match script\n"};

		if($rank eq "infraclass"){$which_rank_value{$pairID}=7};
		if($rank eq "superorder"){$which_rank_value{$pairID}=7};
		if($rank eq "order"){$which_rank_value{$pairID}=6};
		if($rank eq "infraorder"){$which_rank_value{$pairID}=5};
		if($rank eq "parvorder"){$which_rank_value{$pairID}=5};
		if($rank eq "suborder"){$which_rank_value{$pairID}=5};
		if($rank eq "superfamily"){$which_rank_value{$pairID}=4};
		if($rank eq "family"){$which_rank_value{$pairID}=3};
		if($rank eq "subfamily"){$which_rank_value{$pairID}=3};
		if($rank eq "tribe"){$which_rank_value{$pairID}=2};
		if($rank eq "genus"){$which_rank_value{$pairID}=1};
		if($rank eq "subgenus"){$which_rank_value{$pairID}=1};
# nk:NA 5
# rank:clade 115
# rank:cohort 5
# rank:no_rankAschiza 1
# rank:no_rankTermitoidae 1

# rank:infraclass 2
# rank:superorder 11
# rank:order 81
# rank:suborder 119
# rank:infraorder 128
# rank:no_rankSchizophora 26
# rank:parvorder 8

# rank:superfamily 272
# rank:family 523
# rank:subfamily 244
# rank:no_rankStaphylininae_group 2
# rank:no_rankTachyporinae_group 3
# rank:tribe 48
# rank:genus 26
# rank:subgenus 17




		}else{
		print "not parsed:$line\n";
		};
	};
print "\tmatch_counts_parsed:$match_counts_parsed\n";
my @all_ranks_array = keys %all_ranks;@all_ranks_array=sort@all_ranks_array;
foreach my $rank(@all_ranks_array){print "rank:$rank $all_ranks{$rank}\n"};
###############################################################################################################
# TREE CLASSES
open(IN3, $in3) || die "\ncant open file $in3\n";
print "opened $in3\n";
while (my $line = <IN3>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "line:$line\n";
# HeckenhauerFrandsen2022Trichopte.nwk	O	Trichoptera	2000	0	OME
# ThomasFrandsen2020Trichoptera.nwk	G	Trichoptera	6	0	MG
# MalmJohanson2013Trichoptera.nwk	G	Trichoptera	4	0	MG
# GePeng2023Trichoptera.nwk     	M	Trichoptera	13	0	MITO
# LinardArribas2017Trichoptera.nwk	M	Trichoptera	13	0	MITO
# FrandsenHolzenthal2024Trichoptera.nwk	O	Trichoptera	3206	0	TR
# JohansonEspeland2010Ecnomidae.nwk	G	Trichoptera	4	0	MG
# PengGe2022Philopotamidae.nwk  	M	Trichoptera	13	0	MITO

# note final column is molecular subclass, might be restricted to omics only due to most interest in those
# thus may get a bunch of NAs for everything else
# this variable prob not included in model anyways

	if($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/)
		{
		my $newick = $1; my $class = $2; my $order = $3; my $charsMol = $4; my $charsMorph = $5; my $classMol = $6;
		if($order =~ /Coleoptera|Diptera|Hymenoptera|Lepidoptera/)
			{
			}elsif($order =~ /Dermaptera|Plecoptera|Orthoptera|Grylloblattodea|Embioptera|Phasmatodea|Mantodea|Blattodea/)
			{
			$order="Polyneoptera";
			}elsif($order =~ /Ephemeroptera|Odonata/)
			{
			$order="Palaeoptera"
			}elsif($order =~ /Thysanoptera|Hemiptera/)
			{
			$order = "Condylognatha"
			}else{
			 $order="Other";
			};



		$orders_present{$order}++;
		$which_class{$newick}=$class;$which_order{$newick}=$order;$charactersMol{$newick}=$charsMol;$charactersMorph{$newick}=$charsMorph;
		unless($classMol eq "NA"){$which_classMol{$newick}=$classMol;$mol_subclasses{$classMol}=1};
		$classes_parsed++;
		if($newick =~ /^[A-Za-zł]+([12][90]\d\d)/)
			{
			my $year = $1;$which_year{$newick}=$year
			}else{die "\nerror. cant parse year from $newick.\n"};
		}else{
		print "not parsed:$line\n";
		};
	};
print "\tclasses_parsed:$classes_parsed\n";
my @orders_present1 = keys %orders_present;@orders_present1=sort@orders_present1;
foreach my $order(@orders_present1)
	{
	my $pair_count=$orders_present{$order};print "$order\t$pair_count\n";
	};

my @molecular_subclasses = keys %mol_subclasses;@molecular_subclasses=sort @molecular_subclasses;
print "molecular_subclasses:@molecular_subclasses\n";
# die "\nERROR\n";

###############################################################################################################
@pairs_array = keys %pair_PD; @pairs_array = sort @pairs_array;

if($supplementary_figures == 0){open(RESULTS, ">integrated_result") || die "\nerror 62.\n"};
print RESULTS "pairID\tnwk1\tnwk2\tclass1\tclass2\tclass_pair\tsum_overlap\tTreeDistance\tNyeSimilarity\tMatchingSplitDistance\tKendallColijn\tMASTSize",
	"\tyear\tyear2\torder\trank\tsum_loci\tmorph_chars\trank_numeric\tloci_diff\texact_prop\tprop_overlap\tclass.pair2\n";
# this as was previously written later in r as column header:
# nwk1,nwk2,class1,class2,class.pair,sum_overlap,pairPD1,pairPD2,pairPD3,pairPD4,pairPD5,year,year2,order,rank,sum_loci,
#  morph_chars,rank_numeric,loci_diff,exact_prop,prop_overlap,class.pair2



foreach my $pair(@pairs_array)
	{
	my $pairPD = $pair_PD{$pair}; # print "pair:$pair pairPD:$pairPD\n";
# ZhouXiong2021Chalcidoidea.nwk.pairID1796.pruned____Peters2018Chalcidoidea.nwk.pairID1796.pruned pairPD:14
	if($pair =~ /^(\S+nwk)\.(pairID\d+)\.pruned____(\S+\.nwk)\.pairID\d+\.pruned/)
		{
		my $nwk1 = $1; my $pairID = $2; my $nwk2 = $3;
		#########################################################
		if($which_class{$nwk1} =~ /\w/ && $which_class{$nwk2} =~ /\w/ && $sum_overlap{$pairID} =~ /\d/)
			{
			my $year = "NA";if($which_year{$nwk1}=~ /\d/){$year=$which_year{$nwk1}};
			my $year2 = "NA";if($which_year{$nwk2}=~ /\d/){$year2=$which_year{$nwk2}};
			my $rank_numeric = "NA"; if($which_rank_value{$pairID} =~ /\d/){$rank_numeric = $which_rank_value{$pairID}};
			my $overlap_index_current = "NA"; if($overlap_index{$pairID} =~ /\d/){$overlap_index_current = $overlap_index{$pairID}};
			my $proporion_overlap = "NA"; if($proporion_overlaps{$pairID} =~ /\d/){$proporion_overlap = $proporion_overlaps{$pairID}};


			my @classes = ($which_class{$nwk1},$which_class{$nwk2});@classes = sort @classes;my $classes_string = join '',@classes;
			my $current_order = "NA";if($which_order{$nwk1} =~ /\w/){$current_order = $which_order{$nwk1}};
			my $charsMol ="NA"; if($charactersMol{$nwk1} =~ /\d/ && $charactersMol{$nwk2} =~ /\d/){$charsMol = $charactersMol{$nwk1}+$charactersMol{$nwk2}};
			my $charsMol_diff ="NA"; if($charactersMol{$nwk1} =~ /\d/ && $charactersMol{$nwk2} =~ /\d/){$charsMol_diff = $charactersMol{$nwk1}-$charactersMol{$nwk2}};
			my $charsMorph = "NA";if($charactersMorph{$nwk1} =~ /\d/ && $charactersMorph{$nwk2} =~ /\d/){$charsMorph = $charactersMorph{$nwk1}+$charactersMorph{$nwk2}};

			my $classMol = "NA";
			if($which_classMol{$nwk1} =~ /\w/ && $which_classMol{$nwk2} =~ /\w/)
				{
				my @classes2 = ($which_classMol{$nwk1},$which_classMol{$nwk2});@classes2 = sort @classes2;$classMol = join '.',@classes2;

				# insufficient obs for only within AHE,TR,UCE
				# if($classMol =~ /AHE.AHE|TR.TR|UCE.UCE|MITO.MITO|MG.MG/){}else{$classMol="NA"};
				if($classMol =~ /AHE.AHE|AHE.TR|AHE.UCE|TR.TR|TR.UCE|UCE.UCE|MITO.MITO|MITO.TR|MITO.UCE/){}else{$classMol="NA"};

# :AHE.AHE AHE.MG AHE.MITO AHE.OME AHE.SG AHE.TR AHE.UCE MG.MG MG.MITO MG.OME MG.SG MG.TR MG.UCE MITO.MITO MITO.OME MITO.SG MITO.TR MITO.UCE
#  OME.OME OME.SG OME.TR OME.UCE SG.SG SG.TR SG.UCE TR.TR TR.UCE UCE.UCE

# NEXT: to exclude classes, better NA them so other aspects contribute to model.




				$all_mol_classes{$classMol}++;
				};
			$unique_class_combinations{$classes_string}++;
		
			#################################################################
			# those with lots of obs:
		#	if($classes_string =~ /CG|GG|GM|GO|GP|MM|MO|OO/)
		#	all important:
		#	if($classes_string =~ /CC|CG|CM|CO|GG|GO|GM|GP|MM|MO|MP|OO|OP/)
		#     within class only:
		#	if($classes_string =~ /GG|OO|MM|CC|PP/)
			# ALL incl barcode
			if($classes_string =~ /\w/)
				{
			#	print "$classes_string:";
			#	if($classes_string =~ /GG|OO|MM|CC|PP/){}else{$classes_string = "NA"};
			#	if($classes_string =~ /CC|GG|MM|OO|PP/){}else{$classes_string = "XX"};

				# those >100, others NULL
				if($classes_string =~ /CG|GG|GM|GO|GP|MM|MO|OO/){}else{$classes_string = "NA"};


			#	if($classes_string =~ /GP|OP|MP/)
			#		{# betwen mol / morph
			#		$classes_string = "AA"
			#		}elsif($classes_string =~ /CG|CM|GM|GO|MO|CO/)
			#		{# between mol
			#		$classes_string = "BB"
			#		}elsif($classes_string =~ /BB|BC|BG|BM|BO|BP/)	
			#		{
			#		$classes_string = "NA"
			#		};
			#	print "$classes_string\t";

				print RESULTS "$pairID\t$nwk1\t$nwk2\t$which_class{$nwk1}\t$which_class{$nwk2}\t$classes_string\t$sum_overlap{$pairID}\t$pairPD\t$year\t$year2\t$current_order\t$which_rank{$pairID}\t$charsMol\t$charsMorph\t$rank_numeric\t$charsMol_diff\t$overlap_index_current\t$proporion_overlap\t$classMol\n";
				$result_printed++;$classes_written{$classes_string}++;
				$paired_newicks{$nwk1}=1;$paired_newicks{$nwk2}=1;
				}else{
				print "classes string ($classes_string) doesnt match user specification\n";
				};
			#################################################################


			$details_retreived++;
			}else{
			$details_not_retreived++;
			};
		#########################################################
		}else{
		print "cant parse pair string $pairPD\n";
		};
	};

my @all_mol_classes_keys = keys %all_mol_classes;@all_mol_classes_keys=sort@all_mol_classes_keys;print "all_mol_classes_keys:@all_mol_classes_keys\n";

print "
details_retreived:$details_retreived
details_not_retreived:$details_not_retreived

wrote integrated_result, result_printed:$result_printed

";
close RESULTS;


my @calsses = keys %unique_class_combinations;@calsses=sort @calsses;
print "class\tnumber_read\tnumber_written\n";
foreach my $class(@calsses)
	{
	print "$class\t$unique_class_combinations{$class}\t$classes_written{$class}\n";
	};

###############################################################################################################

if($supplementary_figures == 0){open(OUT6, ">list_all_paired_newicks") || die ""}; 
@all_paired_newicks = keys %paired_newicks;@all_paired_newicks = sort @all_paired_newicks;
foreach my $newick(@all_paired_newicks)
	{
	print OUT6 "$newick ";
	};
close OUT6;
print "
count pair newicks $#all_paired_newicks
";

###############################################################################################################

if($supplementary_figures == 1)
	{

open(IPD_DETAILS, "/home/douglas/databases/Phylo/P2_structure/data/sourcetrees_table.202511B") || die "\nerror 298.\n";
while(my $line = <IPD_DETAILS>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print $line;
	if($line =~ /^(\S+\.nwk)\t.+	(https:\/\/doi\S+)\t/)
		{
		my $treeID = $1; my $doi = $2;$dois{$treeID} = $doi;$dois_read++;
		};
	};
close IPD_DETAILS;
print "
dois_read:$dois_read
";

my @pdf_merges=();$merge_batch=1;

foreach my $pair(@pairs_array)
	{
	my $pairPD = $pair_PD{$pair};  print "pair:$pair pairPD:$pairPD\n";
# ZhouXiong2021Chalcidoidea.nwk.pairID1796.pruned____Peters2018Chalcidoidea.nwk.pairID1796.pruned pairPD:14

	my $figure = $figure_filename{$pair};
	my $file_found=0;
	if (-e "/home/douglas/databases/Phylo/P2_structure/analysis3/$figure.pdf") {$file_found=1;$trees_plotted++}else{$trees_not_plotted++};

	

	if($pair =~ /^(\S+nwk)\.(pairID\d+)\.pruned____(\S+\.nwk)\.pairID\d+\.pruned/)
		{
		my $nwk1 = $1; my $pairID = $2; my $nwk2 = $3;

		$sourcetree_sampled{$nwk1}++;$sourcetree_sampled{$nwk2}++;
		my $skip_pair=0;
		if($sourcetree_sampled{$nwk1} >= $sourcetree_limit || $sourcetree_sampled{$nwk2} >= $sourcetree_limit){$skip_pair=1};

		
		#########################################################
		if($which_class{$nwk1} =~ /\w/ && $which_class{$nwk2} =~ /\w/ && $sum_overlap{$pairID} >= $sourcetree_size_limit)
			{
			my $year = "NA";if($which_year{$nwk1}=~ /\d/){$year=$which_year{$nwk1}};
			my $year2 = "NA";if($which_year{$nwk2}=~ /\d/){$year2=$which_year{$nwk2}};
			my $rank_numeric = "NA"; if($which_rank_value{$pairID} =~ /\d/){$rank_numeric = $which_rank_value{$pairID}};
			my $overlap_index_current = "NA"; if($overlap_index{$pairID} =~ /\d/){$overlap_index_current = $overlap_index{$pairID}};
			my $proporion_overlap = "NA"; if($proporion_overlaps{$pairID} =~ /\d/){$proporion_overlap = $proporion_overlaps{$pairID}};
			my $nwk1_doi = "NA";if($dois{$nwk1} =~ /\w/){$nwk1_doi = $dois{$nwk1}};my $nwk2_doi = "NA";if($dois{$nwk2} =~ /\w/){$nwk2_doi = $dois{$nwk2}};

	
			if($plot_only_omics_pairs==1)
				{
				if($which_class{$nwk1} eq "O" && $which_class{$nwk2} eq "O"){}else{$skip_pair=1};
				};

			if($dont_plot_identical_pair==1)
				{
				if($pair_TreeDistance{$pair}=~ /\d/)
					{
					my $current_pair_TD = $pair_TreeDistance{$pair}; # print "current_pair_TD:$current_pair_TD\n";
					if($current_pair_TD <= 0.00001){$skip_pair=1};
					}else{
					print "warning, couldnt retrieve tree distance for currnet pair $pair\n";
					};
				};
		#	die "";

			if($file_found == 1 && $skip_pair ==0 )
				{
				print "$store_taxon{$pairID}\n\t$nwk1, $which_class{$nwk1}, $nwk1_doi\n\t$nwk2, $which_class{$nwk2}, $nwk2_doi\n";
				$final_list++;

my $R_script2 = "
pdf(file = \"/home/douglas/databases/Phylo/P2_structure/analysis3/$figure.pair_titlepage.pdf\", 12, 3)
plot.new()
par(usr = c(0, 1, 0, 1))
text(0.5,0.7,labels=\"$store_taxon{$pairID}\", cex=1.2)
text(0.5,0.5,labels=\"$nwk1, $which_class{$nwk1}, $nwk1_doi\", cex=1.2)
text(0.5,0.3,labels=\"$nwk2, $which_class{$nwk2}, $nwk2_doi\", cex=1.2)
dev.off()
";
	system("rm Current_R_Script2");
	open(CURRENT_R2, ">Current_R_Script2") || die "\nerrr 81\n";
	print CURRENT_R2 "$R_script2\n";
	close CURRENT_R2;

	system("R < Current_R_Script2 --vanilla --slave");

push @pdf_merges, "/home/douglas/databases/Phylo/P2_structure/analysis3/$figure.pair_titlepage.pdf";
push @pdf_merges, "/home/douglas/databases/Phylo/P2_structure/analysis3/$figure.pdf";

if($#pdf_merges >= 400)
	{
	my $pdf_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=/home/douglas/databases/Phylo/P2_structure/analysis3/merged.$merge_batch.pdf @pdf_merges";
	$merge_batch++;
	system($pdf_command); # print "pdf_command:$pdf_command\n";
	@pdf_merges=();
	}else{
	};

				};

			};
		};
	};


# any remaining
my $pdf_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=/home/douglas/databases/Phylo/P2_structure/analysis3/merged.$merge_batch.pdf @pdf_merges";
$merge_batch++;
system($pdf_command); # print "pdf_command:$pdf_command\n";
@pdf_merges=();


print "
trees_plotted:$trees_plotted
trees_not_plotted:$trees_not_plotted
final_list:$final_list
	";

	};

###############################################################################################################





print "\nFIN.\n";
exit;




