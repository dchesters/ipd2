
# 
# study specific script used for work in submission, provided for review
# 

# 
# 
# 
# 
# 2025-04-11:	Initiated.
# 		Taxonomically matches terminals between a pair of phylogenies.
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#######################################################################################################################################

$treeA = $ARGV[0];
$treeB = $ARGV[1];
$pairID = $ARGV[2];unless($pairID =~ /\d/){die "\npair_terminals command error 23, pairID != value \n"}
$COL_tax_counts = $ARGV[3];


# /home/douglas/databases/Phylo/P2_structure/analysis1/constraints.MillerBergsten2012Gyrinidae.nwk
# Adelocera_sp	genusAdelocera
# Agriotes_acuminatus	speciesAgriotes_acuminatus

# /home/douglas/databases/Phylo/P2_structure/analysis1/lineage_assignments.KundrataBocakova2013Elateroidea.nwk
# Anoplischius_sp	root:Hexapoda class:Insecta clade:Dicondylia subclass:Pterygota infraclass:Neoptera cohort:Endopterygota order:Coleoptera suborder:Polyphaga infraorder:Elateriformia superfamily:Elateroidea family:Elateridae subfamily:Elaterinae genus:Anoplischius 
# Anostirus_purpureus	root:Hexapoda class:Insecta clade:Dicondylia subclass:Pterygota infraclass:Neoptera cohort:Endopterygota order:Coleoptera suborder:Polyphaga infraorder:Elateriformia superfamily:Elateroidea family:Elateridae subfamily:Prosterninae genus:Anostirus species:Anostirus_purpureus 

$constraintsA = "constraints.$treeA";$constraintsB = "constraints.$treeB";
$lineagesA = "lineage_assignments.$treeA";$lineagesB = "lineage_assignments.$treeB";


%all_terminals; # unique species of both trees
$species_for_single_tree=0;$species_for_both_trees=0;

open(LOG5, ">pair_terminals_LOGFILE.txt") || die "\nerror 42. cant open log file output\n";

#######################################################################################################################################

# Catalogue of life species counts for each taxon

my $COL_lines=0;
#########open(COL, $COL_tax_counts) || die "\nerror 15, cant open, cant open Catalogue of Life taxon counts file:$COL_tax_counts\n";
while(my $line = <COL>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "line:$line\n";
# tribe:Villini	797
	if($line =~ /^(\S+)\t(\S+)/){my $current_ranktax=$1;my $ranktax_count=$2;$COL_ranktaxs{$current_ranktax}=$ranktax_count};
	$COL_lines++;
	};
close COL;
# print "\nread COL_tax_counts, $COL_lines lines.\n";


#######################################################################################################################################

open(IN1, $constraintsA) || die "\npair_terminal.pl, error 29. cant open file:$constraintsA\n";
print "\topened $constraintsA\n";
while(my $line = <IN1>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "line:$line\n";
# ine:Simplocaria_sp	genusSimplocaria
# line:Spectralia_gracilipes	genusSpectralia

	if($line =~ /^(\S+)\t(\S+)/)
		{
		my $columnA = $1;my $columnB = $2; # print "columnA:$columnA columnB:$columnB\n";
		$assignments{$columnA}=$columnB;$assignments_stored++;
		$all_terminals{$columnA}=1;
		}else{
		print "warning, line not parsed:$line\n";
		};
	};
close IN1;
print "\tassignments_stored:$assignments_stored\n";

#########################################################################

open(IN2, $lineagesB) || die "\nerror 37. cant open lineagesB:$lineagesB\n";
print "\topened $lineagesB\n";
while(my $line = <IN2>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "line:$line\n";

# Rhagonycha_lignosa	root:Hexapoda class:Insecta clade:Dicondylia subclass:Pterygota infraclass:Neoptera cohort:Endopterygota order:Coleoptera suborder:Polyphaga infraorder:Elateriformia superfamily:Elateroidea family:Cantharidae subfamily:Cantharinae genus:Rhagonycha species:Rhagonycha_lignosa 

	if($line =~ /^(\S+)\t(.+)/)
		{
		my $species = $1; my $lineage = $2;$lineage =~ s/\://g; # if($lineage =~ /clade/){ print "lineage:$lineage\n"};
		$all_terminals{$species}=1;


# cohortEndopterygota; cladeDicondylia change to superorder
while($lineage =~ s/ cohortEndopterygota/ superorderEndopterygota/){};
while($lineage =~ s/ cladeDicondylia/ superorderDicondylia/){};
while($lineage =~ s/ cladeObtectomera/ suborderObtectomera/){};
while($lineage =~ s/ cladeAnthophila/ superfamilyAnthophila/){};
while($lineage =~ s/ cladeDitrysia/ suborderDitrysia/){};
while($lineage =~ s/ cladeApoditrysia/ suborderApoditrysia/){};
while($lineage =~ s/ cohortPolyneoptera/ superorderPolyneoptera/){};
while($lineage =~ s/ cohortParaneoptera/ superorderParaneoptera/){};
while($lineage =~ s/ cladeCyclorrhapha/ suborderCyclorrhapha/){};
while($lineage =~ s/ no_rank\S+\s/ /){};
# exclude user specified ranks from consideration; "clade", "cohort", "no_rank"
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
		$lineages{$species} = $lineage; $lineages_stored++;
		while($lineage =~ s/\s(\w+)\s$/ /)
			{
			my $tax = $1; # print " tax:$tax";
			$terminals_for_taxon{$tax}.= "$species,"
			};
		};
	};
close IN2;
 print "\tlineages_stored:$lineages_stored\n";
#######################################################################################################################################

my @species_overall = keys %all_terminals;
$species_in_one_sourcetree=0;
foreach my $terminal(@species_overall)
	{
	if($assignments{$terminal} =~ /\w/){$species_in_one_sourcetree++};
	};

my $proportion_overlap = $species_in_one_sourcetree/scalar @species_overall;

my @terminals = keys %assignments;@terminals = sort @terminals;
foreach my $terminal(@terminals)
	{
	my $tax_assigned = $assignments{$terminal};
	# print "terminal:$terminal highest unique taxon:$tax_assigned\n";
	if($terminals_for_taxon{$tax_assigned} =~ /\w/)
		{
		my $matched_terminals = $terminals_for_taxon{$tax_assigned};$matched_terminals =~ s/\,$//;
		# print "\tmatched_terminals:$matched_terminals\n";
		my @split_matches = split /\,/,$matched_terminals;
		my $exact_match=0;
		foreach my $matched(@split_matches)
			{
			if($matched eq $terminal){$exact_match=1};
			};
		if($exact_match==1)
			{
			$exact_matches{$terminal}=1; $count_exact_matches++;# print "\texact match\n\n"
			}else{
			$inexact_matches{$terminal} = $matched_terminals; # print "\tINexact matches, sorting later.\n\n";
			};
		}else{
		# print "\tno matches in paired tree.\n";
		};
	};

print "\tcount_exact_matches:$count_exact_matches\n";

#############################################################################################################
# need to first tag all exact matches, 
# then choose exemplars for inexact matches.

my @inexactly_matched = keys %inexact_matches; # maybe should keep this randomly sorted
# print "sorting inexact matches ($#inexactly_matched)\n";
foreach my $terminal(@inexactly_matched)
	{
	my $matches = $inexact_matches{$terminal}; # print "\tterminal:$terminal matches:$matches\n";
	my @split_matches = split /\,/,$matches;
 	my @candidates = ();
	foreach my $matched(@split_matches)
		{
		if($exact_matches{$matched}==1)
			{
			# print "\t$matched already an exact matched.\n";
			}elsif($candidates_assigned{$matched}==1)
			{
			# print "\t$matched already an INexact match.\n";
			}else{
			push @candidates, $matched;
			};
		};
	if($#candidates >= 0)
		{
		fisher_yates_shuffle( \@candidates );
		# print "\tselected match:$candidates[0]\n\n";
		$inexact_assigned{$terminal} = $candidates[0];
		$candidates_assigned{$candidates[0]}=1;$paired_names{$terminal} = $candidates[0];
		
		print LOG5 "EQUIVELENTS\t$terminal:$assignments{$terminal}:@candidates\n";
		}

	};


#############################################################################################################

open(TREEA, $treeA) || die "\nerror cant open treeA:$treeA\n";while(my $line = <TREEA>){$treeA_newick .= $line};close TREEA;
open(TREEB, $treeB) || die "\nerror cant open treeB $treeB\n";while(my $line = <TREEB>){$treeB_newick .= $line};close TREEB;
open(TREEA_OUT, ">$treeA.pairID$pairID") || die "\nerror cant open OutA\n";
open(TREEB_OUT, ">$treeB.pairID$pairID") || die "\nerror cant open OutB\n";
open(TREEA_LIST, ">$treeA.pairID$pairID.retain") || die "\nerror cant open\n";
open(TREEB_LIST, ">$treeB.pairID$pairID.retain") || die "\nerror cant open\n";
open(TREEA_ORIGINAL, ">$treeA.pairID$pairID.original_names") || die "\nerror cant open\n";
open(TREEB_ORIGINAL, ">$treeB.pairID$pairID.original_names") || die "\nerror cant open\n";


#############################################################################################################



@exact_matches = keys %exact_matches;@exact_matches = sort @exact_matches;
my @overlapping_binoms;
foreach my $terminal(@exact_matches)
	{
	push @overlapping_binoms, $terminal;
	print TREEA_LIST ">$terminal\n";print TREEB_LIST ">$terminal\n";
	print TREEA_ORIGINAL "$terminal\n";print TREEB_ORIGINAL "$terminal\n";
	print LOG5 "EXACT\t$terminal\n";
#	print "exact terminal match between tree pair:$terminal\n";
	};

@equivelents = keys %paired_names;@equivelents = sort @equivelents;

my $print_string;
foreach my $terminal(@equivelents)
	{
	push @overlapping_binoms, $terminal;
	my $matched_to = $paired_names{$terminal};
	my $tax_assigned = $assignments{$terminal};unless($tax_assigned =~ /\w/){die "\nerror 145, no higher tax assignment returned for $terminal\n"};
	# print "taxonomically equivelent match between tree pair:$terminal/$matched_to, both are $tax_assigned\n";
	if($treeA_newick =~ s/$terminal/$tax_assigned/){$treeA_terminal_conversions++}else{die "\nterminal conversion error 157.\n"};
	if($treeB_newick =~ s/$matched_to/$tax_assigned/){$treeB_terminal_conversions++}else{die "\nterminal conversion error 158.\n"};
	$print_string .= "$tax_assigned($terminal,$matched_to);";
	print TREEA_LIST ">$tax_assigned\n";print TREEB_LIST ">$tax_assigned\n";
	print TREEA_ORIGINAL "$terminal\n";print TREEB_ORIGINAL "$matched_to\n";
	};


# print "\ntreeA_terminal_conversions:$treeA_terminal_conversions
# treeB_terminal_conversions:$treeB_terminal_conversions\n";

print TREEA_OUT $treeA_newick;print TREEB_OUT $treeB_newick;
close TREEA_OUT;close TREEB_OUT;
close TREEA_LIST;close TREEB_LIST;
close TREEA_ORIGINAL;close TREEB_ORIGINAL;


# need this in a different format for r plot
open(TREE_MATCH_LOG2, ">>treepair_match_LOG2.txt");
print TREE_MATCH_LOG2 "pairID$pairID\texact\t", scalar @exact_matches, "\n";
print TREE_MATCH_LOG2 "pairID$pairID\tequivelent\t", scalar @equivelents, "\n";
close TREE_MATCH_LOG2;

close LOG5;

########################################

# need to know shared ranktaxon for remaining terminal set:

my $binoms_list = join '	', @overlapping_binoms;
my $shared_taxon = get_shared_taxa2($binoms_list);
print "\tshared_taxon:$shared_taxon\n";

########################################

open(TREE_MATCH_LOG, ">>treepair_match_LOG.txt");
print TREE_MATCH_LOG "pairID$pairID\t$treeA\t$treeB\t", scalar @exact_matches, "\t", scalar @equivelents, "\t$shared_taxon\t$proportion_overlap\t", $print_string, "\n";
close TREE_MATCH_LOG;

print "\nFIN\n";
exit;

#############################################################################################################

sub fisher_yates_shuffle 	# http://perldoc.perl.org/perlfaq4.html
	{
 	my $deck = shift;my $i = @$deck;
	while (--$i) 
		{my $j = int rand ($i+1);
		@$deck[$i,$j] = @$deck[$j,$i]}
	}

#############################################################################################################


# gets lineage (series of taxa) shared by all members. also option to get only lowest shared taxa


sub get_shared_taxa2
{
my $tax_list = $_[0];my $return_item = $_[1]; # if $_[1] == 2, return inclusive lineage rather than single taxon
my @tax_array = split /\t/ , $tax_list;
my $shared_tax_substring = "";

####################################################################################################################
# 2022-08-06: previous implementation took the first terminal of list, inferred lineage to be compared to remainder.
#		this assumes near complete taxonomic inference for terminals 
#		(since if lineage not retreived for first in list, process would break),
#		which is not the case for parsing source trees which are regulaly non-molecular.
#		thus need to instead look for first instance of lineage across terminal list.
#
my $test_taxonomy = "";
my $terminal_list_index =0;
for my $index(0 .. $#tax_array)
	{
	$terminal_list_index = $index;  # print "\nsub get_shared_taxa, $index OF $#tax_array\n";
	my $test_species = $tax_array[$index];

	# hash complete_lineage_for_this_species, has not just species, 
	# but has lineage for every taxon in ncbi, only modification is space->underscore
	$test_taxonomy = $lineages{$test_species};
	# print "line928. test_species:$test_species test_taxonomy:$test_taxonomy\n";

	unless($test_taxonomy =~ /[a-z][a-z][a-z]/)
		{
		my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
		$test_taxonomy = $lineages{$genusname};
		# print "cant find taxonomy for species ($test_species) trying genus ($genusname)\n";
		if($test_taxonomy =~ /\w+/)
			{
			# print "\tfound using genus name ($genusname), continuing ...\n";
			}else{
			if($test_species =~ /optera_[A-Z]/)
				{print "\nlooks like user format was not processed to an established tax name\n"};
		#	 print "\nerror 1701, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
			$skipnode=1;
			}
		}
	if($test_taxonomy =~ /\w+/)
		{
	#	print "found first lineage at terminal index $index\n";
		last
		}else{
	#	print "keep looking\n"
		};
	};
####################################################################################################################

while($test_taxonomy =~ s/(\s[a-z]+)([A-Z][a-z]+\s)/$1:$2/){}; # hack current dont have colons, this sub was written if inclduede
# print "test_taxonomy:$test_taxonomy\n";

my $most_inclusive_name = "NA";
my $most_inclusive_rank = "NA";
my $most_inclusive_lineage = "";
my $most_inclusive_lineage2 = ""; # including ranks


if($terminal_list_index == $#tax_array)
	{
	# only a single lineage will be found for this list of terminals,
	# which cannot really be called a shared taxon, so return nothing
	print  "\t\tonly a single lineage will be found for this list of terminals, which cannot really be called a shared taxon, so return nothing\n";
	return($most_inclusive_name);
	};

my $skipnode=0;

if($skipnode==0)
{

while($test_taxonomy =~ s/^([^:]+):(\w+)//)
	{
	my $current_rank = $1;my $current_taxname = $2;
	# over counts, should not count taxa at parent node:
	#$number_of_nodes_to_which_this_taxon_has_been_assigned{$current_taxname}++;
	# print "\tcurrent_rank:$current_rank current_taxname:$current_taxname\n";
	$current_rank =~ s/^\s+//;$current_rank =~ s/\s+$//;

	if($#tax_array >= 1)
		{
		# this is the regular situation, there is more than one terminal you have got a lineage for,
		# and you are finding the shared taxa for these

		my $all_members_have_this_name =1;

#		foreach my $tax(@tax_array)
		for my $terminal_index($terminal_list_index .. $#tax_array)
			{
			my $tax = $tax_array[$terminal_index];
			my $test_taxonomy2 = $lineages{$tax};
			my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			unless($test_taxonomy2 =~ /\w/){$test_taxonomy2 = $lineages{$genus_name}};
			while($test_taxonomy2 =~ s/(\s[a-z]+)([A-Z][a-z]+\s)/$1:$2/){};

			# print "\t\tspecies:$tax complete taxonomy:$test_taxonomy2\n";

			# for current taxonomic name in lineage list, see if another terminal also has this
			if($test_taxonomy2 =~ /\w/)
			{
			if($test_taxonomy2 =~ /\:$current_taxname\s/)
				{
				print ""
				}else{
				$all_members_have_this_name =0;#print "0";
			#	print GRAFT_LOG "\t\t\tterminal $tax is not $current_taxname, its lineage:$test_taxonomy2\n";
				}
			};
			}

	#	print "\nall_members_have_this_name:$all_members_have_this_name\n";

		if($all_members_have_this_name == 1)
			{
			$most_inclusive_name = $current_taxname;
			$most_inclusive_lineage .= "$current_taxname\t";
			$most_inclusive_lineage2 .= "$current_rank:$current_taxname\t";
			$most_inclusive_rank = $current_rank;
			};


		}else{

		# alternativly, you only have one taxa, so just parse the whole lineage for this		

		$most_inclusive_name = $current_taxname;
		$most_inclusive_lineage .= "$current_taxname\t";
		$most_inclusive_lineage2 .= "$current_rank:$current_taxname\t";
		$most_inclusive_rank = $current_rank;

		}

	#print "\tmost_inclusive_lineage:$most_inclusive_lineage\n";
	}


};

my $return_object;
my $return_ranktax;

if($return_item == 2)
	{
	$return_object = "$most_inclusive_lineage2"
	}else{ # default
	$return_object = "$most_inclusive_rank\t$most_inclusive_name";
	$return_ranktax = "$most_inclusive_rank:$most_inclusive_name";
	};

###################################################

# print "return_ranktax:$return_ranktax\n";
if($COL_ranktaxs{$return_ranktax} =~ /./)
	{
#	print "\t$return_ranktax found in COL species counts\n";
	
	}else{
#	print "\t$return_ranktax NOT found in COL species counts\n";
#	print "\t$most_inclusive_lineage2\n";

	};

###################################################


return($return_object);

}; # sub get_shared_taxa2

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





