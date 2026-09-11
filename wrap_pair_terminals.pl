


$pairs_file 	= $ARGV[0];
$COL_tax_counts = $ARGV[1]; # not used.




#########################################################################
my $COL_lines=0;
# open(COL, $COL_tax_counts) || die "\nerror 15, cant open $COL_tax_counts\n";
while(my $line = <COL>)
	{

	$COL_lines++;
	};
# close COL;
# print "\nread COL_tax_counts, $COL_lines lines.\n";
#########################################################################


open(LOG, ">missing_lineage_files") || die "\nerror 9\n";

# this file is appended to during each invokation of the wrapped script, so needs deleting here at start.
system("rm treepair_match_LOG.txt");
system("rm treepair_match_LOG2.txt");


#########################################################################
$pair_index=0;
open(IN, $pairs_file) || die "\nerror 10, cant open $pairs_file\n";
while(my $line = <IN>)
	{
	if($line =~ /^(\S+)\t(\S+)/)
		{
		my $newickA = $1; my $newickB = $2;
		
		#################################################
		standardize_pair($newickA,$newickB,$pair_index);# 
		#################################################

		$pair_index++;
		}

	};
close IN;
close LOG;

#########################################################################

sub standardize_pair
{
my $newickA = $_[0];my $newickB= $_[1];my $pairID=$_[2]; # $pairID= "pairID$pairID"; 
print "sub standardize_pair, $pairID, $newickA, $newickB\n";

if (-e "lineage_assignments.$newickA" && -e "lineage_assignments.$newickB") 
	{

	# test pair of phylogenies with overlapping taxa. Arguments are treeA treeB pairID
	system( "perl pair_terminals.pl $newickA $newickB $pairID $COL_tax_counts");
	# outputs are 2 newicks, modified with some higher taxon names where inexact matches have been made between tree pair,
	# list of overlapping terminals to be retained in the tree pair, and a log file.

#	die "";

	}else{
	print "lineage files not found for pair $newickA $newickB, wont invoke pair_terminals.\n";
	print LOG "$newickA\t$newickB\n";
	};

};
#########################################################################


