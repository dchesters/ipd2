
# 
# study specific script used for work in submission, provided for review
# 

# 
# 20251109:	multiple phylo distance indices
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

$match_log = $ARGV[0];
$which_job = $ARGV[1];

#########################################################################################################################################
open(IN, $match_log) || die "";
print "opened \$match_log:$match_log\n";
$line_count=0;
while(my $line = <IN>)
	{
	# print $line;
	# 1	Kundrata2016Elateriformia.nwk	KundrataBocakova2013Elateroidea.nwk	74	7	genusDrilonius(Drilonius_striatulus,Drilonius_sp);genusEubrianax(Eubrianax_edwardsii,Eubrianax_sp);subfamilyOtotretinae(Lamellipalpus_pacholatkoi,Flabellotreta_sp);genusLibnetis(Libnetis_sp,Libnetis_granicollis);genusLyponia(Lyponia_n
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/)
		{
		my $treeID = $1; my $newickA = $2; my $newickB = $3;
		my $prune_file_A = "$newickA.$treeID"; my $prune_file_B = "$newickB.$treeID";
		my $retain_list = "$newickA.$treeID.retain"; # this list should be identical for pair.
		if($treeIDs_encountered{$treeID} == 1)
			{}else{
			# use seqfile rather than list file, as former are those retained, and latter are those removed.
			push @prune_commands, "perl prune_tree.pl -enable_root_pruning -remove_branchlengths -seqfile $retain_list -treefile $prune_file_A -output $prune_file_A.pruned\n";
			push @prune_commands, "perl prune_tree.pl -enable_root_pruning -remove_branchlengths -seqfile $retain_list -treefile $prune_file_B -output $prune_file_B.pruned\n";
			};
		$treeIDs_encountered{$treeID} = 1;
		$pruned_tree_pair{$treeID} = "$prune_file_A.pruned\t$prune_file_B.pruned";
		};
	$line_count++;
	};
close IN;
print "line_count:$line_count\n";
#########################################################################################################################################

if($which_job == 1)
{

open(OUT, ">All_prune_commands") || die "\nerror 31.\n";
foreach my $command(@prune_commands)
	{
	print OUT $command;
	};
close OUT;
print "made list of pruning commands, FIN.\n";
};
#########################################################################################################################################

if($which_job == 2)
{
@all_treeIDs = keys %pruned_tree_pair; @all_treeIDs = sort @all_treeIDs;
$index=0;
foreach my $treeID(@all_treeIDs)
	{
	my $current_pair = $pruned_tree_pair{$treeID};
	print "$index of $#all_treeIDs, treeID:$treeID,  pruned newick pair:$current_pair\n";
	my $nwk1;my $nwk2;
	if($current_pair =~ /(\S+)\t(\S+)/){$nwk1=$1;$nwk2=$2}else{die "\nparse error 56.\n"};

my $R_script = "
library(ape)
treeA <- read.tree(\"$nwk1\")
treeB <- read.tree(\"$nwk2\")
overlapping <- intersect(treeA\$tip.label,treeB\$tip.label)
if(length(overlapping) == length(treeA\$tip.label))
 {
print(\"treeA no processing required\")	
 }else{
nonoverlapping <- setdiff(treeA\$tip.label,overlapping)   
treeA <- drop.tip(treeA, nonoverlapping)
 }
if(length(overlapping) == length(treeB\$tip.label))
 {
print(\"treeB no processing required\")	
 }else{
nonoverlapping <- setdiff(treeB\$tip.label,overlapping)   
treeB <- drop.tip(treeB, nonoverlapping)
 }
disttopo<-as.numeric(dist.topo(treeA, treeB))
write(c(\"$nwk1\",\"$nwk2\",disttopo), file = \"disttopo_RESULTS\", ncolumns = 3, append = T, sep = \"\\t\")
";

	system("rm Current_R_Script");
	open(CURRENT_R, ">Current_R_Script") || die "\nerrr 81\n";
	print CURRENT_R "$R_script\n";
	close CURRENT_R;
	system("R < Current_R_Script --vanilla --slave");

############################################################################################

my $R_script2 = "
library(TreeDist)
library(ape)
treeA <- read.tree(\"$nwk1\")
treeB <- read.tree(\"$nwk2\")
overlapping <- intersect(treeA\$tip.label,treeB\$tip.label)
if(length(overlapping) == length(treeA\$tip.label))
 {
print(\"treeA no processing required\")	
 }else{
nonoverlapping <- setdiff(treeA\$tip.label,overlapping)   
treeA <- drop.tip(treeA, nonoverlapping)
 }
if(length(overlapping) == length(treeB\$tip.label))
 {
print(\"treeB no processing required\")	
 }else{
nonoverlapping <- setdiff(treeB\$tip.label,overlapping)   
treeB <- drop.tip(treeB, nonoverlapping)
 }
dist_results <- rep(NA, length.out=5)
dist_results[1] <- TreeDistance(treeA, treeB)
dist_results[2] <- NyeSimilarity(treeA, treeB)
dist_results[3] <- MatchingSplitDistance(treeA, treeB)
dist_results[4] <- KendallColijn(treeA, treeB)
dist_results[5] <- MASTSize(treeA, treeB)
print(dist_results)
write(c(\"$nwk1\",\"$nwk2\",dist_results), file = \"disttopo_RESULTS2\", ncolumns = 7, append = T, sep = \"\\t\")
";
# # NNIDist(tree1, tree2)



	system("rm Current_R_Script2");
	open(CURRENT_R2, ">Current_R_Script2") || die "\nerrr 81\n";
	print CURRENT_R2 "$R_script2\n";
	close CURRENT_R2;
	system("/opt/R/4.2.3/bin/R < Current_R_Script2 --vanilla --slave");

############################################################################################

	$index+=1;

	};# foreach my $treeID(@all_treeIDs)


# regular:
# library(ape)
# treeA <- read.tree("AcostaMorrone2013Phalacropsyllini.nwk.pairID0.pruned")
# treeB <- read.tree("WhitingWhiting2008Siphonaptera.nwk.pairID0.pruned")
# overlapping <- intersect(treeA$tip.label,treeB$tip.label)
# if(length(overlapping) == length(treeA$tip.label))
#  {
# print("treeA no processing required")	
#  }else{
# nonoverlapping <- setdiff(treeA$tip.label,overlapping)   
# treeA <- drop.tip(treeA, nonoverlapping)
#  }
# if(length(overlapping) == length(treeB$tip.label))
#  {
# print("treeB no processing required")	
#  }else{
# nonoverlapping <- setdiff(treeB$tip.label,overlapping)   
# treeB <- drop.tip(treeB, nonoverlapping)
#  }
# disttopo<-as.numeric(dist.topo(treeA, treeB))
# write(c("AcostaMorrone2013Phalacropsyllini.nwk.pairID0.pruned","WhitingWhiting2008Siphonaptera.nwk.pairID0.pruned",disttopo), file = "disttopo_RESULTS", ncolumns = 3, append = T, sep = "\t")


}; # if($which_job == 2)

#########################################################################################################################################


if($which_job == 3)
{
@all_treeIDs = keys %pruned_tree_pair; @all_treeIDs = sort @all_treeIDs;
$index=0;
foreach my $treeID(@all_treeIDs)
	{
	my $current_pair = $pruned_tree_pair{$treeID};
	print "$index of $#all_treeIDs, treeID:$treeID,  pruned newick pair:$current_pair\n";
	my $nwk1;my $nwk2;
	if($current_pair =~ /(\S+)\t(\S+)/){$nwk1=$1;$nwk2=$2}else{die "\nparse error 56.\n"};


my $R_script2 = "
library(TreeDist)
library(ape)

treeA <- read.tree(\"$nwk1\")
treeB <- read.tree(\"$nwk2\")
overlapping <- intersect(treeA\$tip.label,treeB\$tip.label)
if(length(overlapping) == length(treeA\$tip.label))
 {
print(\"treeA no processing required\")	
 }else{
nonoverlapping <- setdiff(treeA\$tip.label,overlapping)   
treeA <- drop.tip(treeA, nonoverlapping)
 }
if(length(overlapping) == length(treeB\$tip.label))
 {
print(\"treeB no processing required\")	
 }else{
nonoverlapping <- setdiff(treeB\$tip.label,overlapping)   
treeB <- drop.tip(treeB, nonoverlapping)
 }

# pdf dimension barely affect file size
pdf(file = \"/home/douglas/databases/Phylo/P2_structure/analysis3/$nwk1.$nwk2.pdf\", 12+(length(overlapping)*0.1), 6*(length(overlapping)*0.1))

VisualizeMatching(ClusteringInfoDistance, treeA, treeB,  setPar = TRUE, edge.cex = 1.2,value.cex = 1.2)
dev.off()

";

# grDevices::pdf()
# par(usr = c(0, 1, 0, 1))
# text(0.60, -0.17, \"$nwk1.$nwk2\", cex = 0.8, font = 1) 

# jpeg(\"/home/douglas/databases/Phylo/P2_structure/analysis3/$nwk1.$nwk2.jpg\" , width = 1000 , height = 600*(length(overlapping)*0.1))

# does nothing:
# par(mar=c(0, 0, 8, 0))

# text(10,11.2,labels=\"$nwk1.$nwk2\", cex=0.8)
# text(9.5,0.7,labels=\"$nwk1.$nwk2\", cex=0.8)





	system("rm Current_R_Script2");
	open(CURRENT_R2, ">Current_R_Script2") || die "\nerrr 81\n";
	print CURRENT_R2 "$R_script2\n";
	close CURRENT_R2;
	system("/opt/R/4.2.3/bin/R < Current_R_Script2 --vanilla --slave");

############################################################################################

	$index+=1;

# die"";

	};# foreach my $treeID(@all_treeIDs)

};
















