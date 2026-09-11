
# 
# study specific script used for work in submission, provided for review
# 

$newick_list	= $ARGV[0];
$synonyms_list	= $ARGV[1];
# $newicks_path = "/home/douglas/databases/Phylo/Processed_Newick/";
$newicks_path = "./Newicks.202605/";


$just_copy_newicks = 0;










print "file listing newicks:$newick_list\n";
print "path to newicks:$newicks_path\n";


if($just_copy_newicks == 1)
	{}else{

##############################################################################
open(IN5, "$synonyms_list") || die "\nerror cant open synonyms_list\n";
print "\nopened synonyms_list:$synonyms_list\n";
while(my $line = <IN5>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "$line\n";
	$line =~ s/ +$//; # $line =~ s/ +/_/;
	if($line =~ /^([^\t]+)\t(\S.+)/)
		{
		my $canonical_species = $1; my $current_synonyms = $2;$current_synonyms =~ s/\s+$//;my @split_synonyms = split /\s+/, $current_synonyms;
		foreach my $syn(@split_synonyms)
			{
			# print "\tif name $syn is found among termials, change to $canonical_species\n";
			if($syn =~ /\w+_\w+/ && $canonical_species =~ /\w+_\w+/)
				{
				$switch_synonym{$syn} = $canonical_species;
				if($syn eq $canonical_species){die "\nerror 26, why identical $canonical_species\n"};
				}else{	
				print "error 25:$line\n";	
				};	
			};
		$sp_with_synonyms++;
		};
	};
close IN5;
print "sp_with_synonyms:$sp_with_synonyms\n\n";
##############################################################################

	};

open(IN1, "$newick_list") || die "\nerror cant open file listing newicks\n";
open(OUT, ">release_LOG") || die "";close OUT;

my $index = 0;
while(my $line = <IN1>)
	{
	$line =~ s/\n//;$line =~ s/\r//; # print "$line\n";
	$index++;if($index =~ /0$/){print "newick list line $index\n"};
	if($line =~ /^(\S+\.nwk)$/)
		{
		my $current_newick = $1; # print "current_newick:$current_newick\n";

	#	my $command = "perl ~/usr_scripts/newick/convert_scientific_notation_branchlengths.pl $newicks_path$current_newick $current_newick 2 >>release_LOG";
		my $command = "perl convert_scientific_notation_branchlengths.pl $newicks_path$current_newick $current_newick 2";

	#	if($current_newick =~ /Blaimer2023Hymenoptera/)
	#		{
		print "command:$command\n";
	#		};
		system($command);
		$newicks_processed++;

		if (-e $current_newick) {
		print "File exists!\n"
		}else{die "\noutfile missing\n"};


		if($just_copy_newicks == 1)
			{}else{
		###############################
		fix_synonyms($current_newick);#
		###############################

		};
		
		# die "";
		}else{
		print "warning, no newick:$line\n";	
		open(OUT, ">>release_LOG") || die "";print OUT "warning, no newick:$line\n"; close OUT;

		};
	};
close IN1;

print "


newicks_processed:$newicks_processed
newicks_with_synonyms:$newicks_with_synonyms
sucessful_synonym_switches:$sucessful_synonym_switches
synonym_not_switched:$synonym_not_switched
synonym_duplicates:$synonym_duplicates
canon_matches_skipped:$canon_matches_skipped

FIN.
";

################################################################################################################################

sub fix_synonyms
{
my $file = shift;

my $newick3="";
open(IN3, $file) || die "error 56.";
while(my $line = <IN3>)
	{
	$newick3 .= $line;
	};
close IN3;
# print "newick3:$newick3\n";
my $newick3copy=$newick3;
my @current_newick_terminals=();
my %current_newick_terminals2=();
my %new_canon_sp=();

while($newick3copy =~ s/([\(\,])([A-Za-z_]+)([\)\,\:])/$1$3/){my $terminal = $2;push @current_newick_terminals, $terminal;$current_newick_terminals2{$terminal}=1};
# print "current_newick_terminals:@current_newick_terminals\n";
my $log_string = "";my $current_newick_corrected=0;

foreach my $terminal(@current_newick_terminals)
	{
	if($switch_synonym{$terminal} =~ /\w/)
		{
		my $switch = $switch_synonym{$terminal}; # print "FOUND synonym $terminal, will switch for $switch\n";
		$log_string .= "in newick $file, switch synonym $terminal for canon $switch\n";

		if($current_newick_terminals2{$switch}==1)# also need to check the canon species names is not contained in the tree,
			{					# otherwise would be creating duplclicate
			$synonym_duplicates++;
			$log_string .= "\tcanon species already contained, note, conversion will duplicate\n";
			}elsif($new_canon_sp{$switch}==1)
			{# another way duplicates could be formed, if 2 different synonyms of the same canon species, check for this
			$canon_matches_skipped++;print "\tcanon duplicate:$switch\n";
			}else{
			if($newick3 =~ s/([\(\,])$terminal([\)\,\:])/$1$switch$2/)
				{
				$sucessful_synonym_switches++;$new_canon_sp{$switch}=1;
				}else{$synonym_not_switched++};
			$current_newick_corrected=1;
			};

		};	
	};

if($current_newick_corrected == 1)
	{
	system("rm $file");
	open(OUT5, ">$file") || die "error 136.";print OUT5 "$newick3";close OUT5;

	};

if(length($log_string) >= 3)
	{
 open(OUT, ">>release_LOG") || die "";print OUT $log_string; close OUT;
$newicks_with_synonyms++;
	};

};

################################################################################################################################




