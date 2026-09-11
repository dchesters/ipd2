use bignum;


# 
# LOG
# 20250510: also use for branchlength removal
# 20250601: option to remove node labels
# 
# 
# 
# 
# 
###########################################################################


$in1 		= $ARGV[0]; # tree
$outfile 	= $ARGV[1]; # tree
$function1	= $ARGV[2]; # 1=convert bls, 2= remove bls, 3=convert bls, remove node labels

if($in1 eq $outfile){die "\nerror, same name ($in1) given for input and output.\n"};
if($in1 =~ /[\w\d]/){}else{die "\nerror, no file name.\n"};
unless($function1 =~ /\d/){die "\nplease specify function.\n"};


###########################################################################
open(IN1, $in1) || die "\nerror 34 , cant open $in1\n";
while(my $line = <IN1>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /\(/)
		{
		$tree = $line;
		$treecount++;
		};
	};
if($treecount >= 1)
	{
	$treelength = length($tree);
	}else{
	die "\nfailed to read tree from file\n";
	};
print "read tree, treecount:$treecount, treelength:$treelength\n";
###########################################################################


while( $tree =~ s/\:(\d\.\d+)[eE](\-\d+)/"\:" .  ($1 * ( 10 ** $2 ))/ge){};
while( $tree =~ s/\:(\d+)[eE](\-\d+)/"\:" .  ($1 * ( 10 ** $2 ))/ge){};


###########################################################################


if($function1 == 2)
	{
	print "removing branchlengths\n";
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	while($tree =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/){}; 
	# remove regular branchlengths: 0.02048
	while($tree =~ s/\:\-*\d+\.\d+//){}; 
	# remove 0 length branchlengths
	while($tree =~ s/\:\d+//){};
	};

if($function1 == 3)
	{
	print "removing node labels\n"; # print "tree:$tree\n";
	while($tree =~ s/\)[0-9\.]{1,40}/)/){};

	};

###########################################################################

print "writing to $outfile\n";
# print "\n$tree\n";
open(OUT, ">$outfile") || die "error 93"; 
print OUT "$tree\n";
close OUT;

###########################################################################


