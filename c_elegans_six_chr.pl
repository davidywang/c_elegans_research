
#!/usr/bin/perl
# project_six_chr.pl by David Wang
use strict; use warnings FATAL => 'all';

# usage statement
die "usage: project_six_chr.pl <text file> <min>\n" unless @ARGV == 2;

my ($FILE, $min) = @ARGV;

# open file
open (my $IN, "<", $FILE) or die "file cannot be opened";

my $list;

my ($chr1, $chr2, $chr3, $chr4, $chr5, $chrx) = ("I", "II", "III", "IV", "V", "X");

# gathers all mutation info and stores in @all_mut
while (<$IN>) {
	$list .= $_;
}

my @strain_pairs = split(/\n{2,}/, $list);
my $count;
foreach my $strain_pair (@strain_pairs) {
	my @mut = split("\n+", $strain_pair);
	my @chr;
	my %chr_count;

	
	foreach my $mut (@mut) {
		my @info = split("\t+", $mut);
		
		my $chr = $info[0];	
		push @chr, $chr;
		$chr_count{$chr}++;
	}
	
	next unless grep(/^$chr1$/, @chr);
	next unless grep(/^$chr2$/, @chr);
	next unless grep(/^$chr3$/, @chr);
	next unless grep(/^$chr4$/, @chr);
	next unless grep(/^$chr5$/, @chr);
	next unless grep(/^$chrx$/, @chr);
	
	next if !defined $chr_count{$chr1};
	next if !defined $chr_count{$chr2};
	next if !defined $chr_count{$chr3};
	next if !defined $chr_count{$chr4};
	next if !defined $chr_count{$chr5};
	next if !defined $chr_count{$chrx};
	
	next unless $chr_count{$chr1} >= $min;
	next unless $chr_count{$chr2} >= $min;
	next unless $chr_count{$chr3} >= $min;
	next unless $chr_count{$chr4} >= $min;
	next unless $chr_count{$chr5} >= $min;
	next unless $chr_count{$chrx} >= $min;
	
	foreach my $mut (@mut) {
		last if !defined $mut;		
		my @info = split("\t+", $mut);		
		my $start1 = $info[4];
		my $stop2 = $info[10];
		
		if (defined $start1 and defined $stop2) {
			my $left = $start1 - 500;
			my $right = $stop2 + 500;
		
			my $chr = $info[0];
		
			my $entropy1 = `get_entropy.pl c_elegans.PRJNA13758.WS254.genomic.fa $chr $left $start1`;
			my $entropy2 = `get_entropy.pl c_elegans.PRJNA13758.WS254.genomic.fa $chr $stop2 $right`;
		
			print "$mut\t$entropy1\t$entropy2\n"
		}
	}
	print "\n";
	$count++;
}
close $IN;

print $count;