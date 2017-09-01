#!/usr/bin/perl
# project_final.pl by David Wang
use strict; use warnings FATAL => 'all';

# usage statement
die "usage: project_final.pl <gff file> <minimum length> <maximum length> <minimum length difference> <max range> <entropy range>\n" unless @ARGV == 6;

my ($FILE, $min, $max, $min_dif, $max_range, $entropy_range) = @ARGV;

# open file
open (my $IN, "<", $FILE) or die "file cannot be opened";

# all mutations
my @all_mut;

# selected mutations
my @sel_mut;

# gathers all mutation info and stores in @all_mut
while (<$IN>) {
	my @info = split("\t+");	
	if ($info[2] eq "insertion_site" or $info[2] eq "deletion") {		
		my $ref = $info[0];
		my $method = $info[2];
		my $start = $info[3];
		my $stop = $info[4];
		my $group = $info[8];
	
		my $length;
		my @strains;
		my $variation;
	
		# gets length of mutation and strains
		if ($method eq "insertion_site") {
			my @group_items = split(";+", $group);							
			for (my $i = 0; $i < scalar(@group_items); $i++) {
				if ($group_items[$i] =~ m/^ Insertion/) {
					my @insert_split = split("\"+", $group_items[$i]);
					my $insertion = $insert_split[1];
					$length = length($insertion);
				} elsif ($group_items[$i] =~ m/^ Strain/) {
					my @strain_split = split("\"+", $group_items[$i]);
					my $strain = $strain_split[1];
					push (@strains, $strain);
				} elsif ($group_items[$i] =~ m/^Variation/) {
					my @variation_split = split("\"+", $group_items[$i]);
					$variation = $variation_split[1];
				}
			}		
		} else {
			$length = $stop - $start + 1;
			
			my @group_items = split(";+", $group);							
			for (my $i = 0; $i < scalar(@group_items); $i++) {
				if ($group_items[$i] =~ m/^ Strain/) {
					my @strain_split = split("\"+", $group_items[$i]);
					my $strain = $strain_split[1];
					push (@strains, $strain);
				} elsif ($group_items[$i] =~ m/^Variation/) {
					my @variation_split = split("\"+", $group_items[$i]);
					$variation = $variation_split[1];
				}
			}		
		}
	
		# hash of mutation elements	
		if (defined $length) {
			push @all_mut, {
				ref => $ref,
				method => $method,
				start => $start,
				stop => $stop,
				len => $length,
				var => $variation,
				strains => \@strains,
			};
		}
	}
}

close $IN;

my @strain_pairs;
my @var_pairs;

# filters pairs of mutations
for (my $i = 0; $i < @all_mut; $i++) {
	for (my $j = $i + 1; $j < @all_mut; $j++) {
		my $mut1 = $all_mut[$i];
		my $mut2 = $all_mut[$j];
		
		my $ref1 = $mut1->{ref};
		my $ref2 = $mut2->{ref};
		
		my $len1 = $mut1->{len};
		my $len2 = $mut2->{len};
		
		my $start1 = $mut1->{start};
		my $stop1 = $mut1->{stop};
		my $start2 = $mut2->{start};		
		my $stop2 = $mut2->{stop};
		
		my $range = $stop2 - $start1;
		
		my $var1 = $mut1->{var};
		my $var2 = $mut2->{var};
		my $var_pair = "$var1\t$var2";
		
		my @strains1 = @{$mut1->{strains}};
		my @strains2 = @{$mut2->{strains}};
		
		last if $ref1 ne $ref2;
		last if $range > $max_range;

		next if $len1 < $min;
		next if $len1 > $max;
		
		next if $len2 < $min;
		
		next if $len2 > $max;
		
		next if abs($len2 - $len1) < $min_dif;
		
		next if $var1 eq $var2;

		next if grep(/^$var_pair$/, @var_pairs);  
		next if scalar(@strains1) == 0;
		next if scalar(@strains2) == 0;
		
		foreach my $strain1 (@strains1) {
			foreach my $strain2 (@strains2) {
				if ($strain1 ne $strain2) {
					my $strain_pair = "$strain1\t$strain2";
					push (@strain_pairs, $strain_pair);
				}	
			}						
		}
				
		push @sel_mut, $mut1;
		push @sel_mut, $mut2;
		push @var_pairs, $var_pair;
	}
}

print scalar(@sel_mut);
die;
# list of strains with mating plugs
my @mp = qw(AB2 AB3 AB4 CB3196 CB3197 CB3198 CB3199 CB4853 CB4854 CB4855 CB4856 CB4857 CB4858 DR1345 DR1350 ED3005 ED3017 ED3040 ED3042 ED3043 ED3048 ED3049 ED3051 ED3054 ED3063 ED3066 ED3072 ED3073 ED3075 ED3077 EG4348 EG4349 EG4350 EG4351 EG4352 JU1088 JU1171 JU1172 JU258 JU262 JU263 JU299 JU300 JU301 JU302 JU303 JU305 JU306 JU307 JU309 JU310 JU322 JU323 JU342 JU345 JU346 JU347 JU360 JU362 JU363 JU364 JU366 JU368 JU369 JU370 JU395 JU396 JU397 JU398 JU401 JU402 JU438 JU531 JU533 JU561 JU563 JU642 KR314 MY1 MY10 MY11 MY12 MY13 MY14 MY15 MY16 MY17 MY18 MY19 MY2 MY20 MY21 MY22 MY23 MY3 MY4 MY5 MY6 MY7 MY8 MY9 PB303 PB306 PS2025 PX174 RC301);

my @mp_strain_pairs;

# gets rid of non mating plugs
for (my $i = 0; $i < @strain_pairs; $i++) {
	my $strain_pair = $strain_pairs[$i];	
	my @split = split("\t", $strain_pair);
	my $first = $split[0];
	my $second = $split[1];
	if (grep(/^$first$/, @mp) and grep(/^$second$/, @mp)) {
		push (@mp_strain_pairs, $strain_pair);
	}	
}

my @sort_strain_pairs;
my %strain_pair_counts;

# gets pairs which appear multiple times
for (my $i = 0; $i < scalar(@mp_strain_pairs); $i++) {
	my $strain_pair = $mp_strain_pairs[$i];
			
	my @reverse_split = split("\t", $strain_pair);
	my $reverse_1 = $reverse_split[1];
	my $reverse_2 = $reverse_split[0];
	my $reverse_strain = "$reverse_1\t$reverse_2";
	
	splice(@mp_strain_pairs, $i, 1);	
	if (grep(/^$strain_pair$/, @mp_strain_pairs) or grep(/^$reverse_strain$/, @mp_strain_pairs)) {
		push (@sort_strain_pairs, $strain_pair);
		$strain_pair_counts{$strain_pair}++;
	}	
	splice(@mp_strain_pairs, $i, 0, $strain_pair);
}

my @reverse_strain_pairs;
my @nr_strain_pairs;
my %rev_count;

# gets rid of reverse strings
for (my $i = 0; $i < scalar(@sort_strain_pairs); $i++) {
	my $strain_pair = $sort_strain_pairs[$i];
	#my $reverse_strain = rev($strain_pair);

	my @reverse_split = split("\t", $strain_pair);
	my $reverse_1 = $reverse_split[1];
	my $reverse_2 = $reverse_split[0];
	my $reverse_strain = "$reverse_1\t$reverse_2";

	if (grep(/^$reverse_strain$/, @mp_strain_pairs)) {
		if (!exists $rev_count{$reverse_strain}) {		
			push (@nr_strain_pairs, $strain_pair);
			$rev_count{$strain_pair} = "used";	
		}	
	} else {
		push (@nr_strain_pairs, $strain_pair);
	}		
}

my @uniq_strain_pairs = uniq(@nr_strain_pairs);

foreach my $strain_pair (@uniq_strain_pairs) {
	my $count = $strain_pair_counts{$strain_pair};
		
	#if ($count >= 6) {
	my @strain_pair = split("\t", $strain_pair);	
	my $strain1 = $strain_pair[0];
	my $strain2 = $strain_pair[1];
	print "$strain1\t$strain2\n";

	for (my $i = 0; $i < @sel_mut - 1; $i += 2) {
		my $mut1 = $sel_mut[$i];
		my $mut2 = $sel_mut[$i + 1];
			
		my $ref = $mut1->{ref};			
		my @chr_split = split("_", $ref);
		my $chr = $chr_split[1];
					
		my $start1 = $mut1->{start};
		my $stop1 = $mut1->{stop};
		my $start2 = $mut2->{start};
		my $stop2 = $mut2->{stop};
			
		my $range = $stop2 - $start1;
			
		my $len1 = $mut1->{len};
		my $len2 = $mut2->{len};
			
		my $var1 = $mut1->{var};
		my $var2 = $mut2->{var};
	
		my @strains1 = @{$mut1->{strains}};
		my @strains2 = @{$mut2->{strains}};
			
		my $left = $start1 - $entropy_range;
		my $right = $stop2 + $entropy_range;
		
		#my $entropy_1 = `get_entropy.pl c_elegans.PRJNA13758.WS254.genomic.fa $num $left $start_1`;
		#my $entropy_2 = `get_entropy.pl c_elegans.PRJNA13758.WS254.genomic.fa $num $stop_2 $right`;
			
		my $line1 = "$chr\t$range\t$strain1\t$var1\t$start1\t$stop1\t$len1\t$strain2\t$var2\t$start2\t$stop2\t$len2\n";
		my $line2 = "$chr\t$range\t$strain2\t$var2\t$start2\t$stop2\t$len2\t$strain1\t$var1\t$start1\t$stop1\t$len1\n";

		if (grep(/^$strain1$/, @strains1) and grep(/^$strain2$/, @strains2) and !grep(/^$strain1$/, @strains2) 
		or grep(/^$strain1$/, @strains1) and grep(/^$strain2$/, @strains2) and !grep(/^$strain2$/, @strains1)) {				
			print "$line1";
		} elsif (grep(/^$strain1$/, @strains2) and grep(/^$strain2$/, @strains1) and !grep(/^$strain1$/, @strains1) 
		or grep(/^$strain1$/, @strains2) and grep(/^$strain2$/, @strains1) and !grep(/^$strain2$/, @strains2)) {
			print "$line2";				
		}
	}
	print "\n";				
	#}
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub rev {
	my @reverse_split = split("\t", $_);
	my $reverse_1 = $reverse_split[1];
	my $reverse_2 = $reverse_split[0];
	return "$reverse_1\t$reverse_2";

}























