#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT
my $headers = $ARGV[0];
my $genotypes=$ARGV[1];
my $table=$ARGV[2];
my $annotations=$ARGV[3];
my $out= $ARGV[4];
my @variants = ();
my @info = ();
my $a = 0;
my @line=();
my @samples=();

open C, ">$out";

open B, "<$headers";
while (<B>){
	chomp;
 	@line = split;
 	last if ($a>0);
	if ($line[0] !~ m/##/) {
	 	print C "$line[0]\t$line[1]\t$line[3]\t$line[4]";	
		for ($a = 9; $a < @line; $a++) {
			print C "\t$line[$a].genotypes";
			push @samples, $line[$a];
#			print "$line[$a]\n";
		}
		print C "\t";
		$a++;
	}
}
close B;
#print @samples;

open D, "<$genotypes";
while (<D>){
	chomp;
 	last if m/^$/;
 	@line = split;
	push @variants, $_;
}
close D;

my @adindex=();
my @dpindex=();

open F, "<$table";
while (<F>){
	chomp;
 	last if m/^$/;
 	@line = split;
	push @info, $_;
	if ($line[0] eq 'EVENTLENGTH') {
		my $e=0;
		my $d=0;
		@adindex=();
		@dpindex=();
		for ($d=0; $d<@samples; $d++) {
			my $col=$samples[$d] . '.AD';
			for ($e = 0; $e < @line; $e++) {
				if ($col eq $line[$e]) {
					push @adindex, $e;
				}
			}
		}
	}
}
close F;

my $b=0;
my $c=0;
my $cov=0;
##need to add minDepth, maxDepth, and AF values to the step below

open E, "<$annotations";
while (<E>){
	chomp;
	@line = split;
 	last if m/^$/;
 	next if ($line[0] =~ m/##/);
 	if ($b == 0) {
 		print C "$_\t$info[$b]";
 		my $h=0;
 		for ($h=0; $h<@samples; $h++) {
 			print C "\t$samples[$h]" . '.AF';
 		}
 		print C "\n";
 		$b++;
 	}
 	else {
 		print C "$variants[$c]\t$_\t$info[$b]";
 		my $f=0;
 		for ($f=0; $f<@adindex; $f++) {
 			my $gtcol=4+$f;
# 			print "$gtcol\n";
 			my $num=$adindex[$f];
 			my @fields=split(' ',$info[$b]);
 			my @gts=split(' ',$variants[$c]);
 			if ($fields[$num] ne 'NA') {
 				my @counts=split(',',$fields[$num]);
 				my $gt=$gts[$gtcol];
 				$cov=0;
 				my $g=0;
 				for ($g=0; $g<@counts; $g++) {
 					$cov+=$counts[$g];
# 					print "$counts[$g]\n";
 				}
 				my $af='';
 				if ($cov>0) {
	 				$af=$counts[$gt]/$cov;
	 			}
	 			else {
	 				$af=0;
	 			}
 				print C "\t$af";
 			}
 			else {
 				
 				print C "\tNA";
 			}
 		}
 		print C "\n";	
 		$c++;
 		$b++;
 	}
 }
 close E;