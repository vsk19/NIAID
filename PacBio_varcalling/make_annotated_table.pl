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
		}
		print C "\t";
		$a++;
	}
}
close B;

open D, "<$genotypes";
while (<D>){
	chomp;
 	last if m/^$/;
 	@line = split;
	push @variants, $_;
}
close D;

open F, "<$table";
while (<F>){
	chomp;
 	last if m/^$/;
 	@line = split;
	push @info, $_;
}
close F;

my $b=0;
my $c=0;

open E, "<$annotations";
while (<E>){
	chomp;
	@line = split;
 	last if m/^$/;
 	next if ($line[0] =~ m/##/);
 	if ($b == 0) {
 		print C "$_\t$info[$b]\n";
 		$b++;
 	}
 	else {
 		print C "$variants[$c]\t$_\t$info[$b]\n";
 		$c++;
 		$b++;
 	}
 }
 close E;