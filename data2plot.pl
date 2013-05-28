#!/usr/bin/perl

use strict;
use warnings;

my $OrgName=$ARGV[0];
my $ContigAnzahl=$ARGV[1];




for(my $i=0;$i<$ContigAnzahl;$i++){

  
	my $j=$i+1;
	system("./plot.R $OrgName $j");
  
}