# Function: distribution fitting for normalized reactivity on training RNAs
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Jan-14, 2015)

#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

die "perl $0 dataF mode(train|test) Base_col structF1 structF2 ...\n" unless (@ARGV >= 4);
my ($dataF, $mode, $Base_col, @structF) = @ARGV;

my $bases_with_data = 0;
my $bases_all = 0;

my $struct = getStruct(@structF);    # %struct:  RNA => index => 'base' => base
								     #			 RNA => index => 'pair' => pair
my ($data, $header) = getData($dataF);   # %data:  RNA => index => 'column_num' => data

print "$header\n";
foreach my $RNA (sort keys %$data) {
	foreach my $Index (sort {$a <=> $b} keys %{$struct->{$RNA}}) {
		$bases_all ++;
		if (exists $data->{$RNA}{$Index}) {  # bases with data
			print "$RNA\t$Index";
			for my $i (2 .. $Base_col-2) {
				print "\t$data->{$RNA}{$Index}{$i}";
			}
			print "\t$struct->{$RNA}{$Index}{'base'}";
			print "\t$struct->{$RNA}{$Index}{'pair'}" if ($mode =~ /train/i);
			print "\n";
		} else {  # bases without data
			print "$RNA\t$Index";
			for my $i (2 .. $Base_col-2) {
				print "\tNA";
			}
			print "\t$struct->{$RNA}{$Index}{'base'}";
			print "\t$struct->{$RNA}{$Index}{'pair'}" if ($mode =~ /train/i);
			print "\n";
		}
	}
}

print STDERR "Bases_with_data: $bases_with_data\nAll_bases: $bases_all\n";
exit;


#######################################################################################
sub getData {
	my $file = shift;
	my %alldata = ();
	my $header = "";

	open (IN, $file) or die;
	while (<IN>) {
		chomp;
		if (/^RNA/) {  # Header line start with 'RNA'
			$header = $_;
		} else {
			my @t = split("\t", $_);
			my ($RNA, $Index, $Base) = ($t[0], $t[1], $t[$Base_col-1]);
			if ($Base eq $struct->{$RNA}{$Index}{'base'}) {
				for my $i (2 .. $Base_col-2) {
					$alldata{$RNA}{$Index}{$i} = $t[$i];
					$bases_with_data ++;
				}
			}
		}
	}
	close (IN);
	return(\%alldata, $header);
}

sub getStruct {
	my @files = @_;
	my %struct = ();

	foreach my $file (@files) {
		chomp(my $pre=`basename $file .ct`);
		my $data = getCTfile($file);
		$struct{$pre} = $data;
	}
	return(\%struct);
}

sub getCTfile {
	my $file = shift;
	my %data;

	open (IN, $file) or die;
	while (<IN>) {
		chomp;
		my @t = split(" ", $_);   # columns of CT files are separated by blanks
		if (@t == 6) {
			my ($index, $base, $pair) = ($t[0], $t[1], $t[4]);
			$data{$index}{'base'} = $base;
			$data{$index}{'pair'} = $pair;
		} 
	}
	close (IN);
	return(\%data);
}
