#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# my ($input,$output);
#GetOptions('input|i:s' => \$input,
#	'output|o:s' => \$output);

my $base_dir = "/01.Short_Read_Mapping"; 
my $output_dir = "./../output/10kb_depth";
system("mkdir -p $output_dir") unless -d $output_dir;

my @folder = glob "$base_dir/batch_*/D*";
foreach my $folder (@folder){
	print "$folder\n";
	my $sample = basename($folder);
	my $input = "$folder/$sample.depth.txt.gz";
	my $input_fh = read_file($input);
	my $output = "$output_dir/$sample.depth.window.txt.gz";
	#print "$sample\n";exit;
	my $output_fh = write_file($output);
	my $win_size = 10001;
	my %depth;
	while (<$input_fh>){
        	chomp;
                my @line = split /\s+/;
                my ($chr,$pos,$depth_val)=@line[0,1,2];
                #print "$chr\t$pos\t$depth\n";exit;
                my $window = int($pos/$win_size);
                $depth{$chr}->{$window}->{depth}+=$depth_val;
                $depth{$chr}->{$window}->{count}++
        }
	for my $chr (sort{$a cmp $b} keys %depth){
        	for my $window (sort{$a <=> $b} keys %{$depth{$chr}}){
                	my $avg_depth = $depth{$chr}->{$window}->{depth}/$depth{$chr}->{$window}->{count};
                	print $output_fh "$sample\t$chr\t$window\t$avg_depth\n"
        	}
	}
	close $input_fh if $input_fh;
    close $output_fh if $output_fh;
}

sub read_file{
	my $file = shift @_;
	my $fh;
	if ($file =~ /\.gz$/){
		open ($fh, "gunzip -c $file |") or die "cannot open pipe to $file";
	}
	else{
		open ($fh, $file) or die "cannot open $file";
	}
	return $fh
}

sub write_file{
        my $file = shift @_;
        my $fh;
        if ($file =~ /\.gz$/){
                open ($fh, "| gzip -c >$file") or die "cannot open pipe to $file";
        }
        else{
                open ($fh, ">$file") or die "cannot open $file";
        }
        return $fh
}
