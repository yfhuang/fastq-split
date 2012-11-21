#!/usr/bin/perl -w

=comment
input file: single fastq, 4 lines for 1 read
output file: paired-end reads - xxx.pe.r1.fastq and xxx.pe.r2.fastq; mates - xxx.mate.r1.fastq and xxx.mate.r2.fastq

=cut

# module included
use strict;
use Getopt::Long;

# variable
# argument
my $fqfile = "";
my $oprefix = ""
my $format = "";
my $help;

# program used variable
my %hash_read = ();
my %ihash_read = ();
my $hash_read_cnt = 0;

my %hashrq = ();

my @pairlabel = ();

# sub function
sub usage{
  printf("%s -fastq input_fastq_file -outprefix prefix_of_ouput_filename\n", $0);
  printf("\n\nversion: 0.1\ncontact for big fix: Yu-Feng Huang, PhD. email: yfh1202\@sinica.edu.tw\n");
  exit(0);
}

sub getReadid4Cmp{
  open(my $fh, "<$fqfile") or die "file ($fqfile) not found!!\n";
  
  while(my $line = <$fh>){
    chomp($line);
    my $rid = $line;
    
    my $seq = <$fh>;
    chomp($seq);
    
    $line = <$fh>;
    chomp($line);
    
    $qual = <$fh>;
    chomp($qual);
    
    $rid = substr($rid, 1);	# get readid
    
    my @rbuf = split("/", $rid);
    my $pairid = $rbuf[0];
    my $mateid = $rbuf[1];
    
    my $read = "@$rid\n$seq\n+\n$qual\n";
    
    $hashrq{$rid} = $read;
    
    if(defined($hash_read{$pairid}) && $hash_read{$pairid}){
      ;
    }else{
      $hash_read{$pairid} = $hash_read_cnt;
      $ihash_read{$hash_read_cnt} = $pairid;
      $pairlabel[$hash_read_cnt][0] = 0;
      $pairlabel[$hash_read_cnt][1] = 0;
      $hash_read_cnt++;
    }
    
    $pairlabel[$hash_read{$pairid}][$mateid] = 1;
  }
  
  close($fh);
}

sub fastqSplit{
  my $r1_uniq = $oprefix . ".r1.uniq.fastq";
  my $r2_uniq = $oprefix . ".r2.uniq.fastq";
  my $r1_pair = $oprefix . ".r1.pair.fastq";
  my $r2_pair = $oprefix . ".r2.pair.fastq";
  
  open(my $r1u_fh, ">$r1_uniq");
  open(my $r2u_fh, ">$r2_uniq");
  open(my $r1p_fh, ">$r1_pair");
  open(my $r2p_fh, ">$r2_pair");

  for(my $i = 0; $i < $hash_read_cnt; $i++){
    my $pair1id = $ihash_read{$i} . "/1";
    my $pair2id = $ihash_read{$i} . "/2";
    
    if($pairlabel[$i][0] && $pairlabel[$i][1]){
      printf $r1p_fh $hashrq{$pair1id};
      printf $r2p_fh $hashrq{$pair2id};
    }else{
      if($pairlabel[$i][0]){
        printf $r1u_fh $hashrq{$pair1id};
      }elsif($pairlabel[$i][1]){
      printf $r2u_fh $hashrq{$pair1id};
      }
    }
  }
  
  close($r1u_fh);
  close($r2u_fh);
  close($r1p_fh);
  close($r2p_fh);
}
=comment
1. get readid from single fastq
2. identify read1 and read2 and then put them in hash and mark them
=cut

my $result = GetOptions("fastq=s" => \$fqfile,
                        "outprefix=s" => \$oprefix,
                        "help" => \$help);
if($help){
  &usage();
}else{
  if(-f $fqfile && length($oprefix) > 0){
    &getReadid4Cmp();
    &fastqSplit();
  }else{
    print "Error: please check the fastq file or output prefix name.\n\n";
    &usage();
  }
}
