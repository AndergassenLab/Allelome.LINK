#!/usr/bin/perl

## note: sam files have 11 non-optional columns
$col=11;
$SNP_i=$col+5;
$start_col=$col+2;

READS: while (<STDIN>){
  $_=~ s/\n|\r//;
  @infos = split (/\t/, $_);
  $SNP_count=$infos[$col];
  $SNP_p_i=$SNP_i+$SNP_count;
  $SNP_p=$infos[$SNP_p_i];
  $CIGAR = $infos[5];
  $read = $infos[9];
  $read_q = $infos[10];
  $esc = "F";
  $c_m = $infos[$start_col];
  # -1 because c_end indicates the last base of the block before
  $c_end = $c_m-1;
  $M_m=0;
  $CIGAR_n="";
  @C=split(/(?<=[A-Z])/, $CIGAR);
  $C_l=@C;

  CIGAR: for ($x=0;$x<$C_l;$x++){
## note: chop removes last character from original string and returns it
    $c_nr=$C[$x];
    $c_s=chop($c_nr);
    if ($c_s eq "M"){
      $c_end=$c_end+$c_nr;
      SNP: for ($y=0;$y<$SNP_count;$y++){
         $SNP_p = $infos[$SNP_p_i+$y];
        if ( $SNP_p > $c_end ) {
        	last;
        }
        elsif($c_m <= $SNP_p){
          $read_t_nr=($SNP_p-$c_m)+$M_m;
          $read_t=substr($read,0,$read_t_nr);
          $read_q_t=substr($read_q,0,$read_t_nr);
          $CIGAR_n=$CIGAR_n.($read_t_nr-$M_m).$c_s;
          $esc = "T";
          $c_m=$c_m+$c_nr;
          $M_m=$M_m+$c_nr;
          last CIGAR;
        }
      }
      # note: this will only be reached if esc!=T
      $c_m=$c_m+$c_nr;
      $M_m=$M_m+$c_nr;
      $CIGAR_n=$CIGAR_n.$C[$x];
    }
    elsif($c_s eq "N" || $c_s eq "D"){
      $c_m=$c_m+$c_nr;
      $c_end=$c_end+$c_nr;
      $CIGAR_n=$CIGAR_n.$C[$x];
    }
    elsif ($c_s eq "I" || $c_s eq "S"){
          $M_m=$M_m+$c_nr;
          $CIGAR_n=$CIGAR_n.$C[$x];
    }
  }
  if ($esc eq "T"){
    print "$infos[0]\t$infos[1]\t$infos[2]\t$infos[3]\t$infos[4]\t$CIGAR_n\t$infos[6]\t$infos[7]\t$infos[8]\t$read_t\t$read_q_t\n";
  }
}
