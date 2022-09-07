#!/usr/bin/perl

use strict;
use warnings;

my $hbv_reads = shift;
my $human_reads = shift;
my $hbv_out_endo = shift;
my $human_out_endo = shift;
my $hbv_out_exo = shift;
my $human_out_exo = shift;

die "perl $0 <hbv sam file> <human sam file> <hbv_out_Endogenous> <human_out_Endogenous> <hbv_out_exogenous> <human_out_exogenous>\n" unless ($hbv_reads && $human_reads && $hbv_out_endo && $human_out_endo && $hbv_out_exo && $hbv_out_exo);

open HBV, $hbv_reads or die $!;

my %mutual;

my %count1;
my %tmp_hash;
my %hbv_type;
my %record_hbv;

while(<HBV>){
	chomp;
	my @a = split;
	if($a[0]=~/@/){next;}
	my $hbv_cigar = $a[5];
	my $mapq = $a[4];
	my $id = $a[0];
	my $hbv_pos = $a[3];
	my $hbv_ref = $a[2];
	my $hbv_fqseq = $a[9];
	if($hbv_cigar =~ /\*/){next;}
	my $type = $hbv_ref;                                     ## test at 10:21 2012-3-27
	my $info = "$hbv_cigar\t$mapq\t$hbv_ref\t$hbv_pos\t$hbv_fqseq";
	if(not exists $mutual{$id}){
		$mutual{$id} = $info;
		$count1{$id} = 1;
	}else{
		if($mapq > (split /\s+/,$mutual{$id})[1]){
			$mutual{$id} = $info;
			if(exists $record_hbv{$id}){
				delete $tmp_hash{$id};
				$count1{$id} = 1;
			}
		}elsif($mapq == (split /\s+/,$mutual{$id})[1]){
				$mutual{"$id\_$count1{$id}"} = $info;
				push @{$tmp_hash{$id}}, "$id\_$count1{$id}";
				$record_hbv{$id} = 1;
				$count1{$id}++;
		}
	}
}

close HBV;


open HUMAN, $human_reads or die $!;


my %print_hbv;
my %print_human;

my %count2_endo;	## modify at 16:45 2012-10-09
my %count2_exo;       ## modify at 14:52 2012-10-09
my %record_human; my %record_endo;

while(<HUMAN>){
	chomp;
	my @a = split;
	if($a[0]=~/@/){next;}
	my $id = $a[0];
	if($a[5] =~ /\*/){next;}
	if (exists $mutual{$id}){
		my @array_inter;
		if (exists $tmp_hash{$id}){
			@array_inter = @{$tmp_hash{$id}};
			push @array_inter, $id;
		}else{
			@array_inter = ($id);
		}
###################################
		for my $key_id(@array_inter){
####################################
		my $human_cigar = $a[5];
		my $hbv_cigar = (split /\s+/, $mutual{$key_id})[0];
		my $hbv_mapq = (split /\s+/, $mutual{$key_id})[1];
		my $hbv_pos = (split /\s+/, $mutual{$key_id})[3];
		my $hbv_ref = (split /\s+/, $mutual{$key_id})[2];
		my $hbv_fqseq = (split /\s+/, $mutual{$key_id})[4];
		my $human_mapq = $a[4];
		my $human_pos = $a[3];
		my $human_ref = $a[2];
		my $human_fqseq = $a[9];
		my $m = 0; my $mhbv = 0 ; my $mhuman = 0;
		$mhbv += $1 while($hbv_cigar =~ /(\d+)M/g);
		$mhuman += $1 while($human_cigar =~ /(\d+)M/g);
		
		next if ($hbv_cigar!~/S/ || $mhbv<30);
		next if ($human_cigar!~/S/ || $mhuman<30);
		
		my @cigar_split; my @hbv_array_pos; my @hbv_array_pos_s;
		if($hbv_fqseq ne $human_fqseq){
			while($hbv_cigar =~ /(\d+[MSID])/g){
				push @cigar_split,$1;
			}
			my $hbv_cigar_reverse = join ("", reverse(@cigar_split));
			@hbv_array_pos = &treat_cigar_M($hbv_cigar_reverse);
			@hbv_array_pos_s = &treat_cigar_S($hbv_cigar_reverse);
		}else{
			@hbv_array_pos = &treat_cigar_M($hbv_cigar);
			@hbv_array_pos_s = &treat_cigar_S($hbv_cigar);
		}	
		my $flag = 0;
		my $flag_test = 0;
		my @human_array_pos = &treat_cigar_M($human_cigar);
		my @human_array_pos_s = &treat_cigar_S($human_cigar);
		for my $i1(@hbv_array_pos){
			my @hbvpos_tmp = split /\t/, $i1;
			for my $i2(@human_array_pos){
				my @humanpos_tmp = split /\t/, $i2;
				my $dvalue1 = $humanpos_tmp[1] - $hbvpos_tmp[0];
				my $dvalue2 = $hbvpos_tmp[1] - $humanpos_tmp[0];
				if(($dvalue1 > -1) && ($dvalue2 > -1)){                     ## actually all endogenous integration
					$flag_test = 1;
				}
				if (($dvalue1 > 10) && ($dvalue2 > 10)){
					$flag = 1;												## all endogenous integration
					push @{$record_endo{$id}}, "$human_mapq\t$hbv_mapq";
					last;
				}
				if(($dvalue1 >= 20) && ($dvalue2 >= 20)){
					$flag = 2;												## definite endogenous integration
					last;
				}elsif(($dvalue1<20 && $dvalue1>5) && ($dvalue2<20 && $dvalue2>5)){
					$flag = 3;												## undefinite endogenous integration
					last;
				}
			}
		}
		my $flags = 0;	
        for my $i1(@hbv_array_pos_s){
            my @hbvpos_tmp = split /\t/, $i1;
            for my $i2(@human_array_pos_s){
                my @humanpos_tmp = split /\t/, $i2;
                my $dvalue1 = $humanpos_tmp[1] - $hbvpos_tmp[0];
                my $dvalue2 = $hbvpos_tmp[1] - $humanpos_tmp[0];
#				print "\$dvalue1 = \$humanpos_tmp[1] - \$hbvpos_tmp[0]; === $humanpos_tmp[1] - $hbvpos_tmp[0] = $dvalue1\n";
#				print "\$dvalue2 = \$hbvpos_tmp[1] - \$humanpos_tmp[0]; === $hbvpos_tmp[1] - $humanpos_tmp[0] = $dvalue2\n";
                if (($dvalue1 > 10) && ($dvalue2 > 10)){
                    $flags = 1;                                              ## all endogenous integration
                    push @{$record_endo{$id}}, "$human_mapq\t$hbv_mapq";
                    last;
                }
                if(($dvalue1 >= 20) && ($dvalue2 >= 20)){
                    $flags = 2;                                              ## definite endogenous integration
                    last;
                }elsif(($dvalue1<20 && $dvalue1>5) && ($dvalue2<20 && $dvalue2>5)){
                    $flags = 3;                                              ## undefinite endogenous integration
                    last;
                }
            }
        }
		
		if($flag_test == 1){
			if(exists $tmp_hash{$id}){
                $print_hbv{"exo"}{$key_id} = "repeat_reads_$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
            }else{
                $print_hbv{"exo"}{$key_id} = "$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
            }

            if(not exists $print_human{"exo"}{$id}){
                $print_human{"exo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
                $count2_exo{$id} = 1;
            }else{
                next if ($print_human{"exo"}{$id} eq "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq");
                if($human_mapq > (split /\s+/, $print_human{"exo"}{$id})[4]){
                    for my $ppp(keys %{$print_human{"exo"}}){
                        if($ppp=~/_\d+$/ && $ppp=~/$id/){
                            delete $print_human{"exo"}{$ppp};
                            delete $record_human{$id};
                            last;
                        }
                    }
                    $print_human{"exo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
                }elsif($human_mapq == (split /\s+/, $print_human{"exo"}{$id})[4]){
                    if($count2_exo{$id} == 1){
                        my @tmp_item = (split /\t/, $print_human{"exo"}{$id});
                        $print_human{"exo"}{"$id\_0"} = "repeat_reads_$id\t$tmp_item[1]\t$tmp_item[2]\t$tmp_item[3]\t$tmp_item[4]";
                    }
                    $print_human{"exo"}{"$id\_$count2_exo{$id}"} = "repeat_reads_$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";       ## modify at 14:52 2011-11-06; at 11:55 2012-01-10
                    $count2_exo{$id}++;
                    $record_human{$id} = 1;
                }
            }
		}
		
		if($flag != 1 && $flags != 1){                                ## all exogenous integration
			if(exists $tmp_hash{$id}){
				$print_hbv{"exo"}{$key_id} = "repeat_reads_$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
			}else{
				$print_hbv{"exo"}{$key_id} = "$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
			}
				
			if(not exists $print_human{"exo"}{$id}){
				$print_human{"exo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
				$count2_exo{$id} = 1;
			}else{
				next if ($print_human{"exo"}{$id} eq "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq");
				if($human_mapq > (split /\s+/, $print_human{"exo"}{$id})[4]){
					for my $ppp(keys %{$print_human{"exo"}}){
						if($ppp=~/_\d+$/ && $ppp=~/$id/){
							delete $print_human{"exo"}{$ppp};
							delete $record_human{$id};
							last;
						}
					}
					$print_human{"exo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
				}elsif($human_mapq == (split /\s+/, $print_human{"exo"}{$id})[4]){
					if($count2_exo{$id} == 1){
						my @tmp_item = (split /\t/, $print_human{"exo"}{$id});
						$print_human{"exo"}{"$id\_0"} = "repeat_reads_$id\t$tmp_item[1]\t$tmp_item[2]\t$tmp_item[3]\t$tmp_item[4]";
					}
					$print_human{"exo"}{"$id\_$count2_exo{$id}"} = "repeat_reads_$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";		## modify at 14:52 2011-11-06; at 11:55 2012-01-10
					$count2_exo{$id}++;
					$record_human{$id} = 1;
				}
			}
		}else{
			if(exists $tmp_hash{$id}){
                $print_hbv{"endo"}{$key_id} = "repeat_reads_$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
            }else{
                $print_hbv{"endo"}{$key_id} = "$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$hbv_mapq";
            }

            if(not exists $print_human{"endo"}{$id}){
                $print_human{"endo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
                $count2_endo{$id} = 1;
            }else{
				next if ($print_human{"endo"}{$id} eq "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq");
                if($human_mapq > (split /\s+/, $print_human{"endo"}{$id})[4]){
                    for my $ppp(keys %{$print_human{"endo"}}){
                        if($ppp=~/_\d+$/ && $ppp=~/$id/){
                            delete $print_human{"endo"}{$ppp};
                            delete $record_human{$id};
							last;
                        }
                    }
                    $print_human{"endo"}{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
                }elsif($human_mapq == (split /\s+/, $print_human{"endo"}{$id})[4]){
                    if($count2_endo{$id} == 1){
                        my @tmp_item = (split /\t/, $print_human{"endo"}{$id});
                        $print_human{"endo"}{"$id\_0"} = "repeat_reads_$id\t$tmp_item[1]\t$tmp_item[2]\t$tmp_item[3]\t$tmp_item[4]";
                    }
                    $print_human{"endo"}{"$id\_$count2_endo{$id}"} = "repeat_reads_$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";       ## modify at 14:52 2011-11-06; at 11:55 2012-01-10
                    $count2_endo{$id}++;
                    $record_human{$id} = 1;
                }
            }
		}
	}
############################
	}
###########################
}

close HUMAN;

for my $i(keys %record_human){
	delete $print_human{$i};
}

#for my $i(keys %record_hbv){
#	delete $print_hbv{$i};
#}

open HBV_OUT_EXO, ">$hbv_out_exo" or die $!;
open HUMAN_OUT_EXO, ">$human_out_exo" or die $!;
open HBV_OUT_ENDO, ">$hbv_out_endo" or die $!;
open HUMAN_OUT_ENDO, ">$human_out_endo" or die $!;

print HBV_OUT_EXO "reads_id\t\t\t\t\tref\tcigar\tpos\n";
print HBV_OUT_ENDO "reads_id\t\t\t\t\tref\tcigar\tpos\n";
print HUMAN_OUT_EXO "reads_id\t\t\t\t\tref\tcigar\tpos\n";
print HUMAN_OUT_ENDO "reads_id\t\t\t\t\tref\tcigar\tpos\n";

my %record_last;
my %record_tmp;
for my $key1(keys %print_human){
	for my $key2(keys %{$print_human{$key1}}){
		my @a = (split /\s+/, $print_human{$key1}{$key2})[0..3];
		my $p = join "\t", @a;
		print HUMAN_OUT_ENDO "$p\n" if($key1 eq "endo");
		next if $key1 eq "endo";		
		my $flag3 = 0;
		my $human_mapq = (split /\t/, $print_human{$key1}{$key2})[-1];
		my $hbv_mapq;
		if(exists $print_hbv{$key1}{$key2}){
			$hbv_mapq = (split /\t/, $print_hbv{$key1}{$key2})[-1];
		}else{
			for my $k1(keys %print_hbv){
				next if $k1 ne $key1;
				for my $k2(keys %{$print_hbv{$k1}}){
					$k2 =~ s/_\d+$//;
					if($key2 =~ /$k2/){
						$hbv_mapq = (split /\t/, $print_hbv{$key1}{$k2})[-1];
						last;
					}
				}
			}
		}
        for my $i(keys %record_endo){
#			next if $key1 eq "endo";
            if($key2 =~ /$i/){
#				print "$key1\t$p\t$human_mapq\n" if $i eq "chr2_unmapFCD14WRACXX:8:2210:9666:33772#TATCCAGA/1";
				for my $kk(@{$record_endo{$i}}){
					my $thuman_mapq = (split /\t/, $kk)[0];
					my $thbv_mapq = (split /\t/, $kk)[1];
#					print "hbv_mapq\t$hbv_mapq\n";
					$flag3 = 1 if (($thuman_mapq+$thbv_mapq) >= ($human_mapq+$hbv_mapq));
					$record_tmp{$key2} = 1;
#					print "$thuman_mapq\t\t$key1\t$p\n" if (($thuman_mapq+$thbv_mapq) >= ($human_mapq+$hbv_mapq));
				}
            }
			last if (exists $record_tmp{$key2} && $record_tmp{$key2} == 1);
        }
        next if $flag3 == 1;
		if($key1 eq "exo"){
			print HUMAN_OUT_EXO "$p\n";
#			print "HUMAN\t$p\t$human_mapq\n" if($key1 eq "exo");
			$record_last{$key2} = 1;
		}
	}
}

for my $key1(keys %print_hbv){
	for my $key2(keys %{$print_hbv{$key1}}){
		my @a = (split /\t/, $print_hbv{$key1}{$key2})[0..3];
		my $p = join "\t", @a;
		print HBV_OUT_ENDO "$p\n" if ($key1 eq "endo");
		my $flag4 = 0;
		next if $key1 eq "endo";
		for my $i(keys %record_last){
			$i =~ s/_\d+$//;
			if($key2 =~ /$i/){
				$flag4 = 1;
				last;
			}
		}
		if($flag4 == 1){
			print HBV_OUT_EXO "$p\n" if ($key1 eq "exo");
		}
	}
}

close HBV_OUT_EXO;
close HUMAN_OUT_EXO;
close HBV_OUT_ENDO;
close HUMAN_OUT_ENDO;

sub treat_cigar_M{
	my ($cigar) = @_;
#	print $cigar,"\n";
	my $index = 0; my $read_pos = 1; my @array_pos;
	while($cigar =~ /((\d+)[MSI])/g){
           my $len = $2;
           if($1 =~ /M/){
                next if $len < 5;
                my $m_spos = $index + $read_pos;
                my $m_epos = $index + $read_pos + $len - 1 ;
                push @array_pos, "$m_spos\t$m_epos";
#				print "$m_spos\t$m_epos\n";
			}
            $index += $len;
#			print "xxxxxxx\n";
     }
	 return @array_pos;
}

sub treat_cigar_S{
    my ($cigar) = @_;
#   print $cigar,"\n";
    my $index = 0; my $read_pos = 1; my @array_pos;
    while($cigar =~ /((\d+)[MSI])/g){
           my $len = $2;
           if($1 =~ /S/){
                next if $len < 5;
                my $m_spos = $index + $read_pos;
                my $m_epos = $index + $read_pos + $len - 1 ;
                push @array_pos, "$m_spos\t$m_epos";
#               print "$m_spos\t$m_epos\n";
            }
            $index += $len;
#           print "\n";
     }
#	 print "\n";
     return @array_pos;
}

