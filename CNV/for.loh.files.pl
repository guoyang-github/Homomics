#!/usr/bin/perl
#for exomeCNV

$in = shift;
$normal = shift;
$tumor = shift;
$n_col=shift;
$t_col=shift;

open (IN, "$in");
open (NORMAL, ">$normal");
open (TUMOR, ">$tumor");

print NORMAL "chr\tposition\tcoverage\tbaf\n";
print TUMOR "chr\tposition\tcoverage\tbaf\n";

while ($line = <IN>)
{
	chomp($line);
	@tmp = split(/\t/,$line);
	$chr = $tmp[0];
	$pos = $tmp[1];
	$snp = $tmp[2];
	$ref = $tmp[3];
	$var = $tmp[4];
	$filter = $tmp[6];
	$n_info = $tmp[$n_col+8];
	$t_info = $tmp[$t_col+8];
	
	@tmp2 = split(/:/,$n_info);
	$n_geno=$tmp2[0];
	@n_cov=split(/,/,$tmp2[1]);
	$n_ref_cov=$n_cov[0];
	$n_var_cov=$n_cov[1];

	@tmp3 = split(/:/,$t_info);
        $t_geno=$tmp3[0];
        @t_cov=split(/,/,$tmp3[1]);
        $t_ref_cov=$t_cov[0];
        $t_var_cov=$t_cov[1];

	if($filter eq "PASS")
	{
		if($n_geno eq "0/1")
		{
			$n_cov = $n_ref_cov + $n_var_cov;
			$t_cov = $t_ref_cov + $t_var_cov;
			print NORMAL "$chr\t$pos\t$n_cov\t$n_var_cov\n";
			print TUMOR "$chr\t$pos\t$t_cov\t$t_var_cov\n";
		}
		else{}
	}
	else {}
}
close (IN);
close (NORMAL);
close (TUMOR);
