package YIW::stat;

use strict;
use YIW::basic;
#use warnings;

BEGIN {
	require Exporter;
# set the version for version checking
	our $VERSION = 1.00;
# Inherit from Exporter to export functions and variables
	our @ISA = qw(Exporter);
# Functions and variables which are exported by default
	our @EXPORT = qw(permute_array array2quant corr_pearson array2rank array2var array2entropy find_index density_norm_s density_norm kernel_smooth_xy);
# Functions and variables which can be optionally exported
	our @EXPORT_OK = qw();
}

############################################################
#	permute_array($rarr)
#	array2quant($rarr,$ndat,$quan)
#	corr_pearson($rarx,$rary,$narr)
#	array2rank($rarr,$ndat,$rran)
#	array2var($rarr,$ndat)
#	array2entropy($rarr,$ndat)
#	find_index($rpro,$xxxx)
#	density_norm_s(x)
#	kernel_smooth_xy($rarx,$rary,$roux,$rouy,$nout,$band)
############################################################

my $Pi =  3.141592653589793;
my $GAUSS = 1.0/sqrt(2*$Pi);

return 1;

############################################################
#	permute_array($rarr)
############################################################
# "proper" permutation accorting to Knuth
sub permute_array
{
 my $rarr = shift;
 my $ndat = @$rarr;
 
 for(my $i=0;$i<$ndat-1;$i++){
  my $r = $i + int(rand($ndat-$i));
  my $tmp = $$rarr[$i];
  $$rarr[$i] = $$rarr[$r];
  $$rarr[$r] = $tmp;
 }
}

############################################################
#	array2quant($rarr,$ndat,$quan)
############################################################
# quantile in a sorted array
# intelligent re integer boundaries
sub array2quant
{
 my $rarr = shift;
 my $ndat = shift;
 my $quan = shift;

 return 0 if($ndat==0);
 return $$rarr[0] if($ndat==1);
 my $indr = $ndat*$quan;
 $indr = 0.5 if($indr<=0);
 $indr = $ndat - 0.5 if($indr>=$ndat);
 my $indi = int($indr);
 return $$rarr[$indi] if($indi<$indr);
 return ($$rarr[$indi-1]+$$rarr[$indi])/2;
}

############################################################
#	array2rank($rarr,$ndat,$rran)
############################################################
sub array2rank
{
 my $rarr = shift;
 my $ndat = shift;
 my $rran = shift;
 
 my %d2c = (); my %d2r = ();
 for(my $i=0;$i<$ndat;$i++){
  $d2c{$$rarr[$i]}++;
 }
 
 my $rr = 1;
 foreach my $xx (sort {$b<=>$a} keys %d2c){
  my $nn = $d2c{$xx};
  my $rc = $rr+($nn-1)/2;
  $d2r{$xx} = $rc;
  $rr += $nn;
 }

 for(my $i=0;$i<$ndat;$i++){
  $$rran[$i] = $d2r{$$rarr[$i]};
 }
}

############################################################
#	corr_pearson($rarx,$rary,$narr)
############################################################
# numerically stable computation of Pearson correlation coefficient
# Ronald A. Thisted.
# Elements of Statistical Computing: Numerical Computation,
# CRC Press, 1988
# pp. 84-91
sub corr_pearson
{
 my $rarx = shift;
 my $rary = shift;
 my $narr = shift;

 my $s2x = 0; my $s2y = 0;
 my $sco = 0;
 my $avx = $$rarx[0]; my $avy = $$rary[0];

 for(my $i=1; $i<$narr;$i++){
  my $swe = $i/($i+1);
  my $dex = $$rarx[$i] - $avx;
  my $dey = $$rary[$i] - $avy;
  $s2x += $dex*$dex*$swe;
  $s2y += $dey*$dey*$swe;
  $sco += $dex*$dey*$swe;
  $avx += $dex/($i+1);
  $avy += $dey/($i+1);
 }
 my $corr = 0;
 $corr = $sco/sqrt($s2x*$s2y) if($s2x*$s2y>0);
 return $corr;
}

############################################################
#	array2var($rarr,$ndat)
############################################################
# returns ($av,$sd)
# numerically stable online algorithm
# Donald E. Knuth (1998).
# The Art of Computer Programming, volume 2: Seminumerical Algorithms
# 3rd edn., p. 232. Boston: Addison-Wesley.
sub array2var
{
 my $rarr = shift;
 my $ndat = shift;

 my $av = 0; my $sd = 0;
 
 for(my $i=0;$i<$ndat;$i++){
  my $dd = $$rarr[$i] - $av;
  $av += $dd/($i+1);
  $sd += $dd*($$rarr[$i] - $av);
 }
 $sd /= $ndat - 1 if($ndat>1);
 $sd = 0 if($sd<0);
 return ($av,sqrt($sd));
}

############################################################
#	array2entropy($rarr,$ndat)
############################################################
# returns Shannon entropy in bits
# assumes arbitrary non-negative weights
# safe re 0s
sub array2entropy
{
 my $rarr = shift;
 my $ndat = shift;

 my $sumw = 0;
 for(my $i=0;$i<$ndat;$i++){ $sumw += $$rarr[$i] if($$rarr[$i]>0);}
 return 0 if($sumw<=0);
 my $entr = 0;
 for(my $i=0;$i<$ndat;$i++){
  my $pp = $$rarr[$i]/$sumw;
  $entr -= $pp*log($pp) if($pp>0);
 }
 return $entr/log(2);
}

############################################################
#	find_index($rpro,$xxxx)
############################################################
# @$rpro is a sorted numerical array; last element >=max({$xxxx})
# returns index, such that $$rpro[$i]<=$xxxx
sub find_index
{
 my $rpro = shift;
 my $xxxx = shift;

 	#printf "Start find_index()\n";
 	#for(my $i=0;$i<@$rpro;$i++){ printf "\t%3d %.5f",$i,$$rpro[$i];}
 	#print "\n";
 my $i0 = 0;
 my $i1 = (scalar @$rpro) - 1;
 my $ii = int(($i0+$i1)/2);
 	
 	#printf "%d\t%d\t%d\t%.5f\n",$i0,$ii,$i1,$xxxx;
 while($i0<$i1){
  $ii = int(($i0+$i1)/2);
  if($xxxx<$$rpro[$ii]){ $i1 = $ii; $i1 = $i0 if($i1<$i0);}
  else{ $i0 = $ii+1; $i0 = $i1 if($i0>$i1);}
 	#printf "%d\t%.5f\t%d\t%d\n",$ii,$$rpro[$ii],$i0,$i1;
 }
 #$ii = int(($i0+$i1)/2);
 	#printf "%d\t%d\t%d\n",$i0,$ii,$i1;
 return $i0;
}

############################################################
#	density_norm_s(x)
############################################################
# density for Normal(x,0,1)
sub density_norm_s
{
 my $x = shift;
 
 return $GAUSS*exp(-($x*$x)/2);
}

############################################################
#	density_norm(x)
############################################################
# density for Normal(x,m,s)
sub density_norm
{
 my $x = shift;
 my $m = shift;
 my $s = shift;
 
 return density_norm_s(($x-$m)/$s)/$s;
}

############################################################
#	kernel_smooth_xy($rarx,$rary,$roux,$rouy,$nout,$band)
############################################################
# kernel smooths $rary over $rarx
# writes smoothed data in ($roux,$rouy)
# makes $nout evenly spaced points or $ndat points for each point if $nout<2
# returns number of points
# $band = ($xmax-$xmin)/20 if $band==0
sub kernel_smooth_xy
{
 my $rarx = shift;
 my $rary = shift;
 my $roux = shift;
 my $rouy = shift;
 my $nout = shift; $nout = int($nout);
 my $band = shift;

 my $xmin = minr($rarx);
 my $xmax = maxr($rarx);
 $band = ($xmax-$xmin)/20 if($band<=0);			# if not specified
 $band = 1 if($band<=0);				# if $xmax==$xmin
 
 if($nout>=2){						# $nout evenly spaced points
  my $step = ($xmax-$xmin)/($nout-1);
  for(my $i=0,my $x=$xmin;$i<$nout;$i++,$x+=$step){
   $x = $xmax if($x>$xmax);
   my $y = find_kernel_smooth_xy($rarx,$rary,$x,$band);
   push @$roux,$x;
   push @$rouy,$y;
  }
 }else{							# point by point
  foreach my $x (@$rarx){
   my $y = find_kernel_smooth_xy($rarx,$rary,$x,$band);
   next if($y<=-1e300);
   push @$roux,$x;
   push @$rouy,$y;
  }
 }
 return scalar @$roux;
}

############################################################
#	find_kernel_smooth_xy($rarx,$rary,$x,$band)
############################################################
# scans the whole array for ponts within +/-3*$band of $x
# returns the averaged value or -2e300 if none
sub find_kernel_smooth_xy
{
 my $rarx = shift;
 my $rary = shift;
 my $xval = shift;
 my $band = shift;

 my $sumy = 0; my $sumw = 0;
 for(my $i=0;$i<@$rarx;$i++){
  my $dd = ($$rarx[$i]-$xval)/$band;
  next if($dd<-3 or $dd>3);
  my $wt = density_norm_s($dd);
  $sumw += $wt;
  $sumy += $wt*$$rary[$i];
 }
 return -2e300 if($sumw<=0);
 return $sumy/$sumw;
}
