package YIW::stat;

use strict;
use YIW::basic;
use feature 'state';
#use warnings;
use Math::Random::MT qw/rand srand/;

BEGIN {
	require Exporter;
# set the version for version checking
	our $VERSION = 1.00;
# Inherit from Exporter to export functions and variables
	our @ISA = qw(Exporter);
# Functions and variables which are exported by default
	our @EXPORT = qw(permute_array array2quant corr_pearson array2rank array2var array2entropy array2simpson find_index density_norm_s density_norm cdf_norm_s z_norm_s kernel_smooth_xy window_smooth_xy rand_normal rand_exp rand_poisson);
# Functions and variables which can be optionally exported
	our @EXPORT_OK = qw(array2wqprep array2wquant array2smoothmax distr_bhattacharyya_coeff distr_bhattacharyya_distance distr_hellinger_distance distr_earthmover_distance contingency_f1 contingency_p4 data2auc match_norm_miq match_norm_qq);
}

############################################################
#	permute_array($rarr)
#	array2quant($rarr,$ndat,$quan)
#	corr_pearson($rarx,$rary,$narr)
#	array2rank($rarr,$ndat,$rran)
#	array2var($rarr,$ndat)
#	array2entropy($rarr,$ndat)
#	array2simpson($rdat,$ndat)
#	find_index($rpro,$xxxx)
#	density_norm_s(x)
#	density_norm(x,m,s)
#	cdf_norm_s($zz)
#	z_norm_s($pp)
#	kernel_smooth_xy($rarx,$rary,$roux,$rouy,$nout,$band)
#	window_smooth_xy($rarx,$rary,$roux,$rouy,$step,$wind,$post)
#	rand_normal($mm,$ss)
#	rand_exp($ll)
#	rand_poisson($ll)
#	YIW::stat::array2wqprep($rarr,$rwei,$ndat,$rval,$rwsu)
#	YIW::stat::array2wquant($rval,$rwsu,$quan)
#	YIW::stat::array2smoothmax($rarr,$ndat,$alpha)
#	YIW::stat::distr_bhattacharyya_coeff($rda1,$rda2,$ndat,$step)
#	YIW::stat::distr_bhattacharyya_distance($bhat,$DMAX)
#	YIW::stat::distr_hellinger_distance($bhat)
#	YIW::stat::distr_earthmover_distance($rda1,$rda2,$ndat,$flag)
#	YIW::stat::contingency_f1($xpp,$xpn,$xnp,$xnn,$beta)
#	YIW::stat::contingency_p4($xpp,$xpn,$xnp,$xnn)
#	YIW::stat::data2auc($rar1,$rar2)
#	YIW::stat::match_norm_miq($rarr)
#	YIW::stat::match_norm_qq($rarr,$qlo,$qhi)
############################################################

my $Pi =  3.141592653589793;
my $GAUSS = 1.0/sqrt(2*$Pi);

return 1;

############################################################
#	permute_array($rarr)
############################################################
# "proper" permutation according to Knuth
# https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
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
#	array2wqprep($rarr,$rwei,$ndat,$rval,$rwsu)
############################################################
# prepare to calculate
# quantile of an arbitrary numerical array with non-negative weights
sub array2wqprep
{
 my $rarr = shift;
 my $rwei = shift;
 my $ndat = shift;
 my $rval = shift;
 my $rwsu = shift;

 return 0 if($ndat==0);
 return $$rarr[0] if($ndat==1);

 my %data = ();
 for(my $i=0;$i<$ndat;$i++){ $data{$$rarr[$i]} += $$rwei[$i] if($$rwei[$i]>0);}
 
 @$rval = sort{$a<=>$b} keys %data;

 for(my $i=0;$i<@$rval;$i++){
  $$rwsu[$i] = $data{$$rval[$i]};
  $$rwsu[$i] += $$rwsu[$i-1] if($i>0);
 }
}

############################################################
#	array2wquant($rval,$rwsu,$quan)
############################################################
# quantile of an arbitrary numerical array with non-negative weights
sub array2wquant
{
 my $rval = shift;
 my $rwsu = shift;
 my $quan = shift;

 my $nval = @$rval;
 return 0 if($nval==0);
 return $$rval[0] if($nval==1);

 my $wsum = $$rwsu[$nval-1];
 my $qval = $wsum*$quan;
 $qval = 0 if($qval<0);
 $qval = $wsum if($qval>$wsum);

 my $indx = find_index($rwsu,$qval);
 return $$rval[$indx];
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
#	array2simpson($rdat,$ndat)
############################################################
# returns Simpson index
# https://en.wikipedia.org/wiki/Diversity_index#Inverse_Simpson_index
# assumes arbitrary non-negative weights
# safe re 0s
sub array2simpson
{
 my $rdat = shift;
 my $ndat = shift;

 my $sn = 0;
 my $s2 = 0;
 for(my $i=0;$i<$ndat;$i++){
  $sn += $$rdat[$i];
  $s2 += $$rdat[$i]*$$rdat[$i];
 }
 return $s2/$sn/$sn if($sn>0);
 return 0;
}

############################################################
#	array2smoothmax($rarr,$ndat,$alpha)
############################################################
# returns smooth maximum https://en.wikipedia.org/wiki/Smooth_maximum
# safe when exp($$rarr[$i]*$alpha) is safe
# alpha>0 - smooth maximum
# alpha=0 - arithmetic mean
# alpha<0 - smooth minimum
sub array2smoothmax
{
 my $rarr = shift;
 my $ndat = shift;
 my $alpha = shift;

 return $$rarr[0] if($ndat<=0);
 my $suma = 0; my $sumb = 0;
 for(my $i=0;$i<$ndat;$i++){
  my $ww = exp($$rarr[$i]*$alpha);
  $suma += $$rarr[$i]*$ww;
  $sumb += $ww;
 }
 return 0 if($sumb<=0);
 return $suma/$sumb;
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
#	density_norm(x,m,s)
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
#	cdf_norm_s($zz)
############################################################
# cumulative distribution function for standard normal distribution
# Numerical Recipes in Fortran 77: The Art of Scientific Computing (ISBN 0-521-43064-X), 1992, page 214, Cambridge University Press.
# https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
sub cdf_norm_s
{
 my $zz = shift;

 return 0 if($zz<-37);
 return 0.5 if($zz==0);
 return 1 if($zz>8);

 my $tt = 1/(1+abs($zz)/2/sqrt(2));
 
 my @qq = (-1.26551223,1.00002368,0.37409196,0.09678418,-0.18628806,0.27886807,-1.13520398,1.48851587,-0.82215223,0.17087277);

 my $gg = $qq[0]+$tt*($qq[1]+$tt*($qq[2]+$tt*($qq[3]+$tt*($qq[4]+$tt*($qq[5]+$tt*($qq[6]+$tt*($qq[7]+$tt*($qq[8]+$tt*$qq[9]))))))));
 
 my $tau = $tt*exp(-$zz*$zz/2+$gg);
 
 if($zz<=0){ return $tau/2;}
 else{ return 1-$tau/2;}
}

############################################################
#	z_norm_s($pp)
############################################################
# inverse for standard normal distribution
# Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse"
# https://www.academia.edu/9730974/A_handy_approximation_for_the_error_function_and_its_inverse
# https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
sub z_norm_s
{
 my $pp = shift;

 return -38 if($pp<=0);
 return 0 if($pp==0.5);
 return 9 if($pp>=1);

 my $lx = log(4)+log($pp)+log(1-$pp);
 
 my $aa = 0.147;
 my $vx = 2/$Pi/$aa + $lx/2;
 my $s1 = $vx*$vx - $lx/$aa; $s1 = 0 if($s1<0);
 my $s2 = sqrt($s1) - $vx; $s2 = 0 if($s2<0);
 my $ee = sqrt($s2); $ee = -$ee if($pp<0.5);
 return $ee*sqrt(2);
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

############################################################
#	window_smooth_xy($rarx,$rary,$roux,$rouy,$step,$wind,$post)
############################################################
# smooths $rary over $rarx using a sliding window of length $wind
# centered window by default; preceding window if $post>0
# writes smoothed data in ($roux,$rouy)
# makes points starting at >=$xmin, ending at $xmax with $step ($step<=1 -> $step=1)
# returns number of points
# $wind = ($xmax-$xmin)/20 if $wind==0
sub window_smooth_xy
{
 my $rarx = shift;
 my $rary = shift;
 my $roux = shift;
 my $rouy = shift;
 my $step = shift;
 my $wind = shift;
 my $post = shift;

 my $xmin = minr($rarx);
 my $xmax = maxr($rarx);
 $wind = ($xmax-$xmin)/20 if($wind<=0);			# if not specified
 #$wind = $xmax - $xmin if($wind>$xmax-$xmin);		# if overspecified
 
 $step = 1 if($step<=0);				# need a positive $step

 my $xlow = $xmax - $step*int(($xmax-$xmin)/$step);	# start position

 	#printf "xmin\t%.1f\n",$xmin;
 	#printf "xmax\t%.1f\n",$xmax;
 	#printf "step\t%.1f\n",$step;
 	#printf "xlow\t%.1f\n",$xlow;
 	#printf "wind\t%.1f\n",$wind;
 	#exit;
 
 for(my $xpos = $xlow;$xpos<=$xmax;$xpos+=$step){
  my $icur = find_index($rarx,$xpos);
  my $xbeg = $xpos - $wind/2; $xbeg = $xpos - $wind if($post>0);
  my $xend = $xpos + $wind/2; $xend = $xpos if($post>0);
 	#printf "icur\t%d\n",$icur;
 	#printf "xpos\t%.1f\n",$xpos;
 	#printf "xbeg\t%.1f\n",$xbeg;
 	#printf "xend\t%.1f\n",$xend;
  my $sumy = 0; my $npts = 0;
  for(my $i=$icur;$i>=0;$i--){
   if($$rarx[$i]>=$xbeg){
 	#printf "\t%d\t%.1f\t%.1f\n",$i,$$rarx[$i],$$rary[$i];
    $sumy += $$rary[$i];
    $npts++;
   }else{
    last;
   }
  }
  for(my $i=$icur+1;$i<@$rarx;$i++){
   if($$rarx[$i]<=$xend){
 	#printf "\t%d\t%.1f\t%.1f\n",$i,$$rarx[$i],$$rary[$i];
    $sumy += $$rary[$i];
    $npts++;
   }else{
    last;
   }
  }
 	#printf "\t%.1f\t%d\n",$sumy,$npts;
  $sumy /= $npts if($npts>0);
  push @$roux,$xpos;
  push @$rouy,$sumy;
 }

 return scalar @$roux;
}

############################################################
#	rand_normal($mm,$ss)
############################################################
# randomly distributed Normal(mm,ss)
# https://en.wikipedia.org/wiki/Marsaglia_polar_method
sub rand_normal
{
 my $mm = shift;
 my $ss = shift;
 $ss = 1 if($ss<=0);

 state $pair = 0;
 state $rn = 0;
 
 if($pair>0){
  $pair = 0;
  return $rn*$ss + $mm;
 }

 my $rr = 0;
 my $uu = 0;
 my $vv = 0;
 while($rr<=0 or $rr>=1){
  $uu = 2*rand() - 1;
  $vv = 2*rand() - 1;
  $rr = $uu*$uu + $vv*$vv;
 }
 my $xx = sqrt(-2*log($rr)/$rr);
 $rn = $vv*$xx;
 $pair = 1;
 return $uu*$xx*$ss + $mm;
}

############################################################
#	rand_exp($ll)
############################################################
# randomly distributed Exp(lambda)
# expectation = 1/lambda (!)
# https://en.wikipedia.org/wiki/Exponential_distribution
sub rand_exp
{
 my $ll = shift;
 $ll = 1 if($ll<=0);

 my $pp = rand();
 return 0 if($pp==0);
 return -log(1-$pp)/$ll;
}

############################################################
#	rand_poisson($ll)
############################################################
# randomly distributed Poisson(lambda)
# expectation = lambda
# https://en.wikipedia.org/wiki/Poisson_distribution
# inverse transform sampling for lambda < 256
# https://doi.org/10.1007%2F978-1-4613-8643-8_10
# normal approximation for lambda >= 256
sub rand_poisson
{
 my $ll = shift;
 $ll = 1 if($ll<=0);

 if($ll>=256){				# use normal approximation for lambda >= 256
  my $rr = -1; 
  while($rr<0){ $rr = rand_normal($ll,sqrt($ll));}
  return int($rr+0.5);
 }

 my $kk = 0;
 my $pp = exp(-$ll); $pp>0 or die "Underflow with Poisson parameter lambda [$ll] in rand_poisson()";
 my $ss = $pp;
 my $rr = rand();

 while($ss<$rr){
  $kk++;
  $pp = $pp*$ll/$kk;
  $ss += $pp;
 }

 return $kk;
}

############################################################
#	YIW::stat::distr_bhattacharyya_coeff($rda1,$rda2,$ndat,$step)
############################################################
# takes two 1:1 arrays of probability densities
# with points, separated by the given step
# set step=1 for arrays of probability masses in 1:1 bins
# https://en.wikipedia.org/wiki/Bhattacharyya_distance
sub distr_bhattacharyya_coeff
{
 my $rda1 = shift;
 my $rda2 = shift;
 my $ndat = shift;
 my $step = shift;

 my $bhat = 0;
 my $sum1 = 0;
 my $sum2 = 0;
 for(my $i=0;$i<$ndat;$i++){
  my $prod = $$rda1[$i]*$$rda2[$i];
  next if($prod<0);
  $sum1 += $$rda1[$i];
  $sum2 += $$rda2[$i];
  $bhat += sqrt($prod);
 }
 $bhat /= $sum1 if($sum1>0);
 $bhat /= $sum2 if($sum2>0);
 $bhat /= $step;
 $bhat = 0 if($bhat<0);
 $bhat = 1 if($bhat>1);
 return $bhat;
}

############################################################
#	YIW::stat::distr_bhattacharyya_distance($bhat,$DMAX)
############################################################
# takes Bhattacharyya coefficient
# range [0..DMAX]
# https://en.wikipedia.org/wiki/Bhattacharyya_distance
sub distr_bhattacharyya_distance
{
 my $bhat = shift;
 my $DMAX = shift;

 return 0 if($bhat>=1);
 $DMAX = 12 if($DMAX<=0);
 return $DMAX if($bhat<exp(-$DMAX));
 return -log($bhat);
}

############################################################
#	YIW::stat::distr_hellinger_distance($bhat)
############################################################
# takes Bhattacharyya coefficient
# range [0..1]
# https://en.wikipedia.org/wiki/Hellinger_distance
sub distr_hellinger_distance
{
 my $bhat = shift;

 return 0 if($bhat>=1);
 return 1 if($bhat<=0);
 return sqrt(1-$bhat);
}

############################################################
#	YIW::stat::distr_earthmover_distance($rda1,$rda2,$ndat,$flag)
############################################################
# takes two arrays of probability masses in 1:1 bins
# flag>0 to force normalization
# range [0..inf] ; -1 if wrong data
# https://en.wikipedia.org/wiki/Earth_mover%27s_distance
sub distr_earthmover_distance
{
 my $rda1 = shift;
 my $rda2 = shift;
 my $ndat = shift;
 my $flag = shift;

 if($flag>0){
  my $sum1 = 0;
  my $sum2 = 0;
  for(my $i=0;$i<$ndat;$i++){
   return -1 if($$rda1[$i-1]<0 or $$rda2[$i-1]<0);
   $sum1 += $$rda1[$i] if($$rda1[$i]>0);
   $sum2 += $$rda2[$i] if($$rda2[$i]>0);
  }
  return -1 if($sum1<=0 or $sum2<=0);
  for(my $i=0;$i<$ndat;$i++){
   $$rda1[$i] /= $sum1;
   $$rda2[$i] /= $sum2;
  }
 }
 my $ecur = 0;
 my $esum = 0;
 for(my $i=1;$i<$ndat;$i++){
  return -1 if($$rda1[$i-1]<0 or $$rda2[$i-1]<0);
  $ecur += $$rda1[$i-1] - $$rda2[$i-1];
  $esum += abs($ecur);
 }
 return $esum;
}

############################################################
#	YIW::stat::contingency_f1($xpp,$xpn,$xnp,$xnn,$beta)
############################################################
# $xpp = TP
# $xpn = FN
# $xnp = FP
# $xnn = TN
# assumes non-negative counts or weights, non-degenerate
# $beta > 0; default 1;
# https://en.wikipedia.org/wiki/F-score
sub contingency_f1
{
 my $xpp = shift;
 my $xpn = shift;
 my $xnp = shift;
 my $xnn = shift;
 my $beta = shift;

 $beta = 1 if($beta<=0);
 
 return 0 if($xpp<0 or $xpn<0 or $xnp<0 or $xnn<0 or $xpp+$xpn+$xnp<=0);
 return (1+$beta*$beta)*$xpp/((1+$beta*$beta)*$xpp+$beta*$beta*$xpn+$xnp);
}

############################################################
#	YIW::stat::contingency_p4($xpp,$xpn,$xnp,$xnn)
############################################################
# $xpp = TP
# $xpn = FN
# $xnp = FP
# $xnn = TN
# assumes non-negative counts or weights, non-degenerate
# https://en.wikipedia.org/wiki/P4-metric
sub contingency_p4
{
 my $xpp = shift;
 my $xpn = shift;
 my $xnp = shift;
 my $xnn = shift;

 my $de = 4*$xpp*$xnn + ($xpp+$xnn)*($xpn+$xnp);
 
 return 0 if($xpp<0 or $xpn<0 or $xnp<0 or $xnn<0 or $de<=0);
 return 4*$xpp*$xnn/$de;
}

############################################################
#	YIW::stat::data2auc($rar1,$rar2,$stat)
############################################################
# 2 arrays of numerical values; "positive" and "negative" scores
# computes Mann–Whitney z-score if stat>0
# returns (AUC,z-score)
# https://en.wikipedia.org/wiki/Mann–Whitney_U_test
sub data2auc
{
 my $rar1 = shift;
 my $rar2 = shift;
 my $stat = shift;

 my $n1 = @$rar1;
 my $n2 = @$rar2;

 my @data = (@$rar1,@$rar2);
 my @datr = ();
 array2rank(\@data,$n1+$n2,\@datr);

 my $r1 = 0;
 for(my $i=0;$i<$n1;$i++){
  $r1 += $datr[$i+0];
 }
 my $r2 = 0;
 for(my $i=0;$i<$n2;$i++){
  $r2 += $datr[$i+$n1];
 }
 my $u1 = $n1*$n2 - $r1 + $n1*($n1+1)/2;
 my $auc1 = 1; $auc1 = 0 if($n1<=0); $auc1 = $u1/$n1/$n2 if($n1>0 and $n2>0);
 my $zz = 0;

 return ($auc1,$zz) if($stat<=0 or $n1+$n2<=1);

 my %dcnt = ();
 for(my $i=0;$i<@data;$i++){ $dcnt{$data[$i]}++;}
 my $tt = 0;
 foreach my $xx (keys %dcnt){
  my $nn = $dcnt{$xx};
  $tt += $nn*($nn*$nn-1) if($nn>1);
 }
 my $sd = $n1*$n2/12*($n1+$n2+1-$tt/($n1+$n2)/($n1+$n2-1));
 if($sd>0){ $sd = sqrt($sd);}
 
 $zz = ($u1 - $n1*$n2/2)/$sd if($sd>0);
 return ($auc1,$zz);
}

############################################################
#	YIW::stat::data2pBruMun($rar1,$rar2,$nper)
############################################################
# UNREASONABLY SLOW
# 2 arrays of numerical values; "positive" and "negative" scores
# uses nper permtations to compute Brunner Munzel Test p-value
# returns (probability,p-value)
# https://en.wikipedia.org/wiki/Brunner_Munzel_Test
sub data2pBruMun
{
 my $rar1 = shift;
 my $rar2 = shift;
 my $nper = shift;

 my $n1 = @$rar1;
 my $n2 = @$rar2;

 return (0.5,0.5) if($n1<1 or $n2<1 or $n1+$n2<=2 or $nper<2);
 
 my $nn = $n1 + $n2;
 my @data = (@$rar1,@$rar2);

 my $pp = prob_array_diff(\@data,$n1,$n2);
 
 my $cnt = 0;
 for(my $i=0;$i<$nper;$i++){
  permute_array(\@data);
  my $px = prob_array_diff(\@data,$n1,$n2);
  $cnt++ if($px>$pp);
  $cnt+= 0.5 if($px==$pp);
 }

 return ($pp,$cnt/$nper);
}

############################################################
#	YIW::stat::prob_array_diff($rarr,$n1,$n2)
############################################################
sub prob_array_diff
{
 my $rarr = shift;
 my $n1 = shift;
 my $n2 = shift;

 return 0.5 if($n1<1 or $n2<1 or $n1+$n2<=2);

 my $s1 = min($n1,1000);
 my $s2 = min($n2,1000);

 my $pp = 0;

 for(my $i=0;$i<$s1;$i++){
  my $idx1 = ($s1==$n1)?$i:int(rand($n1));
  my $v1 = $$rarr[$idx1];
  for(my $j=0;$j<$s2;$j++){
   my $idx2 = ($s2==$n2)?$j:int(rand($n2));
   my $v2 = $$rarr[$idx2+$n1];
   $pp++ if($v1>$v2);
   $pp+= 0.5 if($v1==$v2);
  }
 }
 
 return $pp/($s1*$s2);
}

############################################################
#	YIW::stat::match_norm_miq($rarr)
############################################################
# for a numerical array matches normal distribution using
# median and interquartile distance
# returns ($av,sd)
sub match_norm_miq
{
 my $rarr = shift;

 my $ndat = @$rarr;
 my @ltmp = sort {$a<=>$b} @$rarr;

 my $q25 = array2quant(\@ltmp,$ndat,0.25);
 my $q50 = array2quant(\@ltmp,$ndat,0.50);
 my $q75 = array2quant(\@ltmp,$ndat,0.75);
 my $siq = ($q75-$q25)/(z_norm_s(0.75)-z_norm_s(0.25));
 return ($q50,$siq);
}

############################################################
#	YIW::stat::match_norm_qq($rarr,$plo,$phi)
############################################################
# for a numerical array matches normal distribution using
# 0<plo<0.5 and 0.5<phi<1
# returns ($av,sd)
sub match_norm_qq
{
 my $rarr = shift;
 my $plo = shift;
 my $phi = shift;


 return (0,0) if($plo>=0.5 or $plo<=0 or $phi<=0.5 or $phi>=1);
 my $ndat = @$rarr;
 my @ltmp = sort {$a<=>$b} @$rarr;

 my $qlo = array2quant(\@ltmp,$ndat,$plo);
 my $qhi = array2quant(\@ltmp,$ndat,$phi);
 my $zlo = z_norm_s($plo);
 my $zhi = z_norm_s($phi);
 my $av = ($qlo*$zhi-$qhi*$zlo)/($zhi-$zlo);
 my $sd = ($qlo-$av)/($zlo);
 return ($av,$sd);
}
