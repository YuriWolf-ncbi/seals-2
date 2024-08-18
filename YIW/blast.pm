package YIW::blast;

############################################################
#	Date:	November 21 2019
############################################################
use strict;
#use warnings;
use lib "/home/wolf/perl5/lib/perl5";
use YIW::basic;

BEGIN {
	require Exporter;
# set the version for version checking
	our $VERSION = 1.00;
# Inherit from Exporter to export functions and variables
	our @ISA = qw(Exporter);
# Functions and variables which are exported by default
	our @EXPORT = qw(blast_read_btab blast_process_hits blast_expand_footptint blast_smooth_path blast_reconcile_paths blast_process_nthits);
# Functions and variables which can be optionally exported
	our @EXPORT_OK = qw(calc_overlap calc_path_overlap path_length clear_subject_hits aggregate_ranges);
}

my $delim = "[/|]";
my $LOEVAL = 1e-200;
my $HIEVAL = 1e10;
my $VERYLONG = 1e10;

return 1;

############################################################
#	blast_read_btab($fnam,$ethr,$qwrd,$hwrd,$rdat,$rlen,$q300,$fpl1)
#	blast_process_hits($qq,$hh,$ql,$hl,$rhit,$othr,$wbas,$dquo,$cthr,$gmax,$fraw)
#	blast_expand_footptint($patq,$path,$ql,$hl)
#	blast_smooth_path($path,$sl,$wbas,$dquo)
#	blast_reconcile_paths($rpat,$rnew,$lpro,$wbas,$dquo,$hardexpand,$doexpand)
#	YIW::blast::calc_overlap($xb,$xe,$yb,$ye)
#	YIW::blast::calc_path_overlap($patx,$paty)
#	YIW::blast::path_length($path)
#	YIW::blast::clear_subject_hits($rhit,$path,$othr)
#	blast_process_nthits($ql,$hl,$rhit,$orth,$oath)
#	YIW::blast::aggregate_ranges($rori,$ragg)
############################################################

############################################################
#	blast_read_btab($fnam,$ethr,$qwrd,$hwrd,$rdat,$rlen,$q300,$fpl1)
############################################################
# assumes -outfmt "7 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore"
# puts arrays of hits for query-subject tab-delimited pairs into %$rdat as
#	$ss,$qb,$qe,$hb,$he,$ev
# $fpl1 increments all coordinates (MMSEQS2)
# puts query and subject lengths into %$rlen
sub blast_read_btab
{
 my $fnam = shift;
 my $ethr = shift;
 my $qwrd = shift;
 my $hwrd = shift;
 my $rdat = shift;
 my $rlen = shift;
 my $q300 = shift;
 my $fpl1 = shift;

 open HAND,"<$fnam" or die "Can't read \"$fnam\"";
 while(<HAND>){
  chomp;
  next if(m/^\#/);
  my ($qq,$hh,$ql,$hl,$qb,$qe,$hb,$he,$ev,$ss) = split/\t/;
  next if($ss<=0 or $ql<=0);
  my $elim = $ethr; $elim *= (300/$ql)**$q300 if($q300>0);
  next if($ev>$elim);
  if($fpl1){ $qb++; $qe++; $hb++; $he++;}
  $ev = $LOEVAL if($ev<$LOEVAL);
  my $qid = $qq;
  $qid = (split /$delim/,$qq)[$qwrd-1] if($qwrd>0 and $qq=~m/$delim/);
  my $hid = $hh;
  $hid = (split /$delim/,$hh)[$hwrd-1] if($hwrd>0 and $hh=~m/$delim/);
  my $qqhh = join "\t",($qid,$hid);
  if($$rdat{$qqhh} eq ""){ my @tmp = (); $$rdat{$qqhh} = \@tmp;}
  my $rhit = $$rdat{$qqhh};
  push @$rhit,(join "\t",($ss,$qb,$qe,$hb,$he,$ev));
  $$rlen{$qid} = $ql + 0;
  $$rlen{$hid} = $hl + 0;
 }
 close HAND;
}

############################################################
#	blast_process_hits($qq,$hh,$ql,$hl,$rhit,$othr,$wbas,$dquo,$cthr,$gmax,$fraw)
############################################################
# reads array of hits for query-subject pair as read by blast_read_btab()
#	$ss,$qb,$qe,$hb,$he,$ev
# sorts @$rhit by score
# returns ($scor,$code,$patq,$path,$qple,$hple,$eval) array
# segmented paths are segmented in both query and subject
sub blast_process_hits
{
 my $qq = shift;
 my $hh = shift;
 my $ql = shift;
 my $hl = shift;
 my $rhit = shift;
 my $othr = shift;
 my $wbas = shift;
 my $dquo = shift;
 my $cthr = shift;
 my $gmax = shift;
 my $fraw = shift;
 
#---	sort by score --------------------------------------
 @$rhit = sort {$b<=>$a} @$rhit;

#---	find greedy non-overlapping set of segments --------

 my $plst = find_best_path($rhit,$othr);
 	#printf "#\t%s\t%s\t%s\n",$qq,$hh,$plst;

 my @hlst = ();
 foreach my $hidx (split/,/,$plst){
  my ($ss,$qb,$qe,$hb,$he,$ev) = split/\t/,$$rhit[$hidx-1];
  my $lhit = min($qe-$qb+1,$he-$hb+1);
  my $sden = $ss/$lhit;
  push @hlst,(join "\t",($qb,$qe,$hb,$he,$ss,$sden,$ev));
  $$rhit[$hidx-1] = "";
 }

#---	sort by position -----------------------------------
 @hlst = sort {$a<=>$b} @hlst;

#---	resolve conflicts ----------------------------------
 resolve_overlaps(\@hlst) unless($fraw>0);
 resolve_long_gaps(\@hlst,$gmax);
 #@@# bridge_short_gaps(\@hlst,$ql,$hl,$wbas,$dquo);
 # not here, need a separate procedure to smooth path pairs

#---	calc code ------------------------------------------
 my ($code,$scor,$eval) = calc_match_code(\@hlst,$ql,$hl,$wbas,$dquo,$cthr);
 $code = -1 if($qq eq $hh);

#---	get paths ------------------------------------------
 my ($patq,$path) = assemble_paths(\@hlst);

 my $qple = path_length($patq);
 my $hple = path_length($path);
 
 return ($scor,$code,$patq,$path,$qple,$hple,$eval);
}

############################################################
#	find_best_path($rhit,$othr)
############################################################
sub find_best_path
{
 my $rhit = shift;
 my $othr = shift;

 my @hlst = ();
 
 for(my $i=0;$i<@$rhit;$i++){
  my ($ss,$qb,$qe,$hb,$he,$ev) = split/\t/,$$rhit[$i];
  next if($ss<=0);
  my $lenq = $qe - $qb + 1;
  my $lenp = $he - $hb + 1;
  my $good = 1;
  for(my $j=0;$j<@hlst;$j++){
   my $jidx = $hlst[$j];
   my ($xs,$xb,$xe,$yb,$ye,$xv) = split/\t/,$$rhit[$jidx-1];
   my $xenq = $xe - $xb + 1;
   my $xenp = $ye - $yb + 1;
   my $overq = calc_overlap($qb,$qe,$xb,$xe);
   my $overp = calc_overlap($hb,$he,$yb,$ye);
   if($overq/min($lenq,$xenq)>$othr or $overp/min($lenp,$xenp)>$othr or calc_compatibility($qb,$hb,$xb,$yb)==0){
    $good = 0;
    last;
   }
  }
  push @hlst,$i+1 if($good);
 }
 return (join ",",@hlst);
}

############################################################
#	YIW::blast::calc_overlap($xb,$xe,$yb,$ye)
############################################################
sub calc_overlap
{
 my $xb = shift;
 my $xe = shift;
 my $yb = shift;
 my $ye = shift;
 
 my $over = min($xe,$ye) - max($xb,$yb) + 1;
 return 0 if($over<0);
 return $over;
}

############################################################
#	calc_compatibility($x1,$y1,$x2,$y2)
############################################################
sub calc_compatibility
{
 my $x1 = shift;
 my $y1 = shift;
 my $x2 = shift;
 my $y2 = shift;
 
 my $dx = $x2 - $x1;
 my $dy = $y2 - $y1;

 return 1 if($dx*$dy>0);
 return 0;
}

############################################################
#	YIW::blast::calc_path_overlap($patx,$paty)
############################################################
sub calc_path_overlap
{
 my $patx = shift;
 my $paty = shift;
 
 my $over = 0;

 foreach my $fx (split/=/,$patx){
  my ($xb,$xe) = split/-/,$fx;
  foreach my $fy (split/=/,$paty){
   my ($yb,$ye) = split/-/,$fy;
   $over += calc_overlap($xb,$xe,$yb,$ye);
  }
 }
 return 0 if($over<0);
 return $over;
}

############################################################
#	resolve_overlaps($rhit)
############################################################
sub resolve_overlaps
{
 my $rhit = shift;

 for(my $i=1;$i<@$rhit;$i++){
  my ($q1b,$q1e,$h1b,$h1e,$ss1,$sd1,$ev1) = split/\t/,$$rhit[$i-1];
  my ($q2b,$q2e,$h2b,$h2e,$ss2,$sd2,$ev2) = split/\t/,$$rhit[$i];
  my $over = max($q1e-$q2b+1,$h1e-$h2b+1);
  next if($over<=0);
  my $lor1 = min($q1e-$q1b+1,$h1e-$h1b+1);
  my $lor2 = min($q2e-$q2b+1,$h2e-$h2b+1);
  my $w1 = $sd1/($sd1+$sd2);
  my $c1 = int($over*$w1+0.5);

  $q1e = $q1e - $over + $c1;
  $h1e = $h1e - $over + $c1;
  $q2b = $q2b + $c1;
  $h2b = $h2b + $c1;
  my $lne1 = min($q1e-$q1b+1,$h1e-$h1b+1);
  my $lne2 = min($q2e-$q2b+1,$h2e-$h2b+1);
  $ss1 = sprintf "%.1f",$ss1*$lne1/$lor1;
  $ss2 = sprintf "%.1f",$ss2*$lne2/$lor2;
  $$rhit[$i-1] = join "\t",($q1b,$q1e,$h1b,$h1e,$ss1,$sd1,$ev1);
  $$rhit[$i] = join "\t",($q2b,$q2e,$h2b,$h2e,$ss2,$sd2,$ev2);
 }
}

############################################################
#	resolve_long_gaps($rhit,$gmax)
############################################################
sub resolve_long_gaps
{
 my $rhit = shift;
 my $gmax = shift;

 $gmax = $VERYLONG if($gmax<=0);

 my $long = 0;
 for(my $i=1;$i<@$rhit;$i++){
  my ($q1b,$q1e,$h1b,$h1e,$ss1,$sd1,$ev1) = split/\t/,$$rhit[$i-1];
  my ($q2b,$q2e,$h2b,$h2e,$ss2,$sd2,$ev2) = split/\t/,$$rhit[$i];
  my $lgap = max($q2b-$q1e,$h2b-$h1e);
  if($lgap>$gmax){
   $long = 1;
   last;
  }
 }
 return if($long==0);
 
 my $bests = -1; my $bestf = "";
 for(my $i=0;$i<@$rhit;$i++){
  my ($qxb,$qxe,$hxb,$hxe,$ssx,$sdx,$evx) = split/\t/,$$rhit[$i];
  if($ssx>$bests){
   $bests = $ssx;
   $bestf = $$rhit[$i];
  }
 }
 return if($bests<0);
 @$rhit = ();
 push @$rhit,$bestf;
}

############################################################
#	bridge_short_gaps($rhit,$ql,$hl,$wbas,$dquo)
############################################################
sub bridge_short_gaps
{
 my $rhit = shift;
 my $ql = shift;
 my $hl = shift;
 my $wbas = shift;
 my $dquo = shift;
 
 my $xwin = min($wbas,min($ql,$hl)*$dquo);
 	#printf "#\t%s\t%s\t%s\n",$ql,$hl,$xwin;

 for(my $i=0;$i<@$rhit;$i++){
  my ($qb,$qe,$hb,$he,$ss,$sd,$ev) = split/\t/,$$rhit[$i];
  if($i==0){						# gap at start
   my $gq = $qb - 1;
   my $gh = $hb - 1;
   my $gmin = min($gq,$gh);
   if($gmin>0 and $gmin<=$xwin){				# short gap
    $qb -= $gmin;
    $hb -= $gmin;
    $$rhit[$i] = join "\t",($qb,$qe,$hb,$he,$ss,$sd,$ev);
   }
  }
  if($i>=1){						# gap between i and i-1
   my ($q1b,$q1e,$h1b,$h1e,$ss1,$sd1,$ev1) = split/\t/,$$rhit[$i-1];
   my $gq = $qb - $q1e - 1;
   my $gh = $hb - $h1e - 1;
   my $gmin = min($gq,$gh);
   if($gmin>0 and $gmin<=$xwin){				# short gap
    my $w1 = $sd1/($sd1+$sd);
    my $c1 = int($gmin*$w1+0.5);
    my $c2 = $gmin - $c1;
    $q1e += $c1;
    $h1e += $c1;
    $qb -= $c2;
    $hb -= $c2;
    $$rhit[$i-1] = join "\t",($q1b,$q1e,$h1b,$h1e,$ss1,$sd1,$ev1);
    $$rhit[$i] = join "\t",($qb,$qe,$hb,$he,$ss,$sd,$ev);
   }
  }
  if($i==(@$rhit-1)){					# gap at end
   my $gq = $ql - $qe;
   my $gh = $hl - $he;
   my $gmin = min($gq,$gh);
   if($gmin>0 and $gmin<=$xwin){				# short gap
    $qe += $gmin;
    $he += $gmin;
    $$rhit[$i] = join "\t",($qb,$qe,$hb,$he,$ss,$sd,$ev);
   }
  }
 }
}

############################################################
#	calc_match_code($rhit,$ql,$hl,$wbas,$dquo,$cthr)
############################################################
sub calc_match_code
{
 my $rhit = shift;
 my $ql = shift;
 my $hl = shift;
 my $wbas = shift;
 my $dquo = shift;
 my $cthr = shift;

 my $xwin = min($wbas,min($ql,$hl)*$dquo);
 $xwin = max($ql,$hl) if($xwin==0);			# code <- 0 if not specified

 my $covq = 0; my $gapq = 0;
 my $covh = 0; my $gaph = 0;
 my $scor = 0;
 my $eval = $HIEVAL;
 for(my $i=0;$i<@$rhit;$i++){
  my ($qb,$qe,$hb,$he,$ss,$sd,$ev) = split/\t/,$$rhit[$i];
  $covq += $qe - $qb + 1;
  $covh += $he - $hb + 1;
  $scor += $ss;
  $eval = $ev if($ev<$eval);
  if($i==0){						# gap at start
   my $gq = $qb - 1; $gapq = $gq if($gq>$gapq);
   my $gh = $hb - 1; $gaph = $gh if($gh>$gaph);
  }
  if($i>=1){						# gap between i and i-1
   my ($q1b,$q1e,$h1b,$h1e,$ss1,$sd1) = split/\t/,$$rhit[$i-1];
   my $gq = $qb - $q1e - 1; $gapq = $gq if($gq>$gapq);
   my $gh = $hb - $h1e - 1; $gaph = $gh if($gh>$gaph);
  }
  if($i==(@$rhit-1)){					# gap at end
   my $gq = $ql - $qe; $gapq = $gq if($gq>$gapq);
   my $gh = $hl - $he; $gaph = $gh if($gh>$gaph);
  }
 }
 $covq /= $ql;
 $covh /= $hl;
 my $code = 0;
 $code += 2 if($covq<$cthr or $gapq>$xwin);
 $code += 1 if($covh<$cthr or $gaph>$xwin);

 return ($code,$scor,$eval);
}

############################################################
#	assemble_paths($rhit)
############################################################
sub assemble_paths
{
 my $rhit = shift;

 my $patq = "";
 my $patp = "";
 for(my $i=0;$i<@$rhit;$i++){
  my ($qb,$qe,$hb,$he,$ss,$sd) = split/\t/,$$rhit[$i];
  $patq .= "=" if($patq ne ""); $patq .= $qb."-".$qe;
  $patp .= "=" if($patp ne ""); $patp .= $hb."-".$he;
 }
 return ($patq,$patp);
}

############################################################
#	YIW::blast::path_length($path)
############################################################
sub path_length
{
 my $path = shift;

 my $plen = 0;
 foreach my $fr (split/=/,$path){
  my ($p1,$p2) = split/-/,$fr;
  $plen += $p2 - $p1 + 1;
 }
 return $plen;
}

############################################################
#	blast_expand_footptint($patq,$path,$ql,$hl)
############################################################
# expands query footprint on the hit sequence to the equivalent of the full query
# expects co-segmented paths
sub blast_expand_footptint
{
 my $patq = shift;
 my $path = shift;
 my $ql = shift;
 my $hl = shift;
#my $wbas = shift;
#my $dquo = shift;

#my $xwin = min($wbas,min($ql,$hl)*$dquo);

#---	split into segments --------------------------------
 my @frq = ();
 foreach my $xx (split /=/,$patq){
  push @frq,$xx;
 }
 
 my @frh = ();
 foreach my $yy (split /=/,$path){
  push @frh,$yy;
 }
 
#---	amend missegmented paths ---------------------------
 if((scalar @frq)!=(scalar @frh)){			# mismatched segmentation - only ends!
  @frq = ();
  my ($pq1) = ($patq =~ m/^(\d+)/);
  my ($pq2) = ($patq =~ m/(\d+)$/);
  push @frq,($pq1."-".$pq2);
  @frh = ();
  my ($ph1) = ($path =~ m/^(\d+)/);
  my ($ph2) = ($path =~ m/(\d+)$/);
  push @frh,($ph1."-".$ph2);
 }

#---	close gaps relative to query -----------------------
 for(my $i=0;$i<@frq;$i++){
  my ($pq1,$pq2) = split/-/,$frq[$i];
  my ($ph1,$ph2) = split/-/,$frh[$i];
  if($i==0){						# gap at start
   my $ming = min($pq1-1,$ph1-1);
   if($ming>0){
    $ph1 -= $ming;
    $frh[$i] = $ph1."-".$ph2;
   }
  }
  if($i>=1){						# gap between i and i-1
   my ($qq1,$qq2) = split/-/,$frq[$i-1];
   my ($qh1,$qh2) = split/-/,$frh[$i-1];
   my $ming = min($pq1-$qq2-1,$ph1-$qh2-1);
   if($ming>0){
    my $lsqp = $pq2 - $pq1 + 1;
    my $lsqq = $qq2 - $qq1 + 1;
    my $llef = int($ming*$lsqq/($lsqq+$lsqp)+0.5); $llef = $ming if($llef>$ming);
    my $lrig = $ming - $llef;
    $qh2 += $llef;
    $ph1 -= $lrig;
    $frh[$i]   = $ph1."-".$ph2;
    $frh[$i-1] = $qh1."-".$qh2;
   }
  }
  if($i==(@frq-1)){					# gap at end
   my $ming = min($ql-$pq2,$hl-$ph2);
   if($ming>0){
    $ph2 += $ming;
    $frh[$i] = $ph1."-".$ph2;
   }
  }
 }

#---	close short remaining gaps -------------------------
# doesn't belong here
#for(my $i=0;$i<@frh;$i++){
# my ($ph1,$ph2) = split/-/,$frh[$i];
# if($i==0){						# gap at start
#  my $lgap = $ph1 - 1;
#  if($lgap>0 and $lgap<=$xwin){
#   $ph1 = 1;
#   $frh[$i] = $ph1."-".$ph2;
#  }
# }
# if($i>=1){						# gap between i and i-1
#  my ($qh1,$qh2) = split/-/,$frh[$i-1];
#  my $lgap = $ph1 - $qh2 - 1;
#  if($lgap>0 and $lgap<=$xwin){
#   my $lshp = $ph2 - $ph1 + 1;
#   my $lshq = $qh2 - $qh1 + 1;
#   my $llef = int($lgap*$lshq/($lshq+$lshp)+0.5); $llef = $lgap if($llef>$lgap);
#   my $lrig = $lgap - $llef;
#   $qh2 += $llef;
#   $ph1 -= $lrig;
#   $frh[$i]   = $ph1."-".$ph2;
#   $frh[$i=1] = $qh1."-".$qh2;
#  }
# }
# if($i==(@frh-1)){					# gap at end
#  my $lgap = $lh - $ph2 - 1;
#  if($lgap>0 and $lgap<=$xwin){
#   $ph2 += $lgap;
#   $frh[$i] = $ph1."-".$ph2;
#  }
# }
#}

#---	export hit path ------------------------------------
 my $patn = "";
 foreach my $yy (@frh){
  $patn .= "=" if($patn ne ""); $patn .= $yy;
 }
 return $patn;
}

############################################################
#	blast_smooth_path($path,$sl,$wbas,$dquo)
############################################################
# eliminates short segments uncovered by the path
# might destroy original segmentation
sub blast_smooth_path
{
 my $path = shift;
 my $sl = shift;
 my $wbas = shift;
 my $dquo = shift;

my $xwin = min($wbas,$sl*$dquo);

#---	split into segments --------------------------------
 my @fra = ();
 foreach my $xx (split /=/,$path){
  push @fra,$xx;
 }
 
#---	close short gaps -----------------------------------
 for(my $i=0;$i<@fra;$i++){
  my ($pp1,$pp2) = split/-/,$fra[$i];
  if($i==0){						# gap at start
   my $lgap = $pp1-1;
   if($lgap<=$xwin){						# short gap
    $pp1 = 1;
    $fra[$i] = $pp1."-".$pp2;						# start from 1
   }
  }
  if($i>=1){						# gap between i and i-1
   my ($qq1,$qq2) = split/-/,$fra[$i-1];
   my $lgap = $pp1 - $qq2 - 1;
   if($lgap<=$xwin){						# short gap
    $pp1 = $qq1;
    $fra[$i] = $pp1."-".$pp2;						# subsume i-1 into i
    $fra[$i-1] = "";							# destroy i-1
   }
  }
  if($i==(@fra-1)){					# gap at end
   my $lgap = $sl - $pp2;
   if($lgap<=$xwin){						# short gap
    $pp2 = $sl;
    $fra[$i] = $pp1."-".$pp2;						# end at end
   }
  }
 }

#---	export path ----------------------------------------
 my $patn = "";
 foreach my $xx (@fra){
  next if($xx eq "");
  $patn .= "=" if($patn ne ""); $patn .= $xx;
 }
 return $patn;
}

############################################################
#	YIW::blast::clear_subject_hits($rhit,$path,$othr)
############################################################
# reads array of hits for query-subject pair as read by blast_read_btab()
#	$ss,$qb,$qe,$hb,$he
# removes segments, overlapping the given path
sub clear_subject_hits
{
 my $rhit = shift;
 my $path = shift;
 my $othr = shift;
 
 foreach my $hh (split/=/,$path){
  my ($p1,$p2) = split/-/,$hh;
  my $ll = $p2 - $p1 + 1;
  for(my $i=0;$i<@$rhit;$i++){
   my ($ss,$qb,$qe,$hb,$he) = split/\t/,$$rhit[$i];
   next if($ss<=0);
   my $lx = $he - $hb + 1;
   my $over = calc_overlap($p1,$p2,$hb,$he);
   $$rhit[$i] = "" if($over/min($ll,$lx)>$othr);
  }
 }
}

############################################################
#	blast_reconcile_paths($rpat,$rnew,$lpro,$wbas,$dquo,$hardexpand,$noexpand)
############################################################
# reads list of compatible raw paths for a particular subject
#	$scor,$code,$path,$patq,$hple,$qple,$hlen,$qlen,$hh,$qq,$ev
# produces paths expanded by query and/or smoothed
#	new paths are in $$rnew[$i] corresponding to the $$rpat[$i]
# set $wbas<0 for expand only
# set $hardexpand>0 to expand class 2 and 3 paths
# set $noexpand>0 to never expand paths
sub blast_reconcile_paths
{
 my $rpat = shift;
 my $rnew = shift;
 my $lpro = shift;
 my $wbas = shift;
 my $dquo = shift;
 my $hardexpand = shift;
 my $noexpand = shift;


#---	expand paths; collect segments ---------------------
 my @lfra = ();

 for(my $i=0;$i<@$rpat;$i++){
  my ($sx,$cx,$path,$patq,$plex,$plqx,$hl,$ql) = split/\t/,$$rpat[$i];
   	#printf STDERR "===\t%d\n",$i;
   	#printf STDERR "%d\t%s\t%s\n",$i,$path,$patq;
  my $patn = $path;
  $patn = blast_expand_footptint($patq,$path,$ql,$hl) if(($noexpand<=0) and ($cx<=1 or $hardexpand));
   	#printf STDERR "%d\t%s\t%s\n",$i,$path,$patn;
  my @lfrp = split/=/,$path;
  my @lfrn = split/=/,$patn;
  for(my $j=0;$j<@lfrp;$j++){
   push @lfra,(join "\t",($lfrp[$j],$lfrn[$j],($i+1)));
   	#printf STDERR "%d\t%d\t%s\n",$i,$j,(join "\t",($lfrp[$j],$lfrn[$j],($i+1)));
  }
 }
 
#---	sort segments --------------------------------------
 @lfra = sort {$a<=>$b} @lfra;

	#for(my $i=0;$i<@lfra;$i++){
	# printf STDERR "# %d\t%s\n",$i,$lfra[$i];
	#}

#---	calculate score denstities -------------------------
 my @dpat = ();
 for(my $i=0;$i<@$rpat;$i++){
  my ($sx,$cx,$path,$patq,$plen,$pleq) = split/\t/,$$rpat[$i];
  my $minl = min($plen,$pleq);  $dpat[$i] = $sx/$minl if($minl>0);
 }

#---	gaps and conflicts ---------------------------------
 my $xwin = min($wbas,$lpro*$dquo);

 for(my $i=0;$i<@lfra;$i++){
  my ($xfr,$xfn,$xi) = split/\t/,$lfra[$i];
  my ($xb,$xe) = split/-/,$xfn;
  my $xden = $dpat[$xi-1];
  if($i==0){						# at the start
   if($xb<=$xwin){
    $xb = 1;
    $xfn = join "-",($xb,$xe);
    $lfra[$i] = join "\t",($xfr,$xfn,$xi);
   }
  }
  if($i>=1){						# gap between i and i-1
   my ($yfr,$yfn,$yi) = split/\t/,$lfra[$i-1];
   my ($yb,$ye) = split/-/,$yfn;
   my $yden = $dpat[$yi-1];
   my $lgap = $xb - $ye - 1;
   	#printf STDERR "# %d-%d\t%d\n\t%s\n\t%s\n",$i,$i-1,$lgap,$lfra[$i-1],$lfra[$i];
   if($lgap>0 and $lgap<=$xwin){				# need to close
    my $yadd = int($lgap*$yden/($yden+$xden)+0.5);
    my $xsub = $lgap - $yadd;
    $ye += $yadd;
    $xb -= $xsub;
    $xfn = join "-",($xb,$xe);
    $lfra[$i] = join "\t",($xfr,$xfn,$xi);
    $yfn = join "-",($yb,$ye);
    $lfra[$i-1] = join "\t",($yfr,$yfn,$yi);
   }
   if($lgap<0){							# need to resolve
   	#printf STDERR "# %d-%d\t%d\n\t%s\n\t%s\n",$i,$i-1,$lgap,$lfra[$i-1],$lfra[$i];
    my $ov = calc_overlap($xb,$xe,$yb,$ye);				# test for special case of pathological overlap
   	#printf STDERR "\t%d\n",$ov;
    my $lx = $xe - $xb + 1; my $ly = $ye - $yb + 1;
    my $ox = 1; $ox = $ov/$lx if($lx>0);				# relative overlap on x
    my $oy = 1; $oy = $ov/$ly if($ly>0);				# relative overlap on y
    $xden = 0 if($ox>=0.5 and $oy<0.5);					# relative overlap on x too much, supress
    $yden = 0 if($oy>=0.5 and $ox<0.5);					# relative overlap on y too much, supress
   	#printf STDERR "\t%.3f\t%.4f\n\t%.3f\t%.4f\n",$ox,$xden,$oy,$yden;
    $lgap = $ye - $xb + 1;
    my $ysub = int($lgap*$xden/($yden+$xden)+0.5);
    my $xadd = $lgap - $ysub;
    $ye -= $ysub;
    $xb += $xadd;
    $xfn = join "-",($xb,$xe);
    $lfra[$i] = join "\t",($xfr,$xfn,$xi);
    $yfn = join "-",($yb,$ye);
    $lfra[$i-1] = join "\t",($yfr,$yfn,$yi);
   }
  }
  if($i==@lfra-1){					# at the end
   if($lpro-$xe<=$xwin){
    $xe = $lpro;
    $xfn = join "-",($xb,$xe);
    $lfra[$i] = join "\t",($xfr,$xfn,$xi);
   }
  }
 }

	#for(my $i=0;$i<@lfra;$i++){
	# printf STDERR "# %d\t%s\n",$i,$lfra[$i];
	#}

#---	reconstruct paths ----------------------------------
 for(my $i=0;$i<@lfra;$i++){
  my ($xfr,$xfn,$xi) = split/\t/,$lfra[$i];
  my ($xb,$xe) = split/-/,$xfn;
  next if($xe-$xb<0);
  $$rnew[$xi-1] .= "=" if($$rnew[$xi-1] ne "");
  $$rnew[$xi-1] .= $xfn;
 }
	#for(my $i=0;$i<@$rnew;$i++){
	# printf STDERR "// %d\t%s\n",$i,$$rnew[$i];
	#}

#---	desegment paths ------------------------------------
 for(my $i=0;$i<@$rnew;$i++){
  my $pats = blast_smooth_path($$rnew[$i],$lpro,0,0);
  $$rnew[$i] = $pats;
 }
}

############################################################
#	blast_process_nthits($ql,$hl,$rhit,$orth,$oath)
############################################################
# reads array of hits for query-subject pair as read by blast_read_btab()
#	$ss,$qb,$qe,$hb,$he,$ev
# sorts @$rhit by score
# in @$rhit regularizes ($qb,$qe) and ($hb,$he) pairs; adds ($qdir,$hdir)
# returns ($scor,$patq,$path,$qple,$hple) array
# paths are comma-delimited
# paths are unresolved w.r.t. conflicts, ignore alignable compatibility, possibly inverted
# segmented paths are segmented in both query and subject
sub blast_process_nthits
{
 my $ql = shift;
 my $hl = shift;
 my $rhit = shift;
 my $orth = shift;
 my $oath = shift;

#---	recode strand --------------------------------------
 for(my $i=0;$i<@$rhit;$i++){
  my ($ss,$qb,$qe,$hb,$he,$ev) = split/\t/,$$rhit[$i];
  my $qdir = 1;
  if($qb>$qe){
   my $tmp = $qb; $qb = $qe; $qe = $tmp;
   $qdir = -1;
  }
  my $hdir = 1;
  if($hb>$he){
   my $tmp = $hb; $hb = $he; $he = $tmp;
   $hdir = -1;
  }
  $$rhit[$i] = join "\t",($ss,$qb,$qe,$hb,$he,$ev,$qdir,$hdir);
 }

#---	sort by score --------------------------------------
 @$rhit = sort {$b<=>$a} @$rhit;

#---	greedily process hits ------------------------------
 my @lhit = ();
 for(my $i=0;$i<@$rhit;$i++){
  my ($ssx,$qbx,$qex,$hbx,$hex) = split/\t/,$$rhit[$i];
  my $lqx = $qex - $qbx + 1;
  my $lhx = $hex - $hbx + 1;
  my $good = 1;
  for(my $j=0;$j<@lhit;$j++){
   my ($ssy,$qby,$qey,$hby,$hey) = split/\t/,$lhit[$j];
   my $lqy = $qey - $qby + 1;
   my $lhy = $hey - $hby + 1;
   my $ovq = calc_overlap($qbx,$qex,$qby,$qey);
   if($ovq/min($lqx,$lqy)>$orth or $ovq>$oath){
    $good = 0; last;
   }
   my $ovh = calc_overlap($hbx,$hex,$hby,$hey);
   if($ovh/min($lhx,$lhy)>$orth or $ovh>$oath){
    $good = 0; last;
   }
  }
  next unless($good);
  push @lhit,$$rhit[$i];
 }

#---	sort by position -----------------------------------
 for(my $i=0;$i<@lhit;$i++){
  my ($ss,$qb,$qe,$hb,$he,$ev,$qdir,$sdir) = split/\t/,$lhit[$i];
  $lhit[$i] = join "\t",($qb,$qe,$hb,$he,$ev,$ss,$qdir,$sdir);
 }
 @lhit = sort {$a<=>$b} @lhit;

#---	make paths and score -------------------------------
 my $ssco = 0;
 my $qpat = ""; my $qpal = 0;
 my $hpat = ""; my $hpal = 0;
 for(my $i=0;$i<@lhit;$i++){
  my ($qb,$qe,$hb,$he,$ev,$ss,$qdir,$hdir) = split/\t/,$lhit[$i];
  my $lq = $qe - $qb + 1;
  $qpal += $lq;
  my $qseg = ($qdir>0)?($qb."-".$qe):($qe."-".$qb);
  $qpat .= "," if($qpat ne "");
  $qpat .= $qseg;
  my $lh = $he - $hb + 1;
  $hpal += $lh;
  my $hseg = ($hdir>0)?($hb."-".$he):($he."-".$hb);
  $hpat .= "," if($hpat ne "");
  $hpat .= $hseg;
  $ssco += $ss;
 }
 $qpal = $ql if($qpal>$ql);
 $hpal = $hl if($hpal>$hl);

#---	return ---------------------------------------------
 return ($ssco,$qpat,$hpat,$qpal,$hpal);
}

############################################################
#	YIW::blast::aggregate_ranges($rori,$ragg)
############################################################
# reads an array of tab-delimited footprints p1,p2 (p1<=p2)
# sorts the original array
# merges overlapping footprints
# returns sorted array of non-overlapping footprints
sub aggregate_ranges
{
 my $rori = shift;
 my $ragg = shift;

 @$rori = sort {$a<=>$b} @$rori;			# sort query ranges
 
 push @$ragg,$$rori[0];					# 1st always
 for(my $i=1;$i<@$rori;$i++){				# from 2nd
  my $nmap = @$ragg;
  my ($p1,$p2) = split/\t/,$$ragg[$nmap-1];
  my ($q1,$q2) = split/\t/,$$rori[$i];
  if($q1-$p2<=1){						# overlap or adjacent, expand
   $p2 = max($p2,$q2);
   $$ragg[$nmap-1] = $p1."\t".$p2;
  }else{							# no overlap, new range
   push @$ragg,$$rori[$i];
  }
 }
}

