#!/usr/bin/perl

## 
## Copyright (c) 2014 Hannah Petersen
## 

# 4pi multiplicity
use Math::Trig;
 $pi=3.1415;
 $ptmin=0;
## this sets the pt bin width. lets try doubling it
## originally set at 0.2
 $dpt=0.4;
 $nch=0;
 $noe = 0;
while(<>) {
 @stuff=split(' ');

 if(/UQMD/){
 }
 elsif(/projectile:/){
  $protar=$_;
  $pro=$stuff[3];
 }
 elsif(/equation_of_state:/){
  $ecm=$stuff[5];
  $midy=log($ecm/0.938);
  $elab=$stuff[3];
  $plab=$stuff[7];
  $ee=$elab;
#  if ($ecm>30){goto ENDE}
  $noeold=1;
 }
 elsif($#stuff==1){
  if($stuff[0]>1){
    $noe++;                     # numbers of events: noe
  }
 }
 elsif($#stuff==1){
  $npart=$stuff[0];
  if($npart==2*$pro){
   $noe_el++;
  } 
 }
 else{
     ## loop over every line in the std input (i.e over all particles ever)
  if($#stuff==14||$#stuff==13){

   @x=@stuff[0..3];
   @p=@stuff[4..7];
   @id=@stuff[8..11];

   $mass=$id[0];
   $ityp=$id[1];
   $iz2=$id[2];
   $charge=$id[3];

# ab hier steht folgende info zur verfuegung:

# stuff[0] = r_0    = $x[0]
# stuff[1] = r_x    = $x[1]
# stuff[2] = r_y    .
# stuff[3] = r_z    .
# stuff[4] = p_0    = $p[0]
# stuff[5] = p_x    = $p[1]
# stuff[6] = p_y    .
# stuff[7] = p_z    .
# stuff[8] = mass   = $id[0]
# stuff[9] = ityp   = $id[1]
# stuff[10]= iz2    = $id[2]
# stuff[11]= charge = $id[3]
  
   $n++;
   $eta=99; 

   if($noe>$noeold){
       ## if this is a new event
      
    foreach($j=0;$j<=25;$j++){ ## loop over the pt bins
        ## what do we do per bin?
     if($nop[$j]!=0){
      $nch=0; ## reset the number of charged particles in this bin!
      $noebin[$j]++;  ## increment the number of events that contribute to bin j
      for($i=1;$i<=$nop[$j];$i++){ ## now loop over all the particles in bin j
       if($ch[$i][$j]!=0 && abs($rap[$i][$j])<1){ # check that this particle is safe to contribute
        $nch++; 
        foreach($j2=0;$j2<=25;$j2++){  ## over pt bins *again*
         if($nop[$j2]!=0){  
      
          for($k=1;$k<=$nop[$j2];$k++){  ## now lets take each particle in this pt bin and compute the
           if(abs($rap[$k][$j2])<2){     ## sum of pt*sin(2*phi) etc inside the rap 2 range
            if($j2==$j && $k==$i){ ## where we are explicitly excluding the pt bin from the outer loop

             next;
            } 
           
           $nloop++;
           $sinsum+=$r[$k][$j2]*sin(2*$phipart[$k][$j2]);
           $cossum+=$r[$k][$j2]*cos(2*$phipart[$k][$j2]);

           $trisinsum+=$r[$k][$j2]*sin(3*$phipart[$k][$j2]);
           $tricossum+=$r[$k][$j2]*cos(3*$phipart[$k][$j2]);
           }
          }
         }
        }

        ## so now we have the average vn contribs from the particles outside the pt bin we're looking at $j
        ## and inside the large rapidity window (2)
        ## 
        ## we use these to construct our event plane
        $sinmean=$sinsum/$nloop; 
        $cosmean=$cossum/$nloop;

        $psitwo=(atan2($sinmean,$cosmean))/2;

        $trisinmean=$trisinsum/$nloop;
        $tricosmean=$tricossum/$nloop;

        $psithree=(atan2($trisinmean,$tricosmean))/3;  
 

        ## now we can accumulate the contributions for the particles in our actual bin of interest $j
        $vtwoevsum[$j]+=cos(2*($phipart[$i][$j]-$psitwo)); 
        $vthreeevsum[$j]+=cos(3*($phipart[$i][$j]-$psithree));   
   
        $sinsum=0;
        $cossum=0;
        $trisinsum=0;
        $tricossum=0;
        $nloop=0; 
  
       }
      }
      if($nch!=0){
      $vtwoevmean[$j]+=$vtwoevsum[$j]/$nch;   
      $vthreeevmean[$j]+=$vthreeevsum[$j]/$nch;  
      $vtwoevmeanb[$j]+=($vtwoevsum[$j]/$nch)**2;   
      $vthreeevmeanb[$j]+=($vthreeevsum[$j]/$nch)**2;   
      }
     }
    }

    ## reset the event count
    $noeold=$noe;

    ## clear out the pt and phi bits (up to 100?)
    for($j=0;$j<=100;$j++){
     for ($i=1;$i<=$nop[$j];$i++){
    
      $phipart[$i][$j]=0; 
      $r[$i][$j]=0;
   
     }  
       $nop[$j]=0; 
       $vtwoevsum[$j]=0;
       $vthreeevsum[$j]=0;

    }  
       $ntot=0;
   }
   ## process particles per event
     $pabs=sqrt($p[1]**2+$p[2]**2+$p[3]**2); ## get |p|
 
     if (($pabs-$p[3])>0 && ($pabs+$p[3])>0){
      $eta=0.5*log(($pabs+$p[3])/($pabs-$p[3])); ## compute rap, check its sane
      
      $pt=sqrt($p[1]**2+$p[2]**2);      
      $ipt=int(($pt-$ptmin)/$dpt); ## get pt bin

      $nop[$ipt]++; ## increment count in pt bin
      $r[$nop[$ipt]][$ipt]=$pt;  ## store this particles actual pt
      $ch[$nop[$ipt]][$ipt]=$id[3]; ## store charge 
      $ntot++;      ## increment total 
      $rap[$nop[$ipt]][$ipt]=$eta; ## store rapidity

      if($p[1]>=0){ 
       $phipart[$nop[$ipt]][$ipt]=asin($p[2]/$pt); 
       if($p[1]==0&&$p[2]==0){
        $phipart[$nop[$ipt]][$ipt]=0;
       }
      }
      else{ ## store the angle phi
       $phipart[$nop[$ipt]][$ipt]=-asin($p[2]/$pt)+$pi; 
      } 
     
    } 
  }
 }
}
ENDE:
 

foreach($j=0;$j<=100;$j++){
 if($noebin[$j]>1){
  $meanvtwo[$j]=$vtwoevmean[$j]/$noebin[$j];  
  $meanvthree[$j]=$vthreeevmean[$j]/$noebin[$j];

  $meanvtwob[$j]=$vtwoevmeanb[$j]/$noebin[$j];  
  $meanvthreeb[$j]=$vthreeevmeanb[$j]/$noebin[$j];

 
  if((($meanvtwob[$j]-$meanvtwo[$j]**2)/($noebin[$j]-1))>=0){
  $errvtwo=sqrt(($meanvtwob[$j]-$meanvtwo[$j]**2)/($noebin[$j]-1));
  }
  else{
   $errvtwo=0;
  }
  if((($meanvthreeb[$j]-$meanvthree[$j]**2)/($noebin[$j]-1))>=0){
  $errvthree=sqrt(($meanvthreeb[$j]-$meanvthree[$j]**2)/($noebin[$j]-1));
  }
  else{
   $errvthree=0;
  }

  print($j*$dpt+$ptmin+$dpt/2,"\t",$noebin[$j], "\t", $meanvtwo[$j],"\t",$errvtwo,"\t",$meanvthree[$j],"\t",$errvthree,"\n");
 }
}
