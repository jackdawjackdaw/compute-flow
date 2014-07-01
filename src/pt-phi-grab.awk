## run as awk -f pt-phi-grab.awk < stuff
## match event lines and print something
BEGIN{count = 0;}
/^event#/{ count++;
    print( "# event", count)
}

#!/(^[a-z])/{print $0}

function abs(value)
{
  return (value<0?-value:value);
}


function asin(x) { return atan2(x, sqrt(1-x*x)) }

## on lines which are particles output pt, phi and rap
!/(^[a-zA-Z])/{ ## forgot to screen out upper case chars (dumb)
    if(NF > 10){ ## ok so only go on potential particle lines
        en = $5
        px = $6;
        py = $7;
        pz = $8;
        ch = $12;
        pabs = sqrt(px*px + py*py + pz*pz)
        ##printf("(%lf %lf %lf) %lf %lf\n", px, py, pz, en, pabs)
        if(abs(pabs) > 0 && abs(pz) > 0 && abs(pabs-pz) > 0 && abs(pabs+pz)>0){
            rap = 0.5*log((pabs + pz)/(pabs - pz))
        } else {
            rap = 0;
        }
        pt = sqrt(px*px + py*py);
        ## hp uses asin(p[2]/pt), or minus, depending on px
        #phi = atan2(px, py);
        
        #printf("## %lf %lf %lf\n", px, py, pz)
        if(px == 0 && py == 0 || pt == 0){ 
            phi = 0
        } else {
            if(px > 0){
                phi = asin(py/pt)
            } else {
                phi = -asin(py/pt) + 3.14
            }
        }

        printf("%d %lf %lf %lf\n", ch, pt, phi, rap);
    }
}

