// credits: Robert Gerbicz

long long int gcd(long long int a,long long int b)
{  // lnko(a,b)-vel tér vissza
// az osztások nagy részét elkerülő gyors algoritmus, bár ezt
// most legfeljebb 2-szer számolunk így a sebessége nem lényeges
 if(a<0) a=-a;
 if(b<0) b=-b;
 if(a==0)  return b;
 if(b==0)  return a;

 long long int c;

 while(b>0)  {
    if(a>=b)  {
       a-=b;
       if(a>=b)  {
          a-=b;
          if(a>=b)  {
             a-=b;
             if(a>=b)  {
                a-=b;
                if(a>=b)  {
                   a-=b;
                   if(a>=b)  {
                      a-=b;
                      if(a>=b)  {
                         a-=b;
                         if(a>=b)  {
                            a-=b;
                            if(a>=b)  a%=b;
             }}}}}}}}
    c=a,a=b,b=c;
 }
 return a;
}
