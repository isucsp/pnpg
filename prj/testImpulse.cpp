#include <cmath>
#include <cstdio>
#define EPS (1e-9)

int main(){
    int nr=1024, nc=1024; //The origin is top-left corner
    const int xC=0, yC=0;
    int theta=89;
    double beamWidth=1;

    theta%=180; // test value of theta
    int t1, t2;
    double height;
    const double pi=3.14159265;
    if(theta<45){
        t1=theta+45; t2=45-theta;
        height=1/cos(theta*pi/180);
    }else if(theta<90){
        t1=135-theta; t2=theta-45;
        height=1/sin(theta*pi/180);
    }else if(theta<135){
        t1=theta-45; t2=135-theta;
        height=1/sin(theta*pi/180);
    }else{
        t1=225-theta; t2=theta-135;
        height=-1/cos(theta*pi/180);
    }
    double bias1=sqrt(2)/2*cos(t1*pi/180);
    double bias2=sqrt(2)/2*cos(t2*pi/180);
    double areaT=0.5*height*(bias2-bias1);
    double areaR=height*bias1*2;

    printf("#dist\tweight\n");

    int xIdx,yIdx;
    double hbeamW=beamWidth/2;
    double weight,d,dist;
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++){
            xIdx=j-512; yIdx=511-i;
            dist=xIdx*cos(theta*pi/180)+yIdx*sin(theta*pi/180);
            double temp=floor(dist+bias2+beamWidth/2-EPS);
            for(int k=ceil(dist-bias2-beamWidth/2+EPS); k<=temp; k++){
                d=dist-k;
                if(k-hbeamW<=dist-bias2)
                    if(k+hbeamW<=dist-bias1)
                        weight=areaT*pow((k+hbeamW-(dist-bias2))/(bias2-bias1),2);
                    else if(k+hbeamW<=dist+bias1)
                        weight=areaT+height*(k+hbeamW-(dist-bias1));
                    else
                        weight=areaT*2+areaR-\
                               areaT*pow((dist+bias2-k-hbeamW)/(bias2-bias1),2);
                else if(k-hbeamW<=dist-bias1)
                    if(k+hbeamW<=dist-bias1)
                        weight=height/(bias2-bias1)*(k-dist+bias2)*beamWidth;
                    else if(k+hbeamW<=dist+bias1)
                        weight=areaT*(1-pow((k-hbeamW-(dist-bias2))/(bias2-bias1),2))+\
                               height*(k+hbeamW-(dist-bias1));
                    else if(k+hbeamW<=dist+bias2)
                        weight=areaT*(1-pow((k-hbeamW-(dist-bias2))/(bias2-bias1),2))+\
                               areaR+\
                               areaT*(1-pow((dist+bias2-k-hbeamW)/(bias2-bias1),2));
                    else
                        weight=areaT*(1-pow((k-hbeamW-(dist-bias2))/(bias2-bias1),2))+\
                               areaR+areaT;
                else if(k-hbeamW<=dist+bias1)
                    if(k+hbeamW<=dist+bias1)
                        weight=height*beamWidth;
                    else if(k+hbeamW<=dist+bias2)
                        weight=height*(dist+bias1-k+hbeamW)+\
                               areaT*(1-pow((dist+bias2-k-hbeamW)/(bias2-bias1),2));
                    else
                        weight=height*(dist+bias1-k+hbeamW)+areaT;
                else{
                    if(k+hbeamW<=dist+bias2)
                        weight=(height)/(bias2-bias1)*(dist+bias2-k)*beamWidth;
                    else
                        weight=0.5*height/(bias2-bias1)*pow(dist+bias2-k+hbeamW,2);
                }
                printf("%e\t%e\n",d,weight);
            }
        }
    return 0;
}
