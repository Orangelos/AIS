#include <stdio.h>
#include<math.h>


int main(void){
    double x1,x2,x1x2;
    printf("Введите коэфиценты перед x1^2 x1x2 x2^2\n");
    scanf("%lf%lf%lf",&x1,&x1x2,&x2);
    printf("Уровнение =   %.1lf*x^2 %.1lf*x1x2 %.1lf*x2^2\n",x1,x1x2,x2);
    double e1,e2,m,x0,x01;
    printf("Введите x0 e1 e2 M\n");
    scanf("%lf %lf %lf %lf %lf",&x0,&x01,&e1,&e2,&m);
    
    
    double dx1,dx2;
    dx1=2*x1;
    dx2=2*x2;
    int k=0;
    double t=0.5;
    double x11,x12,x21,x22,x00,x001=0,x5,x6;
    
    x00=x0; //k
    x001=x01;
    x5=x0;   //k-1
    x6=x01;
    m2:
    x11=dx1*x00+x1x2*x001;  //f(x)k
    x12=dx2*x001+x1x2*x00;
    printf("k = %d  ------  gf(x)=[%.5fl;%.5fl]   xk= [%.5fl;%.5fl]   f(xk)=%f  -----------------\n",k,x11,x12,x00,x001,x1*x00*x00+x1x2*x00*x001+x2*x001*x001);
    
    
    if((sqrt(x11*x11+x12*x12)<e1)){
        x0=x00;
        x01=x001;
    }
    else {
        printf("Вычислим ||gf(x)|| > e =   %.5fl > %.5fl\n",sqrt(x11*x11+x12*x12),e1);
        if(k>=m)
        {
            x0=x00;
            x01=x001;
            goto m3;
        }
        
        
        
        
        
        m1:

        
        x21=x00-t*x11;  //k+1
        x22=x001-t*x12;
        printf("Вычислим xk+1= [%.5fl;%fl]   f(xk+1)= %.5fl    t = %.5fl\n",x21,x22,(x1*x21*x21+x1x2*x21*x22+x2*x22*x22),t);
        if((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)<0)//||((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))<-e1*sqrt(x11*x11+x12*x12)*sqrt(x11*x11+x12*x12)){
        {
            
            
          if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)))<e2)&&(sqrt((x00-x5)*(x00-x5)+(x001-x6)*(x001-x6))<e2) &&  (abs(((x1*x00*x00+x1x2*x00*x001+x2*x001*x001)-(x1*x5*x5+x1x2*x5*x6+x2*x6*x6)))<e2))
          {
              printf("Проверим на выполение обоих условий для k и k=k-1 ||xk+1 - xk ||= %.5fl > %.5fl  ------   |f(xk+1) -  f(xk) |=%fl > %.5fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
              x0=x21;
              x01=x22;
              printf("ERRR\n");
          }
            else
            {
                printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %.5fl > %.5fl  ------   |f(xk+1) -  f(xk) |=%.5fl > %.5fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
                k=k+1;
                x00=x21;
                x001=x22;
                x5=x00;
                x6=x001;
                printf("\n");
                printf("\n");
                printf("\n");
                goto m2;
                
            }
                
                }
        
        
        else{
            t=t/2;
            goto m1;
        }
        
    }
    
    m3:
    printf("Ответ x* = [%lf;%lf]     f(x*)=%fl \n",x0,x01,x1*x0*x0+x1x2*x0*x01+x2*x01*x01);
    
    
    
    
}
