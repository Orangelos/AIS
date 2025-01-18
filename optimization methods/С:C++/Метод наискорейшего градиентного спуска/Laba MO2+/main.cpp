#include <stdio.h>
#include<math.h>



double proizv(double x1, double x2, double x3)
{printf("Функция получила x1= %fl x2=%fl  x3= %fl \n",x1,x2,x3);
    double t;
    t=-x2/(2*x1);
    if(x1>0){
        return t;
    }
    else
        return 0.5;
}




double zoloto(double x1, double x2, double x3, double a, double b, double e) {
   // printf("Функция получила x1= %fl x2=%fl  x3= %fl \n",x1,x2,x3);
   // printf("------------------------------------------------------------------------------------------------------------------\n");
    //printf("                                        МЕТОД ЗОЛОТОГО СЕЧЕНИЯ                                        \n");
    double l = 0;
    int k = 0;
    double y = a + ((3 - sqrt(5)) / 2) * (b - a), z = a + b - y;
m1:
    k += 1;
    double fy = x1 * y * y + x2 * y + x3;
    double fz = x1 * z * z + x2 * z + x3;
    
    if (fy <= fz) {
       // printf("f(y)<=f(x)   ");
        b = z;
        z = y;
        y = a + b - y;
    }
    else
    {
        
        a = y;
        y = z;
        z = a + b - z;
    }
    l = abs(a - b);
    if (l <= e) {
      // printf("ОТВЕТ \n");
       printf("x -- [%.0f ; %.0f]  x*=%.0f  f(x)=%.0f k = %d\n ", a, b, ((a+b)/2), x1 * ((a + b) / 2) * ((a + b) / 2) * ((a + b) / 2) + x2 * ((a + b) / 2) * ((a + b) / 2) + x3,k-1);
        return ((a+b)/2);
    }
    else {
       printf("PROM RESULT  //////////x -- [%lf ; %lf]  x*=%lf   f(x)=%lf |    Шаг k = %d     ||      y=%lf z=%lf f(y)=%lf f(z)=%lf////////////\n", a, b, ((a + b) / 2), x1 * ((a + b) / 2) * ((a + b) / 2) * ((a + b) / 2) + x2 * ((a + b) / 2) * ((a + b) / 2) + x3,k,y, z, fy, fz, l);
        goto m1;
    }
    
  //  printf("------------------------------------------------------------------------------------------------------------------\n");
    return 0;
}


int main(void){
    double x1,x2,x1x2;
    printf("Введите коэфиценты перед x1^2 x1x2 x2^2\n");
    scanf("%lf%lf%lf",&x1,&x1x2,&x2);
    printf("Уровнение =   %fl*x^2 %fl*x1x2 %fl*x2^2\n",x1,x1x2,x2);
    double e1,e2,m,x0,x01;
    printf("Введите x0 e1 e2 M\n");
    scanf("%lf %lf %lf %lf %lf",&x0,&x01,&e1,&e2,&m);
    
    double y1,y2,y3;
    double dx1,dx2;
    dx1=2*x1;
    dx2=2*x2;
    int k=0;
    double t=0.5;
    double x11,x12,x21,x22,x00,x001=0;
    x00=x0; //k
    x001=x01;
    m2:
    x11=dx1*x00+x1x2*x001;  //f(x)k
    x12=dx2*x001+x1x2*x00;
    printf("k = %d  ------  f(x)=[%fl;%fl]   xk= [%fl;%fl]   -----------------\n",k,x11,x12,x00,x001);
    
    
    if((sqrt(x11*x11+x12*x12)<e1)){
        x0=x00;
        x01=x001;
    }
    else {
        printf("Вычислим ||f(x)|| > e =   %fl > %fl\n",sqrt(x11*x11+x12*x12),e1);
        if(k>=m)
        {
            x0=x00;
            x01=x001;
            goto m3;
        }
        
        
        
        
        
        
        y3=x1*x00*x00+x1*x2*x00*x001+x2*x001*x001;
        y1=x1x2*x11*x12+x11*x11*x1+x12*x12*x2;
        y2=-2*x11*x00*x1-2*x12*x001-x1x2*x00*x12-x1x2*x001*x11;
        
        t=proizv(y1, y2, y3);
        
        

        printf("t = %fl ---------------- \n",t);
        
        x21=x00-t*x11;  //k+1
        x22=x001-t*x12;
          if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)))<e2))
          {
              printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %fl > %fl  ------   |f(xk+1) -  f(xk) |=%fl > %fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
              x0=x21;
              x01=x22;
              printf("ERRR\n");
          }
            else
            {
                printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %fl > %fl  ------   |f(xk+1) -  f(xk) |=%fl > %fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
                k=k+1;
                x00=x21;
                x001=x22;
                printf("\n");
                printf("\n");
                printf("\n");
                goto m2;
                
            }
                
                
        
        
       
        
    }
    
    m3:
    printf("Ответ x* = [%lf;%lf]     f(x*)=%fl \n",x0,x01,x1*x0*x0+x1x2*x0*x01+x2*x01*x01);
    
    
    
    
}

