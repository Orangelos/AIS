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
    double hess[2][2];
    double dk[2];
    x00=x0; //k
    x001=x01;
    x5=x0;   //k-1
    x6=x01;

    
    m2:
    x11=dx1*x00+x1x2*x001;  //f(x)k
    x12=dx2*x001+x1x2*x00;
    printf("k = %d  ------  f(x)=[%f;%f]   xk= [%f;%f]   -----------------\n",k,x11,x12,x00,x001);
    
    
    if((sqrt(x11*x11+x12*x12)<e1)){
        x0=x00;
        x01=x001;
        printf("Градиент xk < e1 \n");
    }
    else {
        printf("Градиент xk > e1 \n");
        printf("Вычислим ||f(x)|| > e =   %f > %f\n",sqrt(x11*x11+x12*x12),e1);
        if(k>=m)
        {
            x0=x00;
            x01=x001;
           
        }
        else
        {
            
            hess[0][0]=dx1;
            hess[0][1]=x1x2;
            hess[1][0]=x1x2;
            hess[1][1]=dx2;
            
            double opr =hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0];
            double trans[2][2];
            trans[0][0]=hess[1][1];
            trans[0][1]=-hess[1][0];
            trans[1][0]=-hess[0][1];
            trans[1][1]=hess[0][0];
            printf("obr[0][0] = %f obr[0][1]= %f obr[1][0]= %f obr[1][1] = %f\n",trans[0][0],trans[0][1],trans[1][0],trans[1][1]);
            
            hess[0][0]=trans[0][0]*(1/fabs(opr));
            hess[0][1]=trans[1][0]*(1/fabs(opr));
            hess[1][0]=trans[0][1]*(1/fabs(opr));
            hess[1][1]=trans[1][1]*(1/fabs(opr));
            printf("hess[0][0] = %f hess[0][1]= %f hess[1][0]= %f hess[1][1] = %f\n",hess[0][0],hess[0][1],hess[1][0],hess[1][1]);
            if((hess[1][1]>0)&&((hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0])>0)){
                dk[0]=-1*(hess[0][0]*x11+hess[0][1]*x11);
                dk[1]=-1*(hess[1][0]*x12+hess[1][1]*x12);
                t=1;
                x21=x00+dk[0];
                x22=x001+dk[1];
                
            }
            else
            {
                dk[0]=-1*x11;
                dk[1]=-1*x12;
                m1:
                x21=x00-t*dk[0];
                x22=x001-t*dk[1];
                if((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)<0);
                    else
                    {
                        printf("lwkjebvw\n");
                        t=t/2;
                        goto m1;
                    }
                
            }
            
            if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)))<e2)&&(sqrt((x00-x5)*(x00-x5)+(x001-x6)*(x001-x6))<e2) &&  (abs(((x1*x00*x00+x1x2*x00*x001+x2*x001*x001)-(x1*x5*x5+x1x2*x5*x6+x2*x6*x6)))<e2)){
                printf("Проверим на выполение обоих условий для k и k=k-1 ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %.5fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
                printf("Условие выполнено \n");
                x0=x21;
                x01=x22;
            }
                else{
                    printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %f \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
                    printf("Условие не выполнено \n");
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
        
    }
    
    printf("Ответ x* = [%lf;%lf]     f(x*)=%lf \n",x0,x01,x1*x0*x0+x1x2*x0*x01+x2*x01*x01);
}
