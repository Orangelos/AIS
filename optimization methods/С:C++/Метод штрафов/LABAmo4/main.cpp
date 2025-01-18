#include<stdio.h>
#include<math.h>
//double* minfunc(double x222,double y22, double s, double y, double x,double xy, double* xt ){
//    printf("MINFUNC\n");
//    printf("\n");
//    double x1=x222, x2=y, x1x2=xy;
//    double x0=xt[0],x01=xt[1];
//    double e1=0.15,e2=0.2,m=10;
//    
//    
//    
//    double dx1,dx2;
//    dx1=2*x1;
//    dx2=2*x2;
//    int k=0;
//    double t=0.5;
//    double x11,x12,x21,x22,x00,x001=0,x5,x6;
//    double hess[2][2];
//    double dk[2];
//    x00=x0; //k
//    x001=x01;
//    x5=x0;   //k-1
//    x6=x01;
//
//    
//    m2:
//    x11=dx1*x00+x1x2*x001;  //f(x)k
//    x12=dx2*x001+x1x2*x00;
//    printf("k = %d  ------  f(x)=[%f;%f]   xk= [%f;%f]   -----------------\n",k,x11,x12,x00,x001);
//    
//    
//    if((sqrt(x11*x11+x12*x12)<e1)){
//        x0=x00;
//        x01=x001;
//        printf("Градиент xk < e1 \n");
//    }
//    else {
//        printf("Градиент xk > e1 \n");
//        printf("Вычислим ||f(x)|| > e =   %f > %f\n",sqrt(x11*x11+x12*x12),e1);
//        if(k>=m)
//        {
//            x0=x00;
//            x01=x001;
//           
//        }
//        else
//        {
//            
//            hess[0][0]=dx1;
//            hess[0][1]=x1x2;
//            hess[1][0]=x1x2;
//            hess[1][1]=dx2;
//            
//            double opr =hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0];
//            double trans[2][2];
//            trans[0][0]=hess[1][1];
//            trans[0][1]=-hess[1][0];
//            trans[1][0]=-hess[0][1];
//            trans[1][1]=hess[0][0];
//           
//            
//            hess[0][0]=trans[0][0]*(1/fabs(opr));
//            hess[0][1]=trans[1][0]*(1/fabs(opr));
//            hess[1][0]=trans[0][1]*(1/fabs(opr));
//            hess[1][1]=trans[1][1]*(1/fabs(opr));
//            printf("hess[0][0] = %f hess[0][1]= %f hess[1][0]= %f hess[1][1] = %f\n",hess[0][0],hess[0][1],hess[1][0],hess[1][1]);
//            if((hess[1][1]>0)&&((hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0])>0)){
//                dk[0]=-1*(hess[0][0]*x11+hess[0][1]*x11);
//                dk[1]=-1*(hess[1][0]*x12+hess[1][1]*x12);
//                t=1;
//                x21=x00+dk[0];
//                x22=x001+dk[1];
//                
//            }
//            else
//            {
//                dk[0]=-1*x11;
//                dk[1]=-1*x12;
//                m1:
//                x21=x00+t*dk[0];
//                x22=x001+t*dk[1];
//                if((x1*x21*x21+x1x2*x21*x22+x2*x22*x22+x*x21+y*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001+x*x00+y*x001)<0);
//                    else
//                    {
//                        printf("lwkjebvw\n");
//                        t=t/2;
//                        goto m1;
//                    }
//                
//            }
//            
//            if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001)))<e2)&&(sqrt((x00-x5)*(x00-x5)+(x001-x6)*(x001-x6))<e2) &&  (abs(((x1*x00*x00+x1x2*x00*x001+x2*x001*x001)-(x1*x5*x5+x1x2*x5*x6+x2*x6*x6)))<e2)){
//                printf("Проверим на выполение обоих условий для k и k=k-1 ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %.5fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
//                printf("Условие выполнено \n");
//                x0=x21;
//                x01=x22;
//            }
//                else{
//                    printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %f \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
//                    printf("Условие не выполнено \n");
//                    k=k+1;
//                    x00=x21;
//                    x001=x22;
//                    x5=x00;
//                    x6=x001;
//                    printf("\n");
//                    printf("\n");
//                    printf("\n");
//                    goto m2;
//                }
//                
//            
//            
//            
//            
//        }
//        
//    }
//    double dx1,dx2;
//    dx1=2*x1;
//    dx2=2*x2;
//    int k=0;
//    double t=0.5;
//    double x11,x12,x21,x22,x00,x001=0,x5,x6;
//    double hess[2][2];
//    double dk[2];
//    x00=x0; //k
//    x001=x01;
//    x5=x0;   //k-1
//    x6=x01;
//    
//    
//  m2:
//    x11=dx1*x00+x1x2*x001;  //f(x)k
//    x12=dx2*x001+x1x2*x00;
//    printf("k = %d  ------  f(x)=[%f;%f]   xk= [%f;%f]   -----------------\n",k,x11,x12,x00,x001);
//    
//    
//    if((sqrt(x11*x11+x12*x12)<e1)){
//        x0=x00;
//        x01=x001;
//        printf("Градиент xk < e1 \n");
//    }
//    else {
//        printf("Градиент xk > e1 \n");
//        printf("Вычислим ||f(x)|| > e =   %f > %f\n",sqrt(x11*x11+x12*x12),e1);
//        if(k>=m)
//        {
//            x0=x00;
//            x01=x001;
//            
//        }
//        else
//        {
//            printf("hess  dx1=%lf dx2=%lf =x222=%lf   y22=%lf     s= %lf    y=%lf     x=%lf     xy= %lf    xt=[%lf ;%lf ]\n",dx1,dx2,x222,y22,s,y,x,x1x2,xt[0],xt[1]);
//            hess[0][0]=dx1;
//            hess[0][1]=x1x2;
//            hess[1][0]=x1x2;
//            hess[1][1]=dx2;
//            
//            double opr =hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0];
//            double trans[2][2];
//            trans[0][0]=hess[1][1];
//            trans[0][1]=-hess[1][0];
//            trans[1][0]=-hess[0][1];
//            trans[1][1]=hess[0][0];
//            printf("OPR = %lf\n", opr);
//            
//            hess[0][0]=trans[0][0]*(1/fabs(opr));
//            hess[0][1]=trans[1][0]*(1/fabs(opr));
//            hess[1][0]=trans[0][1]*(1/fabs(opr));
//            hess[1][1]=trans[1][1]*(1/fabs(opr));
//            
//            
//            printf("hess[0][0] = %lf hess[0][1]= %lf hess[1][0]= %lf hess[1][1] = %lf\n",hess[0][0],hess[0][1],hess[1][0],hess[1][1]);
//            if((hess[1][1]>0)&&((hess[0][0]* hess[1][1]-hess[0][1]*hess[1][0])>0)){
//                dk[0]=-1*(hess[0][0]*x11+hess[0][1]*x11);
//                dk[1]=-1*(hess[1][0]*x12+hess[1][1]*x12);
//                t=1;
//                x21=x00+dk[0];
//                x22=x001+dk[1];
//                
//            }
//            else
//            {
//                dk[0]=-1*x11;
//                dk[1]=-1*x12;
//            m1:
//                x21=x00-t*dk[0];
//                x22=x001-t*dk[1];
//                if(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))<0);
//                else //+x*x21+y*x22         x*x00+y*x001
//                {
//                   // printf("t\n");
//                    t=t/2;
//                    goto m1;
//                }
//                
//            }
//            
//           // if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22+x*x21+y*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001+x*x00+y*x001)))<e2)&&(sqrt((x00-x5)*(x00-x5)+(x001-x6)*(x001-x6))<e2) &&  (abs(((x1*x00*x00+x1x2*x00*x001+x2*x001*x001)-(x1*x5*x5+x1x2*x5*x6+x2*x6*x6)))<e2)){
//            if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001)))<e2)){
//                printf("Проверим на выполение обоих условий ");//для k и k=k-1 ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %.5fl \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
//                printf("Условие выполнено \n");
//                x0=x21;
//                x01=x22;
//            }
//            else{
//                //printf("Проверим на выполение обоих условий ||xk+1 - xk ||= %f > %f  ------   |f(xk+1) -  f(xk) |=%f > %f \n",sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001)),e2,abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))),e2);
//                printf("Условие не выполнено \n");
//                k=k+1;
//                x00=x21;
//                x001=x22;
//                x5=x00;
//                x6=x001;
//                printf("\n");
//                printf("\n");
//                printf("\n");
//                goto m2;
//            }
////        }
////    }
//    printf("Ответ x* = [%lf;%lf]     f(x*)=%lf \n",x0,x01,x1*x0*x0+x1x2*x0*x01+x2*x01*x01);
//    xt[0]=x0;
//    xt[1]=x01;
//    return xt;
//}
    
  //  printf("Ответ x* = [%lf;%lf]     f(x*)=%lf \n",x0,x01,x1*x0*x0+x1x2*x0*x01+x2*x01*x01);

double* minfunc(double x222,double y22, double s, double y, double x,double xy, double* xt ){


    
    
    
    
        double x1=x222,x2=y22,x1x2=xy;
//        printf("Введите коэфиценты перед x1^2 x1x2 x2^2\n");
//        scanf("%lf%lf%lf",&x1,&x1x2,&x2);
        printf("Уровнение =   %.1lf*x^2 %.1lf*x1x2 %.1lf*x2^2\n",x1,x1x2,x2);
        double e1=0.15,e2=0.2,m=10,x0=xt[0],x01=xt[1];

        
        
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
        x11=dx1*x00+x1x2*x001+x;  //f(x)k
        x12=dx2*x001+x1x2*x00+y;
        printf("k = %d  ------  gf(x)=[%.5fl;%.5fl]   xk= [%.5fl;%.5fl]   f(xk)=%f  -----------------\n",k,x11,x12,x00,x001,x1*x00*x00+x1x2*x00*x001+x2*x001*x001);
        
        
        if((sqrt(x11*x11+x12*x12)<e1)){
            x0=x00;
            x01=x001;
        }
        else {
            //printf("Вычислим ||gf(x)|| > e =   %.5fl > %.5fl\n",sqrt(x11*x11+x12*x12),e1);
            if(k>=m)
            {
                x0=x00;
                x01=x001;
                goto m3;
            }
            
            
            
            
            
            m1:

            
            x21=x00-t*x11;  //k+1
            x22=x001-t*x12;
          //  printf("Вычислим xk+1= [%.5fl;%fl]   f(xk+1)= %.5fl    t = %.5fl\n",x21,x22,(x1*x21*x21+x1x2*x21*x22+x2*x22*x22),t);
            if((x1*x21*x21+x1x2*x21*x22+x2*x22*x22+x*x21+y*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001+x*x00+y*x001)<0)//||((x1*x21*x21+x1x2*x21*x22+x2*x22*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001))<-e1*sqrt(x11*x11+x12*x12)*sqrt(x11*x11+x12*x12)){
            {
                
                
              if(  (sqrt((x21-x00)*(x21-x00)+(x22-x001)*(x22-x001))<e2) &&  (abs(((x1*x21*x21+x1x2*x21*x22+x2*x22*x22+x*x21+y*x22)-(x1*x00*x00+x1x2*x00*x001+x2*x001*x001+x*x00+y*x001)))<e2)&&(sqrt((x00-x5)*(x00-x5)+(x001-x6)*(x001-x6))<e2) &&  (abs(((x1*x00*x00+x1x2*x00*x001+x2*x001*x001+x*x00+y*x001)-(x1*x5*x5+x1x2*x5*x6+x2*x6*x6+x*x5+y*x6)))<e2))
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
    
    xt[0]=x0;
    xt[1]=x01;
    printf("x1=%lf x2=%lf\n",x0,x01);
    return xt;
   
}
    
    
    
double* Func(double x1,double x2,double sv,double g1,double g2,double svg, double* xt,double r){
    double x,y,xy,s,x22,y22;
    x22=g1*g1;
    y22=g2*g2 ;
    s=svg*svg ;
    y=2*g2*svg;
    x=2*g1*svg;
    xy=2*g1*g2;

    
    if((x22*xt[0]+y22*xt[1]+x*xt[0]+y*xt[1]+xy*xt[0]*xt[1]+s)>0){
        x22=x1+g1*g1*r;
        y22=x2+g2*g2*r ;
         s=sv+svg*svg*r ;
         y=2*g2*svg*r;
         x=2*g1*svg*r;
         xy=2*g1*g2*r;
        return minfunc(x22,y22,s,y,x,xy,xt);
    }
    else
    {
         x22=x1+g1*g1*r/2;
         y22=x2+g2*g2*r/2 ;
          s=sv+svg*svg*r/2 ;
          y=2*g2*svg*r/2;
          x=2*g1*svg*r/2;
          xy=2*g1*g2*r/2;
        return minfunc(x22,y22,s,y,x,xy,xt);
        
        
        
    }
  //  double x,y,xy,s,x22,y22;
  // x22=x1+g1*g1*r;
  // y22=x2+g2*g2*r ;
  //  s=sv+svg*svg*r ;
  //  y=2*g2*svg*r;
  //  x=2*g1*svg*r;
  //  xy=2*g1*g2*r;
  //  return minfunc(x22,y22,s,y,x,xy,xt);

}
    
    
    
double Pfunc(double x1,double x2,double sv,double g1,double g2,double svg, double* xt,double r){
    double x,y,xy,s,x22,y22;
    x22=g1*g1;
    y22=g2*g2 ;
    s=svg*svg ;
    y=2*g2*svg;
    x=2*g1*svg;
    xy=2*g1*g2;

    
    if((x22*xt[0]+y22*xt[1]+x*xt[0]+y*xt[1]+xy*xt[0]*xt[1]+s)>0){
        
        return x22*xt[0]*r+y22*xt[1]*r+x*xt[0]*r+y*xt[1]*r+xy*xt[0]*r*xt[1]*r;
    }
    else
        return x22*xt[0]*r/2+y22*xt[1]*r/2+x*xt[0]*r/2+y*xt[1]*r/2+xy*xt[0]*r/2*xt[1]*r/2;
}
    
    
    
    
    
    //-----------------------------------------------------------
    
    
    

int main(){
    double x1, x2, sv, g1, g2, svg;
    printf("Введите коэффициенты перед f(x)= x1^2 x2^2 C\n");
    scanf("%lf%lf%lf", &x1, &x2, &sv);
    
    printf("Введите коэффициенты перед g(x)= x1 x2 C\n");
    scanf("%lf%lf%lf", &g1, &g2, &svg);
    
    printf("Введите начальную точку x1 x2, r, c, e\n");
    double xt[2],r,c,e;
    scanf("%lf %lf %lf %lf %lf",&xt[0],&xt[1],&r,&c,&e);
    double* xtt;
    double fx;
    
    int k=0;
    m1:
    printf("Найдем точку безусловного минимума методом спуска с постоянным шагом\n");
    printf("-----------------------------------------------------------------------\n");
    xtt=Func(x1,x2,sv,g1,g2,svg,xt,r);
    xt[0]=xtt[0];
    xt[1]=xtt[1];
    printf("-----------------------------------------------------------------------\n");
    printf("точка безусловного минимума F = [%lf;%lf]\n",xt[0],xt[1]);
    fx=Pfunc(x1,x2,sv,g1,g2,svg,xt,r);
    printf("f(x*)= %lf  ---|||||||||-\n",fx);
    printf("Проверим P(x(r),r))<=e\n");
    if(fx<=e){
        
    printf("Истина, поиск завершен \n");
        
    }
    else
    {
        printf("Ложь, k++ \n");
        
        r=c*r;
        k=k+1;
        goto m1;
        
        
    }
   
    
    printf("Результат = x*[%lf;%lf]        f(x*)= %lf  за k = %d\n",xt[0],xt[1],g1*xt[0]+g2*xt[1]+svg,k );
    
    
    
    
    
    
    
    
    
    
    
    
}
