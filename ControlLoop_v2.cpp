/*
    Authors: Nardini Claizoni (nardinicq@gmail.com) Gabriel Benevides(gdebenevides@gmail.com)
    2018-1-15
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/



#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Util.h"
#include <wiringPi.h>
//#include <opencv2/opencv.hpp>
#include <ctime>
#include <sstream>
#include <opencv2/opencv.hpp>
#include <time.h>
#include <raspicam/raspicam_cv.h>

/*
g++ camera2.cpp -o  teste_camera2 -I/usr/local/include/ -L/opt/vc/lib -lraspicam -lraspicam_cv -lmmal -lmmal_core -lmmal_util -lopencv_core -lopencv_highgui
* g++ camera2.cpp -0 testecamera2 `pkg-config --libs opencv`
* g++ teste_controlLoop_v2.cpp Util.cpp -O3 -o teste_controlLoop  -I/usr/local/include/ -L/opt/vc/lib -lwiringPi -lraspicam -lraspicam_cv -lmmal -lmmal_core -lmmal_util `pkg-config --libs opencv`

* g++ teste_controlLoop_v2.cpp Util.cpp -o teste_controlLoop   -lwiringPi -lraspicam -lraspicam_cv `pkg-config --libs opencv`
*/

/*
Motor direito-> pinos 24,25
direcao - 24->H, 25->L
direcao + 24->L, 25->H
parado  24->L, 25->L

Motor esquerdo-> pinos 23,27
direcao - 23->H, 27->L
direcao + 23->L, 27->H
parado  23->L, 27->L

Odometro direito  -> 
Odometro esquerdo ->
*/
/*
void DelayMicrosecondsNoSleep (int delay_us)
{
	long int start_time;
	long int time_difference;
	struct timespec gettime_now;

	clock_gettime(CLOCK_REALTIME, &gettime_now);
	start_time = gettime_now.tv_nsec;		//Get nS value
	while (1)
	{
		clock_gettime(CLOCK_REALTIME, &gettime_now);
		time_difference = gettime_now.tv_nsec - start_time;
		if (time_difference < 0)
			time_difference += 1000000000;				//(Rolls over every 1 second)
		if (time_difference > (delay_us * 1000))		//Delay for # nS
			break;
	}
}
*/

void  softPWM(const int duty,int PIN1, int PIN2)
{
	const short unsigned int PWMslices=32;
	static int count;
	
	++count;
	count=count%PWMslices;

	if((count<duty)) {
		
		digitalWrite(PIN1,LOW);
		digitalWrite(PIN2,HIGH);//high
		
	}
	else
	{
		digitalWrite(PIN1,LOW);
		digitalWrite(PIN2,LOW);
	}


}
using namespace std; 

float getCenterLine(cv::Mat &image,cv::Mat &dst, const int y0=40,const int y1=70) {//const int y0=20,const int y1=50
			
			using namespace cv;
			
			short unsigned int counterX=0; 
			const short unsigned int Xmod=2;
			
			short unsigned int counterY=0;
			const short unsigned int Ymod=2;
			
			
			cvtColor(image,dst,CV_BGR2HSV);

			float cmx=0,nx=0,cmy=0,ny=0,cmx1=0,nx1=0,cx,cx1;

			for(int y=y0;y<y1;++y)
				for(int x=0;x<320;++x) {
					if(((counterX++)%Xmod==0)&&((counterY++)%Ymod==0))
					{
						cv::Vec3b &pxRGB=image.at<cv::Vec3b>(y,x);
						cv::Vec3b &pxHSV=dst.at<cv::Vec3b>(y,x);
						if((pxHSV.val[0]<60)&&(pxHSV.val[0]>10)&&(pxHSV.val[1]>40)&&(pxHSV.val[1]<200)) //(pxHSV.val[0]<60)&&(pxHSV.val[0]>10)&&(pxHSV.val[1]>145)&&(pxHSV.val[1]<255)
						{
							pxHSV.val[0]=255;
							pxHSV.val[1]=255;
							pxHSV.val[2]=255;
							
							pxRGB.val[0]=255;
							pxRGB.val[1]=255;
							pxRGB.val[2]=255;
						}
						else
						{
							pxHSV.val[0]=0;
							pxHSV.val[1]=0;
							pxHSV.val[2]=0;
							
							pxRGB.val[0]=0;
							pxRGB.val[1]=0;
							pxRGB.val[2]=0;
						}

					
						if ((pxHSV.val[0] == 255)&&(pxHSV.val[1] == 255)) 
						{
							cmx += x;
							nx++;
					
							
							//cmy += y;
							//ny++;               
						}
					}
            
				}
			
			if(cmx!=0&&nx!=0){
				cx=(cmx/nx);
				//circle( image, Point( cx, y0 ), 1, Scalar( 0, 0, 255 ), 2, 8 );
				return (159-cx)/159.0f;
			}
				//cout<<"cx= "<<cx<<endl;
			//cx=(cmx/nx);
			
			//if(cx==cx) {
			//	circle( image, Point( cx, y0 ), 1, Scalar( 0, 0, 255 ), 2, 8 );
		//	}
			
			//return (159-cx)/159.0f; //Pode retornar NaN -> testar utilizando if(x!=x)
}

using namespace std; 

//#define TESTMOTOR
//#define PRINTLOG

int main(void)
{
//criação de arquivo
	FILE *arq;
	arq = fopen("teste_controlLoop_v2.txt", "wt");
//setup
	wiringPiSetup();
	pinMode(24,OUTPUT);
	pinMode(25,OUTPUT);
	pinMode(23,OUTPUT);
	pinMode(27,OUTPUT);
	
	using namespace cv;

   // VideoCapture Camera;

    time_t timer_begin,timer_end;
	raspicam::RaspiCam_Cv Camera;
    cv::Mat image, dst;
    //int erosion_size = 2;  
    //Mat element = getStructuringElement(cv::MORPH_RECT, cv::Size(2 * erosion_size + 1, 2 * erosion_size + 1), cv::Point(erosion_size, erosion_size) );
    //Open camera
    
    Camera.set( CV_CAP_PROP_FORMAT, CV_8UC3 );//em cores
    Camera.set(CV_CAP_PROP_FRAME_WIDTH, 320);
    Camera.set(CV_CAP_PROP_FRAME_HEIGHT, 240);
    
    cout<<"Opening Camera..."<<endl;
    if (!Camera.open()) {cerr<<"Error opening the camera"<<endl;return -1;}
	cout<<"Camera aberta!" << endl;
    


 //Start capture
  
   cv::namedWindow("My0",CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
   //cv::namedWindow("My1",CV_WINDOW_AUTOSIZE);
 


//loop
	double dt0=0,dt0pwm=0, dt1=0;
	const double SAMPLING_TIME=0.08;
	const float dt=0.08;
	const double PWM_DT=312e-6;
	float e0=0,e1=0,e2=0,eTemp=0;
	float m=0,m1=0,m2=0;
	const float ku=2.0,pu=1.28,ti=pu/1.2,td=0;
	const float kp=0.5*ku,ki=0,kd=0;//kp=1.6,ki=0*dt,kd=0/dt;
	float k1=10,k2=1,k3=10,k4=1;
	float reto=18; //16
	float iae=0;
	const int N=375;
	short unsigned int coutCounter=0;
	float y[N];
	float x[N];
	int j=0;
	
	for(int i=0;i<N;++i) { y[i]=0; x[i]=0;};
	//Camera.grab();
	//Camera.retrieve ( image);
	//cv::imshow("My0", image); //show the frame in "MyVideo" window
	
	//delay(0.5);	
	for(;;) {

		
		dt1=getCurrentRealTimer();
		//cout<<"\ndiff=  "<<dt1-dt0<<std::flush;
		if((dt1-dt0)>=SAMPLING_TIME) 
		{	
			
			//aquisição da variável de controle
			Camera.grab();
			Camera.retrieve ( image);
			//erode(image,dst,element);
        
			//if(nx!=0 && ny!=0){
			//	circle( image, Point( cmx/nx, cmy/ny ), 1, Scalar( 0, 0, 255 ), 2, 8 );
		//	}
			//cv::imshow("My1", dst);
			
			
#ifndef TESTMOTOR
//calculo do erro
			eTemp=getCenterLine(image,dst);
			if(eTemp==eTemp) 
			{
				e0=e1;
				e1=e2;
				e2=eTemp;
				iae+=e2*e2;
			}
			else
			{
				//È um NaN
				
			}
			//if(e2<0) e2=-e2;
//CONTROLE
			m = m + kp*(e2-e1)+ ki*e2+(kd*(e2-2*e1 + e0));
			if(m>1) 
				m=1;
			else
			if(m<-1)
				m=-1;
			
			m1 = -k1*m+k2*reto;
			m2 = k3*m+k4*reto;
						
#else
int c=cv::waitKey(1)%255;

switch(c) {
	case 129: ++m1;break;//q, m1 motor da direita
	case 138:--m1;break;//z
	case 135:++m2;break;//w, m2 motor da esquerda
	case 136:--m2;break;//x
	case 43:
		softPWM((int) 0,24,25);  //possivel erro
		softPWM((int) 0,23,27);  //possivel erro
		return 0;
		break;
}
//printf("* %d",c);

#endif

//ajuste de escala do sinal de controle m;
			if(m1>32) m1=32; else { if(m1<1) m1=1;};
			if(m2>32) m2=32; else { if(m2<1) m2=1;};


//atualização para tempo de amostragem		
			dt0=dt1;
			
			//y[i]=
			//x[i]=
			
			
//imshow("My0", image); //show the frame in "MyVideo" window
//waitKey(1);

			
#ifdef PRINTLOG		
			if((coutCounter++)%10==0)
				cout<<" e2="<<e2<<" m="<<m<<" m1="<<m1<<" m2="<<m2<<endl;	
#endif

//salvar m e o tempo
			fprintf(arq,"%d %f\n",j,e2);
			//y[j]=m;
			//m[j]=j*SAMPLING_TIME;
			j++;
			if(j==250)// 250 tempo de trajeto com reto=16 
			{
				softPWM(0,24,25);  
				softPWM(0,23,27);
				fclose(arq);
				
				cout << "iae= " <<iae<<endl;
				cout << "w=" << image.cols << "h=" << image.rows << endl;
				imwrite("lastImage.jpg",image);
				
				return 0;
			}
			
}	 
	

	if((dt1-dt0pwm)>=PWM_DT) 
	{
		
		//softPWM(0,24,25);  
		//softPWM(0,23,27);
		softPWM((int) m1,24,25);  //possivel erro
		softPWM((int) m2,23,27);  //possivel erro
	
		dt0pwm=dt1;
	}
	
}
	

return 0;
}

