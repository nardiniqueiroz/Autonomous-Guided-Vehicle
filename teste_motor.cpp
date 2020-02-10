#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Util.h"
#include <wiringPi.h>
#include <opencv2/opencv.hpp>
#include <ctime>
#include <sstream>
#include <opencv2/opencv.hpp>
#include <time.h>

/*
#include <stdio.h>
#include <iostream> 
#include <wiringPi.h>
#include <sstream>

//g++ teste_motor.cpp Util.cpp -o teste_motor   -lwiringPi
//g++ teste_motor.cpp Util.cpp -o teste_motor   -lwiringPi `pkg-config --libs opencv`
Motor direito-> pinos 24,25
direcao - 24->H, 25->L
direcao + 24->L, 25->H
parado  24->L, 25->L

Motor esquerdo-> pinos 23,27
direcao - 23->H, 27->L
direcao + 23->L, 27->H
parado  23->L, 27->L

*/


using namespace std; 

int main(void) {
	
	int d1=0,d2=0,count=0,deltaT=0,k=0, amostra[2200], tempo[2200];
	wiringPiSetup();
	pinMode(24,OUTPUT);
	pinMode(25,OUTPUT);
	pinMode(23,OUTPUT);
	pinMode(27,OUTPUT);

	for(int i=0;i<10000;i++) 
	{
		
		digitalWrite(24,LOW); 
		digitalWrite(25,LOW); //HIGH
		digitalWrite(23,LOW); 
		digitalWrite(27,HIGH); 
		

		//cout<<digitalRead(13)<<"  "<<std::flush;
		d1=d2;
		d2=digitalRead(13);
		
		if(d1!=d2) 
		{
			count+=1;
			amostra[k]=count;
			tempo[k]=deltaT;
			k++;
		}
		//if(!(deltaT%5))
	//	{
			//cout<<count<<"  "<<std::flush;
			//count=0;
		//}
		
		deltaT+=1;
	
		delay(0.001);
		//digitalWrite(23,LOW); 
		//digitalWrite(27,LOW); 
	//	digitalWrite(24,LOW);
		//digitalWrite(25,LOW);
		//delay(3.125);	

	}
	
		digitalWrite(24,LOW); 
		digitalWrite(25,LOW); 
		digitalWrite(23,LOW); 
		digitalWrite(27,LOW); 
		
		for(int j=0;j<300;j++)
		cout<<tempo[j]<<" "<<amostra[j]<<endl;
		
		
		
}



/*
g++ camera2.cpp -o  teste_camera2 -I/usr/local/include/ -L/opt/vc/lib -lraspicam -lraspicam_cv -lmmal -lmmal_core -lmmal_util -lopencv_core -lopencv_highgui
* g++ camera2.cpp -0 testecamera2 `pkg-config --libs opencv`
* g++ teste_controlLoop.cpp Util.cpp -o teste_controlLoop   -lwiringPi `pkg-config --libs opencv`
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


void  softPWM(const int duty,int PIN1, int PIN2)
{
	static int count;
	
	++count;
	count=count%64;

	if((count<duty)) {
		
		digitalWrite(PIN1,LOW);
		digitalWrite(PIN2,HIGH);
		
	}
	else
	{
		digitalWrite(PIN1,LOW);
		digitalWrite(PIN2,LOW);
	}


}


using namespace std; 

int main(void)
{
//setup
	wiringPiSetup();
	pinMode(24,OUTPUT);
	pinMode(25,OUTPUT);
	pinMode(23,OUTPUT);
	pinMode(27,OUTPUT);
	
	using namespace cv;


//loop
	double dt0=0,dt0pwm=0, dt1=0;
	const double SAMPLING_TIME=0.040000;
	const double PWM_DT=0.000080, PWM_DT2=10;
	float e0=0,e1=0,e2=0;
	float m=0,m1=20,m2=63;
	float kp=1,ki=0,kd=0;
	float k1=1,k2=1,k3=1,k4=1;
	float reto=0; 
	
	
//	long int start_time;
//	long int time_difference;
//	struct timespec gettime_now;
	
//	clock_gettime(CLOCK_REALTIME, &gettime_now);
//	start_time = gettime_now.tv_nsec;	
		dt0=getCurrentRealTimer()/1000000;
		
	for(;;) {
		//clock_gettime(CLOCK_REALTIME, &gettime_now);
		//time_difference = gettime_now.tv_nsec - start_time;
		
		dt1=getCurrentRealTimer()/1000000;

	
	if((dt1-dt0pwm)>=PWM_DT) 
	//if(time_difference>(PWM_DT))
	{
		//cout<<(dt1-dt0pwm)<<"  "<<std::flush;
		
		//softPWM((int) m1,24,25);  //possivel erro
		//softPWM((int) m2,23,27);  //possivel erro
		digitalWrite(24,LOW);
		digitalWrite(25,LOW);
		digitalWrite(23,LOW);
		digitalWrite(27,LOW);
		
		dt0pwm=dt1;
	}
	
	if((dt1-dt0)>=10) 
	{
		digitalWrite(24,LOW);
		digitalWrite(25,LOW);
		digitalWrite(23,LOW);
		digitalWrite(27,LOW);
		
		}
	
	
}


return 0;
}

*/
