#include <time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#define X_SIZE 1920 
#define STRIDE 5760 
#define Y_SIZE  1080 
#define Y_SIZ2  400 
#define MAX  44100*60*10
#define FRAMES 6597

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;
XImage *x_image;
unsigned char *x_buffer;
unsigned char *y_buffer;
unsigned char *z_buffer;
int dot_line;
long astart,aend;
int bloom;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


void motif (int,int,int); //motif
void earth (int,int,int);
void grass (int,int,int);
void sun   (int,int,int);
void trees (int,int,int);
void river (int,int,int);
void flower(int,int,int);
void disp (int);
void do_dots (int,long,long);
void add_image (int,int,int,int,int);
void plot (int,int,int,int,int);
void bplot (int,int,int,int,int);
void clot (int,int,int,int,int,int,int);
void square (int,int,int,int,int,int,int);
void dlot (int,int,int,int,int,int,int,int,int);
void line (int *, int, int, int, int, int, int);
void lline (int, int, int, int, int, int, int, float,float,float,float);
void bline (int *, int, int, int, int, int, int);
void tree (int, int, int, int ,int, int);
void sree (int, int, int, int ,int, int);
void cloud (int,int,int);
void blurr(int);
void atten(int);
void my_load_image( unsigned char *,char *, int);
void show_image( unsigned char *);
void light(int);
void lightf(int,int,int,int,int);
void rain (int, int, int);
void rainbow (int, int, int);
void split (int,int,int);
void mask (int);
void merge(unsigned char *, unsigned char *, int);
void trip(unsigned char *, unsigned char *, unsigned char *);
void shrink(int,int);
void bass(int,int);
void verb(int);
void arc ( int, int , int , int , int , int , int , int , int );
void num (int,int,int);


/* Jpegs */
int read_JPEG_file (char *, unsigned char *, int *);
int jayit(unsigned char *,int, int, char *);
unsigned char *image2,*image3,*image4,*dots,*spiral,*logo;
int *image5;
short *waveform,*left,*right;
int dpx[30],dpy[30],sz[30];
float frames;
int jay,*seq;
int dims[2];
int pims[2];
int *clud;
int *clue;
int cord[4];

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}


int main(int argc,char *argv[])
{
	int *fhead;
	char top,*in;
	FILE *input,*dat;
        int x_size,y_size,x_stride,x,y,frame_rate,samples_per_frame,samples;
	int loop,chan,sample_rate,size,time,e0p,e1p,e2p,e3p,e4p,e5p,e6p,e7p,e8p,e9p,start;
	unsigned char *images[10];
        clud=(int *)malloc(sizeof (int)*X_SIZE*Y_SIZE);
        clue=(int *)malloc(sizeof (int)*X_SIZE*Y_SIZE);
        image2=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // earth buffer
        image4=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // bplot buffer
        image5=(int *)malloc(sizeof (int)*3*X_SIZE*Y_SIZE); // layer buffer
        dots=(unsigned char *)malloc(sizeof (char)*3*4096*4096);
        spiral=(unsigned char *)malloc(sizeof (char)*3*4096*4096);
        logo=(unsigned char *)malloc(sizeof (char)*3*4096*4096);
        seq=(int *)malloc(sizeof (int)*X_SIZE*Y_SIZE*4);
        in=(char *)malloc(sizeof (char)*100);
	x_image = (XImage *)malloc(sizeof(XImage));
	int dp[2][8];
	int yp[2][10];

        for (loop=0;loop<10;loop++){images[loop]=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE);}

	for (x=0;x<X_SIZE*Y_SIZE*4;x++){ seq[x]=rand();}


	for (x=0;x<30;x++){
		dpx[x]=rand()%X_SIZE; dpy[x]=rand()%Y_SIZE; sz[x]=10+rand()%20;
	}

	frame_rate=30;

	jay=0;
	int year;
	int month;
	int day;
	int now;
	int xao;
	now=2021;
	if (argv[1]){ day=atoi(argv[1]); }else{exit(1);}
	if (argv[2]){ month=atoi(argv[2]); }else{exit(1);}
	if (argv[3]){ year=atoi(argv[3]); }else{exit(1);}
	if (argv[4]){ xao=atoi(argv[4]); }else{exit(1);}
	if (argv[5]){ jay=1; }


	init_x();

	int ao;
	ao=xao;

	xao=150;


	printf("here\n");

	float fram;

	frames=1800;

	for (fram=0;fram<frames;fram+=1)
	{

	for (loop=0;loop<X_SIZE*3*Y_SIZE;loop++){image3[loop]=image2[loop];image2[loop]=255;}
	float loap;
	float mloap;

	mloap=10000;

	for (loap=0;loap<mloap;loap++)
	{	
		int xb,xe,yb,ye;

		xb=(X_SIZE/4)+(loap*(X_SIZE/4)*sin(loap*(4+(sin(4*M_PI*fram/frames)))*M_PI/mloap)/mloap);
		yb=(Y_SIZE/2)+(loap*(Y_SIZE/2)*cos(loap*4*M_PI/mloap)/mloap);
		xe=((3*X_SIZE)/4)+(loap*(X_SIZE/4)*sin(loap*4*M_PI/mloap)/mloap);
		ye=(Y_SIZE/2)+(loap*(Y_SIZE/2)*cos(loap*(4+(sin(2*M_PI*fram/frames)))*M_PI/mloap)/mloap);

//void lline (int xs, int ys, int xe, int ye, int r, int g, int b, float ag,int dup,float t,float bb)
		lline (xb, yb, xe,ye, 
				127+(127*(sin(loap*(2+(2*fram/frames))*M_PI/mloap))), 
				127+(127*(sin(loap*(2+(6*fram/frames))*M_PI/mloap))), 
				127+(127*(sin(loap*(2+(10*fram/(frames)))*M_PI/mloap))), 
				fram/(frames*5),
				5-(5*(cos(4*M_PI*fram/frames))),
				4,
				0.2*sin(6*M_PI*fram/frames));
		 //if ((int)loap%500==0){disp(loap/10);}
	}
	disp(fram);
	}
	char ds;
	scanf("%c",&ds);

	close_x();
}	

void num (int x, int y, int no)
{
	int n,l,off;
	if (x<X_SIZE/2) { off=-10;}else{off=10;}
	for (n=0;n<no;n++)
	{
		for (l=-10;l<11;l++)
		{
			plot( x+((n+1)*off),y+l,255,255,255);
		}
	}
}


void trip (unsigned char *a, unsigned char *b, unsigned char *m)
{
	int x;
	for (x=0;x<X_SIZE*Y_SIZE*3;x++)
	{
		image2[x]=(a[x]*m[x]+b[x]*(255-m[x]))/255;
	}
}


void verb(int clean)
{
	int x;
	for (x=0;x<X_SIZE*Y_SIZE*3;x++)
	{
		image2[x]=(image2[x]+(clean*image4[x]))/(clean+1);
		image4[x]=image2[x];
	}
}


void bass(int loop, int b)
{
	long along;
	int plus,cross,max,one,two,three,four,xp,yp,p,l,or[2],ang,perp,xs[4],ys[4],ye[4],xe[4];
	cross=0;
	one=34; two=50; three=73; four=1000;
	one=34; two=59; three=120; four=1000;
	max=0;
	p=0;
	l=(loop%3)*300;
	//if (loop%90==0){for (along=0;along<X_SIZE*Y_SIZE*3;along++){ image2[along]=255;}}
	for (along=0;along<X_SIZE*Y_SIZE*3;along++){ image2[along]=255;}
	if (left[astart]>0){ plus=1;}else{plus=0;}
	for (along=astart;along<aend;along++)
	{
		if (plus && left[along]<0){ plus=0; cross++;}
		else if ( !plus && left[along]>0  ){ plus=1; cross++;}
		if (left[along]>max){ max=left[along]; }
	}
	printf("Loop %d Cross %d\n",loop,cross);
	if (cross<=one){ printf("one\n");p=3;}
	//for (along=0;along<X_SIZE*Y_SIZE*3;along++){ image2[along]=255;}}
	else if (cross<=two) { printf("two\n");p=2;}
	else if (cross<=three) { printf("three\n");p=1;}
	else if (cross<=four) { printf("four\n");p=0;}
	/*for (along=astart;along<aend;along++)
	{
		int x,d,m;
		x=((along-astart)*X_SIZE)/((aend-astart));
		for (m=0;m<4;m++)
		{
			if (m==p){ 
			d=l+100+(m*80)+(16*left[along]/max); } else { d=l+100+(m*80);}
			plot (x,d,255,0,0);
			plot (x,3+d,0,255,0);
			plot (x,5+d,0,0,255);
			plot (x,d-3,0,0,255);
			plot (x,d-5,0,255,0);
		}
	}*/	
	or[0]=-200;or[1]=3*Y_SIZE/4;
	ang=65;
	perp=ang+90;
	//neck
	line (or,208,208,208,2500,ang,300);

	//frets
	int fret,len;
	for (fret=0;fret<15;fret++)
	{
		len=100+(180*fret)-(130*sqrt(fret));
		or[0]=-200;or[1]=3*Y_SIZE/4;
		line (or,255,255,255,2,perp,1);
		line (or,208,208,208,len,ang,1);
		line (or,148,148,138,296,perp,8);
	}

	//strings
	or[0]=-200;or[1]=3*Y_SIZE/4;
	line (or,255,255,255,30,perp,1);
	xs[0]=or[0];ys[0]=or[1];
	line (or,0,0,0,2500,ang,1);
	xe[0]=or[0];ye[0]=or[1];

	or[0]=-200;or[1]=3*Y_SIZE/4;
	line (or,255,255,255,110,perp,1);
	xs[1]=or[0];ys[1]=or[1];
	line (or,0,0,0,2500,ang,1);
	xe[1]=or[0];ye[1]=or[1];

	or[0]=-200;or[1]=3*Y_SIZE/4;
	line (or,255,255,255,190,perp,1);
	xs[2]=or[0];ys[2]=or[1];
	line (or,0,0,0,2500,ang,1);
	xe[2]=or[0];ye[2]=or[1];

	or[0]=-200;or[1]=3*Y_SIZE/4;
	line (or,255,255,255,270,perp,1);
	xs[3]=or[0];ys[3]=or[1];
	line (or,0,0,0,2500,ang,1);
	xe[3]=or[0];ye[3]=or[1];

	for (along=astart;along<aend;along++)
	{
		int x,d,m;
		x=((along-astart)*X_SIZE)/((aend-astart));
		for (m=0;m<4;m++)
		{
			if (m==p){ 
			d=(27*left[along]/max); } else { d=4;}

			float dx,dy;
			dx=(float)(along-astart)*(float)(xe[m]-xs[m])/(float)(aend-astart);
			dy=(float)(along-astart)*(float)(ye[m]-ys[m])/(float)(aend-astart);

			or[0]=xs[m]+dx;or[1]=ys[m]+dy;

			int v;
			v=m*10;
			if (d>0)
			{
			line (or,180-m,100+m,50+m,d,perp,2);
			}else{
			line (or,160+m,90-m,110+m,-d,perp+180,2);
			}
		}
	}
	verb(b);
}	


void shrink (int xp, int yp)
{
	int xo,yo,x_pad,y_pad; 

	if(xp==0){xp=1;}
	if(yp==0){yp=1;}

	x_pad=(X_SIZE-((X_SIZE*100)/xp))/2;
	y_pad=(Y_SIZE-((Y_SIZE*100)/yp))/2;

	for (xo=0;xo<X_SIZE*Y_SIZE*3;xo++){ image4[xo]=255;}



	if (xp>0)
	{
	for (yo=0;yo<Y_SIZE;yo++)
	{
		int p,q,yn;
		p=yo*STRIDE; //4
		yn=((yo*100)/yp)+y_pad;
		if (yn>0 && yn<Y_SIZE)
		{
			q=yn*STRIDE; //2
			for (xo=0;xo<X_SIZE;xo++)
			{
				int xn;
				xn=((xo*100)/xp)+x_pad;
				if (xn>0 && xn<X_SIZE)
				{
				image4[p+(xo*3)]=image2[q+(xn*3)];
				image4[p+(xo*3)+1]=image2[q+(xn*3)+1];
				image4[p+(xo*3)+2]=image2[q+(xn*3)+2];
				}
			}
		}
	}
	}else{
	xp=-xp;
	x_pad=(X_SIZE-((X_SIZE*100)/xp))/2;
	for (yo=0;yo<Y_SIZE;yo++)
	{
		int p,q,yn;
		p=yo*STRIDE; //4
		yn=((yo*100)/yp)+y_pad;
		if (yn>0 && yn<Y_SIZE)
		{
			q=yn*STRIDE; //2
			for (xo=0;xo<X_SIZE;xo++)
			{
				int xn;
				xn=((xo*100)/xp)+x_pad;
				if (xn>0 && xn<X_SIZE)
				{
				int xor;
				xor=(X_SIZE-1-xo)*3;
				image4[p+(xor)]=image2[q+(xn*3)];
				image4[p+(xor)+1]=image2[q+(xn*3)+1];
				image4[p+(xor)+2]=image2[q+(xn*3)+2];
				}
			}
		}
	}
	}


	memcpy(image2,image4,X_SIZE*Y_SIZE*3);
}


void merge(unsigned char *out, unsigned char *in, int amount)
{
	if (amount<0){ amount=0;}
	if (amount>255){ amount=255;}
	int prime,x;
	prime=255-amount;
	for (x=0;x<X_SIZE*Y_SIZE*3;x++)
	{
		out[x]=((amount*out[x])+(prime*in[x]))/256;
	}
}


void rainbow (int tot,int val, int stage)
{
	int xo,yo,r,g,b,or[2],width,radius,xx,yy,diff,wp,valp;
	float phi,ghi,bhi,feta,x,dfeta,rad;
	width=200;
	for (xo=0;xo<X_SIZE*Y_SIZE*3;xo++){ image5[xo]=1024;}
	xo=X_SIZE/2;
	yo=X_SIZE/2;
        diff=aend-astart;
	if (val>255){val=255;}
	valp=255-val;


	for (x=0;x<width;x+=0.7)
	{
		dfeta=1/(float)X_SIZE;
		phi=M_PI*(float)x/width;
		ghi=(M_PI/2)+(M_PI*(float)x/width);
		bhi=(M_PI/4)+(M_PI*(float)x/width);
		r=((255-((127*x)/width))*cos(phi)*cos(phi));
		g=(255*cos(ghi)*cos(ghi));
		b=(255*cos(bhi)*cos(bhi));
		for (feta=0;feta<M_PI;feta+=dfeta)
		{
                        wp=astart+((feta*diff)/M_PI);
			rad=((float)(X_SIZE/2)-x-(float)30)+(stage*(float)left[wp]/64);
			xx=(rad)*(cos(feta));
			yy=(rad)*(sin(feta));
      			bplot (xo-xx,yo-yy,r,g,b);
		}
	}
	for (xo=0;xo<(X_SIZE*Y_SIZE*3);xo++){ if(image5[xo]!=1024){
		image2[xo]=((valp*image2[xo])+(val*image5[xo]))/256;}}
}


void rain (int tot,int val, int dir)
{
	int depth,or[2];
	float xmove,ymove;
	xmove=-(float)val*sin((float)2*M_PI*dir/360);	
	ymove=(float)val*cos((float)2*M_PI*dir/360);	
	for (depth=0;depth<val*200;depth++)
	{
		float speed;
		speed=0.5-((float)(seq[depth+2]%1000)/4000);
		or[0]=seq[depth]%X_SIZE+(xmove*(float)tot*speed);
		or[1]=seq[depth+1]%Y_SIZE+(ymove*(float)tot*speed);
		or[1]=or[1]%Y_SIZE;
		or[0]=or[0]%X_SIZE;
		if (or[0]<0){or[0]=or[0]+X_SIZE;}
		if (or[1]<0){or[1]=or[1]+Y_SIZE;}
      		line (or,255,255,255,val/4,dir,1+(val/50));
	}
}	


void split (int ang, int pha, int depth)
{
	float xrmove,xgmove,xbmove;
	float yrmove,ygmove,ybmove;
	int x,y;

	xrmove=(float)depth*(sin((2*M_PI*(float)ang)/360));
	yrmove=(float)depth*(cos((2*M_PI*(float)ang)/360));

	xgmove=(float)depth*(sin((2*M_PI*(float)(ang+pha))/360));
	ygmove=(float)depth*(cos((2*M_PI*(float)(ang+pha))/360));

	xbmove=(float)depth*(sin((2*M_PI*(float)(ang+(2*pha)))/360));
	ybmove=(float)depth*(cos((2*M_PI*(float)(ang+(2*pha)))/360));

	for (y=0;y<Y_SIZE;y++)
	{
		int p,r,ym;
		p=y*STRIDE;
		ym=y+yrmove;
		if (ym<0){ym=Y_SIZE+ym;} if (ym>=Y_SIZE){ym=Y_SIZE-ym;}
		r=ym*STRIDE;
		for (x=0;x<STRIDE;x+=3)
		{
			int xm;
			xm=x+(3*(int)xrmove);
			if (xm<0){xm=STRIDE+xm;} if (xm>=STRIDE){xm=xm-STRIDE;}
			image4[p+x]=image2[r+xm];
		}
	}
	for (y=0;y<Y_SIZE;y++)
	{
		int p,r,ym;
		p=y*STRIDE;
		ym=y+ygmove;
		if (ym<0){ym=0;} if (ym>=Y_SIZE){ym=Y_SIZE-1;}
		r=ym*STRIDE;
		for (x=1;x<STRIDE;x+=3)
		{
			int xm;
			xm=x+(3*(int)xgmove);
			if (xm<0){xm=0;} if (xm>=STRIDE){xm=STRIDE-2;}
			image4[p+x]=image2[r+xm];
		}
	}
	for (y=0;y<Y_SIZE;y++)
	{
		int p,r,ym;
		p=y*STRIDE;
		ym=y+ybmove;
		if (ym<0){ym=0;} if (ym>=Y_SIZE){ym=Y_SIZE-1;}
		r=ym*STRIDE;
		for (x=2;x<STRIDE;x+=3)
		{
			int xm;
			xm=x+(3*(int)xbmove);
			if (xm<0){xm=0;} if (xm>=STRIDE){xm=STRIDE-1;}
			image4[p+x]=image2[r+xm];
		}
	}
	memcpy(image2,image4,STRIDE*Y_SIZE);
}


//lightening
void light (int count)
{
	int x_start,y_start,angle,loop,len;

	x_start=seq[count]%X_SIZE;
	y_start=seq[count+1]%Y_SIZE/8;
	len=200+seq[count+2]%300;
	angle=180;
	if (x_start<X_SIZE/2){
		angle-=seq[count]%90;}
	else{
		angle+=seq[count]%90;}

	lightf(x_start,y_start,angle+rand()%5,1,len);
}	

void lightf (int x, int y, int angle, int point, int length)
{
	if (point>10){ return;}
	int or[2],xx,yy;
	or[0]=x;
	or[1]=y;
      	line (or,255-rand()%10,255-rand()%20,255,length/point,angle+rand()%20,11-point);
	point++;
	xx=or[0];
	yy=or[1];
	lightf(xx,yy,angle,point,length);
	if(rand()%100>70){lightf(xx,yy,angle+40+rand()%15,point+1,length);}
	if(rand()%100>70){lightf(xx,yy,angle-40+rand()%15,point+1,length);}
}

//motif
void motif(int loop, int blur, int stage) 
{
        printf ("radius %d %ld %ld\n",loop,astart,aend);
        int x,y;
        int pot,radius;
        long wf,diff;
        float perimeter,p,rad;
        int xh,yh,bl,act;
	act=0;

        xh=dims[0]/2; yh=dims[1]/2;
        diff=aend-astart;
        radius=10+(loop*loop/5);
	bl=30-(loop/10);
	if (stage>0 ){radius=490;bl=0;}
	if (stage>1 ){act=1;}

        perimeter=2*M_PI*radius;

        // this cleans the frame with tartan
        for (y=0;y<Y_SIZE;y++)
        {
                for (x=0;x<X_SIZE;x++)
                {
                        plot(x,y,200+(55*(sin(((float)loop/1000)+((float)(x+y)/100)))),200+(55*(sin(((float)loop/1000)+((float)(x-y))/80))),255);
                }
        }

        for (rad=1;rad<radius;rad+=0.5)
        {
                for (p=0;p<perimeter;p+=0.5)
                {
                        float ghi,phi;
                        long wp;
                        int yp,t,gx,gy,gp;
                        phi=2*M_PI*p/perimeter;
                        ghi=phi+(stage*loop*M_PI/300);
                        wp=astart+((p*diff)/perimeter);
                        yp=rad+((act*rad*(float)left[wp])/32768);
                        x=(X_SIZE/2)+(float)(yp)*(sin(phi));
                        y=(Y_SIZ2)+(float)(yp)*(cos(phi));
                        gx=xh+(rad*xh*sin(ghi)/(radius));
                        gy=yh+(rad*yh*cos(ghi)/(radius));
                        gp=(gy*dims[0]*3)+gx*3;
                        plot (x,y,logo[gp],logo[gp+1],logo[gp+2]);
                }
        }
	blurr(bl);
}

//earth
void earth(int loo, int blur, int stage)
{
        int x,y;
        if (stage>0){
                for (x=0;x<3*X_SIZE*Y_SIZE;x++){ image2[x]=image3[x];}

                } else {

	if (loo<0){loo=0;}

        square (0,0,90,90,255,X_SIZE,Y_SIZ2); //sky
        square (0,Y_SIZ2,120,80,10,X_SIZE,Y_SIZ2); //earth
        int ygrain,xgrain;
        xgrain=1;
        ygrain=1;
        for (x=0;x<X_SIZE;x+=xgrain){
                int xm;
                for (y=0;y<Y_SIZE;y+=ygrain){
                        int g,h,gp,r,ym;
                        gp=x+(y*X_SIZE);
                        r=(seq[gp])%2048;
                        for (g=0;g<xgrain;g++)
                        {
                                for (h=0;h<ygrain;h++)
                                {
                                        int pix;
                                        for (pix=0;pix<3;pix++)
                                        {
                                                int p,val;
                                                p=((x+g)*3)+((y+h)*3*X_SIZE)+pix;
                                                val=image2[p]+(((1024-r)*loo)/8096);
                                                if (val<0){val=0;}
                                                if (val>255){val=255;}
                                                image2[p]=val;
                                        }
                                }
                        }
                        ym=(Y_SIZ2)-y;if(ym<0){ym=-ym;}
                        ygrain=1+(ym/20);
                }
                xm=(X_SIZE/2)-x;if(xm<0){xm=-xm;}
                xgrain=20-(xm/100);
        }
        blurr(blur);
        }
}

//grass
void grass (int tot, int blur, int stage)
{
        int x,y,off,e,d,len;
        printf ("%d %ld %ld\n",tot,astart,aend);

        off=0; len=tot;

        if (tot>150 ){ len=150;}
	if (stage>0 ){ len=150; off=1;}
	if (stage>1 ){ off=((20-(dot_line/100))); }

        int xgrain,ygrain,xm;
        xgrain=1;
        ygrain=1;
        int t;
        t=0;
        for (y=1;y<Y_SIZ2;y+=ygrain)
        {
                int ym;
                for ( x=xgrain/2;x<X_SIZE;x+=xgrain )
                {
                        int xo;
                        int or[2];
                        xo=seq[t++]%10;
                        if(off==0){xo+=rand()%2;}
                        or[0]=x+xo; or[1]=y+(Y_SIZ2);
                        line (or,0,64+((y*384)/Y_SIZE),0,2*y*len/Y_SIZE,0+off,1+((xgrain+ygrain)/20));
                        or[0]=x+xo; or[1]=y+(Y_SIZ2);
                        line (or,0,64+((y*127)/Y_SIZE),0,y*len/(Y_SIZE),12+off,1+((xgrain+ygrain)/24));
                        or[0]=x+xo; or[1]=y+(Y_SIZ2);
                        line (or,0,64+((y*127)/Y_SIZE),0,y*len/(Y_SIZE),-12+off,1+((xgrain+ygrain)/24));
                        xm=(X_SIZE/2)-x;if(xm<0){xm=-xm;}
                        xgrain=ygrain+15-((10*xm)/X_SIZE);
                }
                ygrain=1+(y/50);
        }

}

//sun
void sun (int tot, int blur, int stage )
{
        printf ("radius %d %ld %ld\n",tot,astart,aend);

        float radius,rad;
        int x,y,xo,yo,diff,hold;
        diff=aend-astart;
	hold=1;
	if (stage<2){hold=0;}

        radius=70;
	//xo=(X_SIZE/2)+(((X_SIZE/2)-radius)*(sin(((tot)+90)*2*M_PI/360)));
	xo=(X_SIZE/2)+(((X_SIZE/2))*(sin(((tot)+90)*2*M_PI/360)));
	yo=(Y_SIZ2)+(((Y_SIZ2)-radius)*(cos(((tot)+90)*2*M_PI/360)));

        for (rad=1;rad<radius;rad+=0.3)
        {
                float p,perimeter;
                perimeter=2*M_PI*rad;
                for (p=0;p<perimeter;p+=0.3)
                {
                        float ghi,phi;
                        long wp;
                        int yp,h;
                        phi=2*M_PI*p/perimeter;
                        wp=astart+((p*diff)/perimeter);
                        yp=rad+((hold*rad*(float)left[wp])/32768);
                        x=xo+(float)(yp)*(sin(phi));
                        y=yo+(float)(yp)*(cos(phi));
			h=(Y_SIZ2)-y;
			if(h>0)
			{
			int g;
			g=100+(h/2);if(g>255){g=255;}
                        plot (x,y,255,g,0);
			}
                }
        }
}

//trees 
void trees(int loop, int blur, int stage)
{
	int mul;
	mul=30;
	bloom=11-((loop)/25);
	if (bloom<2){bloom=2;}
        if (stage==1){bloom=2;}
        if (stage==2){mul=30+(dot_line/1276);bloom=2;}
       	tree(X_SIZE/2,Y_SIZE-500,130,0,20,mul);
       	tree((X_SIZE/2)-150,Y_SIZE-500,120,0,20,mul);

       	sree((X_SIZE/2)+200,Y_SIZE-500,150,0,30,mul);
       	sree((X_SIZE/2)-200,Y_SIZE-500,200,0,30,mul);

        sree(512,Y_SIZE-300,600,0,50,mul-10);
        sree(1400,Y_SIZE-300,600,0,50,mul-10);

        tree(X_SIZE/6,Y_SIZE-10,500,0,30,mul);
        tree(5*X_SIZE/6,Y_SIZE-50,600,0,30,mul);
}

//river
void river (int loop, int val, int stage)
{
	int x,y,r,g,b,pix,xo,valp;
	int or[2],off;
	r=32; g=32;b=128; 
	pix=1;
	if (val>255){val=255;}
	valp=255-val;
	for (xo=0;xo<X_SIZE*Y_SIZE*3;xo++){ image5[xo]=1024;}
	for (y=Y_SIZ2;y<Y_SIZE;y+=pix)
	{
		int pix;
		or[0]=(X_SIZE/2)+(200*sin(10*(float)y/Y_SIZE))+(seq[y+loop]%12);
		or[1]=y;

		pix=1+(y-(Y_SIZ2))/(Y_SIZE/64);
		//if (y%pix==0)
		if (1)
		{
		r=4+r-seq[y+loop]%10;
		b=4+b-seq[y+1+loop]%10;
		if (r>128) { r=128;} if (r<32) { r=32;}
		if (b>255) { b=255;} if (b<64) { b=64;}
		g=r;
		}
              	bline (or,r,g,b,((y)-(Y_SIZ2))/2,90+(seq[y+loop]%10),pix);
		if (1)
		{
		r=4+r-seq[y+loop]%10;
		b=4+b-seq[y+1+loop]%10;
		if (r>128) { r=128;} if (r<32) { r=32;}
		if (b>255) { b=255;} if (b<64) { b=64;}
		g=r;}
              	bline (or,r,g,b,((y)-(Y_SIZ2))/2,90-(seq[2+y+loop]%10),pix);
	}	
	for (xo=0;xo<(X_SIZE*Y_SIZE*3);xo++){ if(image5[xo]!=1024){
		image2[xo]=((valp*image2[xo])+(val*image5[xo]))/256;}}

}

/*void mask (int level)
{
	int x,y,r,g,b;
	int or[2];
	for (x=0;x<X_SIZE;x++)
	{
		for (y=0;y<Y_SIZE;y++)
		{
			if (seq[x+(y*X_SIZE)]%70000<level)
			{
				int ang;
				r=(seq[x]%5)*64;
				g=(seq[y]%5)*64;
				b=(seq[x+y]%5)*64;
				if (r==0 && g==0 && b==0){ r=255;g=255;b=255;}
				for (ang=0;ang<360;ang+=10)
				{
				or[0]=x; or[1]=y;
              			line (or,r,g,b,10+(seq[x+ang]%100),ang,2);
				}
              			line (or,b,r,g,150+(seq[x+or[1]]%300),seq[y]%360,5);
			}
		}
	}
}*/


void mask (int level)
{
	int x,y,r,g,b;
	int or[2],ang,os[2];
	r=rand()%256;
	g=rand()%256;
	b=rand()%256;
	x=rand()%X_SIZE;
	y=rand()%Y_SIZE;
	for (ang=0;ang<360;ang+=level)
	{
		or[0]=x; or[1]=y;
		int l;
		l=100+(rand()%300);
       		line (or,r,g,b,l,ang,2);
		int t;
		for (t=0;t<360;t+=level)
		{
			os[0]=or[0];
			os[1]=or[1];
			int ll;
			ll=20+(rand()%60);
       			line (os,b,r,g,ll,t,1);
		}
	}
       //line (or,b,r,g,100+(rand()%100),rand()%360,2);
}

//flowers
void flower(int loop, int blur, int stage)
{
	int l,off;
	off=0;
	if (stage>0){loop=200;}
	if (loop>200){loop=200;}
	if (stage>1 ){ off=((20-(dot_line/100))); }
	for (l=0;l<2*loop;l++)
	{
        	int or[2],xs,ys;
		int tlen,r,g,b;
                xs=seq[2*l]%(X_SIZE);
                ys=(Y_SIZ2)+seq[(2*l)+1]%(Y_SIZ2);
		or[0]=xs; or[1]=ys;
		tlen=((ys-(Y_SIZ2))*loop)/(Y_SIZ2);
                line (or,25,235,25,tlen,-off,tlen/20);
		r=64+seq[l]%192;
		g=64+seq[l+1]%192;
		b=64+seq[l+2]%192;
		dlot( or[0], or[1], r, g, b, tlen/8, 0, 360, 5+seq[l]%10);
		or[0]=xs; or[1]=ys;
                line (or,20,255,20,tlen/2,20-off,tlen/20);
		or[0]=xs; or[1]=ys;
                line (or,20,215,20,tlen/2,-20-off,tlen/20);
	}
}

void disp (int fram)
{
	int x,y;
	char input[100];

       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*STRIDE;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/jmage%04d.jpg",fram);
	if (jay){jayit(image2,X_SIZE, Y_SIZE, input);}
}


void atten(int depth)
{
	if (depth>255){depth=255;}
	if (depth<0){depth=0;}
	int x;
	for (x=0;x<X_SIZE*Y_SIZE*3;x++)
	{
		image2[x]=(image2[x]*depth)/255;
	}
}


void blurr (int depth)
{
	unsigned char *image1;
        image1=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // blur buffer
	int x,y,d,pix,ly,xe;
	ly=STRIDE*(Y_SIZE-1);
	xe=STRIDE-3;

	for (d=0;d<depth;d++)
	{

		// corners
		for (pix=0;pix<3;pix++){ 
			image1[pix]=(image2[pix]+image2[pix+3]+image2[pix+3+STRIDE]+image2[pix+STRIDE])/4;
			image1[xe+pix]=(image2[xe+pix]+image2[xe-3+pix]+image2[xe+pix-3+STRIDE]+image2[xe+pix+STRIDE])/4;
			image1[pix+ly]=(image2[pix+ly]+image2[ly+pix+3]+image2[ly+pix+3-STRIDE]+image2[ly+pix-STRIDE])/4;
			image1[xe+pix+ly]=(image2[xe+pix+ly]+image2[xe-3+pix+ly]+image2[xe+pix-3+STRIDE+ly]+image2[xe+pix+STRIDE+ly])/4;
		}
		// top/bottom
		for (x=3;x<STRIDE-4;x+=3)
		{
			for (pix=0;pix<3;pix++){ 
				image1[pix+x]=(image2[pix+x]+image2[pix+3+x]+image2[pix+x-3]+
						image2[pix+x+STRIDE]+image2[pix+3+x+STRIDE]+image2[pix+x-3+STRIDE])/6;
				image1[pix+x+ly]=(image2[pix+x+ly]+image2[pix+3+x+ly]+image2[pix+x-3+ly]+
						image2[pix+x+ly-STRIDE]+image2[pix+3+x+ly-STRIDE]+image2[pix+x-3+ly-STRIDE])/6;
			}
		}
		// left/right
		for (y=1;y<Y_SIZE-1;y++)
		{
			int yp,yr;
			yp=y*STRIDE;
			yr=yp+STRIDE-3;
			for (pix=0;pix<3;pix++){ 
				image1[yp+pix]=(image2[yp+pix]+image2[yp+pix-STRIDE]+image2[yp+pix+STRIDE]+
					image2[yp+pix+3]+image2[yp+pix-STRIDE+3]+image2[yp+pix+STRIDE+3])/6;
				image1[yr+pix]=(image2[yr+pix]+image2[yr+pix-STRIDE]+image2[yr+pix+STRIDE]+
					image2[yr+pix-3]+image2[yr+pix-STRIDE-3]+image2[yr+pix+STRIDE-3])/6;
			}
		}

		for (y=1;y<Y_SIZE-1;y++)
		{
			for (x=3;x<STRIDE-4;x+=3)
			{
				int u,p,d;
				p=(x)+(y*STRIDE);	
				u=p+STRIDE;	
				d=p-STRIDE;	
				for (pix=0;pix<3;pix++)
				{
					image1[p+pix]=(image2[p-3+pix]+image2[p+pix]+image2[p+3+pix]+
						image2[u-3+pix]+image2[u+pix]+image2[u+3+pix]+
						image2[d-3+pix]+image2[d+pix]+image2[d+3+pix]
					  )/9;
				}
			}
		}

		// corners
		for (pix=0;pix<3;pix++){ 
			image2[pix]=(image1[pix]+image1[pix+3]+image1[pix+3+STRIDE]+image1[pix+STRIDE])/4;
			image2[xe+pix]=(image1[xe+pix]+image1[xe-3+pix]+image1[xe+pix-3+STRIDE]+image1[xe+pix+STRIDE])/4;
			image2[pix+ly]=(image1[pix+ly]+image1[ly+pix+3]+image1[ly+pix+3-STRIDE]+image1[ly+pix-STRIDE])/4;
			image2[xe+pix+ly]=(image1[xe+pix+ly]+image1[xe-3+pix+ly]+image1[xe+pix-3+STRIDE+ly]+image1[xe+pix+STRIDE+ly])/4;
		}
		// top/bottom
		for (x=3;x<STRIDE-4;x+=3)
		{
			for (pix=0;pix<3;pix++){ 
				image2[pix+x]=(image1[pix+x]+image1[pix+3+x]+image1[pix+x-3]+
						image1[pix+x+STRIDE]+image1[pix+3+x+STRIDE]+image1[pix+x-3+STRIDE])/6;
				image2[pix+x+ly]=(image1[pix+x+ly]+image1[pix+3+x+ly]+image1[pix+x-3+ly]+
						image1[pix+x+ly-STRIDE]+image1[pix+3+x+ly-STRIDE]+image1[pix+x-3+ly-STRIDE])/6;
			}
		}

		// left/right
		for (y=1;y<Y_SIZE-1;y++)
		{
			int yp,yr;
			yp=y*STRIDE;
			yr=yp+STRIDE-3;
			for (pix=0;pix<3;pix++){ 
				image2[yp+pix]=(image1[yp+pix]+image1[yp+pix-STRIDE]+image1[yp+pix+STRIDE]+
					image1[yp+pix+3]+image1[yp+pix-STRIDE+3]+image1[yp+pix+STRIDE+3])/6;
				image2[yr+pix]=(image1[yr+pix]+image1[yr+pix-STRIDE]+image1[yr+pix+STRIDE]+
					image1[yr+pix-3]+image1[yr+pix-STRIDE-3]+image1[yr+pix+STRIDE-3])/6;
			}
		}

                for (y=1;y<Y_SIZE-1;y++)
                {
                        for (x=3;x<STRIDE-4;x+=3)
                        {
                                int u,p,d;
                                p=(x)+(y*STRIDE);
                                u=p+STRIDE;
                                d=p-STRIDE;
                                for (pix=0;pix<3;pix++)
                                {
                                        image2[p+pix]=(image1[p-3+pix]+image1[p+pix]+image1[p+3+pix]+
                                                image1[u-3+pix]+image1[u+pix]+image1[u+3+pix]+
                                                image1[d-3+pix]+image1[d+pix]+image1[d+3+pix]
                                          )/9;
                                }
                        }
                }

	}
	free (image1);
}

void line ( int *origin, int r, int g, int b, int length, int angle, int thick)
{
	float phi,dx,dy,along,tx,ty,wide;
	phi=(2*M_PI*(float)angle/360);

	dx=(float)length*sin(phi);
	dy=-(float)length*cos(phi);
	tx=(float)thick*sin(phi+(M_PI/2));
	ty=-(float)thick*cos(phi+(M_PI/2));

	for (along=0;along<length;along+=0.6)
	{
		for (wide=0;wide<thick;wide+=0.6)
		{	
			float xp,yp;
			xp=((along*dx)/(float)length)+((tx*wide)/thick);
			yp=((along*dy)/(float)length)+((ty*wide)/thick);
			plot(origin[0]+xp,origin[1]+yp,r,g,b);
		}
	}
	origin[0]+=dx;
	origin[1]+=dy;
}


void lline (int xs, int ys, int xe, int ye, int r, int g, int b, float ag,float dup,float t,float bb)
{
	float tx,ty,dx,dy,along,len,theta,bend,thick,lh,px,py;

	dx=xe-xs;
	dy=ye-ys;

	len=sqrt((dx*dx)+(dy*dy));

	bend=(float)bb*len;
	

	theta=asin(dy/len);
	tx=sin(theta);
	ty=cos(theta);

	for (along=0;along<len;along+=0.2)
	{
		float bb,wig;
		bb=bend*sin(along*M_PI/len);
		wig=ag*len*sin(along*M_PI*2*dup/len);
		for (thick=-t/2;thick<t/2;thick+=0.2)
		{
			plot(xs+((dx*along)/len)+((bb+thick+wig)*tx),ys+((dy*along)/len)+((bb+thick+wig)*ty),
					((r*along)+(g*(len-along-1)))/len, ((g*along)+(b*(len-along-1)))/len, ((b*along)+(r*(len-along-1)))/len);
					//r,g,b);
		}

	}
}



void arc ( int x, int y , int r, int g, int b, int radius, int angle, int seg, int thick)
{
	float si,phi,t,along,wide,xp,yp;

        phi=(2*M_PI*(float)angle/360);
        si=(2*M_PI*(float)seg/360);

        cord[0]=x+((radius+thick)*(sin(phi)));
        cord[1]=y+((radius+thick)*(cos(phi)));


        for (along=phi;along<phi+si;along+=M_PI/(3*(float)radius))
        {
                for (wide=0;wide<thick;wide+=0.6)
                {
                        xp=x+((radius+wide)*(sin(along)));
                        yp=y+((radius+wide)*(cos(along)));
                        plot(xp,yp,r,g,b);
                }
        }
	cord[2]=xp;
	cord[3]=yp;
}


void bline ( int *origin, int r, int g, int b, int length, int angle, int thick)
{
	float phi,dx,dy,along,tx,ty,wide;
	phi=(2*M_PI*(float)angle/360);

	dx=(float)length*sin(phi);
	dy=-(float)length*cos(phi);
	tx=(float)thick*sin(phi+(M_PI/2));
	ty=-(float)thick*cos(phi+(M_PI/2));

	for (along=0;along<length;along+=0.6)
	{
		for (wide=0;wide<thick;wide+=0.6)
		{	
			float xp,yp;
			xp=((along*dx)/(float)length)+((tx*wide)/thick);
			yp=((along*dy)/(float)length)+((ty*wide)/thick);
			bplot(origin[0]+xp,origin[1]+yp,r,g,b);
		}
	}
}

//circles
void dlot( int x, int y, int r , int g, int b, int sz, int start, int end, int wig)
{
	float rad,per,fi,beta;
	beta=(2*M_PI*start)/360;
	fi=(2*M_PI*(end-start))/360;
	for (rad=1;rad<sz;rad+=0.5)
	{
		per=fi*rad;
		float p;
		for (p=0;p<per;p+=0.5)
		{
			float phi,wigl;
			int xp,yp;
			wigl=(rad/5)*sin(phi*wig);
			phi=((fi*p)/(per))+beta;
			xp=x+((rad+wigl)*sin(phi));
			yp=y+((rad+wigl)*cos(phi));
			plot(xp,yp,r,g,b);
		}
	}
}

//coloured squares
void square (int x,int y,int r,int g,int b,int xs, int ys)
{
	int xx,yy;
	for (xx=0;xx<xs;xx++)
	{
		for (yy=0;yy<ys;yy++)
		{
			plot(x+xx,y+yy,r,g,b);
		}
	}
}	

//galaxies
void clot (int x,int y,int r,int g,int b,int bad,int ang)
{
	int rad;
	rad=bad/40;
	if (rad==0){return;}
	float t,p,perim,turns,radius;
	turns=rad/10;

	perim=rad*2*M_PI;
	radius=0;
	for (p=0;p<perim*turns;p+=1)
	{
		float xo,yo,phi;
		radius+=rad/(turns*perim);
		phi=p*2*M_PI/(float)perim;
		xo=radius*sin(phi);
		yo=ang*radius*cos(phi)/360;
		if (rand()%40<1){plot (x+xo,y+yo,r,g,b);
		clot(x+xo,y+yo,r,g,b,bad/10,ang);
		}
	}
		
}

// plot into the final
void plot (int x,int y,int r,int g,int b)
{
	if (x>=X_SIZE){ return;}
	if (x<0){ return;}
	if (y>=Y_SIZE){ return;}
	if (y<0){ return;}
	int xpoint;
	xpoint=(y*X_SIZE*3)+(x*3);
	image2[xpoint]=r;
	image2[xpoint+1]=g;
	image2[xpoint+2]=b;
}	

// plot into the buffer
void bplot (int x,int y,int r,int g,int b)
{
	if (x>=X_SIZE){ return;}
	if (x<0){ return;}
	if (y>=Y_SIZE){ return;}
	if (y<0){ return;}
	int xpoint;
	xpoint=(y*X_SIZE*3)+(x*3);
	image5[xpoint]=r;
	image5[xpoint+1]=g;
	image5[xpoint+2]=b;
}

void cloud (int tot, int darkness, int weight)
{
	int seed,x,y,cp,cb,cu,cd,p,b,loop;
	for (loop=0;loop<X_SIZE*Y_SIZE;loop++){clud[loop]=0;clue[loop]=0;}

	for (seed=0;seed<weight;seed++)
	{
		int sx,sy,ss;
		ss=seq[seed]%100; 
		sx=(seq[seed+1]%X_SIZE+(tot*(5-(seq[seed]%11))))%X_SIZE;
		sy=seq[seed+2]%Y_SIZ2;
		cp=sx+(X_SIZE*sy);

		float rad,per,sz;
		int b;
		b=127+seq[seed]%128;

		sz=(weight/4)+((Y_SIZ2-sy)/8)+ss;

		for (rad=1;rad<sz;rad+=0.5)
		{
			per=2*M_PI*rad;
			float p;
			for (p=0;p<per;p+=0.5)
			{
			float phi;
			int xp,yp;
			phi=((2*M_PI*p)/(per));
			xp=sx+(2*rad*sin(phi));
			yp=sy+(rad*cos(phi));
			if (xp>0 && yp>0 && xp<X_SIZE && yp<Y_SIZ2)
			{
				clud[xp+(yp*X_SIZE)]=b;
			}
			}
		}
	}

	int blur;
	for (blur=0;blur<10;blur++)
	{

	for (x=1;x<X_SIZE-1;x++)
	{
		clue[x]=(clud[x-1]+clud[x+1]+clud[x]+
				clud[X_SIZE+x-1]+clud[X_SIZE+x+1]+clud[X_SIZE+x])/6;
	}
	for (y=1;y<Y_SIZ2;y++)
	{
		p=(X_SIZE*y);
		for (x=1;x<X_SIZE-1;x++)
		{
			cp=p+x;
			clue[cp]=(clud[cp]+clud[cp+1]+clud[cp-1]+
				clud[cp+X_SIZE]+clud[cp+1+X_SIZE]+clud[cp-1+X_SIZE]+
				clud[cp-X_SIZE]+clud[cp+1-X_SIZE]+clud[cp-1-X_SIZE])/9;
		}
	}


	for (x=1;x<X_SIZE-1;x++)
	{
		clud[x]=(clue[x-1]+clue[x+1]+clue[x]+
				clue[X_SIZE+x-1]+clue[X_SIZE+x+1]+clue[X_SIZE+x])/6;
	}
	for (y=1;y<Y_SIZ2;y++)
	{
		p=(X_SIZE*y);
		for (x=1;x<X_SIZE-1;x++)
		{
                        cp=p+x;
                        clud[cp]=(clue[cp]+clue[cp+1]+clue[cp-1]+
                                clue[cp+X_SIZE]+clue[cp+1+X_SIZE]+clue[cp-1+X_SIZE]+
                                clue[cp-X_SIZE]+clue[cp+1-X_SIZE]+clue[cp-1-X_SIZE])/9;
		}
	}
	}
	for (y=0;y<Y_SIZ2;y++)
	{
		cp=(X_SIZE*y);
		for (x=0;x<X_SIZE;x++)
		{
			int val;
			val=clud[cp+x];
			if (val>60)
			{
				plot(x,y,val,val,val);
			}
		}
	}
}	

        //tree(X_SIZE/4,Y_SIZE-300,loop,0,30,mul);
void tree (int x, int y, int len, int angle, int bright, int offset)
{
	static int fib[] ={ 1,2,3,5,8,13,21,34,55,89,144,233,377,610,987};
	int or[2],rlen,fibc,p;
	or[0]=x;
	or[1]=y;
	rlen=len;
	fibc=0;
	while (rlen>=fib[fibc+1])
	{
		rlen-=fib[fibc];
		fibc++;
	}
	line (or,(40*bright)/100,((255-(20*fibc))*bright)/100,(30*bright)/100,rlen,angle,1+(fib[fibc]/5));

	if( len-rlen >1 ){
		tree (or[0],or[1],len-rlen,angle+offset,bright,offset);
		tree (or[0],or[1],len-rlen,angle-offset,bright,offset);}
	for(p=fibc-1;p>bloom;p--)
	{
		//line (or,40,200-(20*p),30,fib[p],angle,p);
		line (or,(40*bright)/100,((255-(20*fibc))*bright)/100,(30*bright)/100,fib[p],angle,1+(fib[p]/5));
		tree (or[0],or[1],len-rlen,angle+offset,bright,offset);
		tree (or[0],or[1],len-rlen,angle-offset,bright,offset);
		rlen+=fib[p];
	}
}


void sree (int x, int y, int len, int angle, int bright, int offset)
{
        static int fib[] ={ 1,2,3,5,8,13,21,34,55,89,144,233,377,610,987};
        int or[2],rlen,fibc,p;
        or[0]=x;
        or[1]=y;
        rlen=len;
        fibc=0;

        while (rlen>fib[fibc+1])
        {
                rlen-=fib[fibc];
                fibc++;
        }
        line (or,(40*bright)/100,((255-(20*fibc))*bright)/100,(30*bright)/100,rlen,angle,fibc);

        if( len-rlen >1 ){
               if (fibc%2==0){sree (or[0],or[1],len-rlen,angle+offset,bright,offset);}else{
               sree (or[0],or[1],len-rlen,angle-offset,bright,offset);}
	}
        for(p=fibc-1;p>bloom;p--)
        {
                line (or,(40*bright)/100,((255-(20*fibc))*bright)/100,(30*bright)/100,fib[p],angle,p);
		if (p%2==0){
                sree (or[0],or[1],len-rlen,angle+offset,bright,offset);} else {
                sree (or[0],or[1],len-rlen,angle-offset,bright,offset);}
                rlen+=fib[p];
        }
	{
	if (bloom<4){plot(or[0],or[1],((2-bloom)*bright),0,0);}
	if (bloom<3){
	plot(or[0],or[1]+1,2*bright,0,0);
	plot(or[0],or[1]-1,2*bright,0,0);
	plot(or[0]+1,or[1],2*bright,0,0);
	plot(or[0]-1,or[1],2*bright,0,0);}
	}
}


void my_load_image( unsigned char *image, char *kk, int r)
{
        int x_size,y_size,x,y;
        unsigned char *buff;
        buff=(unsigned char *)malloc(sizeof (char)*20*4096*4096); // image loader  buffer
        read_JPEG_file (kk, buff, pims);
        printf ("here \n");

        //streach it to fit both ways.
        float x_scale,y_scale,scale;

        x_size=pims[0];
        y_size=pims[1];

        x_scale=(float)x_size/(float)X_SIZE;
        y_scale=(float)y_size/(float)Y_SIZE;

        if (x_scale<=y_scale){ scale=x_scale;}else{scale=y_scale;}

        for (y=0;y<Y_SIZE;y++)
        {
                for (x=0;x<X_SIZE;x++)
                {
                        int b,i;
                        i=(x*3)+(y*X_SIZE*3);
                        //b=((x*scale)*3)+((y*scale)*x_size*3);
                        b=((int)(x*scale)*3)+((int)(y*scale)*x_size*3);
                        image[i]=buff[b];
                        image[i+1]=buff[b+1];
                        image[i+2]=buff[b+2];
                }
        } 
        free (buff);
        printf ("here\n");
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

void show_image( unsigned char *img)
{
	memcpy(image2,img,X_SIZE*Y_SIZE*3);
}
