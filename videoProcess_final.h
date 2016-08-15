#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/gpu/gpu.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/video/background_segm.hpp>
#include<iostream>
#include<conio.h>
#include<stdlib.h>
#include<vector>
#include<math.h>

using namespace std;
using namespace cv;
class Scene
{

	public:
	Mat frameNow,framePrev;
	
	int frameNum;
	int width,height;
	int N,MaxPoints;
	char locMaxFlag[100];
    int locMax[100];
	double PointX[5000],PointY[5000];
	double X[10][2000],P[2000],Y[10][2000],KalmanGain[2000];
	double theta[100];
	float rho,k,R1_height,R1_width,R2_height,R2_width,angle;
	

	Point2f A,B,C,D,E_mid;
	Point Tip,Neck;
	Point E1,E2,Knee1,knee2;
	Mat matDiff,matDiffPrev,backgroundMask;
	Mat sumBack;
	Mat background,isBack,object;
	Mat avgMat,ones,refineMask,grayFrame,obj;
	Mat mean,variance,sigma;
	vector<vector<Point>> contours;
	VideoCapture cap;
	
	Point elbow2,R1_vertices[4],R2_vertices[4];
		 
	Point2f R1_vertices2f[4],R1_centre,R2_vertices2f[4],R2_centre;
	RotatedRect R1,R2;

	void initScene(const char *videoName);
	void convertToChar(Mat frame,unsigned char* image);
	void copyMatrix(unsigned char *source,unsigned char *dest);
	Mat frameDiffCalc(Mat now,Mat prev);
	void findBackground(Mat matDffNow,Mat matDffNowPrev);
	double sumArr(double[],int,int);
	void findBackground(Mat frame);
	Mat findObject(Mat frame);
	void startProcessing();
	Mat updateAvg(Mat avgMat,Mat frame);
	void updateStats(Mat frame);
	Mat refineObject(Mat objectIn);
	void findcontours(Mat,int);
	int distance_curve(Mat,Point2f,vector<vector<Point>>,int);
	void head_detection(Mat,double,int,Point2f);
	void feet_detection(Mat,int,Point2f);
	void hand_detection(Mat,int,Point2f);
	int find_min(vector<Point2f>,double[]);
	int find_head(vector<Point2f>,double[],int);
	void findmax2(Mat,double[],vector<Point2f>,int,char *,Point2f);
	double cal_distance(Point2f,Point2f);
	Point find_tip_head(Mat,Point,Point,Point2f);
	//void findmax2_hand(Mat,double [],vector<Point2f>,int );
	//int kalmanPredict(int,int,int,int);
	//void kalmanUpdate(int,double ,double, int , int);
	void elbow_shoulder_knee_detection(Mat,int,Point2f,Point2f,vector<vector<Point> >,double);
	Point find_elbow_knee(Point2f,Point2f,Point2f,vector<vector<Point>>);
	void cloth_texture_fitting(vector<vector<Point> >,Mat,Point2f);
	void cloth_fitting();

};

