// background_segmentation.cpp : Defines the entry point for the console application.



#include"videoProcess_final.h"
#define MAX 5
using namespace std;
using namespace cv;



//double threshold;

/*typedef struct ptr_mat
{
	double mat[1000][1000];
	ptr_mat *next;
};*/


int main()
{
	
	Scene obj;
	obj.initScene("video40.avi");
	obj.startProcessing();
	_getch();
	return 0;
}


