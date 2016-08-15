//#include "stdafx.h"
//
//#include"videoProcess.h"
//#include<fstream>
//
//void Scene::initScene(const char *videoName)
//{
//	cap.open(videoName);
//	//cap.open(0);
//	cap.read(frameNow);
//	
//	frameNum=0;
//	width=frameNow.cols;
//	height=frameNow.rows;
//	N=20;
//	
//	matDiff.create(height,width,CV_8UC1);
//	matDiffPrev.create(height,width,CV_8UC1);
//
//	//backgroundMask.create(height,width,CV_8UC3);
//	object.create(height,width,CV_8UC3);
//	background.create(height,width,CV_8UC1);
//
//	background=Mat::zeros(height,width,CV_8UC1);
//	backgroundMask=Mat::zeros(height,width,CV_8UC1);
//	sumBack=Mat::zeros(height,width,CV_8UC1);
//	ones=5*Mat::ones(height,width,CV_8UC1);
//	refineMask=Mat::zeros(height,width,CV_8UC1);
//	framePrev.create(height,width,CV_8UC3);
//	avgMat=frameNow.clone();
//	framePrev=frameNow.clone();
//	isBack=Mat(height,width,CV_8UC1,0);
//	mean=frameNow.clone();
//    drawing=Mat::zeros(height,width,CV_8UC1);
//	obj=Mat::zeros(height,width,CV_8UC1);
//	rho=0.01;
//	k=2.5;
//}
//
//void Scene::convertToChar(Mat frame,unsigned char* image)
//{
//	
//	for(int i = 0; i < frame.rows; i++)
//    {
//		for(int j = 0; j < frame.cols; j++)
//	    {
//              image[frame.cols*i+j] = frame.at<cv::Vec3b>(i,j)[0];
//              image[frame.cols*i+j+1] = frame.at<cv::Vec3b>(i,j)[1];
//              image[frame.cols*i+j+2] = frame.at<cv::Vec3b>(i,j)[2];
//        }
//    }
//}
//
//void Scene::copyMatrix(unsigned char *source,unsigned char *dest)
//{
//	memcpy(dest,source,sizeof(unsigned char)*3*width*height);
//}
//
//Mat Scene::frameDiffCalc(Mat now, Mat prev)
//{
//	int i,j,diff=0;
//	float matDiff=0;
//	Mat res(now.rows,now.cols,CV_8UC1);
//	
//	cvtColor(now,now,COLOR_RGB2GRAY);
//	cvtColor(prev,prev,COLOR_RGB2GRAY);
//	
//	res = cv::abs(now - prev);
//	//GaussianBlur(res,res,Size(51,51),0,0);
//	threshold(res, res,12, 255, cv::THRESH_BINARY);
//	
//	//adaptiveThreshold(res,res,255,ADAPTIVE_THRESH_GAUSSIAN_C,THRESH_BINARY_INV,11,7);
//	//imshow("Diff",res);
//	return res;
//}
//
//void Scene::findBackground(Mat matDffNow,Mat matDffNowPrev)
//{
//	Mat res(matDffNow.rows,matDffNow.cols,CV_8UC1);
//	res = cv::abs(matDffNow - matDffNowPrev);
//
//	Mat isDiff,thresBack;
//	threshold(res,isDiff,1,1,CV_THRESH_BINARY_INV);
//	
//	//int a=(res==1);
//	//Mat grayFrame;
//	
//	//cvtColor(frameNow,grayFrame,CV_RGB2GRAY);
//	
//	sumBack=((isDiff==1 & sumBack<80)/255).mul(sumBack+isDiff);
//	
//	//sumBack=sumBack+isDiff;
//	
//	background+=((background==0)/255).mul(((sumBack==79)/255.0).mul(grayFrame));
//	backgroundMask+=((sumBack==79));
//	
//		
//	sumBack+=(sumBack==79).mul(ones/255.0);
//
//
//	int stagnateThres=100;
//	double maxVal=0,minVal;
//	cv::minMaxLoc(sumBack,&minVal,&maxVal);
//
//	//cout<<"Max="<<maxVal;
//	//int i=0,j=0;
//	//for(i=0;i<res.rows;i++)
//	//{
//	//	for(j=0;j<res.cols;j++)
//	//	{
//	//		if(sumBack[i*res.cols+j]>sta1gnateThres)
//	//		{
//	//			backgroundMask.at<cv::Vec3b>(i,j)[0]=255;
//	//			backgroundMask.at<cv::Vec3b>(i,j)[1]=255;
//	//			backgroundMask.at<cv::Vec3b>(i,j)[2]=255;
//
//	//			continue;
//	//		}
//
//	//		if(sumBack[i*res.cols+j]==stagnateThres)
//	//		{	
//	//			background.at<cv::Vec3b>(i,j)[0]=frameNow.at<cv::Vec3b>(i,j)[0];
//	//			background.at<cv::Vec3b>(i,j)[1]=frameNow.at<cv::Vec3b>(i,j)[1];
//	//			background.at<cv::Vec3b>(i,j)[2]=frameNow.at<cv::Vec3b>(i,j)[2];
//	//			sumBack[i*res.cols+j]+=5;
//	//			continue;		
//	//		}
//
//
//	//		backgroundMask.at<cv::Vec3b>(i,j)[0]=0;
//	//		backgroundMask.at<cv::Vec3b>(i,j)[1]=0;
//	//		backgroundMask.at<cv::Vec3b>(i,j)[2]=0;
//
//	//		background.at<cv::Vec3b>(i,j)[0]=0;
//	//		background.at<cv::Vec3b>(i,j)[1]=0;
//	//		background.at<cv::Vec3b>(i,j)[2]=0;
//
//	//		if(res.at<uchar>(i,j)>0)
//	//			sumBack[i*res.cols+j]=0;
//	//		else
//	//			sumBack[i*res.cols+j]=sumBack[i*res.cols+j]+1;
//	//	
//	//	}
//	//}
//	////GaussianBlur(backgroundMask,backgroundMask,Size(25,25),2,2);
//	////threshold(backgroundMask,backgroundMask,100, 1, cv::THRESH_BINARY);
//	imshow("Mask",backgroundMask);
//
//}
//Mat Scene::findObject(Mat frame)
//{
//	Mat grayBack;
//	
//	grayBack=background.clone();
//	//cvtColor(background,grayBack,CV_RGB2GRAY);
//	//cvtColor(frame,grayFrame,CV_RGB2GRAY);
//	
//
//	Mat backDiff=cv::abs(grayFrame-grayBack);
//	//Mat backDiff=cv::abs(frame-grayBack);
//	
//
//	threshold(backDiff,backDiff,17, 1, cv::THRESH_BINARY);
//	GaussianBlur(backDiff,backDiff,Size(25,25),2,2);
//	threshold(backDiff,backDiff,0.5, 1, cv::THRESH_BINARY);
//	imshow("Object",backDiff.mul(grayFrame));
//	//imshow("Object",backDiff.mul(frame));
//
//	return backDiff;
//}
//Mat Scene::refineObject(Mat objectIn)
//{
//	Mat abc;
//	int size=32,smallSize=4;
//	//blur(objectIn,abc,Size(40,40));
//	cv::Scalar avg,avg2;Mat tmpImg;
//	refineMask=Mat::zeros(height,width,CV_8UC1);
//	for(int i=0;i<height/size;i++)
//	{
//		for(int j=0;j<width/size;j++)
//		{
//			cv::Rect roi(j*size,i*size,size,size);
//			avg=	cv::mean(objectIn(roi));
//			//tmpImg=objectIn(roi).clone();
//			if(avg[0]>0.5)
//			{	
//				//refineMask(roi)=Mat::ones(size,size,CV_8UC1);
//				for(int ii=0;ii<size/smallSize;ii++)
//				{
//					for(int jj=0;jj<size/smallSize;jj++)
//					{
//						cv::Rect roi2(j*size+jj*smallSize,i*size+ii*smallSize,smallSize,smallSize);
//						avg2=cv::mean(objectIn(roi2));
//						if(avg2[0]>0.5)
//							refineMask(roi2)=Mat::ones(smallSize,smallSize,CV_8UC1);
//					}
//				}
//			}
//		}
//	}
//	imshow("Refine Mask",refineMask.mul(grayFrame));
//	return (refineMask);
//}
//
//double Scene::sumArr(double D[5000],int low,int upper)
//{
//    double arr_sum=0;
//    for(int i=low;i<=upper;i++)
//        arr_sum+=D[i];
//
//    return arr_sum;
//}
//void Scene::findcontours(Mat ob,int count)
//{ 
//	double largest_area = 0,area=0.0;
//	Moments mu;
//	Point2f mc;
//    double orientation,major_axis,minor_axis;
//
//	//vector<vector<Point>> contours;
//	//vector<vector<Point>> largest_contours;
//	//Mat drawing=Mat::zeros(height,width,CV_8UC1);
//	//Mat obj=Mat::zeros(height,width,CV_8UC1);
//	findContours(ob,contours,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_SIMPLE);
//	vector<vector<Point> > contours_poly( contours.size() );
//	vector<Rect> boundRect( contours.size() );
//	for( int i= 0; i < contours.size(); i++)
//	{  // get the largest contour
//        area = fabs(contourArea(contours[i]));
//        if(area >= largest_area){
//            largest_area = area;
//            largest_contours.clear(); 
//            largest_contours.push_back(contours[i]);
//        }
//    }
//	//Moment calculation
//    for( int i = 0; i < largest_contours.size(); i++ )
//     { 
//		 mu = moments( largest_contours[i], false );
//	 }
//	
//    //  Get the mass centers:
//     for( int i = 0; i < largest_contours.size(); i++ )
//     { 
//		 mc = Point2f( mu.m10/mu.m00 , mu.m01/mu.m00 ); 
//		 orientation = 0.5*atan2((2*mu.m11),(mu.m20 - mu.m02))*180/3.14 ; 
//		 //major_axis=sqrt(((mu.m20/mu.m00) + (mu.m02/mu.m00))/2 + sqrt(((4*(mu.m11/mu.m00)*(mu.m11/mu.m00)) - ((mu.m20/mu.m00 - mu.m02/mu.m00)*(mu.m20/mu.m00 - mu.m02/mu.m00)))/2));
//		 //minor_axis=sqrt(((mu.m20/mu.m00) + (mu.m02/mu.m00))/2 - sqrt(((4*(mu.m11/mu.m00)*(mu.m11/mu.m00)) - ((mu.m20/mu.m00 - mu.m02/mu.m00)*(mu.m20/mu.m00 - mu.m02/mu.m00)))/2))
//		 approxPolyDP( Mat(largest_contours[i]), contours_poly[i], 3, true );
//         boundRect[i] = boundingRect( Mat(contours_poly[i]) );
//		 //lengthRect=abs(boundRect[i].br().y-boundRect[i].tl().y);
//		 //widthRect=abs(boundRect[i].br().x-boundRect[i].tl().x);
//		 A=Point2f((boundRect[i].tl().x+boundRect[i].br().x)/2,boundRect[i].tl().y);
//         B=Point2f((boundRect[i].tl().x+boundRect[i].br().x)/2,boundRect[i].br().y);
//		 C=Point2f(boundRect[i].tl().x,(boundRect[i].tl().y+boundRect[i].br().y)/2);
//		 D=Point2f(boundRect[i].br().x,(boundRect[i].tl().y+boundRect[i].br().y)/2);
//		 //LW=lengthRect/widthRect;
//		 //DWT=abs((lengthRect/478)-(widthRect/638));
//	 }
//	 //cout << "length " << lengthRect << endl;
//	 //cout << "width " << widthRect << endl;
//	 //cout << "Ratio " << DWT << endl;
//
//	 
//	 //Drawing the contour of the largest area
//	 for(int i=0;i<largest_contours.size();i++)
//	 {
//		 if(largest_area>=4500)
//		 {
//		  drawContours(drawing,contours_poly,-1,Scalar(255,0,0),1);
//		  rectangle( drawing, boundRect[i].tl(), boundRect[i].br(),Scalar(255,0,0), 2, 8, 0 );
//		  circle( drawing, mc, 4,Scalar(255,0,0), -1);
//		  //if(major_axis>=0 && minor_axis>=0)
//		  /*if((boundRect[i].br().y-boundRect[i].tl().y) >= 0 && (boundRect[i].br().y-boundRect[i].tl().y)>=0)
//		  {
//			  //cout<<"AM here"<<endl;
//		      //ellipse(drawing,Point(mc.x,mc.y),Size(minor_axis/2,major_axis/2),orientation,0,360,Scalar(255,0,0),2,8,0);	  
//			 // rectangle( drawing, boundRect[i].tl(), boundRect[i].br(),Scalar(255,0,0), 2, 8, 0 );
//			  ellipse(drawing,RotatedRect( Point2f((boundRect[i].tl().x+boundRect[i].br().x)/2,(boundRect[i].tl().y+boundRect[i].br().y)/2),Size(boundRect[i].br().x-boundRect[i].tl().x,boundRect[i].br().y-boundRect[i].tl().y),orientation),Scalar(255,0,0),1,8);
//		  }*/
//		   
//		 }
//	 }
//	 int locMax_index = distance_curve(mc);
//	 head_detection(orientation,locMax_index);
//	// feet_detection(locMax_index);
//	 //hand_detection(locMax_index,mc);
//	 imshow("boundary",drawing);
//}
//	 //Distance calculation
//int Scene :: distance_curve(Point2f mass_centre)
//{
//	 int sizeD=0;
//	 double Dist[5000],Dlow[5000],grad[5000],gradLow[5000];
//	 int index=0,k=0;
//     for(int i=0;i<largest_contours.size();i++)
//	 {
//		 for(int j=0;j<largest_contours[i].size();j++)
//		 {
//			 Dist[index]=sqrt((largest_contours[i][j].x-mass_centre.x)*(largest_contours[i][j].x-mass_centre.x)+(largest_contours[i][j].y-mass_centre.y)*(largest_contours[i][j].y-mass_centre.y));
//			 PointX[index]=largest_contours[i][j].x;
//			 PointY[index]=largest_contours[i][j].y;
//			 index++;
//		 }
//	 }
//	 //low pass filtering of distance array
//	 int len=40;
//	 for(int i=len+1; i<index-len;i++)
//	 {
//        Dlow[i]=sumArr(Dist,i-len,i+len) /(2*len+1);
//		
//	 }
//		for(int i=len+2;i<index-len;i++)
//            grad[i]=Dlow[i]-Dlow[i-1];
//
//     for(int i=2*len+1;i<index-2*len;i++)
//     {
//		 gradLow[i]= sumArr(grad,i-len,i+len)/(2*len+1);
//	 }
//     int diff=index/25;
//	 //finding local maxima
//     for(int i=1;i<index;i++)
//     {
//        if((gradLow[i]/gradLow[i-1])<0 && i>diff+1)
//        {
//            int flag=0;
//            for(int j=i-diff;j<i-1;j++)
//                if(gradLow[j]<0)
//                    flag=1;
//
//            for(int j=i+1;j<i+diff;j++)
//                if(gradLow[j]>0)
//                    flag=1;
//            if(flag==0)
//            {    
//				locMax[k++]=i;
//				
//			}
//        }
//    }
//	 // drawing the convex points(local maxima)
//	for(int i=0;i<k;i++)
//	{
//		circle(drawing,Point2f(PointX[locMax[i]],PointY[locMax[i]]),4,Scalar(255,0,0), -1);
//		locMaxFlag[i]=0;
//	}
//	
//	//finding curvature of convex points
//    vector<Point2f> U1(k);
//    vector<Point2f> U2(k);
//	for(int i=0;i<k;i++)
//	{
//		U1[i].x=PointX[locMax[i]+10]-PointX[locMax[i]];
//		U1[i].y=PointY[locMax[i]+10]-PointY[locMax[i]];
//		U2[i].x=PointX[locMax[i]-10]-PointX[locMax[i]];
//		U2[i].y=PointY[locMax[i]-10]-PointY[locMax[i]];
//		theta[i]=acos(((U1[i].x * U2[i].x)+(U1[i].y * U2[i].y)) / ((sqrt((U1[i].x*U1[i].x) + (U1[i].y*U1[i].y))) * (sqrt((U2[i].x*U2[i].x) + (U2[i].y*U2[i].y))))) * 180.0/3.14;
//	}
//	return k;
//}
//
//	/*cout << "------Frame Start ------" << endl;
//	for(int i=0;i<k;i++)
//	    cout<<"PointY "<<PointY[locMax[i]] <<"C.y "<< C.y <<"A.y " <<A.y << "PointX "<< PointX[locMax[i]] << "C.x "<<C.x <<"D.x "<<D.x << endl;
//	
////	cout << "------Frame End --------" << endl;*/
////	if(orientation>=(15) && orientation<=(175))
////	{
////		//cout  << "orientation is  : " << orientation <<endl;
////	    for(int i=0;i<k;i++)
////		{
////			if(PointY[locMax[i]]<C.y && PointY[locMax[i]]>A.y && PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x )
////			{	
////				POI[poi_index++] = Point2f(PointX[locMax[i]],PointY[locMax[i]]);
////		        //cout<< "here I am" << endl;
////			}
////		}
////	
////		//cout << "Point of interest above line CD " <<poi_index <<endl;
////		if(poi_index>1)
////		{
////			int head_index=find_head(theta,poi_index);
////			//int theta_index=find_min(theta,poi_index);
////			putText(drawing,"Head",POI[head_index],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////		    locMaxFlag[head_index]=1;
////		}
////		else if(poi_index==1)
////		{	
////			putText(drawing,"Head",POI[poi_index-1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////			locMaxFlag[poi_index-1]=1;
////		}
////		else
////			putText(drawing,"Head",A,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////	}
//////	cout << "------Frame Start ------" <<endl;
////	int feet_index=0,feet_count=0,index_feet=0,hand_count=0,index_hand=0;
////	double feet_distance[10],hand_distance[10];
////	for(int i=0;i<k;i++)
////	{
////		if(PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x && PointY[locMax[i]]>C.y && PointY[locMax[i]]<B.y )
////		    POI_Feet[feet_index++]= Point2f(PointX[locMax[i]],PointY[locMax[i]]);
////
////	    for(int j=0;j<feet_index;j++)
////		{
////			//cout<<"curvature angle "<<theta[j]<<endl;
////			if(theta[j] > 100)
////				POI_Feet[j]= Point2f(0.0,0.0);
////			else
////				feet_count++;
////		}
////	}
////	if(feet_count==1)
////	{
////		putText(drawing,"Feet",POI[feet_index-1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////		locMaxFlag[feet_index-1]=1;
////	}
////	else if(feet_count==2)
////	{
////		for(int j=0;j<feet_index;j++)
////		{
////			if(POI_Feet[j]!=Point2f(0.0,0.0))
////			{
////				putText(drawing,"Feet",POI[j],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////				locMaxFlag[j]=1;
////			}
////		}
////	}
////	else if(feet_count>2)
////	{
////		for(int j=0;j<feet_index;j++)
////		{
////			if(POI_Feet[j]!=Point2f(0.0,0.0))
////				feet_distance[index_feet++]=cal_distance(A,POI_Feet[j]);
////		}
////		findmax2(drawing,feet_distance,POI_Feet,index_feet);
////	}
////	int hand_index=0;
////	int handPoints[10]={0};
////	//cout << "------Frame end -------" << endl;
////	for(int i=0;i<k;i++)
////	{
////		if(locMaxFlag[i]!=1)
////			handPoints[hand_index++]=i;
////	}
////	if(hand_index==2)
////	{
////		POI_Hand[0] = Point2f(PointX[locMax[handPoints[0]]],PointY[locMax[handPoints[0]]]);
////		putText(drawing,"Hand",POI_Hand[0],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////
////		POI_Hand[1] = Point2f(PointX[locMax[handPoints[1]]],PointY[locMax[handPoints[1]]]);
////		putText(drawing,"Hand",POI_Hand[1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////	}
////	else if(hand_count>2)
////	{
////		for(int j=0;j<hand_index;j++)
////		{
////			POI_Hand[j] = Point2f(PointX[locMax[handPoints[j]]],PointY[locMax[handPoints[j]]]);
////			hand_distance[index_hand++]=cal_distance(mc,POI_Hand[j]);
////		}
////		findmax2_hand(drawing,hand_distance,POI_Hand,index_hand);
////	}
////	else
////	{
////		putText(drawing,"Hand",C,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////		putText(drawing,"Hand",D,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
////	}
////	/*for(int j=0;j<k;j++)
////		cout<<theta[j]<<endl;
////	cout<<"-----done----";*/
//
//
//void Scene :: head_detection(double orient,int s)
//{
//	int poi_index=0;
//	vector<Point2f> POI(s);
//	if(orient>=15 && orient<=175)
//	{
//		//cout  << "orientation is  : " << orient <<endl;
//	    for(int i=0;i<s;i++)
//		{
//			if(PointY[locMax[i]]<C.y && PointY[locMax[i]]>A.y && PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x )
//			{	
//				POI[poi_index++] = Point2f(PointX[locMax[i]],PointY[locMax[i]]);
//		        cout<< "here I am" << endl;
//			}
//		}
//	
//		//cout << "Point of interest above line CD " <<poi_index <<endl;
//		if(poi_index>1)
//		{
//			int head_index=find_head(theta,poi_index);
//			//int theta_index=find_min(theta,poi_index);
//			putText(drawing,"Head",POI[head_index],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		    locMaxFlag[head_index]=1;
//		}
//		else if(poi_index==1)
//		{	
//			putText(drawing,"Head",POI[poi_index-1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//			locMaxFlag[poi_index-1]=1;
//		}
//		else
//			putText(drawing,"Head",A,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//	}
//}
//
//void Scene :: feet_detection(int s)
//{
//	int feet_index=0,feet_count=0,index_feet=0;
//	double feet_distance[10];
//	vector<Point2f> POI_Feet(s);
//	for(int i=0;i<s;i++)
//	{
//		if(PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x && PointY[locMax[i]]>C.y && PointY[locMax[i]]<B.y )
//		    POI_Feet[feet_index++]= Point2f(PointX[locMax[i]],PointY[locMax[i]]);
//
//	    for(int j=0;j<feet_index;j++)
//		{
//			//cout<<"curvature angle "<<theta[j]<<endl;
//			if(theta[j] > 100)
//				POI_Feet[j]= Point2f(0.0,0.0);
//			else
//				feet_count++;
//		}
//	}
//	if(feet_count==1)
//	{
//		putText(drawing,"Feet",POI_Feet[feet_index-1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		locMaxFlag[feet_index-1]=1;
//	}
//	else if(feet_count==2)
//	{
//		for(int j=0;j<feet_index;j++)
//		{
//			if(POI_Feet[j]!=Point2f(0.0,0.0))
//			{
//				putText(drawing,"Feet",POI_Feet[j],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//				locMaxFlag[j]=1;
//			}
//		}
//	}
//	else if(feet_count>2)
//	{
//		for(int j=0;j<feet_index;j++)
//		{
//			if(POI_Feet[j]!=Point2f(0.0,0.0))
//				feet_distance[index_feet++]=cal_distance(A,POI_Feet[j]);
//		}
//		findmax2(feet_distance,POI_Feet,index_feet,"Feet");
//	}
//}
//
//void Scene :: hand_detection(int s,Point2f mass_centre)
//{
//    int hand_index=0,hand_count=0,index_hand=0;
//	int handPoints[10]={0};
//	double hand_distance[10];
//	vector<Point2f> POI_Hand(s);
//	for(int i=0;i<s;i++)
//	{
//		if(locMaxFlag[i]!=1)
//			handPoints[hand_index++]=i;
//	}
//	if(hand_index==2)
//	{
//		POI_Hand[0] = Point2f(PointX[locMax[handPoints[0]]],PointY[locMax[handPoints[0]]]);
//		putText(drawing,"Hand",POI_Hand[0],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//
//		POI_Hand[1] = Point2f(PointX[locMax[handPoints[1]]],PointY[locMax[handPoints[1]]]);
//		putText(drawing,"Hand",POI_Hand[1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//	}
//	else if(hand_count>2)
//	{
//		for(int j=0;j<hand_index;j++)
//		{
//			POI_Hand[j] = Point2f(PointX[locMax[handPoints[j]]],PointY[locMax[handPoints[j]]]);
//			hand_distance[index_hand++]=cal_distance(mass_centre,POI_Hand[j]);
//		}
//		findmax2(hand_distance,POI_Hand,index_hand,"Hand");
//	}
//	else
//	{
//		putText(drawing,"Hand",C,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		putText(drawing,"Hand",D,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//	}
//}
//
///*void Scene :: findmax2_hand(Mat draw,double arr[10],vector<Point2f> POI,int index)
//{
//	double temp1=arr[0],temp2;
//	int max1=0,max2=0;
//
//	for(int i=1;i<index;i++)
//	{
//		if(arr[i]>temp1)
//		{
//			temp2=temp1;
//			temp1=arr[i];
//			max2=max1;
//			max1=i;
//		}
//	}
//
//    putText(draw,"Hand",POI[max1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//	putText(draw,"Hand",POI[max2],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//
//}*/
//void Scene:: findmax2(double arr[10],vector<Point2f> P,int index,char *str)
//{
//	double temp1=arr[0],temp2;
//	int max1=0,max2=0;
//
//	for(int i=1;i<index;i++)
//	{
//		if(arr[i]>temp1)
//		{
//			temp2=temp1;
//			temp1=arr[i];
//			max2=max1;
//			max1=i;
//		}
//	}
//	if (str == "Feet")
//	{
//		putText(drawing,str,P[max1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		putText(drawing,str,P[max2],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		locMaxFlag[max1]=1;
//		locMaxFlag[max2]=1;
//	}
//	else
//	{
//		putText(drawing,"Hand",P[max1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//		putText(drawing,"Hand",P[max2],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
//    }
//
//}
//
//double Scene:: cal_distance(Point2f init,Point2f final)
//{
//	return sqrt((init.x-final.x)*(init.x-final.x)+(init.y-final.y)*(init.y-final.y));
//}
//int Scene :: find_min(double arr[10],int index)
//{
//	double temp=0.0;
//	int th_index=0;
//	temp=arr[0];
//	for(int i=1;i<index;i++)
//	{
//		if(arr[i]<temp)
//		{
//			temp=arr[i];
//			th_index=i;
//		}
//	}
//	return th_index;
//}
//
//int Scene :: find_head(double arr[10],int index)
//{
//	int c=0,theta_index=0;
//	for(int i=0;i<index;i++)
//	{
//		if(arr[i]>50 && arr[i]<150)
//		{ 
//			theta_index=i;
//			c++;
//		}
//	}
//	if(c>1)
//		theta_index=find_min(arr,index);
//    return theta_index;
//}
////find back grnd by gaussian pdf
//void Scene::findBackground(Mat frame)
//{
//	Mat val;
//	val=cv::abs(frame-mean);
//	val.convertTo(val,CV_32F);
//	
//	
//	val=val.mul(1/sigma);
//	cvtColor(val,val,CV_RGB2GRAY);
//	threshold(val, backgroundMask,k, 255, cv::THRESH_BINARY);
//	threshold(val, isBack,k, 1, cv::THRESH_BINARY);
//	//cvtColor(backgroundMask,backgroundMask,CV_RGB2GRAY);
//	//backgroundMask.convertTo(backgroundMask,CV_8UC1);
//	imshow("mask",backgroundMask);
//
//	//cout<<isBack;
//
//}
//
//void Scene::startProcessing()
//{
//	int i,j;
//	int count=0;
//	
//	do
//	{
//		
//		cap.read(frameNow);
//		cvtColor(frameNow,grayFrame,CV_RGB2GRAY);
//		Mat image32f;
//		if(count==0)
//		{
//			Mat image32f;
//			frameNow.convertTo(image32f,CV_32F);
//
//			Mat mu;
//			blur(image32f, mu, Size(3, 3));
//
//			Mat mu2;
//			blur(image32f.mul(image32f), mu2, Size(3, 3));
//
//			cv::sqrt(mu2 - mu.mul(mu), sigma);
//
//			variance=sigma.mul(sigma);
//
//			//imshow("sigma",sigma);
//
//			count++;
//		}
//		avgMat=updateAvg(avgMat,frameNow);
//		//updateStats(frameNow);
//		//findBackground(frameNow);
//		matDiff=frameDiffCalc(frameNow,framePrev);
//		//framePrev=frameNow;
//		framePrev=frameNow.clone();
//		frameNum++;
//		findBackground(matDiff,matDiffPrev);
//		matDiffPrev=matDiff.clone();
//		object=findObject(frameNow);
//		Mat res=refineObject(object);
//		findcontours(object,count);
//		count++;
//		//Mat element = getStructuringElement(MORPH_RECT,Size(3,3),Point(0,0));
//		
//		//erode(object,object,kernel);
//		//Mat dst=object-eroded;
//		
//
//		imshow("yahoo",frameNow);
//
//		if (waitKey(30) >= 0)
//				break;
//		  
//		 //imshow("sigma",frameNow);
//	 }while(1);
//}
//
//Mat Scene::updateAvg(Mat avgFrame,Mat frame)
//{
//	Mat res;
//	res=(avgFrame*(N-1)+frame)/N;
//	return res;
//}
//void Scene::updateStats(Mat frame)
//{
//	mean=rho*frame+(1-rho)*mean;
//
//	
//	Mat dist=cv::abs(frame-mean);
//	dist.convertTo(dist,CV_32F);
//	variance=rho*dist.mul(dist)+(1-rho)*variance;
//
//	cv::sqrt(variance,sigma);
//}

//#include "stdafx.h"

#include"videoProcess_final.h"
#include<fstream>
#define FRAME_COUNT 430
void Scene::initScene(const char *videoName)
{
	cap.open(videoName);
	//cap.open(0);
	cap.read(frameNow);
	cap.set(CV_CAP_PROP_FRAME_WIDTH, 640);
    cap.set(CV_CAP_PROP_FRAME_HEIGHT, 360);
	frameNum=0;
	width=frameNow.cols;
	height=frameNow.rows;
	N=20;
	//X[0]=0.0,Y[0]=0.0;
	//P[0]=1.0,P[1]=1.0;
	matDiff.create(height,width,CV_8UC1);
	matDiffPrev.create(height,width,CV_8UC1);

	//backgroundMask.create(height,width,CV_8UC3);
	object.create(height,width,CV_8UC3);
	background.create(height,width,CV_8UC1);

	background=Mat::zeros(height,width,CV_8UC1);
	backgroundMask=Mat::zeros(height,width,CV_8UC1);
	sumBack=Mat::zeros(height,width,CV_8UC1);
	ones=5*Mat::ones(height,width,CV_8UC1);
	refineMask=Mat::zeros(height,width,CV_8UC1);
	framePrev.create(height,width,CV_8UC3);
	avgMat=frameNow.clone();
	framePrev=frameNow.clone();
	isBack=Mat(height,width,CV_8UC1,0);
	mean=frameNow.clone();
    obj=Mat::zeros(height,width,CV_8UC1);
	rho=0.01;
	k=2.5;
}

void Scene::convertToChar(Mat frame,unsigned char* image)
{
	
	for(int i = 0; i < frame.rows; i++)
    {
		for(int j = 0; j < frame.cols; j++)
	    {
              image[frame.cols*i+j] = frame.at<cv::Vec3b>(i,j)[0];
              image[frame.cols*i+j+1] = frame.at<cv::Vec3b>(i,j)[1];
              image[frame.cols*i+j+2] = frame.at<cv::Vec3b>(i,j)[2];
        }
    }
}

void Scene::copyMatrix(unsigned char *source,unsigned char *dest)
{
	memcpy(dest,source,sizeof(unsigned char)*3*width*height);
}

Mat Scene::frameDiffCalc(Mat now, Mat prev)
{
	int i=0,j=0,diff=0;
	float matDiff=0;
	Mat res(now.rows,now.cols,CV_8UC1);
	
	cvtColor(now,now,COLOR_RGB2GRAY);
	cvtColor(prev,prev,COLOR_RGB2GRAY);
	
	res = cv::abs(now - prev);
	//GaussianBlur(res,res,Size(51,51),0,0);
	threshold(res, res,12, 255, cv::THRESH_BINARY);
	
	//adaptiveThreshold(res,res,255,ADAPTIVE_THRESH_GAUSSIAN_C,THRESH_BINARY_INV,11,7);
	imshow("Frame Difference",res);
	return res;
}

void Scene::findBackground(Mat matDffNow,Mat matDffNowPrev)
{
	Mat res(matDffNow.rows,matDffNow.cols,CV_8UC1);
	res = cv::abs(matDffNow - matDffNowPrev);

	Mat isDiff,thresBack;
	threshold(res,isDiff,1,1,CV_THRESH_BINARY_INV);
	
	//int a=(res==1);
	//Mat grayFrame;
	
	//cvtColor(frameNow,grayFrame,CV_RGB2GRAY);
	
	sumBack=((isDiff==1 & sumBack<80)/255).mul(sumBack+isDiff);
	
	//sumBack=sumBack+isDiff;
	
	background+=((background==0)/255).mul(((sumBack==79)/255.0).mul(grayFrame));
	backgroundMask+=((sumBack==79));
	
		
	sumBack+=(sumBack==79).mul(ones/255.0);


	int stagnateThres=100;
	double maxVal=0,minVal;
	cv::minMaxLoc(sumBack,&minVal,&maxVal);

	imshow("Mask",backgroundMask);

}
Mat Scene::findObject(Mat frame)
{
	Mat grayBack;
	
	grayBack=background.clone();
	//cvtColor(background,grayBack,CV_RGB2GRAY);
	//cvtColor(frame,grayFrame,CV_RGB2GRAY);
	

	Mat backDiff=cv::abs(grayFrame-grayBack);
	//Mat backDiff=cv::abs(frame-grayBack);
	

	threshold(backDiff,backDiff,17, 1, cv::THRESH_BINARY);
	GaussianBlur(backDiff,backDiff,Size(25,25),2,2);
	threshold(backDiff,backDiff,0.5, 1, cv::THRESH_BINARY);
	imshow("Object",backDiff.mul(grayFrame));
	//imshow("Object",backDiff.mul(frame));

	return backDiff;
}
Mat Scene::refineObject(Mat objectIn)
{
	Mat abc;
	int size=32,smallSize=4;
	//blur(objectIn,abc,Size(40,40));
	cv::Scalar avg,avg2;Mat tmpImg;
	refineMask=Mat::zeros(height,width,CV_8UC1);
	for(int i=0;i<height/size;i++)
	{
		for(int j=0;j<width/size;j++)
		{
			cv::Rect roi(j*size,i*size,size,size);
			avg=	cv::mean(objectIn(roi));
			//tmpImg=objectIn(roi).clone();
			if(avg[0]>0.5)
			{	
				//refineMask(roi)=Mat::ones(size,size,CV_8UC1);
				for(int ii=0;ii<size/smallSize;ii++)
				{
					for(int jj=0;jj<size/smallSize;jj++)
					{
						cv::Rect roi2(j*size+jj*smallSize,i*size+ii*smallSize,smallSize,smallSize);
						avg2=cv::mean(objectIn(roi2));
						if(avg2[0]>0.5)
							refineMask(roi2)=Mat::ones(smallSize,smallSize,CV_8UC1);
					}
				}
			}
		}
	}
	imshow("Refine Mask",refineMask.mul(grayFrame));
	return (refineMask);
}

double Scene::sumArr(double D[5000],int low,int upper)
{
    double arr_sum=0;
    for(int i=low;i<=upper;i++)
        arr_sum+=D[i];

    return arr_sum;
}
void Scene::findcontours(Mat ob,int count)
{ 
	Mat drawing=Mat::zeros(height,width,CV_8UC1);
	double largest_area = 0,area=0.0;
	Moments mu;
	Point2f mc;
    double orientation=0.0;
	findContours(ob,contours,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_SIMPLE);
	vector<vector<Point> > contours_poly( contours.size() );
	vector<Rect> boundRect( contours.size() );
	vector<vector<Point>> largest_contours;

	for( int i= 0; i < contours.size(); i++)
	{  // get the largest contour
        area = fabs(contourArea(contours[i]));
        if(area >= largest_area){
            largest_area = area;
            largest_contours.clear(); 
            largest_contours.push_back(contours[i]);
        }
    }
	//Moment calculation
    for( int i = 0; i < largest_contours.size(); i++ )
    { 
		 mu = moments( largest_contours[i], false );
	}
	
    //  Get the mass centers:
    for( int i = 0; i < largest_contours.size(); i++ )
    { 
		mc = Point2f( mu.m10/mu.m00 , mu.m01/mu.m00 ); 
		orientation = 0.5*atan2((2*mu.m11),(mu.m20 - mu.m02))*180/3.14 ; 
		approxPolyDP( Mat(largest_contours[i]), contours_poly[i], 3, true );
        boundRect[i] = boundingRect( Mat(contours_poly[i]) );
		A=Point2f((boundRect[i].tl().x+boundRect[i].br().x)/2,boundRect[i].tl().y);
        B=Point2f((boundRect[i].tl().x+boundRect[i].br().x)/2,boundRect[i].br().y);
		C=Point2f(boundRect[i].tl().x,(boundRect[i].tl().y+boundRect[i].br().y)/2);
		D=Point2f(boundRect[i].br().x,(boundRect[i].tl().y+boundRect[i].br().y)/2);
	}
	 //cout<<"A : "<<Point2f(A.x,A.y)<<endl;
	 //cout<<"B : "<<Point2f(B.x,B.y)<<endl;
	 //cout<<"C : "<<Point2f(C.x,C.y)<<endl;
	 //cout<<"D : "<<Point2f(D.x,D.y)<<endl;
    //Drawing the contour of the largest area
	for(int i=0;i<largest_contours.size();i++)
	{
		if(largest_area>=4500)
		{
			drawContours(drawing,contours_poly,-1,Scalar(255,0,0),1);
		    rectangle( drawing, boundRect[i].tl(), boundRect[i].br(),Scalar(255,0,0), 2, 8, 0 );
		    circle( drawing, mc, 4,Scalar(255,0,0), -1);
		}
	}
	//putText(drawing,"A",A,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
	//putText(drawing,"B",B,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
	//putText(drawing,"C",C,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
	//putText(drawing,"D",D,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
	int locMax_index = distance_curve(drawing,mc,largest_contours,count);
	Point2f N;
	if(locMax_index)
	{
		for(int i=0;i<largest_contours.size();i++)
		{
			if(largest_area>= 4500)
				N=find_tip_head(drawing,boundRect[i].tl(),boundRect[i].br(),mc);
		}
	
	    head_detection(drawing,orientation,locMax_index,mc);
	    feet_detection(drawing,locMax_index,mc);
	    hand_detection(drawing,locMax_index,mc);
		elbow_shoulder_knee_detection(drawing,locMax_index,mc,N,contours_poly,largest_area);
		//elbow_shoulder_detection(drawing,locMax_index,mc,N,largest_contours,largest_area);
	}
	/*if(frameNum==FRAME_COUNT)
		cloth_texture_fitting(contours_poly,drawing);*/
		//cloth_fitting();
	
	imshow("boundary",drawing);
	//if(frameNum==FRAME_COUNT+10)
	//	system("pause");
}
	 //Distance calculation
int Scene :: distance_curve(Mat draw,Point2f mass_centre,vector<vector<Point>> largest_contours1,int count)
{
	 int sizeD=0;
	 double Dist[5000],Dlow[5000],grad[5000],gradLow[5000];
	 int index=0,k=0;
	 static int flag=1; 
	 //ofstream fout("plot.txt");
	 //cout << "Frame Number is : " << frameNum <<endl;
	 for(int i=0;i<largest_contours1.size();i++)
	 {
		 for(int j=0;j<largest_contours1[i].size();j++)
		 {
			 Dist[index]=sqrt((largest_contours1[i][j].x-mass_centre.x)*(largest_contours1[i][j].x-mass_centre.x)+(largest_contours1[i][j].y-mass_centre.y)*(largest_contours1[i][j].y-mass_centre.y));
			 //fout << Dist[index] << endl;
			 PointX[index]=largest_contours1[i][j].x;
			 PointY[index]=largest_contours1[i][j].y;
			 index++;
			 
		 }
	 }
	 
	 //low pass filtering of distance array
	 int len=13;
	 MaxPoints= index;
	 for(int i=len+1; i<index-len;i++)
	 {
        Dlow[i]=sumArr(Dist,i-len,i+len) /(2*len+1);
		
	 }
	 for(int i=len+2;i<index-len;i++)
         grad[i]=Dlow[i]-Dlow[i-1];
     //fout.close();
     for(int i=2*len+1;i<index-2*len;i++)
     {
		 gradLow[i]= sumArr(grad,i-len,i+len)/(2*len+1);
	 }
     int diff=index/25;
	 //finding local maxima
     for(int i=1;i<index;i++)
     {
        if((gradLow[i]/gradLow[i-1])<0 && i>diff+1)
        {
            int flag=0;
            for(int j=i-diff;j<i-1;j++)
                if(gradLow[j]<0)
                    flag=1;

            for(int j=i+1;j<i+diff;j++)
                if(gradLow[j]>0)
                    flag=1;
            if(flag==0)
            {    
				locMax[k++]=i;
				
			}
        }
    }
	//if(frameNum==965)
	//	 system("pause");
	 // drawing the convex points(local maxima)
	//cout <<"----------------------------------------------------" << endl;
	for(int i=0;i<k;i++)
	{
		////cout << "K " << k << endl;
		//if(count>=400)
		//{
		//	flag= kalmanPredict(count,flag,i,k);
		//	kalmanUpdate(count,PointX[locMax[i]],PointY[locMax[i]],i,k);
		//	//cout << "---------------------------------------------------------" << endl;
		//	//cout << "count " << count << endl;
		//	//cout << "Original -> " << PointX[locMax[i]] << " Estimated point -> " <<X[count] << endl;
		//	//cout << "Original -> " << PointY[locMax[i]] << " Estimated point -> " <<Y[count] << endl;
		//	//cout << "----------------------------END--------------------------" << endl;
		//	//PointX[locMax[i]] = X[i][count];
		//	//PointY[locMax[i]] = Y[i][count];
		//}
		//
		
		//cout <<"ith point -> " << i << "    X co-ordinate  " <<PointX[locMax[i]] <<"    Y co-ordiante "<<PointY[locMax[i]] <<endl;
		circle(draw,Point2f(PointX[locMax[i]],PointY[locMax[i]]),4,Scalar(255,0,0), -1);
		locMaxFlag[i]=0;
	}
	
	//cout <<"------------------------END------------------------" << endl;
	
	//finding curvature of convex points
    vector<Point2f> U1(k);
    vector<Point2f> U2(k);
	String ch;
	//cout<<"-----------start------------"<<endl;
	for(int i=0;i<k;i++)
	{
		U1[i].x=PointX[locMax[i]+20]-PointX[locMax[i]];
		U1[i].y=PointY[locMax[i]+20]-PointY[locMax[i]];
		U2[i].x=PointX[locMax[i]-20]-PointX[locMax[i]];
		U2[i].y=PointY[locMax[i]-20]-PointY[locMax[i]];
		theta[i]=acos(((U1[i].x * U2[i].x)+(U1[i].y * U2[i].y)) / ((sqrt((U1[i].x*U1[i].x) + (U1[i].y*U1[i].y))) * (sqrt((U2[i].x*U2[i].x) + (U2[i].y*U2[i].y))))) * 180.0/3.14;
		/*cout<<"Curvature: "<<theta[i]<<",";
		ch='P'+(char)i;
		putText(draw,ch,Point2f(PointX[locMax[i]],PointY[locMax[i]]),FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);*/
	}
	//cout<<endl<<"---------------end------------"<<endl;

	return k;
}
Point Scene :: find_tip_head(Mat draw,Point top,Point bottom,Point2f mc)
{
	int flag=0;
	int p;
	for (int i=0;i<bottom.x-top.x;i++)
	{
		for(int j=0;j<MaxPoints;j++)
		{
			if(PointX[j]==top.x+i && PointY[j]==top.y)
			{	
				p=top.x+i;
				flag=1;
				break;
			}
		}
	}
	if(flag)
	{
		circle( draw, Point(p,top.y), 4,Scalar(255,0,0), -1);
		putText(draw,"Tip",Point(p,top.y),FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		circle(draw,Point((2*p+mc.x)/3,(2*top.y+mc.y)/3),4,Scalar(255,0,0),1);
		putText(draw,"N",Point((2*p+mc.x)/3,(2*top.y+mc.y)/3),FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		Tip=Point(p,top.y);
		Neck=Point((2*p+mc.x)/3,(2*top.y+mc.y)/3);
		return Point((2*p+mc.x)/3,(2*top.y+mc.y)/3);
	}
	return 0;
}
Point Scene :: find_elbow_knee(Point2f TempElbow,Point2f poi,Point2f stPoint,vector<vector<Point>> contours_poly)
{
	int match=0;
	int x,y;
	Point C1,C2;
	Point Temp1,Temp2;
	//stPoint=Point(stPoint.x,stPoint.y);
	/*for(int i=0;i<MaxPoints;i++)
	{
		if(TempElbow.x == PointX[i] && TempElbow.y==PointY[i])
		{
			return TempElbow;
	    }
	}
	for(int i=0;i<MaxPoints;i++)
	{
		C1.x+=1;C1.y+=1;
		if(PointX[i]==C1.x && PointY[i]==C1.y)
		    break;
	}
	for(int i=0;i<MaxPoints;i++)
	{
		C2.x-=1;C2.y-=1;
		if(PointX[i]==C2.x && PointY[i]==C2.y)
		    break;
	}*/
	if(pointPolygonTest(contours_poly[0],TempElbow,false)==0)
	    //match=1;
		return TempElbow;
	//Check whether point is inside the contour
	else if(pointPolygonTest(contours_poly[0],TempElbow,false)>0)
	{ 
		if((stPoint.y-poi.y)/(stPoint.x-poi.x)<0)
		{
			C1.x=TempElbow.x+1;C1.y=TempElbow.y+1;
			C2.x=TempElbow.x-1;C2.y=TempElbow.y-1;
	
			while(pointPolygonTest(contours_poly[0],C1,false)>0)
			{
				C1.x+=1;
				C1.y+=1;
			}
			while(pointPolygonTest(contours_poly[0],C2,false)>0 && C2.x >=0 && C2.y>=0)
			{	
				C2.x-=1;
				C2.y-=1;
			}
		}
		else if((stPoint.y-poi.y)/(stPoint.x-poi.x)==0)
		{
			C1.x=TempElbow.x+1;C1.y=TempElbow.y;
			C2.x=TempElbow.x-1;C2.y=TempElbow.y;
	
			while(pointPolygonTest(contours_poly[0],C1,false)>0)
			{
				C1.x+=1;
			}
			while(pointPolygonTest(contours_poly[0],C2,false)>0 && C2.x >=0)
			{	
				C2.x-=1;
			}
		}
		else
		{
			C1.x=TempElbow.x-1;C1.y=TempElbow.y+1;
			C2.x=TempElbow.x+1;C2.y=TempElbow.y-1;
			while(pointPolygonTest(contours_poly[0],C1,false)>0 && C1.x >=0)
			{
				C1.x-=1;
				C1.y+=1;
			}
			while(pointPolygonTest(contours_poly[0],C2,false)>0 && C2.y >=0)
			{
				C2.x+=1;
				C2.y-=1;
			}
		}
	}
	//Check whether the point is outside the contour
	else if(pointPolygonTest(contours_poly[0],TempElbow,false) <0)
	{
		if((stPoint.y-poi.y)/(stPoint.x-poi.x)<0)
		{
			//double dis1=sqrt(((stPoint.x-TempElbow.x+1)*(stPoint.x-TempElbow.x+1))+((stPoint.y-TempElbow.y+1)*(stPoint.y-TempElbow.y+1)));
			//double dis2=sqrt(((stPoint.x-TempElbow.x-1)*(stPoint.x-TempElbow.x-1))+((stPoint.y-TempElbow.y-1)*(stPoint.y-TempElbow.y-1)));
			int count1=0,count2=0;
			Temp1.x=TempElbow.x+1;Temp1.y=TempElbow.y+1;
		    Temp2.x=TempElbow.x-1;Temp2.y=TempElbow.y-1;
			while(pointPolygonTest(contours_poly[0],Temp1,false) <0 )
			{
				Temp1.x+=1;Temp1.y+=1;
				count1++;
			}
			while(pointPolygonTest(contours_poly[0],Temp2,false) <0 && Temp2.x>=0 && Temp2.y>=0)
			{
				Temp2.x-=1;Temp2.y-=1;
				count2++;
			}
			
			if(count1<count2)
			{
				C2.x=TempElbow.x+1;C2.y=TempElbow.y+1;
				while(pointPolygonTest(contours_poly[0],C2,false)<0)
				{
					C2.x+=1;
					C2.y+=1;
				}
				C1.x=C2.x+1;
				C1.y=C2.y+1;
				while(pointPolygonTest(contours_poly[0],C1,false)>0)
				{	
					C1.x+=1;
					C1.y+=1;
				}
			}
			else
			{
				C2.x=TempElbow.x-1;C2.y=TempElbow.y-1;
				while(pointPolygonTest(contours_poly[0],C2,false)<0 && C2.x>=0 && C2.y>=0)
				{
					C2.x-=1;
					C2.y-=1;
				}
				C1.x=C2.x-1;
				C1.y=C2.y-1;
				while(pointPolygonTest(contours_poly[0],C1,false)>0 && C1.x>=0 && C1.y>=0)
				{	
					C1.x-=1;
					C1.y-=1;
				}
			}
		}
		else 
		{
			//double dis1=sqrt((stPoint.x-TempElbow.x-1)*(stPoint.x-TempElbow.x-1)+(stPoint.y-TempElbow.y+1)*(stPoint.y-TempElbow.y+1));
			//double dis2=sqrt((stPoint.x-TempElbow.x+1)*(stPoint.x-TempElbow.x+1)+(stPoint.y-TempElbow.y-1)*(stPoint.y-TempElbow.y-1));
			//double dis1=sqrt((poi.x-TempElbow.x-1)*(poi.x-TempElbow.x-1)+(poi.y-TempElbow.y+1)*(poi.y-TempElbow.y+1));
			//double dis2=sqrt((poi.x-TempElbow.x+1)*(poi.x-TempElbow.x+1)+(poi.y-TempElbow.y-1)*(poi.y-TempElbow.y-1));
			int count1=0,count2=0;
			Temp1.x=TempElbow.x-1;Temp1.y=TempElbow.y+1;
		    Temp2.x=TempElbow.x+1;Temp2.y=TempElbow.y-1;
			while(pointPolygonTest(contours_poly[0],Temp1,false) <0 && Temp1.x>=0)
			{
				Temp1.x-=1;Temp1.y+=1;
				count1++;
			}
			while(pointPolygonTest(contours_poly[0],Temp2,false) <0 && Temp2.y>=0)
			{
				Temp2.x+=1;Temp2.y-=1;
				count2++;
			}
			/*if(frameNum==852)
			{
				cout<<"--------------" << endl;
				cout<< "Dis1  " << dis1 << endl;
				cout<< "Dis2  " << dis2 << endl;
			}*/
			/*if(frameNum==460)
			{
				cout << "count1 : " << count1<<endl;
				cout <<"count2 : "<< count2 << endl;
			}*/
			if(count1<count2)
			{
				//cout << "In Else part " << endl;
				C2.x=TempElbow.x-1;
				C2.y=TempElbow.y+1;
				//cout << "C2.x  " <<  C2.x << "C2.y  " << C2.y << endl; 
				while(pointPolygonTest(contours_poly[0],C2,false) < 0 && C2.x >=0)
				{
					C2.x-=1;
					C2.y+=1;
				}
				//cout << "After the while " << endl;
				C1.x=C2.x-1;
				C1.y=C2.y+1;
				while(pointPolygonTest(contours_poly[0],C1,false)> 0 && C1.x >=0 )
				{
					C1.x-=1;
					C1.y+=1;
				}
			}
			else
			{
				//cout << "In Else part " << endl;
				C2.x=TempElbow.x+1;
				C2.y=TempElbow.y-1;
				//cout << "C2.x  " <<  C2.x << "C2.y  " << C2.y << endl; 
				while(pointPolygonTest(contours_poly[0],C2,false) < 0 && C2.y >=0)
				{
					C2.x+=1;
					C2.y-=1;
				}
				//cout << "After the while " << endl;
				C1.x=C2.x+1;
				C1.y=C2.y-1;
				while(pointPolygonTest(contours_poly[0],C1,false)> 0 && C1.y >=0 )
				{
					C1.x+=1;
					C1.y-=1;
				}
			}
		}
	}
	/*cout <<"Point "<< poi << endl;
	cout << "Mass Centre " <<stPoint << endl;
	cout << "TempElbow" <<TempElbow <<endl;
	cout << "C1 is : " <<C1 <<"  C2 is : "<<C2 << endl;*/
	float a1= atan2(C1.y-poi.y,C1.x-poi.x)*180/3.14;
	float a2= atan2(C2.y-poi.y,C2.x-poi.x)*180/3.14;
	/*float a1= atan2(poi.y-C1.y,poi.x-C1.x)*180/3.14;
	float a2= atan2(poi.y-C2.y,poi.x-C2.x)*180/3.14;
	*/
	E1=C1;
	E2=C2;
	//cout << "Angle1:  "<< a1 << "  Point:  " << C1 << endl;
	//cout << "Angle2:  "<< a2 << "  Point:  " << C2 << endl;

    if(a1<a2)
		return C1;
	else return C2;

	//return TempElbow;
}

void Scene :: elbow_shoulder_knee_detection(Mat draw, int s, Point2f mass_centre,Point2f neck,vector<vector<Point>> contours_poly,double largest_area)
{
	Point2f TempElbow,shoulder,knee,elbow,TempKnee; 
	int ElbowCount=0;
	for(int i=0;i<s;i++)
	{
		if(locMaxFlag[i]=='l')
		{
			TempElbow=Point((PointX[locMax[i]]+neck.x)/2,(PointY[locMax[i]]+neck.y)/2);		
			elbow=find_elbow_knee(TempElbow,Point(PointX[locMax[i]],PointY[locMax[i]]),neck,contours_poly);
			
			//shoulder=Point((neck.x+elbow.x)/2,(neck.y+elbow.y)/2);
			circle(draw,elbow,5,Scalar(255,0,0), -1);
			putText(draw,"E",elbow,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
			
			
			ElbowCount++;
			/*if(frameNum==FRAME_COUNT+20)
			{
				Mat cloth = imread("redTshirt.jpg",CV_LOAD_IMAGE_UNCHANGED);
				imwrite( "cloth_Image.jpg", cloth );
				elbow2 = Point(elbow.x,elbow.y+3);
				if(pointPolygonTest(contours_poly[0],elbow2,false)>=0)
				{
					cout << "In the If condition " << endl;
					while(pointPolygonTest(contours_poly[0],elbow2,false)>0)
						elbow2.y+=1;
					E2=elbow2;
					E1=elbow;
				}
				else
				{
					cout << "Elbow2" << elbow2 << endl;
					cout << "Elbow" << elbow << endl;
					cout << "In the Else condition " << endl;
					elbow2=Point(elbow.x,elbow.y-2);
					
					while(pointPolygonTest(contours_poly[0],elbow2,false)>=0)
					{	
						elbow2.y-=1;
						//cout << "In while" << endl;
					}
					E1=elbow2;
					E2=elbow;
				}

				circle(draw,E1,5,Scalar(255,0,0), 1);
			    putText(draw,"E1",E1,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				
				circle(draw,E2,5,Scalar(255,0,0), 1);
			    putText(draw,"E2",E2,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);


				// RECTANGLE FOR BODY PART
				R1_height = cal_distance(Neck,mass_centre);
				R1_width=0;
				Point temp= mass_centre;
				while(pointPolygonTest(contours_poly[0],temp,false)>0)
				{
					temp.x+=1;
					R1_width++;
				}
				temp=mass_centre;
				while(pointPolygonTest(contours_poly[0],temp,false)>0)
				{
					temp.x-=1;
					R1_width++;
				}
				int tshirt_width=200;   // use functions
				int tshirt_height=200;  // use functions
			
				resize(cloth,cloth,Size((tshirt_height*R1_height)/173,(tshirt_width*R1_width)/114));
				
				cout << R1_height << endl;
				cout << R1_width << endl;
                cout << tshirt_height << endl;
				cout << tshirt_width << endl;



				cout << (tshirt_height*R1_height)/173 << endl;
				cout << (tshirt_width*R1_width)/114 << endl;
				
				Mat channel[3];
				Mat shirt;
	            //namedWindow("Cloth",CV_WINDOW_AUTOSIZE);
	            split(cloth,channel);


				for(int i=0;i<channel[2].rows;i++)
				{
					for(int j=0;j<channel[2].cols;j++)
					{
						if(channel[2].at<uchar>(i,j)==255)
						    channel[2].at<uchar>(i,j)=0;
					}
				}


				findNonZero(channel[2],shirt);
			
				Point reference= Point(102,18);
				/*for(int i=0;i<shirt.rows;i++)
				{
					for(int j=0;j<shirt.cols;j++)
					{
						cout << shirt.at<uchar>(i,j) << endl;
					}
				}*/
				/*for(int i=0;i<shirt.rows;i++)
				{
					shirt.at<Point>(i)=shirt.at<Point>(i)-reference+Neck;
				}


				for(int i=0;i<shirt.rows;i++)
				{
					if(pointPolygonTest(contours_poly[0],shirt.at<Point>(i),false)>=0)
					{
						//frameNow.at<Vec3b>(shirt.at<Point>(i))[0]=channel[0].at<uchar>(shirt.at<Point>(i));
						//frameNow.at<Vec3b>(shirt.at<Point>(i))[1]=channel[1].at<uchar>(shirt.at<Point>(i));
						//frameNow.at<Vec3b>(shirt.at<Point>(i))[2]=channel[2].at<uchar>(shirt.at<Point>(i));
					    frameNow.at<Vec3b>(shirt.at<Point>(i))[2]=255;
					}
					//channel[2].at<uchar>(shirt.at<Point>(i))=255;
					//draw.at<uchar>(shirt.at<Point>(i))=255;
				}

			}*/
			if(frameNum>FRAME_COUNT-10)
			{
				//circle(draw,E1,5,Scalar(255,0,0), 1);
			    //putText(draw,"E1",E1,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				
				//circle(draw,E2,5,Scalar(255,0,0), 1);
			    //putText(draw,"E2",E2,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);


				
				// Find one another ebow point on the other side of the contour

				
				elbow2 = Point(elbow.x,elbow.y+3);
				if(pointPolygonTest(contours_poly[0],elbow2,false)>=0)
				{
					//cout << "In the If condition " << endl;
					while(pointPolygonTest(contours_poly[0],elbow2,false)>0)
						elbow2.y+=1;
					E2=elbow2;
					E1=elbow;
				}
				else
				{
					//cout << "Elbow2" << elbow2 << endl;
					//cout << "Elbow" << elbow << endl;
					//cout << "In the Else condition " << endl;
					elbow2=Point(elbow.x,elbow.y-2);
					
					while(pointPolygonTest(contours_poly[0],elbow2,false)>=0)
					{	
						elbow2.y-=1;
						//cout << "In while" << endl;
					}
					E1=elbow2;
					E2=elbow;
				}

				circle(draw,E1,5,Scalar(255,0,0), 1);
			    putText(draw,"E1",E1,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				
				circle(draw,E2,5,Scalar(255,0,0), 1);
			    putText(draw,"E2",E2,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);


				// RECTANGLE FOR BODY PART
				R1_height = cal_distance(Neck,mass_centre);
				R1_width=0;
				Point temp= mass_centre;
				while(pointPolygonTest(contours_poly[0],temp,false)>0)
				{
					temp.x+=1;
					R1_width++;
				}
				temp=mass_centre;
				while(pointPolygonTest(contours_poly[0],temp,false)>0)
				{
					temp.x-=1;
					R1_width++;
				}
				R1_centre= Point2f ((mass_centre.x+Neck.x)/2,(mass_centre.y+Neck.y)/2);
				R1= RotatedRect(R1_centre,Size(R1_height,R1_width),90);
				
				R1.points(R1_vertices2f);
             
				for(int i = 0; i < 4; i++)
					R1_vertices[i] = R1_vertices2f[i];
			
			
				//fillConvexPoly(frameNow,R1_vertices,4,Scalar(0,0,255));


				// RECTANGLE FOR LIMB
				R2_height = cal_distance(E1,E2);
				R2_width= cal_distance(E1,Neck);
				E_mid=Point2f((E1.x+E2.x)/2,(E1.y+E2.y)/2);
				//Point2f Neck_extend = Point2f(Neck.x+R2_height*.5,Neck.y+R2_height*.5);
				
				float DeltaY=E1.y-Neck.y; 
				float DeltaX=Neck.x-E1.x;
				//float DeltaY= mass_centre.y-PointY[locMax[i]];
				//float DeltaX= mass_centre.x-PointX[locMax[i]];
				angle= atan2(DeltaY,DeltaX)*180/3.14;
				//cout <<"Angle: " <<angle << endl;
				if(mass_centre.x<E1.x)
					angle=95+angle;
				else 
					angle= 90-angle;


				R2_centre = Point2f((E_mid.x+Neck.x)/2,(E_mid.y+Neck.y+(E_mid.y-E1.y))/2);
                //Point2f R2_centre = Point2f(E_mid.x+width*.5,E_mid.y);
				R2= RotatedRect(R2_centre,Size(R2_height,R2_width),angle);
				R2.points(R2_vertices2f);
            
				for(int i = 0; i < 4; i++)
					R2_vertices[i] = R2_vertices2f[i];
			
			
				//fillConvexPoly(frameNow,R2_vertices,4,Scalar(0,0,255));

				cloth_texture_fitting(contours_poly,draw,mass_centre);                 

			}
			/*else if(frameNum>FRAME_COUNT-10)
			{

				//RECTANGLE FOR BODY
				Point2f Temp_centre= Point2f ((mass_centre.x+Neck.x)/2,(mass_centre.y+Neck.y)/2);
				RotatedRect Temp_Rect= RotatedRect(Temp_centre,Size(R1_height,R1_width),90);
				Point2f Temp_Vertices2f[4];
				Point Temp_Vertices[4];
				Temp_Rect.points(Temp_Vertices2f);
            
				for(int i = 0; i < 4; i++)
					Temp_Vertices[i] = Temp_Vertices2f[i];
			    //fillConvexPoly(frameNow,Temp_Vertices,4,Scalar(0,0,255));

				//RECTANGLE FOR LIMB

				float DeltaY=elbow.y-Neck.y; 
				float DeltaX=Neck.x-elbow.x;
				//float DeltaY= mass_centre.y-PointY[locMax[i]];
				//float DeltaX= mass_centre.x-PointX[locMax[i]];
				angle= atan2(DeltaY,DeltaX)*180/3.14;
				cout <<"Angle: " <<angle << endl;
				if(mass_centre.x<elbow.x)
					angle=95+angle;
				else 
					angle= 90-angle;


				Point2f Temp2_centre = Point2f((elbow.x+Neck.x)/2,(elbow.y+Neck.y)/2);
                //Point2f R2_centre = Point2f(E_mid.x+width*.5,E_mid.y);
				RotatedRect Temp2_Rect= RotatedRect(Temp2_centre,Size(R2_height,R2_width),angle);
				Point2f Temp2_Vertices2f[4];
				Point Temp2_Vertices[4];
				Temp2_Rect.points(Temp2_Vertices2f);
				           
				for(int i = 0; i < 4; i++)
					Temp2_Vertices[i] = Temp2_Vertices2f[i];
			
			
				//fillConvexPoly(frameNow,Temp2_Vertices,4,Scalar(0,0,255));


			}*/
			/*else if(frameNum==634 && frameNum <637)
			{
				float DeltaY=mass_centre.y-E2.y; 
				float DeltaX=mass_centre.x-E2.x;
				float angle = atan2(DeltaY,DeltaX)*180/3.14;

				RotatedRect R2= RotatedRect(R2_centre,Size(R2_height,R2_width),angle);
			}*/

		}


		if(locMaxFlag[i]=='f')
		{
			TempKnee= Point((PointX[locMax[i]]+mass_centre.x)/2,(PointY[locMax[i]]+mass_centre.y)/2);
			knee=find_elbow_knee(TempKnee,Point(PointX[locMax[i]],PointY[locMax[i]]),mass_centre,contours_poly);
			
			circle(draw,knee,5,Scalar(255,0,0), 1);
			putText(draw,"K",knee,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		}
		

	}
	if(ElbowCount==0)
	{ 
		TempElbow=Point((C.x+neck.x)/2,(C.y+neck.y)/2);
		elbow=find_elbow_knee(TempElbow,C,neck,contours_poly);
		circle(draw,elbow,5,Scalar(255,0,0),-1);
		putText(draw,"E",elbow,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		



		TempElbow= Point((D.x+neck.x)/2,(D.y+neck.y)/2);
		elbow=find_elbow_knee(TempElbow,D,neck,contours_poly);
		circle(draw,elbow,5,Scalar(255,0,0),-1);
		putText(draw,"E",elbow,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		/*if(frameNum==632)
		{
			float w= cal_distance(E1,E2);
			float h = cal_distance(Neck,E1);
			//float angle= 
			//rectangle(draw,E2,Neck,Scalar(255,0,0),1,8,0);
			Point2f centre= Point2f ((E2.x+Neck.x)/2,(E2.y+Neck.y)/2);
			//Point2f centre= Point2f((neck.x+Neck.x+10)/2,(Neck.y+Neck.y+10)/2);
			RotatedRect rrect= RotatedRect(centre,Size(w,30),145);
			Point2f vertices2f[4];Point vertices[4];
            rrect.points(vertices2f);
            
			for(int i = 0; i < 4; i++){
                vertices[i] = vertices2f[i];
            }
			
			
			fillConvexPoly(draw,vertices,4,Scalar(255,0,0));
		}*/
		
		    if(frameNum>FRAME_COUNT)
			{
		    
				for(int indexNo=0;indexNo<2;indexNo++)
				{   
					if(indexNo==0) 
					{
						TempElbow=Point((C.x+neck.x)/2,(C.y+neck.y)/2);
		                elbow=find_elbow_knee(TempElbow,C,neck,contours_poly);
					}
					else
					{
						TempElbow=Point((D.x+neck.x)/2,(D.y+neck.y)/2);
						elbow=find_elbow_knee(TempElbow,D,neck,contours_poly);
					}
					//RECTANGLE FOR BODY
					
					/*Point2f Temp_centre= Point2f ((mass_centre.x+Neck.x)/2,(mass_centre.y+Neck.y)/2);
					RotatedRect Temp_Rect= RotatedRect(Temp_centre,Size(R1_height,R1_width),90);
					Point2f Temp_Vertices2f[4];
					Point Temp_Vertices[4];
					Temp_Rect.points(Temp_Vertices2f);
            
					for(int i = 0; i < 4; i++)
						Temp_Vertices[i] = Temp_Vertices2f[i];
					fillConvexPoly(frameNow,Temp_Vertices,4,Scalar(0,0,255));*/

					float Temp_height = cal_distance(Neck,mass_centre);
					float Temp_width=0;
					Point temp= mass_centre;
					while(pointPolygonTest(contours_poly[0],temp,false)>0)
					{
						temp.x+=1;
						Temp_width++;
					}
					temp=mass_centre;
					while(pointPolygonTest(contours_poly[0],temp,false)>0)
					{
						temp.x-=1;
						Temp_width++;
					}
					Point2f Temp_centre= Point2f ((mass_centre.x+Neck.x)/2,(mass_centre.y+Neck.y)/2);
					RotatedRect Temp_R= RotatedRect(Temp_centre,Size(Temp_height,Temp_width),90);

                 	Point2f Temp_vertices2f[4];
					Point Temp_vertices[4];

					Temp_R.points(Temp_vertices2f);
            
					for(int i = 0; i < 4; i++)
						Temp_vertices[i] = Temp_vertices2f[i];
			
			
					//fillConvexPoly(frameNow,Temp_vertices,4,Scalar(0,0,255));
	
					
					//RECTANGLE FOR LIMB

					float DeltaY=elbow.y-Neck.y; 
					float DeltaX=Neck.x-elbow.x;
					angle= atan2(DeltaY,DeltaX)*180/3.14;
					//cout <<"Angle: " <<angle << endl;
					//if(mass_centre.x<elbow.x)
					//	angle=95+angle;
					//else 
					//	angle= 90-angle;

					//if(mass_centre.x < elbow.x)
					//	angle=angle-90;
					//else
					//	angle=angle-180;

					Point2f Temp2_centre = Point2f((elbow.x+Neck.x)/2,(elbow.y+Neck.y)/2);
					//Point2f R2_centre = Point2f(E_mid.x+width*.5,E_mid.y);
					RotatedRect Temp2_Rect= RotatedRect(Temp2_centre,Size(R2_height,R2_width),angle);
					Point2f Temp2_Vertices2f[4];
					Point Temp2_Vertices[4];
					Temp2_Rect.points(Temp2_Vertices2f);
				           
					for(int i = 0; i < 4; i++)
						Temp2_Vertices[i] = Temp2_Vertices2f[i];
			
			
					//fillConvexPoly(frameNow,Temp2_Vertices,4,Scalar(0,0,255));
				}
			}
	}
}
void Scene :: head_detection(Mat draw,double orient,int s,Point2f mass_centre)
{
	//static ofstream fout("Head.txt"); 
	int poi_index=0;
	vector<Point2f> POI(10);
	//finding points of interest in the upper contour for head location
	//cout<<"Head orientation: "<<orient<<endl;
	//cout << "Frame Number is : " << frameNum << endl;
    if ( frameNum > 200)
	{
		//if(frameNum==441)
			//fout << "C.x       C.y        D.x       D.y      Points" <<endl;
		//else
		//{
	//	fout << C.x <<"     " << C.y << "    " << D.x << "    " << D.y <<"     ";
	
	if(orient>=10 && orient<=180)
	{
	    for(int i=0;i<s;i++)
		{
			if(PointY[locMax[i]]<C.y && PointY[locMax[i]]>=A.y && PointX[locMax[i]]>=(C.x+(D.x-C.x)*.29) && PointX[locMax[i]]<=(C.x+(D.x-C.x)*0.83))// && (C.x > 20))
			{	
				POI[i] = Point2f(PointX[locMax[i]],PointY[locMax[i]]);
				poi_index++;

			}
			else
			{
				//fout <<frameNum<<"    "<< C.x <<"     " << C.y << "    " << A.y<< "    "<<D.x << "    " << D.y <<"     " << PointX[locMax[i]]<<","<<PointY[locMax[i]] << endl;  //<<POI[i].x <<","<<POI[i].y<<endl;
				POI[i]=Point2f( 0.0,0.0);
				
			}
		}
		// in case of more than one potential head points, it calls the find_head() function
		if(poi_index>1)
		{
			int head_index=find_head(POI,theta,poi_index);
			putText(draw,"Head",POI[head_index],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		    locMaxFlag[head_index]='h';
			//line(draw,POI[head_index],mass_centre,Scalar(255,0,0),5,8,0);
			
		}
		else if(poi_index==1)
		{	for(int i=0;i<s;i++)
		    {
				if(POI[i]!=Point2f(0.0,0.0))
				{
					putText(draw,"Head",POI[i],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
					locMaxFlag[i]='h';
					//line(draw,POI[i],mass_centre,Scalar(255,0,0),5,8,0);
					//cout<<" Head Orientation: "<<orient<<endl;
				}
		    }
		}
		else
		{
			putText(draw,"Head",A,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
			//line(draw,A,mass_centre,Scalar(255,0,0),5,8,0);
		}
	}
	}
}

void Scene :: feet_detection(Mat draw,int s,Point2f mass_centre)
{
	int feet_index=0,feet_count=0,index_feet=0;
	double feet_distance[10];
	vector<Point2f> POI_Feet(10);
	//finding point of interest in the lower contour

	for(int i=0;i<s;i++)
	{
		if(PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x && PointY[locMax[i]]>(C.y * 1.25) && PointY[locMax[i]]<=B.y )
		{
			/*if(theta[i]>20 && theta[i]<=100)
			{
				POI_Feet[i]= Point2f(PointX[locMax[i]],PointY[locMax[i]]);
				feet_count++;
			}
		    else
				POI_Feet[i]=Point2f(0.0,0.0);*/
			POI_Feet[i]= Point2f (PointX[locMax[i]],PointY[locMax[i]]);
			feet_count++;
		}
		else
			POI_Feet[i] = Point2f(0.0,0.0);
		
	    /*for(int j=0;j<feet_index;j++)
		{
			if(theta[j] >= 100.0)
				POI_Feet[j]= Point2f(0.0,0.0);
			else
				feet_count++;
		}*/
	}
	if(feet_count==1)
	{
		for(int i=0;i<POI_Feet.size();i++)
		{
			if(POI_Feet[i]!=Point2f(0.0,0.0) && locMaxFlag[i]!='h')
		    {
				putText(draw,"Feet",POI_Feet[i],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				locMaxFlag[i]='f';
				//line(draw,POI_Feet[i],mass_centre,Scalar(255,0,0),8,8,0);
			}
		}
			//cout<<"Feet Index: "<<feet_index-1<<endl;
		//cout<< "Curvature of feet with one point " <<theta[feet_index-1]<<endl;
		//cout<<"----------------------------------------------"<<endl;
	}

	else if(feet_count==2)
	{
		//cout<<"---------------START----------------------------"<<endl;
		//cout<<"Feet Index "<<feet_index<<endl;
		for(int j=0;j<POI_Feet.size();j++)
		{
			if(POI_Feet[j]!=Point2f(0.0,0.0) && locMaxFlag[j]!='h')
			{
				putText(draw,"Feet",POI_Feet[j],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				//line(draw,POI_Feet[j],mass_centre,Scalar(255,0,0),8,8,0);
				locMaxFlag[j]='f';
				//cout<< "Curvature of feet with two points " <<theta[j]<<endl;
			}
			
		}
		//cout<<"---------------END----------------------------"<<endl;
	}
	else if(feet_count>2)
	{
		for(int j=0;j<POI_Feet.size();j++)
		{
			if(POI_Feet[j]!=Point2f(0.0,0.0) && locMaxFlag[j]!='h')
				feet_distance[j]=cal_distance(A,POI_Feet[j]);
		}
		findmax2(draw,feet_distance,POI_Feet,POI_Feet.size(),"Feet",mass_centre);
	}
}

void Scene :: hand_detection(Mat draw,int s,Point2f mass_centre)
{
    //static ofstream fout("hand.txt");
	int hand_index=0,hand_count=0,index_hand=0;
	//int handPoints[10]={0};
	double hand_distance[10];
	vector<Point2f> POI_Hand(10);
	//finding points of interest for hand location
	//cout << "Frame Number is : " << frameNum << endl;
	/*if(frameNum ==580 )
		system("pause");*/
	for(int i=0;i<s;i++)
	{
		if((locMaxFlag[i]!='h') && (PointX[locMax[i]]>C.x && PointX[locMax[i]]<D.x && PointY[locMax[i]]>(A.y*1.1) && PointY[locMax[i]]< (A.y+(B.y-A.y)*0.80)))
		{
			POI_Hand[i]=Point2f(PointX[locMax[i]],PointY[locMax[i]]);
			hand_count++;
		}
		else
		{
		    //fout <<frameNum<<"    "<< C.x <<"     " << D.x << "    " << A.y << "    " << B.y <<"     " << PointX[locMax[i]]<<","<<PointY[locMax[i]] << endl;
			POI_Hand[i]=Point2f(0.0,0.0);
		}
	}
	/*for(int i=0;i<s;i++)
	{
		if(POI_Hand[i]!=Point2f(0.0,0.0) && locMaxFlag[i]!='h')
			fout<<"Theta : "<<theta[i]<<endl;
	}*/
	if(hand_count==2 || hand_count==1)
	{
		for(int i=0;i<POI_Hand.size();i++)
		{
			if(POI_Hand[i]!=Point2f(0.0,0.0) && (locMaxFlag[i]!='h') )
			{
				//POI_Hand[0] = Point2f(PointX[locMax[handPoints[0]]],PointY[locMax[handPoints[0]]]);
				putText(draw,"Limb",POI_Hand[i],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
				locMaxFlag[i]='l';
				//line(draw,POI_Hand[i],mass_centre,Scalar(255,0,0),6,8,0);
				//POI_Hand[1] = Point2f(PointX[locMax[handPoints[1]]],PointY[locMax[handPoints[1]]]);
				//putText(draw,"Hand",POI_Hand[1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
			}
		}
	}
	else if(hand_count>2)
	{
		for(int j=0;j<POI_Hand.size();j++)
		{ 
			if(POI_Hand[j]!=Point2f(0.0,0.0)&& (locMaxFlag[j]!='h'))
			{
				//POI_Hand[j] = Point2f(PointX[locMax[handPoints[j]]],PointY[locMax[handPoints[j]]]);
				hand_distance[j]=cal_distance(mass_centre,POI_Hand[j]);
			}
		}
		findmax2(draw,hand_distance,POI_Hand,POI_Hand.size(),"Limb",mass_centre);
	}
	else 
	{
		putText(draw,"Limb",C,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		putText(draw,"Limb",D,FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		//line(draw,C,mass_centre,Scalar(255,0,0),6,8,0);
		//line(draw,D,mass_centre,Scalar(255,0,0),6,8,0);
	}
}

//finds the maximum of two distances
void Scene:: findmax2(Mat draw1,double arr[10],vector<Point2f> P,int index,char *str,Point2f mass_centre)
{
	double temp1=arr[0],temp2;
	int max1=0,max2=0;

	for(int i=1;i<index;i++)
	{
		if(arr[i]>temp1 && P[i]!=Point2f(0.0,0.0) && locMaxFlag[i]!='h')
		{
			temp2=temp1;
			temp1=arr[i];
			max2=max1;
			max1=i;
		}
	}
	if (str == "Feet")
	{
		putText(draw1,str,P[max1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		putText(draw1,str,P[max2],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		locMaxFlag[max1]='f';
		locMaxFlag[max2]='f';
		//line(draw1,P[max1],mass_centre,Scalar(255,0,0),8,8,0);
		//line(draw1,P[max2],mass_centre,Scalar(255,0,0),8,8,0);
	}
	else if(str == "Limb")
	{
		putText(draw1,"Limb",P[max1],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		putText(draw1,"Limb",P[max2],FONT_HERSHEY_SIMPLEX,1,Scalar(255,0,0),4);
		locMaxFlag[max1]='l';
		locMaxFlag[max2]='l';
		//line(draw1,P[max1],mass_centre,Scalar(255,0,0),6,8,0);
		//line(draw1,P[max2],mass_centre,Scalar(255,0,0),6,8,0);
    }

}
//calculates distance between two points
double Scene:: cal_distance(Point2f init,Point2f final)
{
	return sqrt((init.x-final.x)*(init.x-final.x)+(init.y-final.y)*(init.y-final.y));
}
//finds the returns the index f point with minimum curvature
int Scene :: find_min(vector<Point2f> head,double arr[10])
{
	double temp=0.0;
	int th_index=0;
	temp=arr[0];
	for(int i=1;i<head.size();i++)
	{
		if(arr[i]<temp && head[i]!=Point2f(0.0,0.0))
		{
			temp=arr[i];
			th_index=i;
		}
	}
	return th_index;
}

int Scene :: find_head(vector<Point2f> POI_head,double arr[10],int index)
{
	int c=0,theta_index=0;
	//cout<<"----start---"<<endl;
	for(int i=0;i<POI_head.size();i++)
	{ 
		if(POI_head[i] !=Point2f(0.0,0.0))
		{
			if(arr[i]>=80)
			{ 
				theta_index=i;
				c++;
				//cout<<" Angle : "<<arr[i]<<",";
			}
		}
	}
	//cout<<"----end-----------"<<endl;
	if(c>1)
		theta_index=find_min(POI_head,arr);
    return theta_index;
}
//find back grnd by gaussian pdf
void Scene::findBackground(Mat frame)
{
	Mat val;
	val=cv::abs(frame-mean);
	val.convertTo(val,CV_32F);
	
	
	val=val.mul(1/sigma);
	cvtColor(val,val,CV_RGB2GRAY);
	threshold(val, backgroundMask,k, 255, cv::THRESH_BINARY);
	threshold(val, isBack,k, 1, cv::THRESH_BINARY);
	//cvtColor(backgroundMask,backgroundMask,CV_RGB2GRAY);
	//backgroundMask.convertTo(backgroundMask,CV_8UC1);
	imshow("mask",backgroundMask);

	//cout<<isBack;

}
//int Scene :: kalmanPredict(int t,int flag,int cur,int k)
//{
//	if(flag==1 ||t==400)
//	{
//		X[cur][t-1]=0.0;
//		Y[cur][t-1]=0.0;
//		P[t-1]=1.0;
//		flag=0;
//	}	
//
//	X[cur][t] = X[cur][t-1];
//	Y[cur][t] = Y[cur][t-1];
//	P[t] = P[t-1];
//	
//	return flag;
//}
//void Scene :: kalmanUpdate(int t,double track_x,double track_y,int cur,int k)
//{
//	KalmanGain[t]= P[t]/(P[t]+0.1);
//	X[cur][t] = X[cur][t] + KalmanGain[t]*(track_x-X[cur][t]);
//	Y[cur][t] = Y[cur][t] + KalmanGain[t]*(track_y-Y[cur][t]);
//	
//	P[t] = (1-KalmanGain[t])*P[t];
//	
//	cout << "---------------------------------------------------------" << endl;
//	cout << "count " << t << endl;
//	cout << "Original -> " << track_x << " Estimated point -> " <<X[cur][t] << endl;
//	cout << "Original -> " << track_y << " Estimated point -> " <<Y[cur][t] << endl;
//	cout << "valur of P -> " << P[t] << endl; 
//	cout << "kalman gain -> " << KalmanGain[t] << endl;
//	cout << "----------------------------END--------------------------" << endl;
//}

void Scene :: cloth_fitting()
{
	Point reference;
	//ofstream fout("cloth.txt");
	Mat cloth = imread("redTshirt.jpg",CV_LOAD_IMAGE_UNCHANGED);
	resize(cloth,cloth,Size(200,200));
	imwrite( "cloth_Image.jpg", cloth );
	Mat shirt,temp;
	//Mat final=Mat::zeros(height,width,CV_8UC3);
	Mat channel[3];
	vector<Mat> final;
	//namedWindow("Cloth",CV_WINDOW_AUTOSIZE);
	split(cloth,channel);

	//channel[2]=channel[2]-255;
	for(int i=0;i<channel[2].rows;i++)
	{
		for(int j=0;j<channel[2].cols;j++)
		{
			if(channel[2].at<uchar>(i,j)==255)
			    channel[2].at<uchar>(i,j)=0;
		}
	}


	findNonZero(channel[2],shirt);
	//fout<< shirt <<endl;
	reference= Point(102,17);
	/*for(int i=0;i<shirt.rows;i++)
	{
		for(int j=0;j<shirt.cols;j++)
		{
			cout << shirt.at<uchar>(i,j) << endl;
		}
	}*/
	for(int i=0;i<shirt.rows;i++)
	{
		shirt.at<Point>(i)=shirt.at<Point>(i)-reference+Neck;
	}


	/*for(int i=0;i<shirt.rows;i++)
	{
		if(pointPolygonTest(R1_vertices,shirt.at<Point>(i),false)>=0)
			//frameNow.at<Vec3b>(shirt.at<Point>(i))[2]=255;
			frameNow.at<Vec3b>(shirt.at<Point>(i))=255;
		    //channel[2].at<uchar>(shirt.at<Point>(i))=255;
			//draw.at<uchar>(shirt.at<Point>(i))=255;
	}*/
	Point vertices[10000000];
	for(int i=0;i<shirt.rows;i++)
	{
		vertices[i]=shirt.at<Point>(i);
	}

	fillConvexPoly(frameNow,vertices,4,Scalar(0,0,255));
	
	
	//channel[0]=draw;
	//channel[1]=draw;
	
	//channel[0]=frameNow;
	//channel[1]=frameNow;


	/*final.push_back(draw);
	final.push_back(draw);
	final.push_back(channel[2]);*/

	//cout << shirt.at<uchar>(1,1) <<endl;
	//merge(final,temp);
	//imshow("cloth",cloth);
}



void Scene :: cloth_texture_fitting(vector<vector<Point>> contours_poly,Mat draw,Point2f mass_centre)
{
	Point reference;
	//ofstream fout("cloth.txt");
	Mat cloth = imread("redTshirt.jpg",CV_LOAD_IMAGE_UNCHANGED);
	resize(cloth,cloth,Size(200,200));
	imwrite( "cloth_Image.jpg", cloth );
	Mat shirt,temp;
	Mat Temp_image=Mat::zeros(height,width,CV_8UC1);
	//Mat final=Mat::zeros(height,width,CV_8UC3);
	Mat channel[3];
	vector<Mat> final;
	//namedWindow("Cloth",CV_WINDOW_AUTOSIZE);
	RotatedRect BigR= RotatedRect(R1_centre,Size(R1_height*1.1,2.2*R1_width),90);

    Point2f BigR_vertices2f[4];
	Point BigR_vertices[4];

	BigR.points(BigR_vertices2f);
          
	for(int i = 0; i < 4; i++)
		BigR_vertices[i] = BigR_vertices2f[i];
			
			
	fillConvexPoly(Temp_image,BigR_vertices,4,Scalar(255,255,255));
	findNonZero(Temp_image,shirt);


	float DeltaY=E1.y-Neck.y; 
	float DeltaX=Neck.x-E1.x;

	angle= atan2(DeltaY,DeltaX)*180/3.14;
	//cout <<"Angle: " <<angle << endl;
	if(mass_centre.x<E1.x)
	    angle=95+angle;
	else 
	    angle= 90-angle;


	RotatedRect SmallR= RotatedRect(R2_centre,Size(R2_height,R2_width),angle);

	Point2f SmallR_vertices2f[4];
	Point SmallR_vertices[4];

	SmallR.points(SmallR_vertices2f);
          
	for(int i = 0; i < 4; i++)
		SmallR_vertices[i] = SmallR_vertices2f[i];
			
			
	fillConvexPoly(Temp_image,SmallR_vertices,4,Scalar(255,255,255));
	findNonZero(Temp_image,shirt);



	/*split(cloth,channel);
	int ch=2;
	//channel[2]=channel[2]-255;
	for(int i=0;i<channel[ch].rows;i++)
	{
		for(int j=0;j<channel[ch].cols;j++)
		{
			if(channel[ch].at<uchar>(i,j)==255)
			    channel[ch].at<uchar>(i,j)=0;
		}
	}*/

	
	//findNonZero(channel[ch],shirt);
	//fout<< shirt <<endl;
	//reference= Point(102,17);
	/*for(int i=0;i<shirt.rows;i++)
	{
		for(int j=0;j<shirt.cols;j++)
		{
			cout << shirt.at<uchar>(i,j) << endl;
		}
	}
	for(int i=0;i<shirt.rows;i++)
	{
		shirt.at<Point>(i)=shirt.at<Point>(i)-reference+Neck;
	}*/
	for(int i=0;i<shirt.rows;i++)
	{
		if(pointPolygonTest(contours_poly[0],shirt.at<Point>(i),false)>=0)
		{	//frameNow.at<Vec3b>(shirt.at<Point>(i))[2]=255;
			frameNow.at<Vec3b>(shirt.at<Point>(i))[0]=0;    // B
 			frameNow.at<Vec3b>(shirt.at<Point>(i))[1]=0;    // G
			frameNow.at<Vec3b>(shirt.at<Point>(i))[2]=0;    // R
		}
			//channel[2].at<uchar>(shirt.at<Point>(i))=255;
			//draw.at<uchar>(shirt.at<Point>(i))=255;
	}
	//channel[0]=draw;
	//channel[1]=draw;
	
	//channel[0]=frameNow;
	//channel[1]=frameNow;


	/*final.push_back(draw);
	final.push_back(draw);
	final.push_back(channel[2]);*/

	//cout << shirt.at<uchar>(1,1) <<endl;
	//merge(final,temp);
	//imshow("cloth",cloth);
	
	
	/*imshow("cloth",draw);
	waitKey(0);
	destroyWindow("cloth");*/

}

void Scene::startProcessing()
{
	int i,j;
	int count=0;
	
	do
	{
		
		cap.read(frameNow);
		cvtColor(frameNow,grayFrame,CV_RGB2GRAY);
		Mat image32f;
		if(count==0)
		{
			Mat image32f;
			frameNow.convertTo(image32f,CV_32F);

			Mat mu;
			blur(image32f, mu, Size(3, 3));

			Mat mu2;
			blur(image32f.mul(image32f), mu2, Size(3, 3));

			cv::sqrt(mu2 - mu.mul(mu), sigma);

			variance=sigma.mul(sigma);

			//imshow("sigma",sigma);

			count++;
		}
		avgMat=updateAvg(avgMat,frameNow);
		//updateStats(frameNow);
		//findBackground(frameNow);
		matDiff=frameDiffCalc(frameNow,framePrev);
		//framePrev=frameNow;
		framePrev=frameNow.clone();
		frameNum++;
		findBackground(matDiff,matDiffPrev);
		matDiffPrev=matDiff.clone();
		object=findObject(frameNow);
		Mat res=refineObject(object);
		findcontours(object,count);
		count++;
		//Mat element = getStructuringElement(MORPH_RECT,Size(3,3),Point(0,0));
		
		//erode(object,object,kernel);
		//Mat dst=object-eroded;
		

		imshow("yahoo",frameNow);
		if(frameNum==FRAME_COUNT+1)
		    system("pause");
		if (waitKey(30) >= 0)
				break;
		  
		 //imshow("sigma",frameNow);
	 }while(1);

	
}

Mat Scene::updateAvg(Mat avgFrame,Mat frame)
{
	Mat res;
	res=(avgFrame*(N-1)+frame)/N;
	return res;
}
void Scene::updateStats(Mat frame)
{
	mean=rho*frame+(1-rho)*mean;

	
	Mat dist=cv::abs(frame-mean);
	dist.convertTo(dist,CV_32F);
	variance=rho*dist.mul(dist)+(1-rho)*variance;

	cv::sqrt(variance,sigma);
}




