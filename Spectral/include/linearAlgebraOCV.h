#ifndef LINEARALGEBRAOCV_H
#define LINEARALGEBRAOCV_H

#include <iostream>
#include <opencv2/core.hpp>

void testLinearAlgebra() {

  // Set points
  cv::Point2f pt1, pt2, pt3;

  pt1 = cv::Point2f(1, 0);
  pt2 = cv::Point2f(0, 2);
  pt3 = cv::Point2f(3, 2);

  // Set Matrix
  cv::Mat A = (cv::Mat_<float>(3, 3) << std::pow(pt1.x, 2), pt1.x, 1,
	       std::pow(pt2.x, 2), pt2.x, 1, std::pow(pt3.x, 2), pt3.x, 1);

  cv::Mat b = (cv::Mat_<float>(3, 1) << pt1.y, pt2.y, pt3.y);

  cv::Mat x;

  cv::solve(A, b, x);
  std::cout << x.at<float>(0) << std::endl
	    << x.at<float>(1) << std::endl
	    << x.at<float>(2) << std::endl;

  return;
}

#endif
