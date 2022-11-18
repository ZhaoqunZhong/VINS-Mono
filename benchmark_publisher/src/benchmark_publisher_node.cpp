#include <cstdio>
#include <vector>
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <tf/transform_broadcaster.h>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "alignment/AlignTrajectory.h"
#include <sys/time.h>

using namespace std;
using namespace Eigen;
using namespace SlamTester;

const int SKIP = 50;
string benchmark_output_path;
string estimate_output_path;
template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
    T ans;
    if (n.getParam(name, ans))
    {
        ROS_INFO_STREAM("Loaded " << name << ": " << ans);
    }
    else
    {
        ROS_ERROR_STREAM("Failed to load " << name);
        n.shutdown();
    }
    return ans;
}

struct Data
{
    Data(FILE *f)
    {
        if (fscanf(f, " %lf,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &t,
               &px, &py, &pz,
               &qw, &qx, &qy, &qz,
               &vx, &vy, &vz,
               &wx, &wy, &wz,
               &ax, &ay, &az) != EOF)
        {
            t /= 1e9;
        }
    }
    double t;
    float px, py, pz;
    float qw, qx, qy, qz;
    float vx, vy, vz;
    float wx, wy, wz;
    float ax, ay, az;
};
int idx = 1;
vector<Data> benchmark;

ros::Publisher pub_odom;
ros::Publisher pub_path;
nav_msgs::Path path;
vector<Eigen::Matrix<double, 7, 1>> imu_poses;
vector<double> imu_ts_s;
vector<Eigen::Matrix<double, 7, 1>> gt_poses;
vector<double> gt_ts_s;

int init = 0;
Quaterniond baseRgt;
Vector3d baseTgt;
tf::Transform trans;

// void odom_callback(const nav_msgs::OdometryConstPtr &odom_msg)
// {
//     //ROS_INFO("odom callback!");
//     if (odom_msg->header.stamp.toSec() > benchmark.back().t)
//       return;
  
//     for (; idx < static_cast<int>(benchmark.size()) && benchmark[idx].t <= odom_msg->header.stamp.toSec(); idx++)
//         ;


//     if (init++ < SKIP)
//     {
//         baseRgt = Quaterniond(odom_msg->pose.pose.orientation.w,
//                               odom_msg->pose.pose.orientation.x,
//                               odom_msg->pose.pose.orientation.y,
//                               odom_msg->pose.pose.orientation.z) *
//                   Quaterniond(benchmark[idx - 1].qw,
//                               benchmark[idx - 1].qx,
//                               benchmark[idx - 1].qy,
//                               benchmark[idx - 1].qz).inverse();
//         baseTgt = Vector3d{odom_msg->pose.pose.position.x,
//                            odom_msg->pose.pose.position.y,
//                            odom_msg->pose.pose.position.z} -
//                   baseRgt * Vector3d{benchmark[idx - 1].px, benchmark[idx - 1].py, benchmark[idx - 1].pz};
//         return;
//     }

//     nav_msgs::Odometry odometry;
//     odometry.header.stamp = ros::Time(benchmark[idx - 1].t);
//     odometry.header.frame_id = "world";
//     odometry.child_frame_id = "world";

//     Vector3d tmp_T = baseTgt + baseRgt * Vector3d{benchmark[idx - 1].px, benchmark[idx - 1].py, benchmark[idx - 1].pz};
//     odometry.pose.pose.position.x = tmp_T.x();
//     odometry.pose.pose.position.y = tmp_T.y();
//     odometry.pose.pose.position.z = tmp_T.z();

//     Quaterniond tmp_R = baseRgt * Quaterniond{benchmark[idx - 1].qw,
//                                               benchmark[idx - 1].qx,
//                                               benchmark[idx - 1].qy,
//                                               benchmark[idx - 1].qz};
//     odometry.pose.pose.orientation.w = tmp_R.w();
//     odometry.pose.pose.orientation.x = tmp_R.x();
//     odometry.pose.pose.orientation.y = tmp_R.y();
//     odometry.pose.pose.orientation.z = tmp_R.z();

//     Vector3d tmp_V = baseRgt * Vector3d{benchmark[idx - 1].vx,
//                                         benchmark[idx - 1].vy,
//                                         benchmark[idx - 1].vz};
//     odometry.twist.twist.linear.x = tmp_V.x();
//     odometry.twist.twist.linear.y = tmp_V.y();
//     odometry.twist.twist.linear.z = tmp_V.z();
//     pub_odom.publish(odometry);

//     geometry_msgs::PoseStamped pose_stamped;
//     pose_stamped.header = odometry.header;
//     pose_stamped.pose = odometry.pose.pose;
//     path.header = odometry.header;
//     path.poses.push_back(pose_stamped);
//     pub_path.publish(path);
// }

// Precise gt alignment 
void alignTrajToGt(double t_offset_traj, double max_t_diff, std::vector<Eigen::Matrix<double, 7, 1>> &traj,
 std::vector<double> &traj_times, std::vector<Eigen::Matrix<double, 7, 1>> &gt,
 std::vector<double> &gt_times, Eigen::Matrix4d &transform) {
    auto traj_copy = traj;
    auto traj_time_copy = traj_times;
    auto gt_pose_copy = gt;
    auto gt_time_copy = gt_times;

    AlignUtils::perform_association(t_offset_traj, max_t_diff,
        traj_time_copy, gt_time_copy, traj_copy, gt_pose_copy);

    if (traj_copy.size() < 3)
        return;

    Eigen::Matrix3d R_gtToTraj;
    Eigen::Vector3d t_gtToTraj;
    double s_trajToGt;

    AlignTrajectory::align_trajectory(traj_copy, gt_pose_copy, R_gtToTraj, t_gtToTraj, s_trajToGt, "sim3");

    transform.setIdentity();
    transform.block<3,3>(0,0) = R_gtToTraj * s_trajToGt;
    transform.block<3,1>(0,3) = t_gtToTraj;
}

struct timeval last_align_time;
Eigen::Matrix4d gtToTraj;



void odom_callback(const nav_msgs::OdometryConstPtr &odom_msg) {
    imu_ts_s.push_back(odom_msg->header.stamp.toSec());
    Eigen::Matrix<double, 7, 1> imu_pose;
    imu_pose << odom_msg->pose.pose.position.x, odom_msg->pose.pose.position.y, odom_msg->pose.pose.position.z, 
    odom_msg->pose.pose.orientation.x, odom_msg->pose.pose.orientation.y, odom_msg->pose.pose.orientation.z, odom_msg->pose.pose.orientation.w;
    imu_poses.push_back(imu_pose);

    struct timeval time_now;
    gettimeofday(&time_now, nullptr);
    if ((time_now.tv_sec - last_align_time.tv_sec) * 1.0f + (time_now.tv_usec - last_align_time.tv_usec) / 1000000.f > 2) {
        alignTrajToGt(0, 0.02, imu_poses, imu_ts_s, gt_poses, gt_ts_s, gtToTraj);
        last_align_time = time_now;
    }

    // std::cout << "gtToTraj: \n"<< gtToTraj << "\n";

    nav_msgs::Path aligned_gt;
    for (auto &pt: gt_poses) {
        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = odom_msg->header;
        Vector4d align_to_body = gtToTraj.inverse() * Vector4d{pt[0], pt[1], pt[2], 1};
        pose_stamped.pose.position.x = align_to_body.x();
        pose_stamped.pose.position.y = align_to_body.y();
        pose_stamped.pose.position.z = align_to_body.z();
        aligned_gt.header = odom_msg->header;
        aligned_gt.poses.push_back(pose_stamped);
    }


    pub_path.publish(aligned_gt);
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "benchmark_publisher");
    ros::NodeHandle n("~");

    string csv_file = readParam<string>(n, "data_name");
    std::cout << "load ground truth " << csv_file << std::endl;
    FILE *f = fopen(csv_file.c_str(), "r");
    if (f==NULL)
    {
      ROS_WARN("can't load ground truth; wrong path");
      //std::cerr << "can't load ground truth; wrong path " << csv_file << std::endl;
      return 0;
    }
    char tmp[10000];
    if (fgets(tmp, 10000, f) == NULL)
    {
        ROS_WARN("can't load ground truth; no data available");
    }
    while (!feof(f))
        benchmark.emplace_back(f);
    fclose(f);
    benchmark.pop_back();
    ROS_INFO("Data loaded: %d", (int)benchmark.size());

    for (auto &pt : benchmark) {
        gt_ts_s.push_back(pt.t);
        Eigen::Matrix<double, 7, 1> gt_pose;
        gt_pose << pt.px, pt.py, pt.pz, pt.qx, pt.qy, pt.qz, pt.qw;
        gt_poses.push_back(gt_pose);
    }

    pub_odom = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);

    ros::Subscriber sub_odom = n.subscribe("estimated_odometry", 1000, odom_callback);
    
    ros::Rate r(20);
    ros::spin();
}
