#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <opencv2/core/core.hpp>
#include <nav_msgs/Odometry.h>

#include <chrono>

const double DEG2RAD = M_PI / 180.0;
#define THETA_THRESHOLD 0.0001 //sin a = a 

class Utils{
    // 重力常数
    public:
        static Eigen::Vector3d gravity;
    // 数学函数
    public:
    static Eigen::Matrix3d skew(const Eigen::Vector3d &v); // 反对称矩阵

    static Eigen::Vector3d RToSO3(const Eigen::Matrix3d &R); // 旋转矩阵R -> SO3
    static Eigen::Quaterniond RToq(const Eigen::Matrix3d &R); //  R-> q(四元数)
    static Eigen::Matrix3d R_normalize(const Eigen::Matrix3d& R_in);  // R的标准化矩阵

    static Eigen::Matrix3d qToR(const Eigen::Quaterniond &q); // q-> R
    static Eigen::Vector3d qToSO3(const Eigen::Quaterniond &q); // q->SO3

    static Eigen::Quaterniond SO3Toq(const Eigen::Vector3d &SO3); //SO3->q
    static Eigen::Matrix3d SO3ToR(const Eigen::Vector3d &SO3); // SO3->R

    static Eigen::Matrix3d Jl_SO3(const Eigen::Vector3d &SO3); // SO3 左扰动
    static Eigen::Matrix3d Jr_SO3(const Eigen::Vector3d &SO3); // SO3 右扰动
    static Eigen::Matrix3d Jl_SO3_inv(const Eigen::Vector3d &SO3); // SO3 左扰动逆
    static Eigen::Matrix3d Jr_SO3_inv(const Eigen::Vector3d &SO3); // SO3 右扰动逆

    //[q]_L [p]_R = [p]_R [q]_L 四元数相乘也可以写成矩阵和向量相乘的形式
    static Eigen::Matrix4d ql_mat(const Eigen::Quaterniond &q); //[q]_L
    static Eigen::Matrix4d qr_mat(const Eigen::Quaterniond &q); //[q]_R
    static Eigen::Quaterniond deltaQ(const Eigen::Vector3d &theta); // delta_q
    static Eigen::Isometry3d odomToEigen(const nav_msgs::OdometryConstPtr &odom_msg); // Eigen::Isometry3d 欧式变换矩阵 4x4
    static Eigen::Isometry3d SO3ToT(const Eigen::Matrix3d &rotation, const Eigen::Vector3d &translation);

};

class TicToc{
public:
    TicToc(std::string name_in);
    void tic();
    void toc(int freq);

private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
    unsigned long total_frame;
    double total_time;
    std::string name;
};