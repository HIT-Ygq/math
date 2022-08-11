#include "utils.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

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
    static Eigen::Isometry3d SO3ToT(const Eigen::Vector3d &rotation, const Eigen::Vector3d &translation);

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

//define global variables
Eigen::Vector3d Utils::gravity = Eigen::Vector3d(0,0,9.81); // 重力加速度

Eigen::Matrix3d Utils::skew(const Eigen::Vector3d &v) {  // 反对称矩阵
    Eigen::Matrix3d w;
    w << 0., -v(2), v(1),
        v(2), 0., -v(0),
        -v(1), v(0), 0.;
    return w;
}

// 旋转矩阵R的单位化
Eigen::Matrix3d Utils::R_normalize(const Eigen::Matrix3d &R) {
    Eigen::Quaterniond q(R);
    q.normalize();
    return q.toRotationMatrix(); 
}

Eigen::Vector3d Utils::RToSO3(const Eigen::Matrix3d &R_in) {
    Eigen::Matrix3d R = R_normalize(R_in);
    double theta = acos((R(0,0)+R(1,1)+R(2,2)-1.0)/2.0); // 转角 theta = [trace(R)-1] / 2
    if(theta < THETA_THRESHOLD) {
        return Eigen::Vector3d(R(2,1)-R(1,2), R(0,2)-R(2,0),R(1,0)-R(0,1)) / 2.0;
    }else{
        return theta*Eigen::Vector3d(R(2,1)- R(1,2),R(0,2)-R(2,0),R(1,0)-R(0,1)) / (2.0*sin(theta));
    }
}

Eigen::Quaterniond Utils::RToq(const Eigen::Matrix3d &R) {
    return Eigen::Quaterniond(R);
}

Eigen::Matrix3d Utils::qToR(const Eigen::Quaterniond &q) {
    return q.toRotationMatrix();
}

Eigen::Matrix3d Utils::SO3ToR(const Eigen::Vector3d &SO3) {
    double theta = SO3.norm(); // 向量SO3的范数
    if(theta < THETA_THRESHOLD) {
        Eigen::Matrix3d u_x =  skew(SO3); 
         // u_x 的二阶 泰勒展开式
        return Eigen::Matrix3d::Identity() + u_x + 0.5 * u_x * u_x;
    }else{
        Eigen::Matrix3d u_x = skew(SO3.normalized());
        return Eigen::Matrix3d::Identity() + sin(theta)* u_x + (1-cos(theta)) * u_x * u_x;
    }
}

Eigen::Vector3d Utils::qToSO3(const Eigen::Quaterniond &q){
    return RToSO3(q.toRotationMatrix());
}

Eigen::Matrix3d  Utils::Jl_SO3(const Eigen::Vector3d &SO3) {
    double theta = SO3.norm();
    if(theta < THETA_THRESHOLD) {
        return Eigen::Matrix3d::Identity() + skew(SO3)/2;
    }else{
        Eigen::Vector3d u = SO3.normalized();
        return sin(theta)/theta * Eigen::Matrix3d::Identity() + (1 - sin(theta)/theta)*u*u.transpose() + ((1-cos(theta))/theta)*skew(u);

    }
}

Eigen::Matrix3d Utils::Jr_SO3(const Eigen::Vector3d &SO3) {
    double theta = SO3.norm();
    if(theta < THETA_THRESHOLD) {
        return Eigen::Matrix3d::Identity() - skew(SO3)/2;
    }else {
        Eigen::Vector3d u=SO3.normalized();
        return sin(theta)/theta * Eigen::Matrix3d::Identity() + (1 - sin(theta)/theta)*u*u.transpose() - ((1-cos(theta))/theta)*skew(u);
    }
}

Eigen::Matrix3d Utils::Jl_SO3_inv(const Eigen::Vector3d &SO3) {
    double theta = SO3.norm();
    if(theta < THETA_THRESHOLD) {
        return  cos(theta/2) * Eigen::Matrix3d::Identity() + 0.125 * SO3 * SO3.transpose() - 0.5 * skew(SO3); 
    }else{
        Eigen::Vector3d u = SO3.normalized();
        return 0.5 * theta / tan(theta/2) * Eigen::Matrix3d::Identity() + (1 - 0.5 * theta/tan(theta/2))*u*u.transpose() - 0.5 * skew(SO3); 
    }
}

Eigen::Matrix3d Utils::Jr_SO3_inv(const Eigen::Vector3d &SO3) {
    double theta = SO3.norm();
    if(theta < THETA_THRESHOLD) {
        return cos(theta/2) * Eigen::Matrix3d::Identity() + 0.125 * SO3 * SO3.transpose() - 0.5 * skew(SO3);
    }else{
        Eigen::Vector3d u = SO3.normalized();
        return 0.5 * theta / tan(theta/2) * Eigen::Matrix3d::Identity() + (1 - 0.5 * theta/tan(theta/2))*u*u.transpose() + 0.5 * skew(SO3);
    }
}

// 四元数 q 转换为 q_L矩阵
/*                    | 0 -[q_v]^T|
    [q]_L = q_w * I + |           |
                      |q_v [q_v]^ |
*/

Eigen::Matrix4d Utils::ql_mat(const Eigen::Quaterniond &q) {
    Eigen::Matrix4d m4 = Eigen::Matrix4d::Zero();
    m4.block<3,1>(0,1) = -q.vec();
    m4.block<1,3>(1,0) = q.vec();
    m4.block<3,3>(1,1) = skew(q.vec());
    m4 += Eigen::Matrix4d::Identity() * q.w(); 
    return m4;
}

// 四元数 q 转换为 q_R矩阵
/*                    | 0  -[q_v]^T|
    [q]_L = q_w * I + |           |
                      |q_v -[q_v]^ |
*/
Eigen::Matrix4d Utils::qr_mat(const Eigen::Quaterniond &q) {
    Eigen::Matrix4d m4 = Eigen::Matrix4d::Zero();
    m4.block<3,1>(0,1) = -q.vec();
    m4.block<1,3>(1,0) = q.vec();
    m4.block<3,3>(1,1) = -skew(q.vec());
    m4 += Eigen::Matrix4d::Identity() * q.w();
    return m4;
}

Eigen::Quaterniond Utils::deltaQ(const Eigen::Vector3d &theta) {
    Eigen::Quaterniond dq;
    Eigen::Vector3d half_theta = theta;
    half_theta /= 2.0;
    dq.w() = 1.0; //small value can take as 1
    dq.x() = half_theta.x();
    dq.y() = half_theta.y();
    dq.z() = half_theta.z();
    return dq;
}

// ros nav_msgs::Odometry获取四元数与平移,eigen类型的变量
// Isometry Transform)欧式变换也称为等距变换
Eigen::Isometry3d Utils::odomToEigen(const nav_msgs::OdometryConstPtr& odom_msg){
    Eigen::Isometry3d odom = Eigen::Isometry3d::Identity();
    Eigen::Quaterniond q_temp(odom_msg->pose.pose.orientation.w,odom_msg->pose.pose.orientation.x,
                              odom_msg->pose.pose.orientation.y,odom_msg->pose.pose.orientation.z);

    // .linear()&.linearExt():返回变换的线性部分，对于Isometry而言就是旋转对应的旋转矩阵，Eigen::Block类型
    odom.linear() = q_temp.toRotationMatrix();
    odom.translation().x() = odom_msg->pose.pose.position.x;
    odom.translation().y() = odom_msg->pose.pose.position.y;
    odom.translation().z() = odom_msg->pose.pose.position.z;
    return odom;
}

Eigen::Isometry3d Utils::SO3ToT(const Eigen::Vector3d &rotation, const Eigen::Vector3d &translation){
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();//虽然称为3d，实质上是4x4的矩阵(旋转R+平移t)
    T.linear() = SO3ToR(rotation);//旋转部分赋值
    T.translation() = translation;//平移部分赋值
    return T;
}

TicToc::TicToc(std::string name_in){
    start = std::chrono::system_clock::now();
    total_frame = 0;
    total_time = 0.0;
    name = name_in;
}

void TicToc::tic(){
    start = std::chrono::system_clock::now();
}

void TicToc::toc(int freq){
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start; // 单位 ms
    total_time += elapsed_seconds.count() * 1000;
    total_frame++;
    if(total_frame%freq==0)
        std::cout<<"the average time of " << name <<" is "<< elapsed_seconds.count() * 1000 << "ms"<<std::endl;
}