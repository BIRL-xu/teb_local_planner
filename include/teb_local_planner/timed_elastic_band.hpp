/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2016,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the institute nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Christoph Rösmann
 *********************************************************************/

#include <teb_local_planner/timed_elastic_band.h>

namespace teb_local_planner
{

  
/**
 * 通过参考路径的迭代器来对路径进行轨迹参数化，主要是为每一个点设置姿态角，除起点以外的点设置时间差。
 * 时间差分为两种:以最大速度作匀速运动通过两点间距离时所需时间和以最大加速度作匀加速运动通过两点间距离时所需时间，最终取较大者。
 * 需要注意的是，前进运动和后退运动在计算姿态角时，角度相差180°； 起点和终点被标记为固定点。 
 */
template<typename BidirIter, typename Fun>
bool TimedElasticBand::initTrajectoryToGoal(BidirIter path_start, BidirIter path_end, Fun fun_position, double max_vel_x, double max_vel_theta,
                                     boost::optional<double> max_acc_x, boost::optional<double> max_acc_theta,
                                     boost::optional<double> start_orientation, boost::optional<double> goal_orientation, int min_samples, bool guess_backwards_motion) 
{
    Eigen::Vector2d start_position = fun_position( *path_start );
    Eigen::Vector2d goal_position = fun_position( *boost::prior(path_end) );
    
    bool backwards = false;
    
    double start_orient, goal_orient;
    if (start_orientation)
    {
      start_orient = *start_orientation; // 使用传入的起始姿态角
      
      // check if the goal is behind the start pose (w.r.t. start orientation)
      if (guess_backwards_motion && (goal_position-start_position).dot(Eigen::Vector2d(std::cos(start_orient), std::sin(start_orient))) < 0) // 通过矢量点乘判断向量方向
        backwards = true;
    }
    else
    {
      // 起始姿态角为起点指向目标点。
      Eigen::Vector2d start2goal =  goal_position - start_position;
      start_orient = atan2(start2goal[1],start2goal[0]);
    }

    double timestep = 1; // TODO: time

    if (goal_orientation)
    {
      // 使用传入的目标姿态角
      goal_orient = *goal_orientation;
    }
    else
    {
      // 目标姿态角和起点姿态角一样。
      goal_orient = start_orient;
    }
    
    if (!isInit())
    {
      // 添加起点，标记为固定点，起点没有时间差
      addPose(start_position, start_orient, true); // add starting point and mark it as fixed for optimization		

      // we insert middle points now (increase start by 1 and decrease goal by 1)
      std::advance(path_start,1); // 起始迭代器前进1
      std::advance(path_end,-1);  // 终止迭代器后退1
      int idx=0;
      for (; path_start != path_end; ++path_start) // insert middle-points
      {
            //Eigen::Vector2d point_to_goal = path.back()-*it;
            //double dir_to_goal = atan2(point_to_goal[1],point_to_goal[0]); // direction to goal
            // Alternative: Direction from last path
            Eigen::Vector2d curr_point = fun_position(*path_start); // 获取当前迭代器所指向的点
            // 有前一个点指向当前点的向量。
            Eigen::Vector2d diff_last = curr_point - Pose(idx).position(); // we do not use boost::prior(*path_start) for those cases,
                                                                        // where fun_position() does not return a reference or is expensive.
                                                                        
            double diff_norm = diff_last.norm(); // 向量的模，即两个相邻点的距离。
            
            double timestep_vel = diff_norm/max_vel_x; // constant velocity，计算该段距离以最大速度通过所需的时间，即t=s/v.
            double timestep_acc;

            if (max_acc_x)
            {
              timestep_acc = sqrt(2*diff_norm/(*max_acc_x)); // constant acceleration, 以最大加速度加速通过该段距离所需的时间，即s=(a*t^2)/2 --> t = sqrt(2*s/a).
               // 使用较大的时间值
              if (timestep_vel < timestep_acc && max_acc_x) timestep = timestep_acc;
              else timestep = timestep_vel;
            }
            else timestep = timestep_vel;
            
            if (timestep<0) timestep=0.2; // TODO: this is an assumption
            
            // 当前点的姿态角为由前一个点指向当前点的方向角，如果是后退，则反向。
            double yaw = atan2(diff_last[1],diff_last[0]);
            if (backwards)
                yaw = g2o::normalize_theta(yaw + M_PI); // 反向
            addPoseAndTimeDiff(curr_point, yaw ,timestep); // 添加该点和对应的姿态角，时间差
            
            /*
            // TODO: the following code does not seem to hot-start the optimizer. Instead it recudes convergence time.

            Eigen::Vector2d diff_next = fun_position(*boost::next(path_start))-curr_point; // TODO maybe store the boost::next for the following iteration
            double ang_diff = std::abs( g2o::normalize_theta( atan2(diff_next[1],diff_next[0])
                                                            -atan2(diff_last[1],diff_last[0]) ) );
            
            timestep_vel = ang_diff/max_vel_theta; // constant velocity
            if (max_acc_theta)
            {
                    timestep_acc = sqrt(2*ang_diff/(*max_acc_theta)); // constant acceleration
                    if (timestep_vel < timestep_acc) timestep = timestep_acc;
                    else timestep = timestep_vel;
            }
            else timestep = timestep_vel;
            
            if (timestep<0) timestep=0.2; // TODO: this is an assumption
            
            yaw = atan2(diff_last[1],diff_last[0]); // TODO redundant right now, not yet finished
            if (backwards)
                yaw = g2o::normalize_theta(yaw + M_PI);
            addPoseAndTimeDiff(curr_point, yaw ,timestep);

            */
            
            ++idx;
      }
      // 处理最后一段，即目标点和它的前一个点
      Eigen::Vector2d diff = goal_position-Pose(idx).position();
      double diff_norm = diff.norm();
      double timestep_vel = diff_norm/max_vel_x; // constant velocity
      if (max_acc_x)
      {
            double timestep_acc = sqrt(2*diff_norm/(*max_acc_x)); // constant acceleration
            if (timestep_vel < timestep_acc) timestep = timestep_acc;
            else timestep = timestep_vel;
      }
      else timestep = timestep_vel;

      
      PoseSE2 goal(goal_position, goal_orient);
      
      // 如果已添加的顶点数少于参数指定的最小数减1（目标点还没有添加），则需要在目标点和已添加的最后一个顶点之间添加
      // 更多的点以满足数量。添加规则是: 每次插入一个已添加的最后顶点和目标点的中间点(位置平均值和姿态平均值)，时间差也是减半。
      // Baul:这种插入方法会导致越靠近目标点，插入点越密集，时间差越小，为什么不用需要插入的点数量来等分插入呢？
      // if number of samples is not larger than min_samples, insert manually
      if ( sizePoses() < min_samples-1 )
      {
        ROS_DEBUG("initTEBtoGoal(): number of generated samples is less than specified by min_samples. Forcing the insertion of more samples...");
        while (sizePoses() < min_samples-1) // subtract goal point that will be added later
        {
          // Each inserted point bisects the remaining distance. Thus the timestep is also bisected.
          timestep /= 2;
          // simple strategy: interpolate between the current pose and the goal
          addPoseAndTimeDiff( PoseSE2::average(BackPose(), goal), timestep ); // let the optimier correct the timestep (TODO: better initialization	
        }
      }
      
      // 最后添加目标点，设为固定点，时间差为上一个点到它的时间差。
      // now add goal
      addPoseAndTimeDiff(goal, timestep); // add goal point
      setPoseVertexFixed(sizePoses()-1,true); // GoalConf is a fixed constraint during optimization
    }
    else // size!=0
    {
      ROS_WARN("Cannot init TEB between given configuration and goal, because TEB vectors are not empty or TEB is already initialized (call this function before adding states yourself)!");
      ROS_WARN("Number of TEB configurations: %d, Number of TEB timediffs: %d", sizePoses(), sizeTimeDiffs());
      return false;
    }
    return true;
}  
  

} // namespace teb_local_planner



