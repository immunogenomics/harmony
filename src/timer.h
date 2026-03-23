#pragma once

#include <chrono>
#include <iostream>
#include <map>

using namespace std::chrono;
typedef high_resolution_clock::time_point Timepoint;


struct TimerMeta{
  double elapsed = 0;
  double last = 0;
  unsigned iter = 0;
};

static std::map<std::string, TimerMeta> timers;



void print_timers();

class Timer{
 public:
  Timepoint start;
  std::string task_name;
  double& timer;
  Timer(std::string _task_name, double& _timer);
  Timer(TimerMeta& t);
  double getLapse();
  ~Timer();
};
