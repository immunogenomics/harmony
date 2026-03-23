#include "timer.h"

#include <Rcpp.h>


using namespace Rcpp;

Timer::Timer(std::string _task_name, double& _timer): timer(_timer){
    this->task_name = _task_name;
    this->start = high_resolution_clock::now();    
}

Timer::Timer(TimerMeta& t) : timer(t.elapsed){
    t.last = this->timer;
    this->start = high_resolution_clock::now();
    t.iter++;
}
  
double Timer::getLapse() {
    return duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - this->start).count();
}
  
Timer::~Timer(){
    auto time_span = duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - this->start);
    timer += time_span.count();
    // Rcout << "Task: " << this->task_name << " took " << time_span.count() << " seconds" <<std::endl;
}



