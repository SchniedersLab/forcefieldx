//
// File:       procs.h
//

#define test_start()
#define log_perf(_number, _higherBetter, _numType, _format, ...) printf("Performance Number " _format " (in %s, %s): %g\n",##__VA_ARGS__, _numType, _higherBetter?"higher is better":"lower is better" , _number)
#define log_info printf
#define log_error printf
#define test_finish()
