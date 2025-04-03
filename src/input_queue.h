#ifndef INPUT_QUEUE_H
#define INPUT_QUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>
#include <Rcpp.h>

// Thread-safe blocking queue for Rcpp::List inputs
template<typename T>
class InputQueue {
private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_;
    bool finished_ = false;

public:
    void push(T item) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(item);
        cond_.notify_one();
    }

    T pop() {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [&] { return !queue_.empty() || finished_; });

        if (queue_.empty()) {
            return T();  // Return default (empty) if finished
        }

        T item = queue_.front();
        queue_.pop();
        return item;
    }

    void set_finished() {
        std::lock_guard<std::mutex> lock(mutex_);
        finished_ = true;
        cond_.notify_all();
    }
};

#endif
