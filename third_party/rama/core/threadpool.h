// Bounded per-call pool. Unlike the upstream discarded futures, errors are propagated.
#pragma once
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>
#include <stdexcept>
namespace ramag::rama_source {
class ThreadPool {
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable ready, drained;
    std::size_t pending{};
    bool stop{};
    std::exception_ptr error;
public:
    explicit ThreadPool(std::size_t n) {
        try {
            for (std::size_t i=0;i<n;++i) workers.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock lock(mutex);
                        ready.wait(lock,[this]{ return stop || !tasks.empty(); });
                        if (tasks.empty() && stop) return;
                        task=std::move(tasks.front()); tasks.pop();
                    }
                    try { task(); }
                    catch (...) { std::lock_guard lock(mutex); if (!error) error=std::current_exception(); }
                    {
                        std::lock_guard lock(mutex);
                        if (--pending==0) drained.notify_all();
                    }
                }
            });
        } catch (...) {
            { std::lock_guard lock(mutex); stop=true; }
            ready.notify_all();
            for (auto& w:workers) w.join();
            throw;
        }
    }
    template<class F> void enqueue(F&& f) {
        std::lock_guard lock(mutex);
        if (stop) throw std::runtime_error("enqueue on stopped RaMA pool");
        tasks.emplace(std::forward<F>(f)); ++pending; ready.notify_one();
    }
    void waitAllTasksDone() {
        std::unique_lock lock(mutex);
        drained.wait(lock,[this]{return pending==0;});
        if (error) std::rethrow_exception(error);
    }
    ~ThreadPool() {
        { std::unique_lock lock(mutex);
          drained.wait(lock,[this]{return pending==0;}); stop=true; }
        ready.notify_all();
        for (auto& w:workers) w.join();
    }
    ThreadPool(const ThreadPool&)=delete;
};
}

