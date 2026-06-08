#pragma once
#include <chrono>

class Chrono {
private:
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> end_time;
    bool is_running;

public:
    Chrono() : is_running(false) {}

    void start() {
        start_time = std::chrono::steady_clock::now();
        is_running = true;
    }

    void stop() {
        end_time = std::chrono::steady_clock::now();
        is_running = false;
    }

    double elapsed_ms() const {
        auto current_end = is_running ? std::chrono::steady_clock::now() : end_time;
        std::chrono::duration<double, std::milli> elapsed = current_end - start_time;
        return elapsed.count();
    }
};