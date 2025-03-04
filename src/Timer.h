#pragma once

#include <iostream>
#include <chrono>
#include <string>

class Timer {
public:
    static void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    static void stop(const std::string& message = "Elapsed time") {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        std::cout << message << ": " << elapsed.count() << " s\n";
    }

private:
    static inline std::chrono::high_resolution_clock::time_point start_time;
};
