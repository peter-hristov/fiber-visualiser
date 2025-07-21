#pragma once

#include <iostream>
#include <chrono>
#include <memory>
#include <mutex>
#include <string>

#include "./utility/indicators.hpp"

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

    static std::unique_ptr<indicators::ProgressBar> getLoadingBar()
    {
        using namespace indicators;

        auto bar = std::make_unique<ProgressBar>(
                option::BarWidth{50},
                option::Start{"["},
                option::Fill{"■"},
                option::Lead{"■"},
                option::Remainder{"-"},
                option::End{" ]"},
                option::PostfixText{"Computing Reeb space."},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true}
                );

        return bar;
    }

private:
    static inline std::chrono::high_resolution_clock::time_point start_time;
};
