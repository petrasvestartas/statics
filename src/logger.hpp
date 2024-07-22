#pragma once

#include <string>
#include <fstream>

namespace geo {
    enum class LogLevel {
        INFO = 34, // Blue
        WARNING = 33, // Yellow
        ERROR = 31, // Red
        DEFAULT = 0 // Default terminal color
    };

    void log(const std::string& message, LogLevel color = LogLevel::INFO);

    std::string readFullLog(const std::string& logFilePath);

}
