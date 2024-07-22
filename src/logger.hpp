#pragma once

#include <string>
#include <fstream>

namespace geo {
    enum class LogLevel {
        INFO = 34, // Default terminal color
        WARNING = 33, // Yellow
        ERROR = 31, // Red
        DEFAULT = 0 // Blue
    };

    void log(const std::string& message, LogLevel color = LogLevel::DEFAULT);

    std::string readFullLog(const std::string& logFilePath);

}
