#pragma once

#include <iostream>
#include <sstream>
#include <string>

enum class LogLevel {
    DEBUG = 0,
    INFO,
    WARNING,
    IMPORTANT,
    ERROR
};

// Global log configuration (no locking)
class LogConfig {
public:
    static void set_min_level(LogLevel lvl) { min_level() = lvl; }

    static void set_min_level(const std::string& level_str) {

        std::string lc_level_str = level_str;
        std::transform(lc_level_str.begin(), lc_level_str.end(), lc_level_str.begin(),
                                 [](unsigned char c){ return std::tolower(c); });

        if (lc_level_str == "debug")
            set_min_level(LogLevel::DEBUG);
        else if (lc_level_str == "info")
            set_min_level(LogLevel::INFO);
        else if (lc_level_str == "warning")
            set_min_level(LogLevel::WARNING);
        else if (lc_level_str == "error")
            set_min_level(LogLevel::ERROR);
        else if (lc_level_str == "important")
            set_min_level(LogLevel::IMPORTANT);
        else
            throw std::invalid_argument("Invalid log level: " + lc_level_str);
    }

    static LogLevel get_min_level() { return min_level(); }

private:
    static LogLevel& min_level() {
        static LogLevel lvl = LogLevel::INFO; // default: log everything
        return lvl;
    }
};

// Logger class
class Logger {
public:
    explicit Logger(LogLevel lvl) : level(lvl) {}

    template<typename T>
    Logger& operator<<(const T& value) {
        if (should_log())
            buffer << value;
        return *this;
    }

    // Handle std::endl and similar manipulators
    Logger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (should_log()) {
            buffer << manip;
            //flush();
        }
        return *this;
    }

    ~Logger() {
        if (!flushed && should_log())
            flush();
    }

private:
    LogLevel level;
    std::ostringstream buffer;
    bool flushed = false;

    bool should_log() const {
        return level >= LogConfig::get_min_level();
    }

    void flush() {
        //std::cout << "[" << level_name(level) << "] " << buffer.str();
        std::cout << buffer.str();
        std::cout.flush();
        flushed = true;
    }

    /*static const char* level_name(LogLevel lvl) {
        switch (lvl) {
            case LogLevel::DEBUG:   return "DEBUG";
            case LogLevel::INFO:    return "INFO";
            case LogLevel::WARNING: return "WARNING";
            case LogLevel::ERROR:   return "ERROR";
        }
        return "?";
    }*/
};

// Factory function for convenience
inline Logger log(LogLevel lvl) {
    return Logger(lvl);
}
