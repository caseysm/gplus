#include "utils/logger.h"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

namespace gangsta {
namespace utils {

Logger::Logger() 
    : logLevel(INFO), consoleOutput(true), initialized(false) 
{
}

Logger::~Logger() 
{
    close();
}

Logger& Logger::getInstance() 
{
    static Logger instance;
    return instance;
}

void Logger::initialize(const std::string& logFileName, Level level, bool consoleOutput) 
{
    this->logLevel = level;
    this->consoleOutput = consoleOutput;
    
    // Close any open log file
    close();
    
    // Open a new log file if a name is provided
    if (!logFileName.empty()) {
        logFile.open(logFileName, std::ios::out | std::ios::app);
        if (!logFile.is_open()) {
            std::cerr << "Error: Could not open log file " << logFileName << std::endl;
        }
    }
    
    initialized = true;
    
    // Log start message
    log(INFO, "Logging started");
}

void Logger::setLevel(Level level) 
{
    logLevel = level;
}

void Logger::debug(const std::string& message) 
{
    log(DEBUG, message);
}

void Logger::info(const std::string& message) 
{
    log(INFO, message);
}

void Logger::warning(const std::string& message) 
{
    log(WARNING, message);
}

void Logger::error(const std::string& message) 
{
    log(ERROR, message);
}

void Logger::close() 
{
    if (logFile.is_open()) {
        log(INFO, "Logging stopped");
        logFile.close();
    }
    initialized = false;
}

void Logger::log(Level level, const std::string& message) 
{
    if (level < logLevel) {
        return;
    }
    
    // Get current time
    auto now = std::chrono::system_clock::now();
    auto timeT = std::chrono::system_clock::to_time_t(now);
    std::tm localTime = *std::localtime(&timeT);
    
    // Format the log message
    std::stringstream ss;
    ss << std::setfill('0') 
       << std::setw(4) << localTime.tm_year + 1900 << "-"
       << std::setw(2) << localTime.tm_mon + 1 << "-"
       << std::setw(2) << localTime.tm_mday << " "
       << std::setw(2) << localTime.tm_hour << ":"
       << std::setw(2) << localTime.tm_min << ":"
       << std::setw(2) << localTime.tm_sec 
       << " [" << getLevelString(level) << "] " 
       << message;
    
    std::string formattedMessage = ss.str();
    
    // Write to file if open
    if (logFile.is_open()) {
        logFile << formattedMessage << std::endl;
        logFile.flush();
    }
    
    // Write to console if enabled
    if (consoleOutput) {
        if (level == ERROR) {
            std::cerr << formattedMessage << std::endl;
        } else {
            std::cout << formattedMessage << std::endl;
        }
    }
}

std::string Logger::getLevelString(Level level) 
{
    switch (level) {
        case DEBUG:
            return "DEBUG";
        case INFO:
            return "INFO";
        case WARNING:
            return "WARNING";
        case ERROR:
            return "ERROR";
        default:
            return "UNKNOWN";
    }
}

} // namespace utils
} // namespace gangsta