#ifndef GANGSTA_LOGGER_H
#define GANGSTA_LOGGER_H

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

namespace gangsta {
namespace utils {

/**
 * @brief Simple logging utility
 */
class Logger 
{
public:
    /**
     * @brief Log levels
     */
    enum Level {
        DEBUG = 0,
        INFO = 1,
        WARNING = 2,
        ERROR = 3,
        NONE = 4
    };
    
    /**
     * @brief Get the singleton instance
     * @return Reference to the logger
     */
    static Logger& getInstance();
    
    /**
     * @brief Initialize the logger
     * @param logFile Path to log file (empty for console only)
     * @param level Minimum log level
     * @param consoleOutput Whether to also output to console
     */
    void initialize(const std::string& logFile = "", Level level = INFO, bool consoleOutput = true);
    
    /**
     * @brief Set the log level
     * @param level Minimum log level
     */
    void setLevel(Level level);
    
    /**
     * @brief Log a debug message
     * @param message Message to log
     */
    void debug(const std::string& message);
    
    /**
     * @brief Log an info message
     * @param message Message to log
     */
    void info(const std::string& message);
    
    /**
     * @brief Log a warning message
     * @param message Message to log
     */
    void warning(const std::string& message);
    
    /**
     * @brief Log an error message
     * @param message Message to log
     */
    void error(const std::string& message);
    
    /**
     * @brief Close the log file
     */
    void close();
    
private:
    /**
     * @brief Private constructor (singleton)
     */
    Logger();
    
    /**
     * @brief Private destructor
     */
    ~Logger();
    
    /**
     * @brief Log a message with the given level
     * @param level Message level
     * @param message Message text
     */
    void log(Level level, const std::string& message);
    
    /**
     * @brief Get string representation of log level
     * @param level Log level
     * @return String representation
     */
    std::string getLevelString(Level level);
    
    std::ofstream logFile;       ///< Log file stream
    Level logLevel;              ///< Current log level
    bool consoleOutput;          ///< Whether to output to console
    bool initialized;            ///< Whether the logger is initialized
    
    // Make the logger a singleton
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
};

} // namespace utils
} // namespace gangsta

#endif // GANGSTA_LOGGER_H