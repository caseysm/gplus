#ifndef TEST_FILE_UTILS_H
#define TEST_FILE_UTILS_H

#include <string>
#include <sys/stat.h>
#include <sys/types.h>

namespace test {
namespace utils {

/**
 * @brief Check if a file exists
 * @param filePath Path to the file
 * @return True if file exists
 */
inline bool fileExists(const std::string& filePath) {
    struct stat buffer;
    return (stat(filePath.c_str(), &buffer) == 0);
}

/**
 * @brief Check if a directory exists
 * @param dirPath Path to the directory
 * @return True if directory exists
 */
inline bool directoryExists(const std::string& dirPath) {
    struct stat buffer;
    return (stat(dirPath.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode));
}

/**
 * @brief Create a directory if it doesn't exist
 * @param dirPath Path to the directory
 * @return True if directory exists or was created successfully
 */
inline bool createDirectory(const std::string& dirPath) {
    if (directoryExists(dirPath)) {
        return true;
    }
    
    #ifdef _WIN32
    return _mkdir(dirPath.c_str()) == 0;
    #else
    return mkdir(dirPath.c_str(), 0755) == 0;
    #endif
}

/**
 * @brief Remove a file
 * @param filePath Path to the file
 * @return True if file was removed successfully
 */
inline bool removeFile(const std::string& filePath) {
    return remove(filePath.c_str()) == 0;
}

} // namespace utils
} // namespace test

#endif // TEST_FILE_UTILS_H