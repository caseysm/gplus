#ifndef GANGSTA_FILE_H
#define GANGSTA_FILE_H

#include <string>
#include <fstream>
#include <vector>

namespace gangsta {
namespace utils {

/**
 * @brief File handling utility class
 */
class File 
{
public:
    /**
     * @brief Default constructor
     */
    File();
    
    /**
     * @brief Constructor with file path
     * @param filePath Path to the file
     */
    explicit File(const std::string& filePath);
    
    /**
     * @brief Destructor, ensures file is closed
     */
    ~File();
    
    /**
     * @brief Open a file
     * @param filePath Path to the file
     * @param mode Open mode (default is read)
     * @return True if file opened successfully
     */
    bool open(const std::string& filePath, std::ios_base::openmode mode = std::ios_base::in);
    
    /**
     * @brief Close the file
     */
    void close();
    
    /**
     * @brief Check if file is open
     * @return True if file is open
     */
    bool isOpen() const;
    
    /**
     * @brief Check if file exists
     * @param filePath Path to the file
     * @return True if file exists
     */
    static bool exists(const std::string& filePath);
    
    /**
     * @brief Read entire file into a string
     * @return File contents
     */
    std::string readAll();
    
    /**
     * @brief Read entire file into a vector of lines
     * @return Vector of file lines
     */
    std::vector<std::string> readLines();
    
    /**
     * @brief Read next line from file
     * @param line String to store the line
     * @return True if a line was read
     */
    bool readLine(std::string& line);
    
    /**
     * @brief Write string to file
     * @param text Text to write
     * @return True if write was successful
     */
    bool write(const std::string& text);
    
    /**
     * @brief Write line to file (adds newline)
     * @param line Line to write
     * @return True if write was successful
     */
    bool writeLine(const std::string& line);
    
    /**
     * @brief Get the file path
     * @return File path
     */
    const std::string& getPath() const;
    
    /**
     * @brief Get file size
     * @return File size in bytes
     */
    size_t size();
    
    /**
     * @brief Get filename from path
     * @param filePath Path to the file
     * @return Filename without directory
     */
    static std::string getFileName(const std::string& filePath);
    
    /**
     * @brief Get file extension
     * @param filePath Path to the file
     * @return File extension
     */
    static std::string getExtension(const std::string& filePath);
    
    /**
     * @brief Get directory from path
     * @param filePath Path to the file
     * @return Directory path
     */
    static std::string getDirectory(const std::string& filePath);
    
private:
    std::string filePath;   ///< Path to the file
    std::fstream fileStream; ///< File stream
};

} // namespace utils
} // namespace gangsta

#endif // GANGSTA_FILE_H