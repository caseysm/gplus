#include "utils/file.h"
#include <fstream>
#include <sys/stat.h>

namespace gangsta {
namespace utils {

File::File() 
{
}

File::File(const std::string& filePath) 
{
    open(filePath);
}

File::~File() 
{
    close();
}

bool File::open(const std::string& filePath, std::ios_base::openmode mode) 
{
    // Close any previously opened file
    close();
    
    this->filePath = filePath;
    fileStream.open(filePath, mode);
    
    return fileStream.is_open();
}

void File::close() 
{
    if (fileStream.is_open()) {
        fileStream.close();
    }
}

bool File::isOpen() const 
{
    return fileStream.is_open();
}

bool File::exists(const std::string& filePath) 
{
    struct stat buffer;
    return (stat(filePath.c_str(), &buffer) == 0);
}

std::string File::readAll() 
{
    if (!isOpen()) {
        return "";
    }
    
    // Save current position
    std::streampos currentPos = fileStream.tellg();
    
    // Go to beginning
    fileStream.seekg(0, std::ios::beg);
    
    // Read the whole file
    std::string contents((std::istreambuf_iterator<char>(fileStream)),
                          std::istreambuf_iterator<char>());
    
    // Restore position
    fileStream.seekg(currentPos);
    
    return contents;
}

std::vector<std::string> File::readLines() 
{
    std::vector<std::string> lines;
    std::string line;
    
    if (!isOpen()) {
        return lines;
    }
    
    // Save current position
    std::streampos currentPos = fileStream.tellg();
    
    // Go to beginning
    fileStream.seekg(0, std::ios::beg);
    
    // Read lines
    while (std::getline(fileStream, line)) {
        lines.push_back(line);
    }
    
    // Restore position
    fileStream.clear(); // Clear EOF flag
    fileStream.seekg(currentPos);
    
    return lines;
}

bool File::readLine(std::string& line) 
{
    if (!isOpen()) {
        return false;
    }
    
    return static_cast<bool>(std::getline(fileStream, line));
}

bool File::write(const std::string& text) 
{
    if (!isOpen()) {
        return false;
    }
    
    fileStream << text;
    return !fileStream.fail();
}

bool File::writeLine(const std::string& line) 
{
    return write(line + "\n");
}

const std::string& File::getPath() const 
{
    return filePath;
}

size_t File::size() 
{
    if (!isOpen()) {
        return 0;
    }
    
    // Save current position
    std::streampos currentPos = fileStream.tellg();
    
    // Get file size
    fileStream.seekg(0, std::ios::end);
    size_t fileSize = fileStream.tellg();
    
    // Restore position
    fileStream.seekg(currentPos);
    
    return fileSize;
}

std::string File::getFileName(const std::string& filePath) 
{
    size_t pos = filePath.find_last_of("/\\");
    if (pos == std::string::npos) {
        return filePath;
    }
    
    return filePath.substr(pos + 1);
}

std::string File::getExtension(const std::string& filePath) 
{
    std::string fileName = getFileName(filePath);
    size_t pos = fileName.find_last_of('.');
    if (pos == std::string::npos) {
        return "";
    }
    
    return fileName.substr(pos + 1);
}

std::string File::getDirectory(const std::string& filePath) 
{
    size_t pos = filePath.find_last_of("/\\");
    if (pos == std::string::npos) {
        return "";
    }
    
    return filePath.substr(0, pos);
}

} // namespace utils
} // namespace gangsta