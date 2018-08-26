#include <fstream>

class Logger : public std::ofstream {
friend Logger& operator <<(Logger& logger, const char* str);
friend Logger& operator <<(Logger& logger, int str);

public:
    Logger(std::string log_fname) : std::ofstream(log_fname, std::ios::app) {}
    ~Logger() {
        this->write(buf.c_str(), buf.size());
        this->close();
    }
private:
    std::string buf;
};

Logger& operator <<(Logger& logger, const char* str) {
    logger.buf = logger.buf + std::string(str);
    return logger;
}

Logger& operator <<(Logger& logger, const int str) {
    logger.buf = logger.buf + std::to_string(str);
    return logger;
}
